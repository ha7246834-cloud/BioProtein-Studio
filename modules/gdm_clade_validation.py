from __future__ import annotations

import io
from collections import Counter, defaultdict
from typing import Iterable

import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO

from .gdm_circular import auto_topology_clades, support_text


CLOUD_ML_MAX_SEQUENCES = 300


def _parse_alignment(alignment_text: str) -> dict[str, str]:
    text = str(alignment_text or '').strip()
    if not text:
        return {}
    out = {}
    for rec in SeqIO.parse(io.StringIO(text), 'fasta'):
        name = str(rec.id).split()[0]
        seq = str(rec.seq).upper()
        if name and seq:
            out[name] = seq
    return out


def _pruned_tree(tree_text: str, keep: Iterable[str]):
    tree = Phylo.read(io.StringIO(str(tree_text or '').strip()), 'newick')
    wanted = set(str(x) for x in keep)
    for terminal in list(tree.get_terminals()):
        if terminal.name not in wanted:
            try:
                tree.prune(terminal)
            except Exception:
                pass
    return tree


def alignment_distance_matrix(alignment_text: str, order: Iterable[str] | None = None):
    """Return (names, p-distance matrix) from a protein multiple alignment.

    Only columns with residues in both sequences are compared. Gaps are ignored
    pairwise so the matrix represents sequence disagreement rather than gap
    padding introduced by the aligner.
    """
    seqs = _parse_alignment(alignment_text)
    if order is None:
        names = list(seqs)
    else:
        names = [str(x) for x in order if str(x) in seqs]
    n = len(names)
    if n < 2:
        raise ValueError('At least two aligned protein sequences are required for clade concordance.')

    lengths = {len(seqs[name]) for name in names}
    if len(lengths) != 1:
        raise ValueError('Aligned protein sequences must have equal alignment length.')

    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        a = seqs[names[i]]
        for j in range(i):
            b = seqs[names[j]]
            compared = 0
            mismatches = 0
            for ca, cb in zip(a, b):
                if ca in '-.' or cb in '-.':
                    continue
                compared += 1
                if ca != cb:
                    mismatches += 1
            d = 1.0 if compared == 0 else mismatches / compared
            matrix[i, j] = matrix[j, i] = float(d)
    return names, matrix


def kmedoids(distance: np.ndarray, k: int, max_iter: int = 100):
    """Deterministic PAM-like k-medoids clustering for a precomputed distance matrix."""
    d = np.asarray(distance, dtype=float)
    n = int(d.shape[0])
    if d.ndim != 2 or d.shape[1] != n:
        raise ValueError('Distance matrix must be square.')
    if not (2 <= int(k) <= n):
        raise ValueError('k must be between 2 and the number of sequences.')
    k = int(k)

    # Deterministic farthest-first initialization: central first medoid, then
    # repeatedly choose the point farthest from its closest selected medoid.
    medoids = [int(np.argmin(d.mean(axis=1)))]
    while len(medoids) < k:
        candidates = [i for i in range(n) if i not in medoids]
        nxt = max(candidates, key=lambda i: (float(np.min(d[i, medoids])), -i))
        medoids.append(int(nxt))

    labels = np.argmin(d[:, medoids], axis=1)
    for _ in range(max(1, int(max_iter))):
        new_medoids = list(medoids)
        for cluster in range(k):
            members = np.where(labels == cluster)[0]
            if len(members) == 0:
                candidates = [i for i in range(n) if i not in new_medoids]
                if candidates:
                    new_medoids[cluster] = max(
                        candidates,
                        key=lambda i: float(np.min(d[i, new_medoids])),
                    )
                continue
            sub = d[np.ix_(members, members)]
            costs = sub.sum(axis=1)
            new_medoids[cluster] = int(members[int(np.argmin(costs))])
        if new_medoids == medoids:
            break
        medoids = new_medoids
        labels = np.argmin(d[:, medoids], axis=1)

    return labels.astype(int), medoids


def silhouette_score_precomputed(distance: np.ndarray, labels) -> float:
    d = np.asarray(distance, dtype=float)
    labels = np.asarray(labels)
    n = len(labels)
    if n < 3 or len(set(labels.tolist())) < 2:
        return 0.0
    scores = []
    clusters = sorted(set(labels.tolist()))
    for i in range(n):
        same = np.where(labels == labels[i])[0]
        same = same[same != i]
        if len(same) == 0:
            scores.append(0.0)
            continue
        a = float(np.mean(d[i, same]))
        b_vals = []
        for c in clusters:
            if c == labels[i]:
                continue
            other = np.where(labels == c)[0]
            if len(other):
                b_vals.append(float(np.mean(d[i, other])))
        if not b_vals:
            scores.append(0.0)
            continue
        b = min(b_vals)
        denom = max(a, b)
        scores.append(0.0 if denom <= 0 else (b - a) / denom)
    return float(np.mean(scores))


def adjusted_rand_index(labels_a, labels_b) -> float:
    a = list(labels_a)
    b = list(labels_b)
    if len(a) != len(b):
        raise ValueError('Label vectors must have the same length.')
    n = len(a)
    if n < 2:
        return 1.0

    def c2(x):
        x = int(x)
        return x * (x - 1) / 2.0

    rows = Counter(a)
    cols = Counter(b)
    cells = Counter(zip(a, b))
    sum_cells = sum(c2(v) for v in cells.values())
    sum_rows = sum(c2(v) for v in rows.values())
    sum_cols = sum(c2(v) for v in cols.values())
    total = c2(n)
    if total == 0:
        return 1.0
    expected = (sum_rows * sum_cols) / total
    maximum = 0.5 * (sum_rows + sum_cols)
    denom = maximum - expected
    if abs(denom) < 1e-12:
        return 1.0 if a == b else 0.0
    return float((sum_cells - expected) / denom)


def _labels_from_mapping(names, mapping: dict[str, str]):
    label_ids = {}
    labels = []
    for name in names:
        lab = mapping.get(name)
        if lab is None:
            return None
        if lab not in label_ids:
            label_ids[lab] = len(label_ids)
        labels.append(label_ids[lab])
    return np.asarray(labels, dtype=int)


def clade_sequence_concordance(
    tree_text: str,
    alignment_text: str,
    target_groups: int | None = None,
    max_groups: int = 8,
):
    """Compare topology clades with independent unsupervised sequence clustering.

    The second line of evidence is deterministic k-medoids clustering of the
    MAFFT alignment p-distance matrix. It is intentionally described as
    *computational concordance*, not biological validation: true biological
    subgroup assignment still requires reference genes, support values and
    domain/functional evidence.
    """
    if not str(tree_text or '').strip():
        raise ValueError('Newick tree is required for clade concordance.')
    seqs = _parse_alignment(alignment_text)
    if not seqs:
        raise ValueError('Protein alignment is required for sequence-cluster concordance.')

    raw_tree = Phylo.read(io.StringIO(str(tree_text).strip()), 'newick')
    tree_order = [str(t.name) for t in raw_tree.get_terminals() if t.name]
    names = [name for name in tree_order if name in seqs]
    if len(names) < 4:
        raise ValueError('At least four tree-matched aligned proteins are required for ML concordance.')
    if len(names) > CLOUD_ML_MAX_SEQUENCES:
        raise ValueError(
            f'Unsupervised clade concordance is limited to {CLOUD_ML_MAX_SEQUENCES} proteins in the app. '
            'The circular tree can still be rendered for larger families without claiming ML validation.'
        )

    tree = _pruned_tree(tree_text, names)
    names, distance = alignment_distance_matrix(alignment_text, names)

    if target_groups is not None:
        candidates = [max(2, min(int(target_groups), len(names) - 1, 12))]
    else:
        candidates = list(range(2, min(int(max_groups), len(names) - 1, 12) + 1))

    evaluated = []
    seen_counts = set()
    for requested_k in candidates:
        group_defs, mapping = auto_topology_clades(tree, requested_k)
        actual_k = len(group_defs)
        if actual_k < 2 or actual_k in seen_counts:
            continue
        seen_counts.add(actual_k)
        top_labels = _labels_from_mapping(names, mapping)
        if top_labels is None or len(set(top_labels.tolist())) < 2:
            continue
        ml_labels, medoids = kmedoids(distance, actual_k)
        ari = adjusted_rand_index(top_labels, ml_labels)
        sil = silhouette_score_precomputed(distance, top_labels)
        ml_sil = silhouette_score_precomputed(distance, ml_labels)
        # Prefer topology partitions that are compact in sequence space and are
        # independently recovered by k-medoids. Negative ARI is not rewarded.
        quality = sil + 0.50 * max(0.0, ari) + 0.15 * max(0.0, ml_sil)
        evaluated.append({
            'requested_groups': requested_k,
            'groups': actual_k,
            'group_defs': group_defs,
            'mapping': mapping,
            'top_labels': top_labels,
            'ml_labels': ml_labels,
            'medoids': medoids,
            'ari': float(ari),
            'silhouette': float(sil),
            'ml_silhouette': float(ml_sil),
            'quality': float(quality),
        })

    if not evaluated:
        raise ValueError('Could not derive at least two comparable topology groups from this tree.')

    best = max(evaluated, key=lambda x: (x['quality'], x['ari'], x['silhouette'], -x['groups']))
    ml_labels = best['ml_labels']
    mapping = best['mapping']

    ml_by_name = {name: int(ml_labels[i]) for i, name in enumerate(names)}
    rows = []
    member_rows = []
    for label, clade in best['group_defs']:
        members = [str(t.name) for t in clade.get_terminals() if str(t.name) in ml_by_name]
        counts = Counter(ml_by_name[m] for m in members)
        dominant_cluster, dominant_n = counts.most_common(1)[0] if counts else (-1, 0)
        purity = dominant_n / max(1, len(members))

        idx = [names.index(m) for m in members]
        intra = []
        for ii in range(len(idx)):
            for jj in range(ii):
                intra.append(float(distance[idx[ii], idx[jj]]))
        mean_intra = float(np.mean(intra)) if intra else 0.0

        if purity >= 0.80 and best['ari'] >= 0.50:
            status = 'ML-CONCORDANT'
        elif purity >= 0.65:
            status = 'PARTIAL-CONCORDANCE'
        else:
            status = 'REVIEW'

        rows.append({
            'auto_clade': label,
            'members': len(members),
            'root_support': support_text(clade),
            'dominant_ml_cluster': f'Cluster {dominant_cluster + 1}' if dominant_cluster >= 0 else '',
            'sequence_cluster_purity': round(float(purity), 3),
            'mean_intra_p_distance': round(mean_intra, 4),
            'evidence_status': status,
            'member_ids': ';'.join(members),
        })
        for m in members:
            member_rows.append({
                'gene': m,
                'auto_clade': label,
                'ml_cluster': f'Cluster {ml_by_name[m] + 1}',
                'clade_cluster_match': ml_by_name[m] == dominant_cluster,
            })

    ari = best['ari']
    sil = best['silhouette']
    if ari >= 0.75 and sil >= 0.25:
        overall = 'STRONG COMPUTATIONAL CONCORDANCE'
    elif ari >= 0.50 and sil >= 0.10:
        overall = 'MODERATE COMPUTATIONAL CONCORDANCE'
    else:
        overall = 'REVIEW: TOPOLOGY AND SEQUENCE CLUSTERING ARE NOT STRONGLY CONCORDANT'

    summary = pd.DataFrame([{
        'recommended_clades': int(best['groups']),
        'adjusted_rand_index': round(float(ari), 3),
        'topology_silhouette': round(float(sil), 3),
        'ml_silhouette': round(float(best['ml_silhouette']), 3),
        'status': overall,
        'interpretation': 'Unsupervised k-medoids sequence-distance concordance; not biological subgroup validation.',
    }])

    return {
        'recommended_groups': int(best['groups']),
        'summary': summary,
        'clade_table': pd.DataFrame(rows),
        'member_table': pd.DataFrame(member_rows),
        'ari': float(ari),
        'silhouette': float(sil),
        'status': overall,
    }
