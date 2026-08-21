from __future__ import annotations
import io
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict

import pandas as pd
from Bio import Phylo
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from .gdm_common import fasta_text


def _which_fasttree():
    return shutil.which('FastTree') or shutil.which('fasttree') or shutil.which('fasttreeMP')


def _which_iqtree():
    return shutil.which('iqtree3') or shutil.which('iqtree2') or shutil.which('iqtree')


def external_phylogeny_ready():
    return shutil.which('mafft') is not None and _which_fasttree() is not None


def publication_phylogeny_ready():
    return shutil.which('mafft') is not None and _which_iqtree() is not None


def phylogeny_tool_status():
    return {
        'mafft': shutil.which('mafft'),
        'fasttree': _which_fasttree(),
        'iqtree': _which_iqtree(),
        'external_ready': external_phylogeny_ready(),
        'publication_ready': publication_phylogeny_ready(),
    }


def _clean_proteins(proteins: Dict[str, str]):
    out = {}
    for gene, seq in proteins.items():
        s = str(seq).upper().replace('*', '').replace('-', '')
        if s:
            out[str(gene)] = s
    return out


def _safe_records(proteins):
    mapping = {}
    safe = {}
    for i, (gene, seq) in enumerate(proteins.items(), 1):
        sid = f'SEQ{i:06d}'
        safe[sid] = seq
        mapping[sid] = gene
    return safe, mapping


def _restore_tree_names(tree_text: str, mapping: dict[str, str]):
    tree = Phylo.read(io.StringIO(tree_text), 'newick')
    for terminal in tree.get_terminals():
        if terminal.name in mapping:
            terminal.name = mapping[terminal.name]
    buf = io.StringIO()
    Phylo.write(tree, buf, 'newick')
    return buf.getvalue().strip()


def _restore_alignment_headers(alignment_text: str, mapping: dict[str, str]):
    lines = []
    for line in alignment_text.splitlines():
        if line.startswith('>'):
            key = line[1:].split()[0]
            lines.append('>' + mapping.get(key, key))
        else:
            lines.append(line)
    return '\n'.join(lines).strip() + '\n'


def _run_mafft(proteins: Dict[str, str], td: Path, timeout=900):
    mafft = shutil.which('mafft')
    if not mafft:
        raise RuntimeError('MAFFT is not installed in the active environment.')
    safe, mapping = _safe_records(proteins)
    fasta = td / 'proteins.faa'
    fasta.write_text(fasta_text(safe))
    p = subprocess.run([mafft, '--auto', '--quiet', str(fasta)], capture_output=True, text=True, timeout=timeout)
    if p.returncode != 0 or not p.stdout.strip():
        raise RuntimeError(p.stderr.strip() or 'MAFFT failed to produce an alignment.')
    return p.stdout, mapping


def run_mafft_fasttree(proteins: Dict[str, str], timeout=900):
    proteins = _clean_proteins(proteins)
    if len(proteins) < 3:
        raise ValueError('At least 3 protein sequences are recommended for phylogenetic inference.')
    fasttree = _which_fasttree()
    if not fasttree:
        raise RuntimeError('FastTree is not installed in the active environment.')

    with tempfile.TemporaryDirectory(prefix='bps_tree_') as td0:
        td = Path(td0)
        aln_safe, mapping = _run_mafft(proteins, td, timeout=timeout)
        p2 = subprocess.run([fasttree, '-wag', '-gamma'], input=aln_safe, capture_output=True, text=True, timeout=timeout)
        if p2.returncode != 0 or '(' not in p2.stdout:
            raise RuntimeError(p2.stderr.strip() or 'FastTree failed to produce a Newick tree.')

    tree_text = _restore_tree_names(p2.stdout.strip(), mapping)
    alignment_text = _restore_alignment_headers(aln_safe, mapping)
    qc = pd.DataFrame([{
        'method': 'MAFFT + FastTree',
        'sequences': len(proteins),
        'model': 'WAG + CAT/Gamma rescaling',
        'support': 'FastTree local support where available',
        'status': 'PASS: approximate maximum-likelihood screening tree',
    }])
    return {
        'tree_text': tree_text,
        'alignment_text': alignment_text,
        'method': 'MAFFT + FastTree',
        'qc': qc,
        'warning': 'FastTree is retained as a fast screening method. Use Publication mode (IQ-TREE) for final inference when available.',
        'log_text': p2.stderr or '',
    }


def _parse_iqtree_model(text: str):
    pats = [
        r'Best-fit model according to BIC:\s*([^\s]+)',
        r'Best-fit model:\s*([^\s]+)',
        r'Model of substitution:\s*([^\n]+)',
    ]
    for pat in pats:
        m = re.search(pat, text, flags=re.I)
        if m:
            return m.group(1).strip()
    return 'ModelFinder best-fit model'


def run_mafft_iqtree(proteins: Dict[str, str], bootstrap=1000, alrt=1000, threads='AUTO', timeout=7200):
    proteins = _clean_proteins(proteins)
    if len(proteins) < 3:
        raise ValueError('At least 3 protein sequences are recommended for phylogenetic inference.')
    iqtree = _which_iqtree()
    if not iqtree:
        raise RuntimeError('IQ-TREE is not installed in the active environment.')

    with tempfile.TemporaryDirectory(prefix='bps_iqtree_') as td0:
        td = Path(td0)
        aln_safe, mapping = _run_mafft(proteins, td, timeout=min(timeout, 1800))
        aln = td / 'alignment.fasta'
        aln.write_text(aln_safe)
        prefix = td / 'bps_iqtree'
        cmd = [iqtree, '-s', str(aln), '-m', 'MFP', '-B', str(int(bootstrap)), '-alrt', str(int(alrt)), '-T', str(threads), '--prefix', str(prefix), '-redo']
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        treefile = Path(str(prefix) + '.treefile')
        report = Path(str(prefix) + '.iqtree')
        logfile = Path(str(prefix) + '.log')
        if p.returncode != 0 or not treefile.exists():
            raise RuntimeError((p.stderr or p.stdout or 'IQ-TREE failed to produce a tree.').strip())
        tree_safe = treefile.read_text(errors='replace').strip()
        report_text = report.read_text(errors='replace') if report.exists() else ''
        log_text = logfile.read_text(errors='replace') if logfile.exists() else (p.stdout + '\n' + p.stderr)

    tree_text = _restore_tree_names(tree_safe, mapping)
    alignment_text = _restore_alignment_headers(aln_safe, mapping)
    model = _parse_iqtree_model(report_text + '\n' + log_text)
    qc = pd.DataFrame([{
        'method': 'MAFFT + IQ-TREE',
        'sequences': len(proteins),
        'model': model,
        'support': f'SH-aLRT {int(alrt)} + ultrafast bootstrap {int(bootstrap)}',
        'status': 'PASS: publication-oriented maximum-likelihood inference; biological interpretation still requires support review',
    }])
    return {
        'tree_text': tree_text,
        'alignment_text': alignment_text,
        'method': 'MAFFT + IQ-TREE',
        'qc': qc,
        'warning': '',
        'log_text': log_text,
        'iqtree_report': report_text,
        'command': ' '.join(cmd),
    }


def _pairwise_distance(seq_a: str, seq_b: str, aligner: PairwiseAligner):
    aln = aligner.align(seq_a, seq_b)[0]
    idx = aln.indices
    cols = idx.shape[1]
    if cols == 0:
        return 1.0
    matches = 0
    for c in range(cols):
        ia, ib = int(idx[0, c]), int(idx[1, c])
        if ia >= 0 and ib >= 0 and seq_a[ia] == seq_b[ib]:
            matches += 1
    return max(0.0, min(1.0, 1.0 - matches / cols))


def run_nj_fallback(proteins: Dict[str, str], max_sequences=100):
    proteins = _clean_proteins(proteins)
    n = len(proteins)
    if n < 3:
        raise ValueError('At least 3 protein sequences are recommended for phylogenetic inference.')
    if n > max_sequences:
        raise ValueError(f'Internal NJ fallback is limited to {max_sequences} sequences. Install MAFFT + FastTree/IQ-TREE for larger families.')

    names = list(proteins)
    aligner = PairwiseAligner()
    try:
        aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
        aligner.open_gap_score = -10.0
        aligner.extend_gap_score = -0.5
    except Exception:
        aligner.match_score = 1.0
        aligner.mismatch_score = 0.0
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = -0.1

    lower = []
    for i, a in enumerate(names):
        row = []
        for j in range(i + 1):
            row.append(0.0 if i == j else _pairwise_distance(proteins[a], proteins[names[j]], aligner))
        lower.append(row)

    dm = DistanceMatrix(names=names, matrix=lower)
    tree = DistanceTreeConstructor().nj(dm)
    tree.rooted = False
    buf = io.StringIO()
    Phylo.write(tree, buf, 'newick')
    qc = pd.DataFrame([{
        'method': 'Internal NJ fallback',
        'sequences': n,
        'model': 'BLOSUM62 global pairwise distance',
        'support': 'No resampling support',
        'status': 'REVIEW: screening tree; use MAFFT + IQ-TREE for publication inference',
    }])
    return {
        'tree_text': buf.getvalue().strip(),
        'alignment_text': '',
        'method': 'Internal NJ fallback',
        'qc': qc,
        'warning': 'Phylogeny used the internal NJ fallback because external tools were unavailable or not selected.',
        'log_text': '',
    }


def build_phylogeny(proteins: Dict[str, str], mode='auto', bootstrap=1000, alrt=1000, threads='AUTO'):
    n = len(_clean_proteins(proteins))
    if n < 3:
        return {
            'tree_text': '', 'alignment_text': '', 'method': 'Not applicable',
            'qc': pd.DataFrame([{'method': 'Not applicable', 'sequences': n, 'model': '', 'support': '', 'status': 'REVIEW: at least 3 sequences are required for family phylogeny'}]),
            'warning': 'Phylogeny was skipped because fewer than 3 protein sequences were supplied.', 'log_text': ''
        }

    aliases = {'external': 'fasttree'}
    mode = aliases.get(mode, mode)
    if mode not in {'auto', 'publication', 'fasttree', 'nj'}:
        raise ValueError('Unknown phylogeny mode.')

    if mode == 'publication':
        if not publication_phylogeny_ready():
            raise RuntimeError('Publication mode requires MAFFT + IQ-TREE.')
        return run_mafft_iqtree(proteins, bootstrap=bootstrap, alrt=alrt, threads=threads)

    if mode == 'fasttree':
        if not external_phylogeny_ready():
            raise RuntimeError('FastTree mode requires MAFFT + FastTree.')
        return run_mafft_fasttree(proteins)

    if mode == 'auto':
        if publication_phylogeny_ready():
            return run_mafft_iqtree(proteins, bootstrap=bootstrap, alrt=alrt, threads=threads)
        if external_phylogeny_ready():
            return run_mafft_fasttree(proteins)
    return run_nj_fallback(proteins)
