from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from .gdm_common import fasta_text, norm_id


def datasets_ready() -> bool:
    return shutil.which('datasets') is not None


def miniprot_ready() -> bool:
    return shutil.which('miniprot') is not None


def auto_reference_ready() -> bool:
    return datasets_ready() and miniprot_ready()


def reference_tool_status() -> dict:
    return {
        'NCBI Datasets': shutil.which('datasets') or '',
        'miniprot': shutil.which('miniprot') or '',
    }


def _run(cmd, timeout=900, cwd=None):
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, cwd=cwd)
    if p.returncode:
        msg = (p.stderr or p.stdout or '').strip()
        raise RuntimeError(msg or f'Command failed: {" ".join(map(str, cmd))}')
    return p.stdout


def _flatten_records(text: str):
    records = []
    for line in (text or '').splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            obj = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(obj, dict) and isinstance(obj.get('reports'), list):
            records.extend(x for x in obj['reports'] if isinstance(x, dict))
        elif isinstance(obj, dict):
            records.append(obj)
    return records


def _assembly_rank(rec: dict):
    acc = str(rec.get('accession', ''))
    ai = rec.get('assemblyInfo') or rec.get('assembly_info') or {}
    refcat = str(ai.get('refseqCategory', ai.get('refseq_category', ''))).lower()
    level = str(ai.get('assemblyLevel', ai.get('assembly_level', ''))).lower()
    ref_score = 2 if 'reference' in refcat else 1 if 'representative' in refcat else 0
    level_score = {'complete genome': 4, 'complete': 4, 'chromosome': 3, 'scaffold': 2, 'contig': 1}.get(level, 0)
    source_score = 1 if acc.startswith('GCF_') else 0
    annotated = 1 if (rec.get('annotationInfo') or rec.get('annotation_info')) else 0
    return ref_score, level_score, annotated, source_score


def select_reference_assembly(taxon: str, assembly_accession: str = '', api_key: str = '') -> dict:
    if not datasets_ready():
        raise RuntimeError('NCBI Datasets CLI not found. Install ncbi-datasets-cli in the BioProtein Studio environment.')
    assembly_accession = (assembly_accession or '').strip()
    if assembly_accession:
        if not re.fullmatch(r'GC[AF]_\d+(?:\.\d+)?', assembly_accession, re.I):
            raise ValueError('Assembly accession should look like GCF_... or GCA_....')
        cmd = ['datasets', 'summary', 'genome', 'accession', assembly_accession, '--as-json-lines']
        if api_key:
            cmd += ['--api-key', api_key]
        recs = _flatten_records(_run(cmd, timeout=120))
        if not recs:
            raise RuntimeError(f'NCBI Datasets could not resolve assembly {assembly_accession}.')
        rec = recs[0]
        return _assembly_summary(rec, requested_taxon=taxon)

    taxon = (taxon or '').strip()
    if not taxon:
        raise ValueError('Species/taxon is required when no assembly accession is supplied.')

    common = [
        'datasets', 'summary', 'genome', 'taxon', taxon,
        '--annotated', '--exclude-atypical', '--limit', '25', '--as-json-lines'
    ]
    if api_key:
        common += ['--api-key', api_key]

    # Prefer an NCBI reference genome if one exists; otherwise rank annotated assemblies.
    try:
        recs = _flatten_records(_run(common + ['--reference'], timeout=180))
    except Exception:
        recs = []
    if not recs:
        recs = _flatten_records(_run(common, timeout=180))
    if not recs:
        raise RuntimeError(f'No annotated NCBI genome assembly found for taxon: {taxon}')
    rec = sorted(recs, key=_assembly_rank, reverse=True)[0]
    return _assembly_summary(rec, requested_taxon=taxon)


def _assembly_summary(rec: dict, requested_taxon: str = '') -> dict:
    ai = rec.get('assemblyInfo') or rec.get('assembly_info') or {}
    org = rec.get('organism') or {}
    annot = rec.get('annotationInfo') or rec.get('annotation_info') or {}
    return {
        'accession': rec.get('accession', ''),
        'organism': org.get('organismName', org.get('organism_name', requested_taxon)),
        'taxid': org.get('taxId', org.get('tax_id', '')),
        'assembly_name': ai.get('assemblyName', ai.get('assembly_name', '')),
        'assembly_level': ai.get('assemblyLevel', ai.get('assembly_level', '')),
        'refseq_category': ai.get('refseqCategory', ai.get('refseq_category', '')),
        'annotation_name': annot.get('name', ''),
        'annotation_release_date': annot.get('releaseDate', annot.get('release_date', '')),
    }


def _find_reference_files(root: Path):
    # NCBI Datasets genome packages place files under ncbi_dataset/data/<accession>/.
    all_files = [p for p in root.rglob('*') if p.is_file()]
    genome_candidates = []
    for p in all_files:
        low = p.name.lower()
        if low.endswith(('.fna', '.fa', '.fasta', '.fna.gz', '.fa.gz', '.fasta.gz')):
            score = 0
            if 'genomic' in low or low == 'genomic.fna':
                score += 20
            if 'cds' in low or 'rna' in low or 'protein' in low:
                score -= 30
            genome_candidates.append((score, p))
    if not genome_candidates:
        raise RuntimeError('Downloaded NCBI package did not contain a genome FASTA.')
    genome = sorted(genome_candidates, key=lambda x: (x[0], -len(str(x[1]))), reverse=True)[0][1]
    gffs = [p for p in all_files if p.name.lower().endswith(('.gff', '.gff3'))]
    gff = gffs[0] if gffs else None
    return genome, gff


def prepare_reference_genome(meta: dict, cache_dir: str = '', api_key: str = '') -> dict:
    if not datasets_ready():
        raise RuntimeError('NCBI Datasets CLI not found.')
    acc = str(meta.get('accession', '')).strip()
    if not acc:
        raise ValueError('Reference assembly metadata is missing an accession.')
    base = Path(cache_dir).expanduser() if cache_dir else Path.home() / '.cache' / 'bioprotein-studio' / 'references'
    outdir = base / acc
    outdir.mkdir(parents=True, exist_ok=True)
    marker = outdir / 'reference_manifest.json'

    try:
        genome, gff = _find_reference_files(outdir)
        if genome.exists():
            return {'genome_fasta': str(genome), 'gff3': str(gff) if gff else '', 'cache_dir': str(outdir), 'reused_cache': True, **meta}
    except Exception:
        pass

    with tempfile.TemporaryDirectory(prefix='bps_ncbi_ref_') as td:
        td = Path(td)
        zpath = td / f'{acc}.zip'
        cmd = [
            'datasets', 'download', 'genome', 'accession', acc,
            '--include', 'genome,gff3', '--filename', str(zpath), '--no-progressbar'
        ]
        if api_key:
            cmd += ['--api-key', api_key]
        _run(cmd, timeout=1800)
        if not zpath.exists():
            raise RuntimeError('NCBI Datasets finished without producing the genome package.')
        with zipfile.ZipFile(zpath) as z:
            z.extractall(outdir)

    genome, gff = _find_reference_files(outdir)
    marker.write_text(json.dumps({**meta, 'genome_fasta': str(genome), 'gff3': str(gff) if gff else ''}, indent=2))
    return {'genome_fasta': str(genome), 'gff3': str(gff) if gff else '', 'cache_dir': str(outdir), 'reused_cache': False, **meta}


def _attrs(text: str):
    out = {}
    for token in str(text).split(';'):
        token = token.strip()
        if not token:
            continue
        if '=' in token:
            k, v = token.split('=', 1)
            out[k.strip()] = v.strip()
    return out


def parse_miniprot_gff(text: str, proteins: Dict[str, str]):
    mrnas = []
    cds_rows = []
    for line in (text or '').splitlines():
        if not line or line.startswith('#'):
            continue
        x = line.rstrip().split('\t')
        if len(x) < 9:
            continue
        seqid, source, ftype, start, end, score, strand, phase, attr = x[:9]
        if source != 'miniprot' or ftype not in {'mRNA', 'CDS'}:
            continue
        try:
            start_i, end_i = int(start), int(end)
            score_f = float(score) if score not in {'.', ''} else np.nan
        except ValueError:
            continue
        a = _attrs(attr)
        if ftype == 'mRNA':
            target = a.get('Target', '').split()
            query = norm_id(target[0]) if target else ''
            qstart = int(target[1]) if len(target) >= 3 and target[1].isdigit() else np.nan
            qend = int(target[2]) if len(target) >= 3 and target[2].isdigit() else np.nan
            try:
                identity = 100.0 * float(a.get('Identity', 'nan'))
            except ValueError:
                identity = np.nan
            try:
                rank = int(a.get('Rank', '999'))
            except ValueError:
                rank = 999
            mrnas.append({
                'alignment_id': a.get('ID', ''), 'gene': query, 'contig': seqid,
                'start': min(start_i, end_i), 'end': max(start_i, end_i), 'strand': strand,
                'score': score_f, 'rank': rank, 'identity_pct': identity,
                'query_start_aa': qstart, 'query_end_aa': qend,
                'frameshifts': int(a.get('Frameshift', '0') or 0),
                'stop_codons': int(a.get('StopCodon', '0') or 0),
            })
        else:
            target = a.get('Target', '').split()
            query = norm_id(target[0]) if target else ''
            try:
                identity = 100.0 * float(a.get('Identity', 'nan'))
            except ValueError:
                identity = np.nan
            cds_rows.append({
                'alignment_id': a.get('Parent', ''), 'gene': query, 'contig': seqid,
                'start': min(start_i, end_i), 'end': max(start_i, end_i), 'strand': strand,
                'phase': phase, 'identity_pct': identity,
                'frameshifts': int(a.get('Frameshift', '0') or 0),
                'stop_codons': int(a.get('StopCodon', '0') or 0),
                'noncanonical_donor': a.get('Donor', ''), 'noncanonical_acceptor': a.get('Acceptor', ''),
            })
    mdf = pd.DataFrame(mrnas)
    cdf = pd.DataFrame(cds_rows)
    if mdf.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    # Some miniprot builds place Target only on mRNA. Back-fill CDS gene names from Parent.
    if not cdf.empty:
        id_to_gene = dict(zip(mdf.alignment_id, mdf.gene))
        cdf.loc[cdf.gene.eq(''), 'gene'] = cdf.loc[cdf.gene.eq(''), 'alignment_id'].map(id_to_gene).fillna('')

    qrows = []
    for gene, prot in proteins.items():
        sub = mdf[mdf.gene == gene].sort_values(['rank', 'score'], ascending=[True, False])
        if sub.empty:
            qrows.append({'gene': gene, 'mapping_status': 'UNMAPPED'})
            continue
        best = sub.iloc[0]
        qcov = np.nan
        if not pd.isna(best.query_start_aa) and not pd.isna(best.query_end_aa):
            qcov = 100.0 * (int(best.query_end_aa) - int(best.query_start_aa) + 1) / max(1, len(prot.rstrip('*')))
        secondary_ratio = np.nan
        if len(sub) > 1 and not pd.isna(best.score) and best.score != 0 and not pd.isna(sub.iloc[1].score):
            secondary_ratio = float(sub.iloc[1].score) / float(best.score)
        qrows.append({
            'gene': gene, 'mapping_status': 'MAPPED', 'alignment_id': best.alignment_id,
            'contig': best.contig, 'genomic_start': int(best.start), 'genomic_end': int(best.end),
            'strand': best.strand, 'alignment_score': best.score, 'identity_pct': best.identity_pct,
            'protein_coverage_pct': round(float(qcov), 2) if not pd.isna(qcov) else np.nan,
            'secondary_score_ratio': round(float(secondary_ratio), 4) if not pd.isna(secondary_ratio) else np.nan,
            'frameshifts': int(best.frameshifts), 'stop_codons': int(best.stop_codons),
        })
    return mdf, cdf, pd.DataFrame(qrows)


def _protein_identity(query: str, target: str) -> float:
    q = (query or '').rstrip('*').upper()
    t = (target or '').rstrip('*').upper()
    if not q or not t:
        return 0.0
    n = min(len(q), len(t))
    matches = 0
    compared = 0
    for a, b in zip(q[:n], t[:n]):
        if a in {'X', 'B', 'Z', 'J', 'U', 'O'} or b in {'X', 'B', 'Z', 'J', 'U', 'O'}:
            continue
        compared += 1
        matches += int(a == b)
    denom = max(len(q), len(t), compared, 1)
    return 100.0 * matches / denom


def _protein_coverage(query: str, translated: str) -> float:
    q = (query or '').rstrip('*')
    t = (translated or '').rstrip('*')
    if not q:
        return 0.0
    return min(100.0, 100.0 * min(len(q), len(t)) / max(1, len(q)))


def _sequence_index(genome_fasta: str):
    # SeqIO.index avoids loading an entire reference genome into RAM.
    return SeqIO.index(genome_fasta, 'fasta')


def _gff_attrs_any(text: str):
    """Parse GFF3 key=value and a minimal GTF key "value" attribute syntax."""
    out = _attrs(text)
    if out:
        return out
    for token in str(text).split(';'):
        token = token.strip()
        if not token or ' ' not in token:
            continue
        k, v = token.split(None, 1)
        out[k.strip()] = v.strip().strip('"')
    return out


def _overlap_len(a1, a2, b1, b2):
    return max(0, min(int(a2), int(b2)) - max(int(a1), int(b1)) + 1)


def _annotation_cds_candidates(gff_path: str, map_qc: pd.DataFrame, pad: int = 5000):
    """Stream an NCBI GFF3 once and retain CDS features near mapped loci."""
    if not gff_path or not Path(gff_path).exists() or map_qc is None or map_qc.empty:
        return pd.DataFrame()
    windows = {}
    for _, r in map_qc.iterrows():
        if r.get('mapping_status') != 'MAPPED':
            continue
        contig = str(r.get('contig', ''))
        windows.setdefault(contig, []).append((
            max(1, int(r.get('genomic_start', 1)) - int(pad)),
            int(r.get('genomic_end', 1)) + int(pad),
            str(r.get('strand', '.')),
        ))
    rows = []
    with open(gff_path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            x = line.rstrip('\n').split('\t')
            if len(x) < 9 or x[2] != 'CDS' or x[0] not in windows:
                continue
            try:
                s, e = int(x[3]), int(x[4])
            except ValueError:
                continue
            strand = x[6]
            if not any(strand == ws and _overlap_len(s, e, w1, w2) > 0 for w1, w2, ws in windows[x[0]]):
                continue
            a = _gff_attrs_any(x[8])
            parents = a.get('Parent') or a.get('transcript_id') or a.get('ID') or ''
            parent_list = [p.strip() for p in str(parents).split(',') if p.strip()]
            for parent in parent_list:
                rows.append({
                    'transcript': parent, 'contig': x[0], 'start': min(s, e), 'end': max(s, e),
                    'strand': strand, 'phase': x[7], 'source': x[1],
                    'protein_id': a.get('protein_id', a.get('Name', '')),
                    'gene_attr': a.get('gene', a.get('gene_id', a.get('locus_tag', ''))),
                })
    return pd.DataFrame(rows)


def _extract_cds_from_segments(record, sub: pd.DataFrame, strand: str):
    ordered = sub.sort_values('start', ascending=(strand != '-'))
    parts = []
    for _, r in ordered.iterrows():
        frag = record.seq[int(r.start) - 1:int(r.end)]
        if strand == '-':
            frag = frag.reverse_complement()
        parts.append(str(frag).upper())
    cds = ''.join(parts)
    trim = cds[:len(cds) - (len(cds) % 3)] if len(cds) % 3 else cds
    translated = str(Seq(trim).translate(to_stop=False)).rstrip('*') if trim else ''
    return cds, translated, ordered


def _annotation_reconcile_gene(gene: str, protein: str, mq: pd.Series, ann: pd.DataFrame, record):
    if ann is None or ann.empty:
        return None
    contig = str(mq['contig'])
    strand = str(mq['strand'])
    gstart, gend = int(mq['genomic_start']), int(mq['genomic_end'])
    candidates = []
    for tid, sub in ann[(ann.contig == contig) & (ann.strand == strand)].groupby('transcript'):
        cstart, cend = int(sub.start.min()), int(sub.end.max())
        ov = _overlap_len(gstart, gend, cstart, cend)
        if ov <= 0:
            continue
        reciprocal = 100.0 * ov / max(1, min(gend - gstart + 1, cend - cstart + 1))
        if reciprocal < 50:
            continue
        cds, translated, ordered = _extract_cds_from_segments(record, sub, strand)
        ident = _protein_identity(protein, translated)
        cov = _protein_coverage(protein, translated)
        len_delta = abs(len((protein or '').rstrip('*')) - len(translated))
        score = ident * 2.0 + cov + min(100.0, reciprocal) - min(50.0, len_delta)
        candidates.append({
            'transcript': tid, 'segments': ordered, 'cds': cds, 'translated': translated,
            'translation_identity_pct': ident, 'protein_coverage_pct': cov,
            'overlap_pct': reciprocal, 'score': score, 'span_start': cstart, 'span_end': cend,
        })
    if not candidates:
        return None
    candidates.sort(key=lambda x: (x['score'], x['translation_identity_pct'], x['protein_coverage_pct']), reverse=True)
    best = candidates[0]
    second = candidates[1] if len(candidates) > 1 else None
    best['ambiguous'] = bool(
        second and second['translation_identity_pct'] >= 95 and
        abs(best['score'] - second['score']) < 2.0
    )
    if best['translation_identity_pct'] < 95 or best['protein_coverage_pct'] < 95:
        return None
    return best


def _canonical_one_intron_rescue(protein: str, genomic_oriented: str, max_intron: int = 5000):
    """Conservative rescue for a single canonical GT-AG intron when annotation cannot confirm it.

    This is never promoted to publication PASS by itself; it is returned as REVIEW evidence.
    """
    seq = str(genomic_oriented or '').upper()
    q = (protein or '').rstrip('*').upper()
    if len(seq) < 30 or not q or len(seq) > 50000:
        return None
    expected = len(q) * 3
    donors = [i for i in range(1, len(seq) - 3) if seq[i:i+2] == 'GT']
    acceptors = [j for j in range(3, len(seq)) if seq[j-2:j] == 'AG']
    candidates = []
    for d in donors:
        for a in acceptors:
            intron_len = a - d
            if intron_len < 20 or intron_len > max_intron:
                continue
            spliced = seq[:d] + seq[a:]
            # A complete CDS may include or omit the terminal stop codon.
            if len(spliced) % 3 != 0:
                continue
            if abs(len(spliced) - expected) > max(90, int(expected * 0.15)) and abs(len(spliced) - (expected + 3)) > max(90, int(expected * 0.15)):
                continue
            translated = str(Seq(spliced).translate(to_stop=False)).rstrip('*')
            ident = _protein_identity(q, translated)
            cov = _protein_coverage(q, translated)
            internal_stops = translated[:-1].count('*') if translated else 0
            score = ident * 2 + cov - 20 * internal_stops - abs(len(translated) - len(q))
            if ident >= 95 and cov >= 95 and internal_stops == 0:
                candidates.append({
                    'donor0': d, 'acceptor0': a, 'intron_length': intron_len, 'cds': spliced,
                    'translated': translated, 'translation_identity_pct': ident,
                    'protein_coverage_pct': cov, 'score': score,
                })
    if not candidates:
        return None
    candidates.sort(key=lambda x: (x['score'], x['translation_identity_pct']), reverse=True)
    best = candidates[0]
    second = candidates[1] if len(candidates) > 1 else None
    best['ambiguous'] = bool(second and abs(best['score'] - second['score']) < 2.0)
    return best


def _canonical_splice_pct(rows_for_gene, genomic_oriented: str):
    if not rows_for_gene or len(rows_for_gene) < 2:
        return np.nan
    rows = sorted((int(a), int(b)) for a, b, *_ in rows_for_gene)
    seq = str(genomic_oriented or '').upper()
    checks = []
    for (_, e1), (s2, _) in zip(rows, rows[1:]):
        intron = seq[e1:s2 - 1]
        if len(intron) < 2:
            checks.append(False)
        else:
            checks.append(intron.startswith('GT') and intron.endswith('AG'))
    if not checks:
        return np.nan
    return 100.0 * sum(checks) / len(checks)


def build_structures_from_mapping(proteins: Dict[str, str], genome_fasta: str, gff_text: str, annotation_gff_path: str = ''):
    mdf, cdf, map_qc = parse_miniprot_gff(gff_text, proteins)
    if map_qc.empty:
        raise RuntimeError('miniprot produced no parseable mappings.')
    # Clamp a miniprot Target-derived coverage at 100%; terminal stop coordinates can otherwise yield 100.x%.
    if 'protein_coverage_pct' in map_qc.columns:
        map_qc['protein_coverage_pct'] = pd.to_numeric(map_qc['protein_coverage_pct'], errors='coerce').clip(upper=100.0)
    ann = _annotation_cds_candidates(annotation_gff_path, map_qc) if annotation_gff_path else pd.DataFrame()
    seqs = _sequence_index(genome_fasta)
    structure_rows = []
    qc_rows = []
    bundles = {}
    reconciliation_rows = []
    try:
        for _, mq in map_qc.iterrows():
            gene = mq['gene']
            if mq.get('mapping_status') != 'MAPPED':
                qc_rows.append({'gene': gene, 'method': 'miniprot protein-to-genome', 'status': 'REVIEW: unmapped'})
                continue
            aid = mq['alignment_id']
            sub = cdf[cdf.alignment_id == aid].copy() if not cdf.empty else pd.DataFrame()
            if sub.empty:
                qc_rows.append({'gene': gene, 'method': 'miniprot protein-to-genome', 'status': 'REVIEW: no CDS segments'})
                continue
            contig = str(mq['contig'])
            if contig not in seqs:
                qc_rows.append({'gene': gene, 'method': 'miniprot protein-to-genome', 'status': 'REVIEW: contig missing from FASTA'})
                continue
            strand = str(mq['strand'])
            gstart, gend = int(mq['genomic_start']), int(mq['genomic_end'])
            record = seqs[contig]
            miniprot_genomic = record.seq[gstart - 1:gend]
            if strand == '-':
                miniprot_genomic = miniprot_genomic.reverse_complement()

            # First reconstruct the raw miniprot CDS so we can decide whether annotation improves it.
            raw_cds, raw_translated, raw_ordered = _extract_cds_from_segments(record, sub, strand)
            raw_trans_ident = _protein_identity(proteins.get(gene, ''), raw_translated)
            raw_cov = _protein_coverage(proteins.get(gene, ''), raw_translated)

            chosen_source = 'miniprot'
            chosen_segments = raw_ordered
            chosen_cds = raw_cds
            chosen_translated = raw_translated
            chosen_span_start, chosen_span_end = gstart, gend
            annotation_transcript = ''
            annotation_overlap = np.nan
            computational_rescue = False

            # Reconcile against the reference annotation. This is especially important when miniprot
            # aligns through a short intron as one low-identity CDS block.
            recon = _annotation_reconcile_gene(gene, proteins.get(gene, ''), mq, ann, record)
            if recon is not None and (
                raw_trans_ident < 95 or len(raw_ordered) == 1 or recon['translation_identity_pct'] >= raw_trans_ident - 0.25
            ):
                chosen_source = 'NCBI annotation reconciliation'
                chosen_segments = recon['segments']
                chosen_cds = recon['cds']
                chosen_translated = recon['translated']
                chosen_span_start, chosen_span_end = int(recon['span_start']), int(recon['span_end'])
                annotation_transcript = str(recon['transcript'])
                annotation_overlap = float(recon['overlap_pct'])
                reconciliation_rows.append({
                    'gene': gene, 'reference_contig': contig, 'transcript': annotation_transcript,
                    'miniprot_translation_identity_pct': round(raw_trans_ident, 2),
                    'annotation_translation_identity_pct': round(float(recon['translation_identity_pct']), 2),
                    'annotation_protein_coverage_pct': round(float(recon['protein_coverage_pct']), 2),
                    'annotation_overlap_pct': round(annotation_overlap, 2),
                    'annotation_exons': len(chosen_segments),
                    'ambiguous_annotation': bool(recon.get('ambiguous', False)),
                    'decision': 'USED' if not recon.get('ambiguous', False) else 'USED_WITH_REVIEW',
                })
            elif raw_trans_ident < 95 and len(raw_ordered) == 1:
                rescue = _canonical_one_intron_rescue(proteins.get(gene, ''), str(miniprot_genomic))
                if rescue is not None:
                    computational_rescue = True
                    chosen_source = 'canonical one-intron rescue'
                    chosen_cds = rescue['cds']
                    chosen_translated = rescue['translated']
                    # Relative, transcript-oriented segments; these are appropriate for structure plotting.
                    L = len(str(miniprot_genomic))
                    d, a = int(rescue['donor0']), int(rescue['acceptor0'])
                    chosen_segments = pd.DataFrame([
                        {'start': 1, 'end': d, 'identity_pct': np.nan},
                        {'start': a + 1, 'end': L, 'identity_pct': np.nan},
                    ])
                    chosen_span_start, chosen_span_end = gstart, gend
                    reconciliation_rows.append({
                        'gene': gene, 'reference_contig': contig, 'transcript': '',
                        'miniprot_translation_identity_pct': round(raw_trans_ident, 2),
                        'annotation_translation_identity_pct': np.nan,
                        'annotation_protein_coverage_pct': np.nan,
                        'annotation_overlap_pct': np.nan,
                        'annotation_exons': 2,
                        'ambiguous_annotation': bool(rescue.get('ambiguous', False)),
                        'decision': 'COMPUTATIONAL_RESCUE_REVIEW',
                    })

            # Build relative exon/CDS positions. NCBI/miniprot genomic segments are converted to
            # transcript-oriented coordinates; rescue segments are already relative.
            rows_for_gene = []
            if computational_rescue:
                for _, r in chosen_segments.iterrows():
                    rows_for_gene.append((int(r.start), int(r.end), np.nan))
                genomic_oriented = str(miniprot_genomic).upper()
                ref_start, ref_end = gstart, gend
            else:
                ref_start, ref_end = chosen_span_start, chosen_span_end
                genomic_seq = record.seq[ref_start - 1:ref_end]
                if strand == '-':
                    genomic_seq = genomic_seq.reverse_complement()
                genomic_oriented = str(genomic_seq).upper()
                for _, r in chosen_segments.iterrows():
                    if strand == '-':
                        rel_start = ref_end - int(r.end) + 1
                        rel_end = ref_end - int(r.start) + 1
                    else:
                        rel_start = int(r.start) - ref_start + 1
                        rel_end = int(r.end) - ref_start + 1
                    rows_for_gene.append((rel_start, rel_end, r.get('identity_pct', np.nan)))

            trans_ident = _protein_identity(proteins.get(gene, ''), chosen_translated)
            protein_cov = _protein_coverage(proteins.get(gene, ''), chosen_translated)
            exons = len(rows_for_gene)
            introns = max(0, exons - 1)
            ambiguous_locus = not pd.isna(mq.secondary_score_ratio) and float(mq.secondary_score_ratio) >= 0.95
            ambiguous_annotation = bool(recon.get('ambiguous', False)) if recon is not None and chosen_source.startswith('NCBI') else False
            frame_ok = len(chosen_cds) % 3 == 0

            if chosen_source.startswith('NCBI'):
                pass_status = protein_cov >= 95 and trans_ident >= 95 and not ambiguous_locus and not ambiguous_annotation and frame_ok
                status = 'PASS' if pass_status else 'REVIEW'
                method = 'NCBI GFF annotation reconciliation after miniprot locus mapping'
            elif computational_rescue:
                status = 'REVIEW'
                method = 'miniprot locus + canonical one-intron computational rescue'
            else:
                pass_status = (
                    float(mq.get('protein_coverage_pct', 0) or 0) >= 95
                    and float(mq.get('identity_pct', 0) or 0) >= 95
                    and trans_ident >= 95
                    and int(mq.get('frameshifts', 0) or 0) == 0
                    and not ambiguous_locus
                    and frame_ok
                )
                status = 'PASS' if pass_status else 'REVIEW'
                method = 'miniprot protein-to-genome'

            reasons = []
            if protein_cov < 95: reasons.append('coverage<95%')
            if chosen_source == 'miniprot' and float(mq.get('identity_pct', 0) or 0) < 95: reasons.append('alignment identity<95%')
            if trans_ident < 95: reasons.append('translated CDS identity<95%')
            if int(mq.get('frameshifts', 0) or 0) > 0: reasons.append('frameshift')
            if ambiguous_locus: reasons.append('ambiguous locus')
            if ambiguous_annotation: reasons.append('ambiguous annotation transcript')
            if not frame_ok: reasons.append('CDS frame')
            if computational_rescue: reasons.append('computational splice rescue requires review')
            status_text = status if not reasons else f'{status}: {", ".join(reasons)}'

            for rel_start, rel_end, ident in rows_for_gene:
                structure_rows.append({
                    'gene': gene, 'feature': 'CDS', 'start': int(rel_start), 'end': int(rel_end),
                    'identity_pct': ident, 'strand': strand,
                    'source': chosen_source, 'validation_status': status,
                    'annotation_transcript': annotation_transcript,
                })

            qc_rows.append({
                'gene': gene, 'method': method, 'reference_contig': contig,
                'reference_start': int(ref_start), 'reference_end': int(ref_end), 'strand': strand,
                'exons': exons, 'introns': introns,
                'protein_coverage_pct': round(float(protein_cov), 2),
                'alignment_identity_pct': min(100.0, float(mq.get('identity_pct', np.nan))) if not pd.isna(mq.get('identity_pct', np.nan)) else np.nan,
                'translation_identity_pct': round(float(trans_ident), 2),
                'secondary_score_ratio': mq.get('secondary_score_ratio', np.nan),
                'frameshifts': mq.get('frameshifts', 0),
                'canonical_splice_pct': round(_canonical_splice_pct(rows_for_gene, genomic_oriented), 2) if introns else np.nan,
                'frame_multiple_of_3': frame_ok,
                'annotation_transcript': annotation_transcript,
                'annotation_overlap_pct': round(annotation_overlap, 2) if not pd.isna(annotation_overlap) else np.nan,
                'structure_source': chosen_source,
                'status': status_text,
            })
            bundles[gene] = {
                'cds': chosen_cds, 'genomic': genomic_oriented, 'translated_protein': chosen_translated,
                'reference_record': contig, 'reference_start': int(ref_start), 'reference_end': int(ref_end),
                'strand': strand, 'structure_source': chosen_source, 'annotation_transcript': annotation_transcript,
            }
    finally:
        seqs.close()

    structures = pd.DataFrame(structure_rows)
    if not structures.empty:
        structures = structures.sort_values(['gene', 'start', 'end']).reset_index(drop=True)
    return structures, pd.DataFrame(qc_rows), bundles, map_qc, mdf, cdf, pd.DataFrame(reconciliation_rows)


def auto_resolve_gene_structure(
    proteins: Dict[str, str], taxon: str = '', assembly_accession: str = '',
    api_key: str = '', cache_dir: str = '', threads: int = 4, timeout: int = 3600
):
    if not proteins:
        raise ValueError('Protein FASTA is required for automatic reference mapping.')
    if not miniprot_ready():
        raise RuntimeError('miniprot not found. Install miniprot in the BioProtein Studio environment.')
    meta = select_reference_assembly(taxon, assembly_accession, api_key)
    ref = prepare_reference_genome(meta, cache_dir, api_key)

    with tempfile.TemporaryDirectory(prefix='bps_miniprot_') as td:
        td = Path(td)
        protein_fa = td / 'queries.faa'
        protein_fa.write_text(fasta_text(proteins))
        cmd = [
            'miniprot', '-I', '-t', str(max(1, int(threads))), '--gff-only',
            '--outn=2', '--outs=0.85', '--outc=0.50',
            ref['genome_fasta'], str(protein_fa)
        ]
        gff = _run(cmd, timeout=timeout)
    structures, qc, bundles, map_qc, mrna_hits, cds_hits, reconciliation = build_structures_from_mapping(
        proteins, ref['genome_fasta'], gff, ref.get('gff3', '')
    )
    ref_table = pd.DataFrame([{k: v for k, v in ref.items() if k not in {'genome_fasta'}}])
    return {
        'structures': structures,
        'qc': qc,
        'bundles': bundles,
        'mapping_qc': map_qc,
        'mrna_hits': mrna_hits,
        'cds_hits': cds_hits,
        'annotation_reconciliation': reconciliation,
        'reference': ref,
        'reference_table': ref_table,
        'raw_gff': gff,
        'command': ' '.join(cmd),
    }
