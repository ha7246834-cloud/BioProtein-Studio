from pathlib import Path

# gdm_reference.py: cap miniprot threads on shared cloud.
p = Path('modules/gdm_reference.py')
s = p.read_text()

anchor = """from .gdm_common import fasta_text, norm_id


def datasets_ready() -> bool:
"""
replacement = """from .gdm_common import fasta_text, norm_id


def _is_shared_cloud() -> bool:
    return bool(os.environ.get('STREAMLIT_SHARING_MODE')) or Path('/mount/src').exists()


def _safe_miniprot_threads(threads: int) -> int:
    requested = max(1, int(threads))
    return min(requested, 2) if _is_shared_cloud() else requested


def datasets_ready() -> bool:
"""
if anchor in s:
    s = s.replace(anchor, replacement, 1)
elif "def _safe_miniprot_threads" not in s:
    raise SystemExit('reference helper anchor not found')

old = """        cmd = [
            'miniprot', '-I', '-t', str(max(1, int(threads))), '--gff-only',
            '--outn=2', '--outs=0.85', '--outc=0.50',
            ref['genome_fasta'], str(protein_fa)
        ]
        gff = _run(cmd, timeout=timeout)
"""
new = """        safe_threads = _safe_miniprot_threads(threads)
        cmd = [
            'miniprot', '-I', '-t', str(safe_threads), '--gff-only',
            '--outn=2', '--outs=0.85', '--outc=0.50',
            ref['genome_fasta'], str(protein_fa)
        ]
        effective_timeout = min(int(timeout), 1800) if _is_shared_cloud() else int(timeout)
        gff = _run(cmd, timeout=effective_timeout)
"""
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit('miniprot command block not found')
p.write_text(s)

# gdm_structure.py: exact contiguous CDS should not spawn est2genome.
p = Path('modules/gdm_structure.py')
s = p.read_text()
old = """def est2genome_pair(gene:str,cds:str,genomic:str,timeout=180):
    if not est2genome_ready():return _exact(gene,cds,genomic)
    with tempfile.TemporaryDirectory(prefix='bps_est2genome_') as td:
"""
new = """def est2genome_pair(gene:str,cds:str,genomic:str,timeout=180):
    # Fast exact-contiguous route is scientifically definitive for intronless
    # CDS and avoids launching one EMBOSS process per gene unnecessarily.
    try:
        return _exact(gene, cds, genomic)
    except RuntimeError:
        if not est2genome_ready():
            raise
    with tempfile.TemporaryDirectory(prefix='bps_est2genome_') as td:
"""
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit('est2genome_pair block not found')
p.write_text(s)

print('GENE_STRUCTURE_PERFORMANCE_PATCH_COMPLETE')
