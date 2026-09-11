from pathlib import Path

# Add an optional progress callback to CDD batching.
p = Path('modules/gdm_cdd_meme.py')
s = p.read_text()
old = "def run_cdd(proteins, evalue=0.01, timeout=360, batch_size=CDD_BATCH_SIZE):"
new = "def run_cdd(proteins, evalue=0.01, timeout=360, batch_size=CDD_BATCH_SIZE, progress=None):"
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit('run_cdd signature not found')

old = """    for i, batch in enumerate(batches, 1):
        if time.time() >= deadline:
"""
new = """    for i, batch in enumerate(batches, 1):
        if callable(progress):
            try:
                progress(i - 1, len(batches), f'Submitting NCBI CDD batch {i}/{len(batches)} ({len(batch)} proteins)')
            except Exception:
                pass
        if time.time() >= deadline:
"""
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit('CDD loop start not found')

old = """        raw_parts.append(
            f'# BioProtein Studio CDD batch {i}/{len(batches)}; sequences={len(batch)}; Search-ID={rid}\\n{raw}'
        )
        if i < len(batches):
"""
new = """        raw_parts.append(
            f'# BioProtein Studio CDD batch {i}/{len(batches)}; sequences={len(batch)}; Search-ID={rid}\\n{raw}'
        )
        if callable(progress):
            try:
                progress(i, len(batches), f'Completed NCBI CDD batch {i}/{len(batches)}')
            except Exception:
                pass
        if i < len(batches):
"""
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit('CDD loop completion block not found')
p.write_text(s)

# Surface an explicit run plan in the Streamlit progress bar.
p = Path('pages/2_Gene_Structure_Domain_Motif.py')
s = p.read_text()
old = """        r['proteins'] = proteins
        r['sequence_qc'] = protein_qc(proteins)
        r['tree_text'] = tree_text
        protein_txt = fasta_text(proteins)
        if (r['sequence_qc'].status == 'REVIEW').any():
            r['warnings'].append('Sequence-integrity QC flagged one or more proteins.')

        prog = st.progress(5, text='Phylogeny...')
"""
new = """        r['proteins'] = proteins
        r['sequence_qc'] = protein_qc(proteins)
        r['tree_text'] = tree_text
        protein_txt = fasta_text(proteins)
        family_nseq = len(proteins)
        family_residues = sum(len(str(x).replace('-', '').replace('*', '')) for x in proteins.values())
        if (r['sequence_qc'].status == 'REVIEW').any():
            r['warnings'].append('Sequence-integrity QC flagged one or more proteins.')

        if uploaded_tree_text:
            phylo_progress_text = f'Phylogeny: validating uploaded Newick for {family_nseq} proteins...'
        elif not auto_tree:
            phylo_progress_text = 'Phylogeny: automatic inference disabled.'
        elif phylo_mode == 'auto' and (family_nseq > 60 or family_residues > 40000):
            phylo_progress_text = (
                f'Phylogeny: large family ({family_nseq} proteins, {family_residues:,} aa) → '
                'MAFFT + FastTree cloud-safe screening...'
            )
        elif phylo_mode == 'publication':
            phylo_progress_text = f'Phylogeny: MAFFT + IQ-TREE publication inference for {family_nseq} proteins...'
        else:
            phylo_progress_text = f'Phylogeny: analysing {family_nseq} proteins...'

        prog = st.progress(5, text=phylo_progress_text)
"""
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit('phylogeny progress block not found')

old = """        prog.progress(15, text='Gene structure...')
"""
new = """        if auto_reference and (reference_taxon.strip() or reference_accession.strip()) and not (structure_up or (cds and genomic)):
            gene_progress_text = f'Gene structure: mapping {family_nseq} proteins to the selected NCBI reference...'
        elif cds and genomic:
            gene_progress_text = f'Gene structure: validating matching CDS/genomic inputs for {family_nseq} proteins...'
        else:
            gene_progress_text = 'Gene structure: resolving available reference evidence...'
        prog.progress(15, text=gene_progress_text)
"""
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit('gene progress block not found')

old = """        prog.progress(42, text='Conserved domains...')
        if do_cdd:
"""
new = """        cdd_batches = max(1, (family_nseq + 199) // 200)
        prog.progress(42, text=f'Conserved domains: NCBI CDD for {family_nseq} proteins ({cdd_batches} batch(es))...')
        if do_cdd:
"""
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit('CDD progress heading not found')

old = """            try:
                r['domains_raw'], r['cdd_rid'], r['cdd_raw'] = run_cdd(proteins, float(cdd_e))
                r['domains'] = collapse_domains(r['domains_raw'])
"""
new = """            try:
                def _cdd_progress(done, total, message):
                    pct = 42 + int(24 * done / max(1, total))
                    prog.progress(min(66, pct), text=f'CDD {done}/{total}: {message}')

                r['domains_raw'], r['cdd_rid'], r['cdd_raw'] = run_cdd(
                    proteins, float(cdd_e), progress=_cdd_progress
                )
                r['domains'] = collapse_domains(r['domains_raw'])
"""
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit('CDD run call not found')

old = """        prog.progress(68, text='MEME motifs...')
"""
new = """        prog.progress(
            68,
            text=f'MEME motifs: {family_nseq} proteins, up to {int(nm)} motifs (cloud safety limits apply)...'
        )
"""
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit('MEME progress block not found')

p.write_text(s)
print('LARGE_FAMILY_PROGRESS_PATCH_COMPLETE')
