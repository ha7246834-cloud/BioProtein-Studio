from pathlib import Path

p = Path('pages/2_Gene_Structure_Domain_Motif.py')
s = p.read_text()

old = "from modules.gdm_cdd_meme import run_cdd, collapse_domains, domain_qc, meme_ready, run_meme, motif_qc, parse_cdd, parse_meme_xml"
new = "from modules.gdm_cdd_meme import run_cdd, collapse_domains, domain_qc, meme_ready, run_meme, motif_qc, parse_cdd, parse_meme_xml, CDD_BATCH_SIZE, CDD_MAX_SEQUENCES, CLOUD_MEME_MAX_SEQUENCES, CLOUD_MEME_MAX_RESIDUES"
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit('CDD/MEME import line not found')

old = "from modules.gdm_phylogeny import build_phylogeny, external_phylogeny_ready, publication_phylogeny_ready, phylogeny_tool_status"
new = "from modules.gdm_phylogeny import build_phylogeny, external_phylogeny_ready, publication_phylogeny_ready, phylogeny_tool_status, AUTO_IQTREE_MAX_SEQUENCES, AUTO_IQTREE_MAX_RESIDUES, CLOUD_PUBLICATION_MAX_SEQUENCES, CLOUD_PUBLICATION_MAX_RESIDUES"
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit('phylogeny import line not found')

anchor = """    with s4:
        model = st.selectbox('MEME model', ['zoops', 'oops', 'anr'], format_func=lambda x: {'zoops': 'ZOOPS — zero/one', 'oops': 'OOPS — exactly one', 'anr': 'ANR — any number'}[x])

    if st.button('🚀 Run complete analysis', type='primary', width='stretch'):
"""

block = """    with s4:
        model = st.selectbox('MEME model', ['zoops', 'oops', 'anr'], format_func=lambda x: {'zoops': 'ZOOPS — zero/one', 'oops': 'OOPS — exactly one', 'anr': 'ANR — any number'}[x])

    # Lightweight preflight: show the execution plan before a potentially long run.
    preview_text = file_text(up) if up else paste.strip()
    if preview_text:
        try:
            if mode.startswith('Protein'):
                preview_proteins = parse_fasta(preview_text, 'protein')
            else:
                preview_cds = parse_fasta(preview_text, 'dna')
                preview_proteins, _preview_qc = translate_cds(preview_cds)

            preview_n = len(preview_proteins)
            preview_residues = sum(
                len(str(seq).replace('-', '').replace('*', ''))
                for seq in preview_proteins.values()
            )
            preview_batches = max(1, (preview_n + CDD_BATCH_SIZE - 1) // CDD_BATCH_SIZE)
            shared_cloud = Path('/mount/src').exists()

            with st.expander('🧪 Run preflight & resource plan', expanded=True):
                m1, m2, m3 = st.columns(3)
                m1.metric('Proteins', preview_n)
                m2.metric('Total residues', f'{preview_residues:,}')
                m3.metric('CDD batches', preview_batches if do_cdd else 'Off')

                if tree_up:
                    tree_plan = 'Uploaded Newick — no tree inference required'
                    tree_level = 'success'
                elif not auto_tree:
                    tree_plan = 'Phylogeny disabled'
                    tree_level = 'info'
                elif phylo_mode == 'auto':
                    if preview_n <= AUTO_IQTREE_MAX_SEQUENCES and preview_residues <= AUTO_IQTREE_MAX_RESIDUES:
                        tree_plan = 'Auto → MAFFT + IQ-TREE publication-oriented inference'
                    else:
                        tree_plan = 'Auto → MAFFT + FastTree cloud-safe screening'
                    tree_level = 'success'
                elif phylo_mode == 'publication' and shared_cloud and (
                    preview_n > CLOUD_PUBLICATION_MAX_SEQUENCES
                    or preview_residues > CLOUD_PUBLICATION_MAX_RESIDUES
                ):
                    tree_plan = (
                        'Publication IQ-TREE is too large for shared Cloud and will be blocked. '
                        'Use Auto here, or run Publication locally/HPC and upload Newick.'
                    )
                    tree_level = 'warning'
                else:
                    tree_plan = {
                        'publication': 'Publication → MAFFT + IQ-TREE',
                        'fasttree': 'Fast screening → MAFFT + FastTree',
                        'nj': 'Internal Neighbor-Joining screening',
                    }.get(phylo_mode, str(phylo_mode))
                    tree_level = 'success'

                getattr(st, tree_level)(f'**Phylogeny plan:** {tree_plan}')

                if structure_up:
                    st.info('**Gene structure plan:** use supplied exon/CDS annotation.')
                elif cds_up and gen_up:
                    st.info('**Gene structure plan:** exact CDS/genomic validation first; EMBOSS est2genome only where spliced alignment is needed.')
                elif auto_reference and (reference_taxon.strip() or reference_accession.strip()):
                    thread_note = ' (Cloud capped at 2 threads)' if shared_cloud else ''
                    st.info(f'**Gene structure plan:** NCBI reference + miniprot mapping{thread_note}.')
                else:
                    st.warning('**Gene structure plan:** no complete genomic-reference route is currently specified; gene structure may remain unavailable.')

                if do_cdd:
                    if preview_n > CDD_MAX_SEQUENCES:
                        st.warning(
                            f'**CDD:** {preview_n} proteins exceeds the in-app limit of {CDD_MAX_SEQUENCES}. '
                            'Split the family into smaller CDD jobs.'
                        )
                    else:
                        st.info(
                            f'**CDD:** all {preview_n} proteins will be analysed in {preview_batches} '
                            f'NCBI batch(es) of up to {CDD_BATCH_SIZE}; no subsampling.'
                        )

                if do_meme:
                    meme_too_large = (
                        preview_n > CLOUD_MEME_MAX_SEQUENCES
                        or preview_residues > CLOUD_MEME_MAX_RESIDUES
                    )
                    if shared_cloud and meme_too_large:
                        st.warning(
                            f'**MEME:** shared-Cloud safety guard will skip this full-family MEME run '
                            f'({preview_n} proteins; {preview_residues:,} aa). Nothing will be subsampled. '
                            'Run MEME locally/HPC and import MEME XML for full-family motifs.'
                        )
                    else:
                        st.success('**MEME:** input is within the configured run-safety envelope.')

                st.caption(
                    'Preflight does not alter or subsample your sequences. It only predicts the route and resource safeguards.'
                )
        except Exception as e:
            st.warning(f'Preflight could not validate the current FASTA yet: {e}')

    if st.button('🚀 Run complete analysis', type='primary', width='stretch'):
"""

if anchor not in s:
    if "🧪 Run preflight & resource plan" in s:
        print('already patched')
    else:
        raise SystemExit('preflight anchor not found')
else:
    s = s.replace(anchor, block, 1)

p.write_text(s)
print('CLOUD_PREFLIGHT_PATCH_COMPLETE')
