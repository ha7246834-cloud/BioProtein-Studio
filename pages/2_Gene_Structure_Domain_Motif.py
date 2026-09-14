import json, platform
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from modules.gdm_common import parse_fasta, fasta_text, protein_qc, translate_cds, looks_ncbi_accession, newick_order, choose_order, zip_files
from modules.gdm_structure import est2genome_ready, gene_structure_batch, ncbi_structures, parse_gene_structure_annotation
from modules.gdm_cdd_meme import run_cdd, collapse_domains, domain_qc, meme_ready, run_meme, motif_qc, parse_cdd, parse_meme_xml, CDD_BATCH_SIZE, CDD_MAX_SEQUENCES, CLOUD_MEME_MAX_SEQUENCES, CLOUD_MEME_MAX_RESIDUES
from modules.gdm_phylogeny import build_phylogeny, external_phylogeny_ready, publication_phylogeny_ready, phylogeny_tool_status, AUTO_IQTREE_MAX_SEQUENCES, AUTO_IQTREE_MAX_RESIDUES, CLOUD_AUTO_IQTREE_MAX_SEQUENCES, CLOUD_AUTO_IQTREE_MAX_RESIDUES, CLOUD_PUBLICATION_MAX_SEQUENCES, CLOUD_PUBLICATION_MAX_RESIDUES
from modules.gdm_plot import gene_structure, missing_structure_figure, architecture, combined, phylogeny_figure, fig_bytes
from modules.gdm_style import STYLE_PRESETS, style_from_preset, assign_colors
from modules.gdm_reference import auto_resolve_gene_structure, auto_reference_ready, datasets_ready, miniprot_ready, reference_tool_status

st.set_page_config(page_title='Gene Structure, Domains & Motifs | BioProtein Studio', page_icon='🧬', layout='wide')
st.title('🧬 Gene Structure, Conserved Domains & Motifs')
st.caption('GSDS-like exon–intron mapping + NCBI CDD + MEME motif discovery + transparent scientific QC.')


def file_text(f):
    return f.getvalue().decode('utf-8', 'replace') if f else ''


def empty_result():
    return dict(
        sequence_qc=pd.DataFrame(), cds_qc=pd.DataFrame(), gene_structures=pd.DataFrame(), gene_qc=pd.DataFrame(), domains_raw=pd.DataFrame(),
        domains=pd.DataFrame(), domain_qc=pd.DataFrame(), motifs=pd.DataFrame(), motif_summary=pd.DataFrame(), motif_qc=pd.DataFrame(),
        errors=[], warnings=[], figures={}, raw_est={}, ncbi_bundles={}, auto_reference_bundles={}, tree_text='', alignment_text='', phylogeny_method='', phylogeny_qc=pd.DataFrame(), reference_qc=pd.DataFrame(), reference_meta=pd.DataFrame(), reference_mapping=pd.DataFrame(), annotation_reconciliation=pd.DataFrame(), miniprot_raw='', phylogeny_log='', iqtree_report='', phylogeny_command=''
    )


def validation_summary(r):
    rows = []
    if not r.get('phylogeny_qc', pd.DataFrame()).empty:
        q = r['phylogeny_qc'].iloc[0]
        p = str(q.get('status', '')).startswith('PASS')
        rows.append(dict(analysis='Phylogeny', criterion=str(q.get('method','')) + '; ' + str(q.get('support','')), PASS=int(p), REVIEW=int(not p)))
    if not r['gene_qc'].empty:
        rows.append(dict(analysis='Gene structure', criterion='CDS coverage >=95% and identity >=95%; splice/frame reviewed', PASS=int(r['gene_qc'].status.astype(str).str.startswith('PASS').sum()), REVIEW=int((~r['gene_qc'].status.astype(str).str.startswith('PASS')).sum())))
    if not r['domain_qc'].empty:
        p = r['domain_qc'].scientific_status.str.startswith('PASS', na=False)
        rows.append(dict(analysis='CDD domains', criterion='Specific hits prioritized; E<=1e-5 strong; raw full output retained', PASS=int(p.sum()), REVIEW=int((~p).sum())))
    if not r['motif_qc'].empty:
        p = r['motif_qc'].status == 'PASS'
        rows.append(dict(analysis='MEME motifs', criterion='MEME E<0.05 + >=50% family prevalence for PASS', PASS=int(p.sum()), REVIEW=int((~p).sum())))
    return pd.DataFrame(rows)


def filter_domains(domains, mode='specific'):
    if domains is None or domains.empty:
        return pd.DataFrame()
    d = domains.copy()
    if mode == 'specific':
        s = d[d.hit_type.str.contains('specific', case=False, na=False) & ~d.hit_type.str.contains('non', case=False, na=False)]
        return s if not s.empty else d
    if mode == 'pass':
        s = domain_qc(d)
        keep = s[s.scientific_status.str.startswith('PASS', na=False)][['gene', 'start', 'end', 'domain']]
        return d.merge(keep, on=['gene', 'start', 'end', 'domain']) if not keep.empty else d
    return d


def filter_motifs(motifs, motif_qc_df, mode='pass'):
    if motifs is None or motifs.empty:
        return pd.DataFrame()
    m = motifs.copy()
    if mode == 'pass' and motif_qc_df is not None and not motif_qc_df.empty:
        keep = set(motif_qc_df[motif_qc_df.status == 'PASS'].motif)
        f = m[m.motif.isin(keep)]
        return f if not f.empty else m
    return m


def package(r, protein_txt, cds_txt, gen_txt, tree_txt, params, structure_txt=''):
    files = {
        'inputs/proteins.fasta': protein_txt,
        'METHODS_AND_QC.txt': 'Gene structure: supplied annotation, NCBI reference CDS feature, EMBOSS est2genome, or automatic protein-to-reference-genome mapping with miniprot.\nAutomatic reference genomes are resolved/downloaded with NCBI Datasets CLI.\nDomains: NCBI Batch CD-Search full mode.\nMotifs: MEME Suite protein mode.\nPASS/REVIEW rules are visible in the app and must be biologically reviewed.\n'
    }
    if cds_txt:
        files['inputs/cds.fasta'] = cds_txt
    if gen_txt:
        files['inputs/genomic.fasta'] = gen_txt
    if tree_txt:
        files['phylogeny/tree.nwk'] = tree_txt
    if structure_txt:
        files['inputs/gene_structure_annotation.txt'] = structure_txt
    if r.get('alignment_text'):
        files['phylogeny/alignment.fasta'] = r['alignment_text']
    if isinstance(r.get('phylogeny_qc'), pd.DataFrame) and not r['phylogeny_qc'].empty:
        files['tables/phylogeny_qc.csv'] = r['phylogeny_qc'].to_csv(index=False)
    if r.get('phylogeny_log'):
        files['phylogeny/phylogeny.log.txt'] = r['phylogeny_log']
    if r.get('iqtree_report'):
        files['phylogeny/iqtree_report.txt'] = r['iqtree_report']
    if r.get('phylogeny_command'):
        files['phylogeny/command.txt'] = r['phylogeny_command']
    for k, n in [
        ('sequence_qc', 'tables/protein_sequence_qc.csv'), ('cds_qc', 'tables/cds_translation_qc.csv'), ('gene_structures', 'tables/gene_structures.csv'), ('gene_qc', 'tables/gene_structure_qc.csv'),
        ('domains_raw', 'tables/cdd_hits_full.csv'), ('domains', 'tables/cdd_hits_nonredundant.csv'), ('domain_qc', 'tables/domain_validation.csv'), ('motifs', 'tables/motif_sites.csv'), ('motif_summary', 'tables/motif_summary.csv'), ('motif_qc', 'tables/motif_validation.csv')
    ]:
        d = r.get(k)
        if isinstance(d, pd.DataFrame) and not d.empty:
            files[n] = d.to_csv(index=False)
    for g, t in r.get('raw_est', {}).items():
        files[f'raw/est2genome/{g}.txt'] = t
    if r.get('cdd_raw'):
        files['raw/cdd_output.txt'] = r['cdd_raw']
    if r.get('cdd_rid'):
        files['raw/cdd_search_id.txt'] = r['cdd_rid']
    if r.get('meme_xml'):
        files['raw/meme.xml'] = r['meme_xml']
    if r.get('meme_zip'):
        files['raw/meme_complete_output.zip'] = r['meme_zip']
    for g, b in r.get('ncbi_bundles', {}).items():
        files[f'reference/{g}.cds.fasta'] = f'>{g}\n{b["cds"]}\n'
        files[f'reference/{g}.genomic.fasta'] = f'>{g}|{b["reference_record"]}\n{b["genomic"]}\n'
    if isinstance(r.get('reference_meta'), pd.DataFrame) and not r['reference_meta'].empty:
        files['reference/auto_reference_metadata.csv'] = r['reference_meta'].to_csv(index=False)
    if isinstance(r.get('reference_mapping'), pd.DataFrame) and not r['reference_mapping'].empty:
        files['reference/auto_reference_mapping_qc.csv'] = r['reference_mapping'].to_csv(index=False)
    if isinstance(r.get('annotation_reconciliation'), pd.DataFrame) and not r['annotation_reconciliation'].empty:
        files['reference/annotation_reconciliation.csv'] = r['annotation_reconciliation'].to_csv(index=False)
    if r.get('miniprot_raw'):
        files['raw/miniprot.gff3'] = r['miniprot_raw']
    auto_bundles = r.get('auto_reference_bundles', {})
    if auto_bundles:
        files['reference/auto_cds.fasta'] = ''.join(f'>{g}\n{b.get("cds", "")}\n' for g,b in auto_bundles.items())
        files['reference/auto_genomic.fasta'] = ''.join(f'>{g}|{b.get("reference_record", "")}:{b.get("reference_start", "")}-{b.get("reference_end", "")}({b.get("strand", "")})\n{b.get("genomic", "")}\n' for g,b in auto_bundles.items())
    for n, b in r.get('figures', {}).items():
        files['figures/' + n] = b
    files['run_manifest.json'] = json.dumps(dict(software='BioProtein Studio', module_version='5.5.0-experimental', python=platform.python_version(), parameters=params, warnings=r['errors'] + r['warnings']), indent=2)
    return zip_files(files)


shared_cloud_runtime = bool(phylogeny_tool_status().get('shared_cloud', False))

c1, c2, c3, c4, c5 = st.columns(5)
if shared_cloud_runtime:
    c1.metric('Phylogeny', 'MAFFT + FastTree cloud mode' if external_phylogeny_ready() else 'NJ fallback available')
else:
    c1.metric('Phylogeny', 'IQ-TREE publication mode' if publication_phylogeny_ready() else ('MAFFT + FastTree' if external_phylogeny_ready() else 'NJ fallback available'))
c2.metric('MEME local', 'Ready' if meme_ready() else 'Not installed')
c3.metric('Gene mapping', 'miniprot Ready' if miniprot_ready() else 'miniprot missing')
c4.metric('Reference fetch', 'NCBI Datasets Ready' if datasets_ready() else 'datasets missing')
c5.metric('CDD', 'NCBI remote')
if not meme_ready() or not est2genome_ready():
    st.info('For full one-click mode use the supplied Linux/WSL Conda environment. The app does not replace MEME or spliced alignment with a homemade predictor.')
if shared_cloud_runtime:
    st.info('Shared Streamlit Cloud uses MAFFT + FastTree for automatic phylogeny. IQ-TREE publication inference is intentionally local/HPC-only to prevent shared-process crashes; a validated Newick tree can be uploaded here.')
elif not publication_phylogeny_ready():
    st.info('Publication phylogeny needs MAFFT + IQ-TREE. If IQ-TREE is unavailable, Auto mode falls back to MAFFT + FastTree and then NJ screening.')
if not auto_reference_ready():
    st.info('Automatic protein → reference genome → exon/intron mapping needs both miniprot and NCBI Datasets CLI. Manual CDS+genomic and annotation routes remain available.')
if 'gdm_result' not in st.session_state:
    st.session_state.gdm_result = {}
auto, imp, qctab = st.tabs(['🚀 One-click Auto Analysis', '📥 Import / Re-validate', '✅ Scientific QC Rules'])

with auto:
    mode = st.radio('Primary input', ['Protein FASTA', 'CDS FASTA (auto-translate)'], horizontal=True)
    a, b = st.columns(2)
    with a:
        up = st.file_uploader('Upload primary FASTA', type=['fa', 'fasta', 'faa', 'fna', 'txt'], key='gdm_primary')
        paste = st.text_area('Or paste FASTA', height=180, key='gdm_paste')
    with b:
        cds_up = st.file_uploader('Optional matching CDS FASTA', type=['fa', 'fasta', 'fna', 'txt'], key='gdm_cds')
        gen_up = st.file_uploader('Optional matching genomic FASTA', type=['fa', 'fasta', 'fna', 'txt'], key='gdm_gen')
        structure_up = st.file_uploader('Optional exon/CDS annotation table (CSV/TSV/GFF3/GTF)', type=['csv','tsv','txt','gff','gff3','gtf'], key='gdm_structure')
        tree_up = st.file_uploader('Optional Newick tree', type=['nwk', 'newick', 'tree', 'txt'], key='gdm_tree')
    with st.expander('NCBI reference auto-fetch for gene structure'):
        auto_ncbi = st.checkbox('Try NCBI when CDS + genomic FASTA are absent', True)
        email = st.text_input('Entrez email')
        api = st.text_input('NCBI API key (optional)', type='password')
        st.caption('Only FASTA IDs that look like real NCBI protein accessions are queried. Custom labels are never guessed by BLAST.')
    with st.expander('🧭 Automatic gene structure from Protein + Species', expanded=True):
        auto_reference = st.checkbox('Automatically resolve a reference genome and map proteins when gene structure inputs are absent', True)
        rc1, rc2 = st.columns(2)
        with rc1:
            reference_taxon = st.text_input('Species / NCBI taxon', placeholder='e.g. Carica papaya')
        with rc2:
            reference_accession = st.text_input('Optional assembly accession', placeholder='e.g. GCF_... or GCA_...')
        reference_threads = st.slider('Reference mapping threads', 1, 16, 2 if shared_cloud_runtime else 4)
        st.caption('Scientific route: NCBI Datasets selects/downloads an annotated reference assembly → miniprot maps each protein splice-aware to the genome → CDS/genomic FASTA and exon coordinates are extracted → translated CDS is checked against the input protein. Ambiguous/low-confidence loci are marked REVIEW, not forced.')
    with st.expander('Phylogeny settings', expanded=True):
        auto_tree = st.checkbox('Automatically build phylogeny when Newick is not uploaded', True)
        phylo_options = ['auto', 'fasttree', 'nj'] if shared_cloud_runtime else ['auto', 'publication', 'fasttree', 'nj']
        phylo_mode = st.selectbox(
            'Automatic tree method', phylo_options,
            format_func=lambda x: {
                'auto': 'Auto — MAFFT + FastTree cloud-safe' if shared_cloud_runtime else 'Auto — adaptive: IQ-TREE for small families, FastTree for larger families',
                'publication': 'Publication — MAFFT + IQ-TREE + ModelFinder + SH-aLRT + ultrafast bootstrap',
                'fasttree': 'Fast screening — MAFFT + FastTree',
                'nj': 'Internal Neighbor-Joining — quick screening only'
            }[x]
        )
        if shared_cloud_runtime:
            phylo_bootstrap, phylo_alrt, phylo_threads = 1000, 1000, '2'
            st.caption('Cloud stability mode: IQ-TREE controls are hidden because publication inference is local/HPC-only. Auto uses MAFFT + FastTree without subsampling.')
        else:
            pc1, pc2, pc3 = st.columns(3)
            with pc1:
                phylo_bootstrap = st.selectbox('Ultrafast bootstrap replicates', [1000, 2000, 5000], index=0)
            with pc2:
                phylo_alrt = st.selectbox('SH-aLRT replicates', [1000, 2000, 5000], index=0)
            with pc3:
                phylo_threads = st.selectbox('IQ-TREE threads', ['2', '4', '8', 'AUTO'], index=0)
            st.caption('Publication mode uses MAFFT alignment followed by IQ-TREE maximum-likelihood inference, automatic ModelFinder selection, SH-aLRT support and ultrafast bootstrap. FastTree remains a screening option.')
    s1, s2, s3, s4 = st.columns(4)
    with s1:
        do_cdd = st.checkbox('Run NCBI CDD', True)
        cdd_e = st.selectbox('CDD E-value', [0.01, 0.001, 1e-5])
    with s2:
        do_meme = st.checkbox('Run MEME', True)
        nm = st.slider('Max motifs', 3, 20, 10)
    with s3:
        minw = st.number_input('Min motif width', 3, 50, 6)
        maxw = st.number_input('Max motif width', 6, 100, 50)
    with s4:
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
            shared_cloud = shared_cloud_runtime

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
                    if shared_cloud:
                        tree_plan = 'Auto → MAFFT + FastTree cloud-safe screening'
                    elif preview_n <= AUTO_IQTREE_MAX_SEQUENCES and preview_residues <= AUTO_IQTREE_MAX_RESIDUES:
                        tree_plan = 'Auto → MAFFT + IQ-TREE publication-oriented inference'
                    else:
                        tree_plan = 'Auto → MAFFT + FastTree screening'
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
        r = empty_result()
        primary = file_text(up) if up else paste.strip()
        uploaded_tree_text = file_text(tree_up)
        tree_text = uploaded_tree_text
        if not primary:
            st.error('Provide FASTA input.')
            st.stop()
        try:
            if mode.startswith('Protein'):
                proteins = parse_fasta(primary, 'protein')
                cds = parse_fasta(file_text(cds_up), 'dna') if cds_up else {}
            else:
                cds = parse_fasta(primary, 'dna')
                proteins, r['cds_qc'] = translate_cds(cds)
            genomic = parse_fasta(file_text(gen_up), 'dna') if gen_up else {}
        except Exception as e:
            st.error(f'Sequence parsing failed: {e}')
            st.stop()

        r['proteins'] = proteins
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
        elif phylo_mode == 'auto' and shared_cloud_runtime:
            phylo_progress_text = (
                f'Phylogeny: shared-Cloud safe route ({family_nseq} proteins, {family_residues:,} aa) → '
                'MAFFT + FastTree screening...'
            )
        elif phylo_mode == 'auto' and (
            family_nseq > AUTO_IQTREE_MAX_SEQUENCES or family_residues > AUTO_IQTREE_MAX_RESIDUES
        ):
            phylo_progress_text = (
                f'Phylogeny: large-family route ({family_nseq} proteins, {family_residues:,} aa) → '
                'MAFFT + FastTree screening...'
            )
        elif phylo_mode == 'publication':
            phylo_progress_text = f'Phylogeny: MAFFT + IQ-TREE publication inference for {family_nseq} proteins...'
        else:
            phylo_progress_text = f'Phylogeny: analysing {family_nseq} proteins...'

        prog = st.progress(5, text=phylo_progress_text)
        if uploaded_tree_text:
            tree_text = uploaded_tree_text
            r['phylogeny_method'] = 'Uploaded Newick'
            r['phylogeny_qc'] = pd.DataFrame([{
                'method': 'Uploaded Newick', 'sequences': len(proteins), 'model': 'User-supplied',
                'support': 'Preserved from uploaded tree where present',
                'status': 'REVIEW: verify the alignment/model/support used to create the uploaded tree'
            }])
        elif auto_tree:
            try:
                phy = build_phylogeny(proteins, phylo_mode, bootstrap=int(phylo_bootstrap), alrt=int(phylo_alrt), threads=phylo_threads)
                tree_text = phy['tree_text']
                r['alignment_text'] = phy.get('alignment_text', '')
                r['phylogeny_method'] = phy.get('method', '')
                r['phylogeny_qc'] = phy.get('qc', pd.DataFrame())
                r['phylogeny_log'] = phy.get('log_text', '')
                r['iqtree_report'] = phy.get('iqtree_report', '')
                r['phylogeny_command'] = phy.get('command', '')
                if phy.get('warning'):
                    r['warnings'].append(phy['warning'])
            except Exception as e:
                r['errors'].append('Phylogeny: ' + str(e))
                tree_text = ''
        r['tree_text'] = tree_text

        if auto_reference and (reference_taxon.strip() or reference_accession.strip()) and not (structure_up or (cds and genomic)):
            gene_progress_text = f'Gene structure: mapping {family_nseq} proteins to the selected NCBI reference...'
        elif cds and genomic:
            gene_progress_text = f'Gene structure: validating matching CDS/genomic inputs for {family_nseq} proteins...'
        else:
            gene_progress_text = 'Gene structure: resolving available reference evidence...'
        prog.progress(15, text=gene_progress_text)
        structure_txt = file_text(structure_up) if structure_up else ''
        if structure_txt:
            try:
                r['gene_structures'], r['gene_qc'] = parse_gene_structure_annotation(structure_txt, proteins.keys())
                r['warnings'].append('Gene structure loaded from imported annotation; coordinates were not re-aligned to genomic sequence.')
            except Exception as e:
                r['errors'].append('Gene structure annotation: ' + str(e))
        elif cds and genomic:
            try:
                r['gene_structures'], r['gene_qc'], r['raw_est'] = gene_structure_batch(cds, genomic)
            except Exception as e:
                r['errors'].append('Gene structure: ' + str(e))
        elif auto_ncbi and '@' in email and any(looks_ncbi_accession(x) for x in proteins):
            ids = [x for x in proteins if looks_ncbi_accession(x)]
            try:
                r['gene_structures'], r['gene_qc'], r['ncbi_bundles'], missing = ncbi_structures(ids, email, api or None)
                if missing:
                    r['warnings'].append('NCBI direct-accession structure unresolved: ' + ', '.join(missing[:20]))
            except Exception as e:
                r['errors'].append('NCBI direct-accession structure: ' + str(e))
        elif auto_reference and (reference_taxon.strip() or reference_accession.strip()):
            try:
                rr = auto_resolve_gene_structure(
                    proteins, taxon=reference_taxon.strip(), assembly_accession=reference_accession.strip(),
                    api_key=api or '', threads=int(reference_threads)
                )
                r['gene_structures'] = rr['structures']
                r['gene_qc'] = rr['qc']
                r['auto_reference_bundles'] = rr['bundles']
                r['reference_mapping'] = rr['mapping_qc']
                r['reference_meta'] = rr['reference_table']
                r['annotation_reconciliation'] = rr.get('annotation_reconciliation', pd.DataFrame())
                r['miniprot_raw'] = rr['raw_gff']
                ref = rr.get('reference', {})
                r['warnings'].append('Automatic gene structure used NCBI reference assembly ' + str(ref.get('accession','')) + '. Review mapping QC before publication.')
                if isinstance(r.get('annotation_reconciliation'), pd.DataFrame) and not r['annotation_reconciliation'].empty:
                    used = int(r['annotation_reconciliation'].decision.astype(str).str.startswith('USED').sum())
                    rescued = int((r['annotation_reconciliation'].decision == 'COMPUTATIONAL_RESCUE_REVIEW').sum())
                    if used:
                        r['warnings'].append(f'Reference annotation reconciled {used} mapped gene structure(s) against NCBI GFF3.')
                    if rescued:
                        r['warnings'].append(f'{rescued} gene structure(s) required computational splice rescue and remain REVIEW until annotation/manual confirmation.')
            except Exception as e:
                r['errors'].append('Automatic reference gene structure: ' + str(e))
        if r['gene_structures'].empty:
            if auto_reference and not (reference_taxon.strip() or reference_accession.strip()):
                r['warnings'].append('Gene structure unavailable: enter Species/NCBI taxon (or an assembly accession) to let the app retrieve a reference genome automatically.')
            else:
                r['warnings'].append('Gene structure unavailable for this run. No validated genomic reference mapping was produced.')

        cdd_batches = max(1, (family_nseq + 199) // 200)
        prog.progress(42, text=f'Conserved domains: NCBI CDD for {family_nseq} proteins ({cdd_batches} batch(es))...')
        if do_cdd:
            if len(proteins) > 200:
                r['warnings'].append(
                    f'CDD large-family mode: all {len(proteins)} proteins are analysed in transparent NCBI batches; '
                    'no sequences are subsampled.'
                )
            try:
                def _cdd_progress(done, total, message):
                    pct = 42 + int(24 * done / max(1, total))
                    prog.progress(min(66, pct), text=f'CDD {done}/{total}: {message}')

                r['domains_raw'], r['cdd_rid'], r['cdd_raw'] = run_cdd(
                    proteins, float(cdd_e), progress=_cdd_progress
                )
                r['domains'] = collapse_domains(r['domains_raw'])
                r['domain_qc'] = domain_qc(r['domains'])
            except Exception as e:
                r['errors'].append('CDD: ' + str(e))

        prog.progress(
            68,
            text=f'MEME motifs: {family_nseq} proteins, up to {int(nm)} motifs (cloud safety limits apply)...'
        )
        if do_meme:
            try:
                m = run_meme(proteins, int(nm), int(minw), int(maxw), model)
                r['motifs'] = m['sites']
                r['motif_summary'] = m['summary']
                r['meme_xml'] = m['xml']
                r['meme_zip'] = m['zip']
                r['logos'] = m['logos']
            except Exception as e:
                msg = str(e)
                if msg.startswith('CLOUD_RESOURCE_LIMIT:'):
                    r['warnings'].append('MEME cloud safety: ' + msg.split(':', 1)[1].strip())
                else:
                    r['errors'].append('MEME: ' + msg)

        r['motif_qc'] = motif_qc(r['motifs'], r['motif_summary'], len(proteins), r['domains'])
        order = choose_order(
            set(proteins)
            | set(r['gene_structures'].gene if not r['gene_structures'].empty else [])
            | set(r['domains'].gene if not r['domains'].empty else [])
            | set(r['motifs'].gene if not r['motifs'].empty else []),
            newick_order(tree_text) if tree_text else []
        )
        r['order'] = order

        prog.progress(88, text='Preparing reproducibility package...')
        # Figures are deliberately rendered on demand below. Generating every
        # panel in four high-resolution formats during analysis is expensive
        # and can exhaust shared-cloud memory for large protein families.
        r['figures'] = {}
        if len(order) > 60:
            r['warnings'].append(
                f'Large-family visualization mode enabled for {len(order)} sequences. '
                'Figures are paginated and rendered on demand; the full tree/alignment remain downloadable.'
            )

        params = dict(cdd_evalue=float(cdd_e), meme_nmotifs=int(nm), meme_min_width=int(minw), meme_max_width=int(maxw), meme_model=model, plot_domains='specific', plot_motifs='pass', tree_provided=bool(uploaded_tree_text), phylogeny_method=r.get('phylogeny_method',''), phylogeny_mode=phylo_mode, phylogeny_bootstrap=int(phylo_bootstrap), phylogeny_alrt=int(phylo_alrt), auto_phylogeny=bool(auto_tree), auto_reference=bool(auto_reference), reference_taxon=reference_taxon.strip(), reference_accession=reference_accession.strip())
        r['package'] = package(r, protein_txt, file_text(cds_up) if cds_up else (primary if mode.startswith('CDS') else ''), file_text(gen_up), tree_text, params, structure_txt)
        # Streamlit Community Cloud uses ephemeral storage. Avoid writing duplicate
        # ZIP archives on every run there; local/WSL users still receive autosave.
        if not Path('/mount/src').exists():
            try:
                outdir = Path.cwd() / 'results'
                outdir.mkdir(parents=True, exist_ok=True)
                latest = outdir / 'BioProtein_Studio_GDM_Results_v5_5_latest.zip'
                latest.write_bytes(r['package'])
                stamped = outdir / f'BioProtein_Studio_GDM_Results_v5_5_{datetime.now().strftime("%Y%m%d_%H%M%S")}.zip'
                stamped.write_bytes(r['package'])
                r['autosave_path'] = str(latest)
            except Exception as e:
                r['warnings'].append('Automatic result-package save failed: ' + str(e))
        st.session_state.gdm_result = r
        prog.progress(100, text='Complete')

    r = st.session_state.gdm_result
    if r:
        for x in r.get('errors', []):
            st.error(x)
        for x in r.get('warnings', []):
            st.warning(x)

        st.subheader('Validation dashboard')
        st.dataframe(validation_summary(r), width='stretch', hide_index=True)
        with st.expander('Sequence QC'):
            st.dataframe(r['sequence_qc'], width='stretch', hide_index=True)
            if not r['cds_qc'].empty:
                st.dataframe(r['cds_qc'], width='stretch', hide_index=True)

        if isinstance(r.get('reference_meta'), pd.DataFrame) and not r['reference_meta'].empty:
            with st.expander('Automatic reference genome & mapping QC', expanded=True):
                st.dataframe(r['reference_meta'], width='stretch', hide_index=True)
                if isinstance(r.get('reference_mapping'), pd.DataFrame) and not r['reference_mapping'].empty:
                    st.dataframe(r['reference_mapping'], width='stretch', hide_index=True)
                if isinstance(r.get('annotation_reconciliation'), pd.DataFrame) and not r['annotation_reconciliation'].empty:
                    st.markdown('**Reference annotation reconciliation**')
                    st.dataframe(r['annotation_reconciliation'], width='stretch', hide_index=True)

                auto_bundles = r.get('auto_reference_bundles', {})
                if auto_bundles:
                    auto_cds_txt = ''.join(f'>{g}\n{b.get("cds", "")}\n' for g,b in auto_bundles.items())
                    auto_gen_txt = ''.join(f'>{g}|{b.get("reference_record", "")}:{b.get("reference_start", "")}-{b.get("reference_end", "")}({b.get("strand", "")})\n{b.get("genomic", "")}\n' for g,b in auto_bundles.items())
                    dl1, dl2 = st.columns(2)
                    with dl1:
                        st.download_button('Download auto-generated CDS FASTA', auto_cds_txt, 'auto_cds.fasta', 'text/plain', width='stretch')
                    with dl2:
                        st.download_button('Download auto-generated genomic FASTA', auto_gen_txt, 'auto_genomic.fasta', 'text/plain', width='stretch')

        st.subheader('Visualization options')
        vc1, vc2 = st.columns(2)
        with vc1:
            domain_mode = st.radio('Domain display', ['specific', 'all'], horizontal=True, format_func=lambda x: {'specific': 'Specific hits only', 'all': 'All retained hits'}[x])
        with vc2:
            motif_mode = st.radio('Motif display', ['pass', 'all'], horizontal=True, format_func=lambda x: {'pass': 'PASS motifs only', 'all': 'All motifs'}[x])
        dom_plot = filter_domains(r['domains'], domain_mode)
        mot_plot = filter_motifs(r['motifs'], r['motif_qc'], motif_mode)

        full_order = list(r['order'])
        n_family = len(full_order)
        display_order = full_order
        page_num = 1

        if n_family > 60:
            st.info(
                f'Large-family viewer: {n_family} sequences detected. Analysis uses all sequences, '
                'while browser figures are paginated to protect CPU and memory.'
            )
            pg1, pg2 = st.columns(2)
            with pg1:
                page_size = st.selectbox(
                    'Sequences shown per figure', [40, 60, 100], index=1, key='gdm_page_size'
                )
            page_count = max(1, (n_family + page_size - 1) // page_size)
            with pg2:
                page_num = int(st.number_input(
                    'Figure page', min_value=1, max_value=page_count, value=1, step=1,
                    key='gdm_page_number'
                ))
            lo = (page_num - 1) * page_size
            hi = min(n_family, lo + page_size)
            display_order = full_order[lo:hi]
            st.caption(
                f'Displaying sequences {lo + 1}–{hi} of {n_family}. '
                'Full Newick/alignment downloads still contain the complete family.'
            )

        with st.expander('🎨 Graph Studio — auto styles, colors & layout', expanded=True):
            st.caption('Click Generate alternatives to preview several publication-ready color versions. Select any preset, then optionally fine-tune individual colors. All figure tabs and downloads update to the selected style.')
            if 'gdm_show_style_variants' not in st.session_state:
                st.session_state.gdm_show_style_variants = False
            if n_family > 60:
                st.session_state.gdm_show_style_variants = False
                st.caption('Alternative multi-panel previews are disabled for large families to protect cloud memory.')
            if st.button(
                '✨ Generate alternative color versions',
                key='gdm_generate_styles',
                width='stretch',
                disabled=(n_family > 60),
            ):
                st.session_state.gdm_show_style_variants = True
            if st.session_state.gdm_show_style_variants:
                pv_names = list(STYLE_PRESETS.keys())[:4]
                pv_cols = st.columns(2)
                for pi, pname in enumerate(pv_names):
                    with pv_cols[pi % 2]:
                        pf = combined(r.get('tree_text',''), r['gene_structures'], dom_plot, mot_plot, display_order, style_from_preset(pname))
                        if pf:
                            st.caption(pname)
                            st.pyplot(pf, width='stretch')
                            plt.close(pf)
            preset_name = st.radio('Use figure style', list(STYLE_PRESETS.keys()), horizontal=True, key='gdm_style_preset')
            graph_style = style_from_preset(preset_name)
            fine_tune = st.checkbox('Fine-tune colors and figure sizing', False, key='gdm_fine_tune')
            if fine_tune:
                gc1, gc2, gc3, gc4 = st.columns(4)
                with gc1:
                    graph_style['tree_color'] = st.color_picker('Tree branches', graph_style['tree_color'])
                    graph_style['support_color'] = st.color_picker('Support labels', graph_style['support_color'])
                with gc2:
                    graph_style['exon_color'] = st.color_picker('Validated exons', graph_style['exon_color'])
                    graph_style['review_exon_color'] = st.color_picker('REVIEW exons', graph_style['review_exon_color'])
                with gc3:
                    graph_style['intron_color'] = st.color_picker('Introns', graph_style['intron_color'])
                    graph_style['backbone_color'] = st.color_picker('Protein backbone', graph_style['backbone_color'])
                with gc4:
                    graph_style['edge_color'] = st.color_picker('Feature borders', graph_style['edge_color'])
                    graph_style['title_size'] = st.slider('Title size', 10, 20, int(graph_style['title_size']))
                graph_style['label_size'] = st.slider('Gene label size', 6, 14, int(graph_style['label_size']))
                graph_style['line_width'] = st.slider('Line width', 0.5, 2.5, float(graph_style['line_width']), 0.1)

                motif_labels = list(dict.fromkeys(mot_plot['motif'].astype(str))) if not mot_plot.empty else []
                domain_labels = list(dict.fromkeys(dom_plot['domain'].astype(str))) if not dom_plot.empty else []
                default_motif = assign_colors(motif_labels, graph_style['motif_palette'])
                default_domain = assign_colors(domain_labels, graph_style['domain_palette'])
                if motif_labels:
                    st.markdown('**Motif colors**')
                    mcols = st.columns(min(5, len(motif_labels)))
                    for mi, lab in enumerate(motif_labels):
                        with mcols[mi % len(mcols)]:
                            graph_style['motif_colors'][lab] = st.color_picker(lab, default_motif[lab], key='motif_color_'+lab)
                if domain_labels:
                    st.markdown('**Domain colors**')
                    dcols = st.columns(min(4, len(domain_labels)))
                    for di, lab in enumerate(domain_labels):
                        with dcols[di % len(dcols)]:
                            graph_style['domain_colors'][lab] = st.color_picker(lab, default_domain[lab], key='domain_color_'+lab)

        st.subheader('Figure viewer')

        if r.get('tree_text'):
            td1, td2 = st.columns(2)
            with td1:
                st.download_button(
                    'Download full Newick tree', r['tree_text'], 'phylogenetic_tree_full.nwk',
                    'text/plain', width='stretch'
                )
            if r.get('alignment_text'):
                with td2:
                    st.download_button(
                        'Download full protein alignment', r['alignment_text'], 'protein_alignment_full.fasta',
                        'text/plain', width='stretch'
                    )

        figure_choice = st.radio(
            'Figure to display',
            ['Phylogeny', 'Gene Structure', 'Domains', 'Motifs', 'Integrated Figure'],
            horizontal=True,
            key='gdm_figure_choice',
        )

        selected_fig = None
        selected_table = pd.DataFrame()
        stem = 'figure'

        if figure_choice == 'Phylogeny':
            selected_fig = phylogeny_figure(r.get('tree_text', ''), display_order, graph_style)
            selected_table = r.get('phylogeny_qc', pd.DataFrame())
            stem = 'phylogenetic_tree'
        elif figure_choice == 'Gene Structure':
            selected_fig = (
                gene_structure(r['gene_structures'], display_order, graph_style)
                if not r['gene_structures'].empty
                else missing_structure_figure(display_order, graph_style)
            )
            selected_table = r['gene_qc']
            stem = 'gene_structure'
        elif figure_choice == 'Domains':
            selected_fig = architecture(
                dom_plot, display_order, 'domain', 'Conserved Domain Architecture', style=graph_style
            )
            selected_table = r['domain_qc']
            stem = 'domain_architecture'
        elif figure_choice == 'Motifs':
            selected_fig = architecture(
                mot_plot, display_order, 'motif', 'Conserved Motif Architecture', style=graph_style
            )
            selected_table = r['motif_qc']
            stem = 'motif_architecture'
        elif figure_choice == 'Integrated Figure':
            selected_fig = combined(
                r.get('tree_text', ''), r['gene_structures'], dom_plot, mot_plot,
                display_order, graph_style
            )
            stem = 'integrated_architecture'

        if selected_fig:
            st.pyplot(selected_fig, width='stretch')
            if isinstance(selected_table, pd.DataFrame) and not selected_table.empty:
                st.dataframe(selected_table, width='stretch', hide_index=True)

            suffix = f'_page{page_num}' if n_family > 60 else ''
            key_base = f'{stem}_{domain_mode}_{motif_mode}_{page_num}'

            export_fmt = st.selectbox(
                'Figure download format', ['SVG', 'PDF', 'PNG', 'TIFF'],
                index=0, key='gdm_export_format'
            )
            if st.button('Prepare selected figure download', key=key_base + '_prepare', width='stretch'):
                fmt = export_fmt.lower()
                dpi = 300 if fmt in {'svg', 'pdf', 'tiff'} else 220
                data = fig_bytes(selected_fig, fmt, dpi)
                mime = {
                    'svg': 'image/svg+xml',
                    'pdf': 'application/pdf',
                    'png': 'image/png',
                    'tiff': 'image/tiff',
                }[fmt]
                st.download_button(
                    f'Download {export_fmt}', data, f'{stem}{suffix}.{fmt}', mime,
                    key=key_base + '_download_' + fmt, width='stretch'
                )
            plt.close(selected_fig)
        else:
            st.info('No result available for this figure.')

        st.caption(
            'Cloud-safe rendering: one figure is generated at a time. Large families are paginated for display, '
            'while the complete analysis, Newick tree, and alignment remain preserved.'
        )

        st.download_button('📦 Download complete analysis package', r['package'], 'BioProtein_Studio_GDM_Results_v5_5.zip', 'application/zip', type='primary', width='stretch')
        if r.get('autosave_path'):
            st.caption('Auto-saved locally: ' + r['autosave_path'])

with imp:
    st.write('Import your previous manual NCBI CDD and MEME results and re-validate them.')
    c = st.file_uploader('CDD hit file', type=['txt', 'tsv'], key='imp_cdd')
    m = st.file_uploader('MEME XML', type=['xml'], key='imp_meme')
    if c:
        try:
            st.dataframe(domain_qc(collapse_domains(parse_cdd(file_text(c)))), width='stretch', hide_index=True)
        except Exception as e:
            st.error(str(e))
    if m:
        try:
            s, u = parse_meme_xml(file_text(m))
            st.dataframe(u, width='stretch', hide_index=True)
            st.dataframe(s, width='stretch', hide_index=True)
        except Exception as e:
            st.error(str(e))

with qctab:
    st.markdown('''**Gene structure:** reference CDS feature or EMBOSS `est2genome`; PASS requires >=95% CDS coverage and >=95% weighted exon identity in the alignment route. Non-canonical splice evidence and frame problems are flagged.\n\n**CDD domains:** NCBI specific hits are strongest. Other hits are ranked by E-value; raw full results and Search-ID are kept. Redundant highly overlapping footprints are collapsed only for visualization.\n\n**MEME motifs:** real MEME Suite is run in protein mode. ZOOPS is the default for gene families. Motif E < 0.05 is required, and the app adds a >=50% family-prevalence rule for PASS; lower-prevalence significant motifs remain biologically reviewable.\n\n**No fabricated result:** protein sequence alone does not define introns. With custom gene IDs and no genomic/CDS reference, gene structure is deliberately left unavailable.''')

st.caption('BioProtein Studio • Phylogeny / Gene Structure / CDD / MEME module • v5.5 experimental')
