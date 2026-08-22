import json, platform
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from modules.gdm_common import parse_fasta, fasta_text, protein_qc, translate_cds, looks_ncbi_accession, newick_order, choose_order, zip_files
from modules.gdm_structure import est2genome_ready, gene_structure_batch, ncbi_structures, parse_gene_structure_annotation
from modules.gdm_cdd_meme import run_cdd, collapse_domains, domain_qc, meme_ready, run_meme, motif_qc, parse_cdd, parse_meme_xml
from modules.gdm_phylogeny import build_phylogeny, external_phylogeny_ready, publication_phylogeny_ready, phylogeny_tool_status
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


c1, c2, c3, c4, c5 = st.columns(5)
c1.metric('Phylogeny', 'IQ-TREE publication mode' if publication_phylogeny_ready() else ('MAFFT + FastTree' if external_phylogeny_ready() else 'NJ fallback available'))
c2.metric('MEME local', 'Ready' if meme_ready() else 'Not installed')
c3.metric('Gene mapping', 'miniprot Ready' if miniprot_ready() else 'miniprot missing')
c4.metric('Reference fetch', 'NCBI Datasets Ready' if datasets_ready() else 'datasets missing')
c5.metric('CDD', 'NCBI remote')
if not meme_ready() or not est2genome_ready():
    st.info('For full one-click mode use the supplied Linux/WSL Conda environment. The app does not replace MEME or spliced alignment with a homemade predictor.')
if not publication_phylogeny_ready():
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
        reference_threads = st.slider('Reference mapping threads', 1, 16, 4)
        st.caption('Scientific route: NCBI Datasets selects/downloads an annotated reference assembly → miniprot maps each protein splice-aware to the genome → CDS/genomic FASTA and exon coordinates are extracted → translated CDS is checked against the input protein. Ambiguous/low-confidence loci are marked REVIEW, not forced.')
    with st.expander('Phylogeny settings', expanded=True):
        auto_tree = st.checkbox('Automatically build phylogeny when Newick is not uploaded', True)
        phylo_mode = st.selectbox(
            'Automatic tree method', ['auto', 'publication', 'fasttree', 'nj'],
            format_func=lambda x: {
                'auto': 'Auto — IQ-TREE publication mode if installed, otherwise FastTree, then NJ fallback',
                'publication': 'Publication — MAFFT + IQ-TREE + ModelFinder + SH-aLRT + ultrafast bootstrap',
                'fasttree': 'Fast screening — MAFFT + FastTree',
                'nj': 'Internal Neighbor-Joining — quick screening only'
            }[x]
        )
        pc1, pc2, pc3 = st.columns(3)
        with pc1:
            phylo_bootstrap = st.selectbox('Ultrafast bootstrap replicates', [1000, 2000, 5000], index=0)
        with pc2:
            phylo_alrt = st.selectbox('SH-aLRT replicates', [1000, 2000, 5000], index=0)
        with pc3:
            phylo_threads = st.selectbox('IQ-TREE threads', ['AUTO', '2', '4', '8'], index=0)
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
        if (r['sequence_qc'].status == 'REVIEW').any():
            r['warnings'].append('Sequence-integrity QC flagged one or more proteins.')

        prog = st.progress(5, text='Phylogeny...')
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

        prog.progress(15, text='Gene structure...')
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

        prog.progress(42, text='Conserved domains...')
        if do_cdd:
            try:
                r['domains_raw'], r['cdd_rid'], r['cdd_raw'] = run_cdd(proteins, float(cdd_e))
                r['domains'] = collapse_domains(r['domains_raw'])
                r['domain_qc'] = domain_qc(r['domains'])
            except Exception as e:
                r['errors'].append('CDD: ' + str(e))

        prog.progress(68, text='MEME motifs...')
        if do_meme:
            try:
                m = run_meme(proteins, int(nm), int(minw), int(maxw), model)
                r['motifs'] = m['sites']
                r['motif_summary'] = m['summary']
                r['meme_xml'] = m['xml']
                r['meme_zip'] = m['zip']
                r['logos'] = m['logos']
            except Exception as e:
                r['errors'].append('MEME: ' + str(e))

        r['motif_qc'] = motif_qc(r['motifs'], r['motif_summary'], len(proteins), r['domains'])
        order = choose_order(
            set(proteins)
            | set(r['gene_structures'].gene if not r['gene_structures'].empty else [])
            | set(r['domains'].gene if not r['domains'].empty else [])
            | set(r['motifs'].gene if not r['motifs'].empty else []),
            newick_order(tree_text) if tree_text else []
        )
        r['order'] = order

        prog.progress(88, text='Figures and reproducibility package...')
        plot_domains = filter_domains(r['domains'], 'specific')
        plot_motifs = filter_motifs(r['motifs'], r['motif_qc'], 'pass')
        figs = {
            'phylogenetic_tree': phylogeny_figure(tree_text, order),
            'gene_structure': gene_structure(r['gene_structures'], order),
            'domain_architecture': architecture(plot_domains, order, 'domain', 'Conserved Domain Architecture'),
            'motif_architecture': architecture(plot_motifs, order, 'motif', 'Conserved Motif Architecture'),
            'integrated_architecture': combined(tree_text, r['gene_structures'], plot_domains, plot_motifs, order)
        }
        for stem, fig in figs.items():
            if fig:
                for fmt in ['png', 'svg', 'pdf', 'tiff']:
                    r['figures'][f'{stem}.{fmt}'] = fig_bytes(fig, fmt, 600)
                plt.close(fig)

        params = dict(cdd_evalue=float(cdd_e), meme_nmotifs=int(nm), meme_min_width=int(minw), meme_max_width=int(maxw), meme_model=model, plot_domains='specific', plot_motifs='pass', tree_provided=bool(uploaded_tree_text), phylogeny_method=r.get('phylogeny_method',''), phylogeny_mode=phylo_mode, phylogeny_bootstrap=int(phylo_bootstrap), phylogeny_alrt=int(phylo_alrt), auto_phylogeny=bool(auto_tree), auto_reference=bool(auto_reference), reference_taxon=reference_taxon.strip(), reference_accession=reference_accession.strip())
        r['package'] = package(r, protein_txt, file_text(cds_up) if cds_up else (primary if mode.startswith('CDS') else ''), file_text(gen_up), tree_text, params, structure_txt)
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

        with st.expander('🎨 Graph Studio — auto styles, colors & layout', expanded=True):
            st.caption('Click Generate alternatives to preview several publication-ready color versions. Select any preset, then optionally fine-tune individual colors. All figure tabs and downloads update to the selected style.')
            if 'gdm_show_style_variants' not in st.session_state:
                st.session_state.gdm_show_style_variants = False
            if st.button('✨ Generate alternative color versions', key='gdm_generate_styles', width='stretch'):
                st.session_state.gdm_show_style_variants = True
            if st.session_state.gdm_show_style_variants:
                pv_names = list(STYLE_PRESETS.keys())[:4]
                pv_cols = st.columns(2)
                for pi, pname in enumerate(pv_names):
                    with pv_cols[pi % 2]:
                        pf = combined(r.get('tree_text',''), r['gene_structures'], dom_plot, mot_plot, r['order'], style_from_preset(pname))
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

        tabs = st.tabs(['Phylogeny', 'Gene Structure', 'Domains', 'Motifs', 'Integrated Figure'])
        items = [
            (tabs[0], phylogeny_figure(r.get('tree_text',''), r['order'], graph_style), r.get('phylogeny_qc', pd.DataFrame())),
            (tabs[1], gene_structure(r['gene_structures'], r['order'], graph_style) if not r['gene_structures'].empty else missing_structure_figure(r['order'], graph_style), r['gene_qc']),
            (tabs[2], architecture(dom_plot, r['order'], 'domain', 'Conserved Domain Architecture', style=graph_style), r['domain_qc']),
            (tabs[3], architecture(mot_plot, r['order'], 'motif', 'Conserved Motif Architecture', style=graph_style), r['motif_qc']),
            (tabs[4], combined(r.get('tree_text', ''), r['gene_structures'], dom_plot, mot_plot, r['order'], graph_style), pd.DataFrame())
        ]
        stems = ['phylogenetic_tree', 'gene_structure', 'domain_architecture', 'motif_architecture', 'integrated_architecture']
        for i, (tab, fig, table) in enumerate(items):
            with tab:
                if fig:
                    st.pyplot(fig)
                    fig_bytes_now = {fmt: fig_bytes(fig, fmt, 600) for fmt in ['png', 'svg', 'pdf', 'tiff']}
                    plt.close(fig)
                    if not table.empty:
                        st.dataframe(table, width='stretch', hide_index=True)
                    for fmt in ['png', 'svg', 'pdf', 'tiff']:
                        st.download_button(f'Download {fmt.upper()}', fig_bytes_now[fmt], f'{stems[i]}.{fmt}', key=f'{stems[i]}_{fmt}_{domain_mode}_{motif_mode}')
                else:
                    st.info('No result available for this analysis.')

        styled_files = {}
        styled_specs = {
            'phylogenetic_tree': phylogeny_figure(r.get('tree_text',''), r['order'], graph_style),
            'gene_structure': gene_structure(r['gene_structures'], r['order'], graph_style) if not r['gene_structures'].empty else missing_structure_figure(r['order'], graph_style),
            'domain_architecture': architecture(dom_plot, r['order'], 'domain', 'Conserved Domain Architecture', style=graph_style),
            'motif_architecture': architecture(mot_plot, r['order'], 'motif', 'Conserved Motif Architecture', style=graph_style),
            'integrated_architecture': combined(r.get('tree_text',''), r['gene_structures'], dom_plot, mot_plot, r['order'], graph_style),
        }
        for stem_name, sf in styled_specs.items():
            if sf:
                for fmt in ['png','svg','pdf','tiff']:
                    styled_files[f'{stem_name}.{fmt}'] = fig_bytes(sf, fmt, 600)
                plt.close(sf)
        styled_files['GRAPH_STYLE.json'] = json.dumps(graph_style, indent=2)
        st.download_button('🎨 Download customized figure set', zip_files(styled_files), 'BioProtein_Studio_Custom_Figures.zip', 'application/zip', width='stretch')

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
