import io

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from modules.gdm_circular import circular_phylogeny_figure
from modules.gdm_clade_validation import clade_sequence_concordance
from modules.gdm_common import parse_fasta
from modules.gdm_phylogeny import build_phylogeny
from modules.gdm_plot import fig_bytes
from modules.gdm_style import STYLE_PRESETS, style_from_preset


st.set_page_config(page_title='Circular Phylogeny | BioProtein Studio', page_icon='🌳', layout='wide')
st.title('🌳 Circular Phylogeny & Sequence-Supported Clade View')
st.caption(
    'Build a circular tree directly from protein FASTA or re-visualize a validated Newick tree. '
    'Automatic clades can be cross-checked against independent unsupervised sequence-distance clustering.'
)

r = st.session_state.get('gdm_result') or {}
session_tree = str(r.get('tree_text', '') or '').strip()
session_alignment = str(r.get('alignment_text', '') or '').strip()
session_order = list(r.get('order') or [])
session_method = str(r.get('phylogeny_method', '') or '').strip()
session_qc = r.get('phylogeny_qc')

st.info(
    'Clade evidence is kept method-correct. FastTree local support is not relabelled as IQ-TREE support. '
    'The optional ML concordance is an independent unsupervised k-medoids check of MAFFT sequence distances. '
    'It can strengthen or challenge an automatic topology partition, but it is not a substitute for '
    'reference-gene, functional, taxonomic, or experimental validation of biological subgroups.'
)

source_options = []
if session_tree:
    source_options.append('Latest full GDM analysis')
source_options += ['Protein FASTA → automatic tree', 'Upload / paste Newick']
source = st.radio('Tree source', source_options, horizontal=True)

uploaded_tree = None
pasted_tree = ''
iqtree_report = ''
method_label = ''
order = None
alignment_text = ''
tree_text = ''
phylo_qc = pd.DataFrame()

if source == 'Latest full GDM analysis':
    tree_text = session_tree
    alignment_text = session_alignment
    order = session_order or None
    method_label = session_method
    phylo_qc = session_qc if isinstance(session_qc, pd.DataFrame) else pd.DataFrame()
    if not phylo_qc.empty:
        st.markdown('**Inference evidence from the latest GDM run**')
        st.dataframe(phylo_qc, width='stretch', hide_index=True)
    if 'FastTree' in method_label:
        st.warning(
            'Current tree is a Cloud screening tree (MAFFT + FastTree). Branch lengths and FastTree local '
            'support are preserved. ModelFinder, SH-aLRT and ultrafast bootstrap were not computed.'
        )
    elif 'IQ-TREE' in method_label:
        st.success('Current tree contains IQ-TREE publication-oriented inference evidence from the GDM run.')

elif source == 'Protein FASTA → automatic tree':
    st.markdown('### Protein input')
    protein_up = st.file_uploader(
        'Upload protein FASTA', type=['fa', 'fasta', 'faa', 'txt'], key='circular_protein_upload'
    )
    protein_paste = st.text_area(
        'Or paste protein FASTA', height=160, key='circular_protein_paste',
        placeholder='>Gene1\nM...\n>Gene2\nM...'
    )
    protein_text = (
        protein_up.getvalue().decode('utf-8', 'replace') if protein_up else protein_paste
    )
    if st.button('Build tree from protein sequences', type='primary', width='stretch'):
        try:
            proteins = parse_fasta(protein_text, 'protein')
            if len(proteins) < 3:
                raise ValueError('At least three protein sequences are required for phylogenetic inference.')
            with st.spinner(f'Aligning and inferring phylogeny for {len(proteins)} proteins...'):
                auto_result = build_phylogeny(proteins, mode='auto')
            auto_result['input_ids'] = list(proteins)
            st.session_state.circular_auto_result = auto_result
        except Exception as exc:
            st.error(f'Automatic protein phylogeny failed: {exc}')

    auto_result = st.session_state.get('circular_auto_result') or {}
    tree_text = str(auto_result.get('tree_text', '') or '').strip()
    alignment_text = str(auto_result.get('alignment_text', '') or '').strip()
    order = list(auto_result.get('input_ids') or []) or None
    method_label = str(auto_result.get('method', '') or '')
    phylo_qc = auto_result.get('qc') if isinstance(auto_result.get('qc'), pd.DataFrame) else pd.DataFrame()
    if tree_text:
        if not phylo_qc.empty:
            st.dataframe(phylo_qc, width='stretch', hide_index=True)
        warning = str(auto_result.get('warning', '') or '').strip()
        if warning:
            st.warning(warning)
        st.success(f'Automatic tree ready: {method_label or "phylogeny inferred"}.')

else:
    c1, c2 = st.columns(2)
    with c1:
        uploaded_tree = st.file_uploader(
            'Newick / IQ-TREE .treefile',
            type=['nwk', 'newick', 'tree', 'treefile', 'txt'],
        )
    with c2:
        iq_report_up = st.file_uploader(
            'Optional IQ-TREE report (.iqtree / .log)',
            type=['iqtree', 'log', 'txt'],
        )
    pasted_tree = st.text_area('Or paste Newick', height=130)
    if uploaded_tree:
        tree_text = uploaded_tree.getvalue().decode('utf-8', 'replace').strip()
    else:
        tree_text = pasted_tree.strip()
    if iq_report_up:
        iqtree_report = iq_report_up.getvalue().decode('utf-8', 'replace')
        method_label = 'Uploaded IQ-TREE tree/report (support semantics preserved from source)'
        with st.expander('IQ-TREE report', expanded=False):
            st.code(iqtree_report[:20000], language='text')
    else:
        method_label = 'Uploaded Newick (support semantics user-supplied)'

    aligned_up = st.file_uploader(
        'Optional aligned protein FASTA for independent sequence-cluster concordance',
        type=['fa', 'fasta', 'faa', 'txt'],
        key='circular_alignment_upload',
        help='Use the MAFFT alignment corresponding to this tree. Without an alignment, the tree can still be plotted but ML concordance is not claimed.',
    )
    if aligned_up:
        alignment_text = aligned_up.getvalue().decode('utf-8', 'replace')

if not tree_text:
    st.warning('Provide protein FASTA, run the full GDM analysis, or supply a Newick tree.')
    st.stop()

st.subheader('Circular tree & clade evidence options')
a, b, c, d = st.columns(4)
with a:
    group_choice = st.selectbox(
        'Automatic clade groups',
        ['Data-driven ML Auto'] + list(range(2, 13)),
        help='Data-driven Auto tests several topology partitions and selects the one most concordant with the independent sequence-distance clustering when an alignment is available.'
    )
    requested_groups = None if group_choice == 'Data-driven ML Auto' else int(group_choice)
with b:
    midpoint = st.checkbox('Midpoint-root for display', value=False)
with c:
    show_support = st.checkbox('Show internal support', value=True)
with d:
    preset = st.selectbox('Figure style', list(STYLE_PRESETS.keys()), index=0)

style = style_from_preset(preset)
e1, e2 = st.columns(2)
with e1:
    style['label_size'] = st.slider('Tip-label size', 5, 14, int(style['label_size']))
with e2:
    style['line_width'] = st.slider('Branch line width', 0.5, 2.5, float(style['line_width']), 0.1)

validation = None
plot_groups = requested_groups
if alignment_text:
    try:
        validation = clade_sequence_concordance(
            tree_text,
            alignment_text,
            target_groups=requested_groups,
        )
        if requested_groups is None:
            plot_groups = int(validation['recommended_groups'])
        st.markdown('### Independent sequence-cluster concordance')
        st.dataframe(validation['summary'], width='stretch', hide_index=True)
        if str(validation['status']).startswith('STRONG'):
            st.success(
                'Tree partition and unsupervised sequence-distance clustering show strong computational concordance.'
            )
        elif str(validation['status']).startswith('MODERATE'):
            st.info(
                'Tree partition and unsupervised sequence-distance clustering show moderate computational concordance.'
            )
        else:
            st.warning(
                'Topology and independent sequence clustering are not strongly concordant. Treat automatic clades as REVIEW.'
            )
    except Exception as exc:
        st.warning(f'Sequence-cluster concordance was not assigned: {exc}')
else:
    st.warning(
        'No matching protein alignment is available, so automatic clades are topology-only. '
        'Provide protein FASTA through the automatic route or an aligned FASTA to enable independent ML concordance.'
    )

try:
    fig, topology_clade_table = circular_phylogeny_figure(
        tree_text,
        order=order or None,
        style=style,
        target_groups=plot_groups,
        midpoint_root=midpoint,
        show_support=show_support,
        method_label=method_label,
    )
except Exception as exc:
    st.error(f'Circular phylogeny could not be rendered: {exc}')
    st.stop()

st.subheader('Automatic circular clade figure')
st.pyplot(fig, width='stretch')

if validation is not None:
    st.markdown('**Clade evidence table**')
    evidence_table = validation['clade_table']
    st.dataframe(evidence_table, width='stretch', hide_index=True)
    st.download_button(
        'Download sequence-supported clade evidence CSV',
        evidence_table.to_csv(index=False),
        'sequence_supported_clade_evidence.csv',
        'text/csv',
        width='stretch',
    )
    with st.expander('Per-protein topology vs ML-cluster assignment'):
        st.dataframe(validation['member_table'], width='stretch', hide_index=True)
else:
    st.markdown('**Topology-derived clade table**')
    st.dataframe(topology_clade_table, width='stretch', hide_index=True)
    st.download_button(
        'Download topology clade membership CSV',
        topology_clade_table.to_csv(index=False),
        'automatic_topology_clades.csv',
        'text/csv',
        width='stretch',
    )

st.caption(
    'ML-CONCORDANT means the monophyletic topology partition is independently recovered by unsupervised '
    'sequence-distance clustering. It does not by itself establish a named biological subgroup. '
    'For publication subgroup claims, combine this evidence with branch support, curated references, '
    'domains/motifs, species context and other relevant biology.'
)

st.subheader('Publication export')
fmt = st.selectbox('Figure format', ['SVG', 'PDF', 'PNG', 'TIFF'], index=0)
if st.button('Prepare circular tree download', type='primary', width='stretch'):
    out_fmt = fmt.lower()
    dpi = 300 if out_fmt in {'svg', 'pdf', 'tiff'} else 240
    data = fig_bytes(fig, out_fmt, dpi)
    mime = {
        'svg': 'image/svg+xml',
        'pdf': 'application/pdf',
        'png': 'image/png',
        'tiff': 'image/tiff',
    }[out_fmt]
    st.download_button(
        f'Download circular phylogeny ({fmt})',
        data,
        f'circular_phylogeny_sequence_supported.{out_fmt}',
        mime,
        width='stretch',
    )

st.download_button('Download Newick used for circular plot', tree_text, 'tree_used_for_circular_plot.nwk', 'text/plain')
if alignment_text:
    st.download_button('Download alignment used for ML concordance', alignment_text, 'alignment_used_for_clade_concordance.fasta', 'text/plain')
if iqtree_report:
    st.download_button('Download supplied IQ-TREE report', iqtree_report, 'iqtree_report.txt', 'text/plain')

plt.close(fig)
