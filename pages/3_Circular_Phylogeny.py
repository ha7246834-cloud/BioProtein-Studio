import io

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from modules.gdm_circular import circular_phylogeny_figure
from modules.gdm_plot import fig_bytes
from modules.gdm_style import STYLE_PRESETS, style_from_preset


st.set_page_config(page_title='Circular Phylogeny | BioProtein Studio', page_icon='🌳', layout='wide')
st.title('🌳 Circular Phylogeny & Automatic Clade View')
st.caption(
    'Gene-family-agnostic circular tree visualization with topology-derived clade colours, '
    'branch lengths, support labels, and publication-ready export.'
)

r = st.session_state.get('gdm_result') or {}
session_tree = str(r.get('tree_text', '') or '').strip()
session_order = list(r.get('order') or [])
session_method = str(r.get('phylogeny_method', '') or '').strip()
session_qc = r.get('phylogeny_qc')

st.info(
    'The circular plot does not invent IQ-TREE statistics. A FastTree tree displays FastTree local support. '
    'An IQ-TREE tree/report can display its SH-aLRT/UFBoot support and model evidence. '
    'Automatic clade colours are topology-derived display groups and must not be presented as validated '
    'functional/evolutionary subgroups without biological review.'
)

source_options = []
if session_tree:
    source_options.append('Latest GDM analysis')
source_options.append('Upload / paste Newick')
source = st.radio('Tree source', source_options, horizontal=True)

uploaded_tree = None
pasted_tree = ''
iqtree_report = ''
method_label = session_method
order = session_order if source == 'Latest GDM analysis' else None

if source == 'Latest GDM analysis':
    tree_text = session_tree
    if isinstance(session_qc, pd.DataFrame) and not session_qc.empty:
        st.markdown('**Inference evidence from the latest GDM run**')
        st.dataframe(session_qc, width='stretch', hide_index=True)
    if 'FastTree' in method_label:
        st.warning(
            'Current tree is a Cloud screening tree (MAFFT + FastTree). Branch lengths and FastTree local '
            'support are preserved, but ModelFinder, SH-aLRT and ultrafast bootstrap were not computed.'
        )
    elif 'IQ-TREE' in method_label:
        st.success('Current tree contains IQ-TREE publication-oriented inference evidence from the GDM run.')
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

if not tree_text:
    st.warning('Run the GDM analysis first or provide a Newick tree to build the circular figure.')
    st.stop()

st.subheader('Circular tree options')
a, b, c, d = st.columns(4)
with a:
    group_choice = st.selectbox('Automatic clade groups', ['Auto'] + list(range(2, 13)))
    target_groups = None if group_choice == 'Auto' else int(group_choice)
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

try:
    fig, clade_table = circular_phylogeny_figure(
        tree_text,
        order=order or None,
        style=style,
        target_groups=target_groups,
        midpoint_root=midpoint,
        show_support=show_support,
        method_label=method_label,
    )
except Exception as exc:
    st.error(f'Circular phylogeny could not be rendered: {exc}')
    st.stop()

st.subheader('Automatic circular clade figure')
st.pyplot(fig, width='stretch')

st.markdown('**Topology-derived clade table**')
st.dataframe(clade_table, width='stretch', hide_index=True)
st.download_button(
    'Download clade membership CSV',
    clade_table.to_csv(index=False),
    'automatic_topology_clades.csv',
    'text/csv',
    width='stretch',
)

st.caption(
    'Auto Clade A/B/C… labels are display partitions obtained only from tree topology. '
    'Rename/interpret clades biologically only after reference-gene and support-based validation.'
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
        f'circular_phylogeny_auto_clades.{out_fmt}',
        mime,
        width='stretch',
    )

if source == 'Upload / paste Newick':
    st.download_button('Download supplied Newick', tree_text, 'tree_used_for_circular_plot.nwk', 'text/plain')
    if iqtree_report:
        st.download_button('Download supplied IQ-TREE report', iqtree_report, 'iqtree_report.txt', 'text/plain')

plt.close(fig)
