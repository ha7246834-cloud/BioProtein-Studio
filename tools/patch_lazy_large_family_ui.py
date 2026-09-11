from pathlib import Path

p = Path('pages/2_Gene_Structure_Domain_Motif.py')
s = p.read_text()


def replace_once(old, new, label):
    global s
    if old in s:
        s = s.replace(old, new, 1)
        print('patched:', label)
        return
    if new in s:
        print('already patched:', label)
        return
    raise SystemExit(f'Expected block not found: {label}')

replace_once(
    "'auto': 'Auto — IQ-TREE publication mode if installed, otherwise FastTree, then NJ fallback',",
    "'auto': 'Auto — adaptive: IQ-TREE for small families, FastTree for larger families',",
    'adaptive Auto label',
)

replace_once(
    "phylo_threads = st.selectbox('IQ-TREE threads', ['AUTO', '2', '4', '8'], index=0)",
    "phylo_threads = st.selectbox('IQ-TREE threads', ['2', '4', '8', 'AUTO'], index=0)",
    'safe default IQ-TREE threads',
)

old_figs = """        prog.progress(88, text='Figures and reproducibility package...')
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
"""
new_figs = """        prog.progress(88, text='Preparing reproducibility package...')
        # Figures are deliberately rendered on demand below. Generating every
        # panel in four high-resolution formats during analysis is expensive
        # and can exhaust shared-cloud memory for large protein families.
        r['figures'] = {}
        if len(order) > 60:
            r['warnings'].append(
                f'Large-family visualization mode enabled for {len(order)} sequences. '
                'Figures are paginated and rendered on demand; the full tree/alignment remain downloadable.'
            )
"""
replace_once(old_figs, new_figs, 'remove eager multi-format rendering')

old_anchor = """        dom_plot = filter_domains(r['domains'], domain_mode)
        mot_plot = filter_motifs(r['motifs'], r['motif_qc'], motif_mode)

        with st.expander('🎨 Graph Studio — auto styles, colors & layout', expanded=True):
"""
new_anchor = """        dom_plot = filter_domains(r['domains'], domain_mode)
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
"""
replace_once(old_anchor, new_anchor, 'large-family pagination')

old_variants = """            if 'gdm_show_style_variants' not in st.session_state:
                st.session_state.gdm_show_style_variants = False
            if st.button('✨ Generate alternative color versions', key='gdm_generate_styles', width='stretch'):
                st.session_state.gdm_show_style_variants = True
"""
new_variants = """            if 'gdm_show_style_variants' not in st.session_state:
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
"""
replace_once(old_variants, new_variants, 'large-family style-preview guard')

s = s.replace(
    "combined(r.get('tree_text',''), r['gene_structures'], dom_plot, mot_plot, r['order'], style_from_preset(pname))",
    "combined(r.get('tree_text',''), r['gene_structures'], dom_plot, mot_plot, display_order, style_from_preset(pname))",
)

start_marker = "        tabs = st.tabs(['Phylogeny', 'Gene Structure', 'Domains', 'Motifs', 'Integrated Figure'])"
end_marker = "\n        st.download_button('📦 Download complete analysis package'"
start = s.find(start_marker)
end = s.find(end_marker, start)

if start >= 0 and end >= 0:
    lazy = """        st.subheader('Figure viewer')

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

            # Prepare only the displayed figure. SVG/PDF are ideal for publication;
            # PNG is capped by gdm_plot.fig_bytes to avoid giant allocations.
            svg_bytes = fig_bytes(selected_fig, 'svg', 300)
            pdf_bytes = fig_bytes(selected_fig, 'pdf', 300)
            png_bytes = fig_bytes(selected_fig, 'png', 220)

            d1, d2, d3 = st.columns(3)
            with d1:
                st.download_button('Download SVG', svg_bytes, f'{stem}{suffix}.svg', 'image/svg+xml', key=key_base + '_svg', width='stretch')
            with d2:
                st.download_button('Download PDF', pdf_bytes, f'{stem}{suffix}.pdf', 'application/pdf', key=key_base + '_pdf', width='stretch')
            with d3:
                st.download_button('Download PNG', png_bytes, f'{stem}{suffix}.png', 'image/png', key=key_base + '_png', width='stretch')

            if st.checkbox('Prepare TIFF download (slower)', False, key='gdm_prepare_tiff'):
                tiff_bytes = fig_bytes(selected_fig, 'tiff', 300)
                st.download_button(
                    'Download TIFF', tiff_bytes, f'{stem}{suffix}.tiff', 'image/tiff',
                    key=key_base + '_tiff'
                )
            plt.close(selected_fig)
        else:
            st.info('No result available for this figure.')

        st.caption(
            'Cloud-safe rendering: one figure is generated at a time. Large families are paginated for display, '
            'while the complete analysis, Newick tree, and alignment remain preserved.'
        )
"""
    s = s[:start] + lazy + s[end:]
elif "st.subheader('Figure viewer')" in s:
    print('already patched: lazy figure viewer')
else:
    raise SystemExit('Expected figure viewer block not found')

p.write_text(s)
print('UI_PATCH_COMPLETE')
