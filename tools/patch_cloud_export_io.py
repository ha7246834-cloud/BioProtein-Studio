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

old_autosave = """        try:
            outdir = Path.cwd() / 'results'
            outdir.mkdir(parents=True, exist_ok=True)
            latest = outdir / 'BioProtein_Studio_GDM_Results_v5_5_latest.zip'
            latest.write_bytes(r['package'])
            stamped = outdir / f'BioProtein_Studio_GDM_Results_v5_5_{datetime.now().strftime(\"%Y%m%d_%H%M%S\")}.zip'
            stamped.write_bytes(r['package'])
            r['autosave_path'] = str(latest)
        except Exception as e:
            r['warnings'].append('Automatic result-package save failed: ' + str(e))
"""
new_autosave = """        # Streamlit Community Cloud uses ephemeral storage. Avoid writing duplicate
        # ZIP archives on every run there; local/WSL users still receive autosave.
        if not Path('/mount/src').exists():
            try:
                outdir = Path.cwd() / 'results'
                outdir.mkdir(parents=True, exist_ok=True)
                latest = outdir / 'BioProtein_Studio_GDM_Results_v5_5_latest.zip'
                latest.write_bytes(r['package'])
                stamped = outdir / f'BioProtein_Studio_GDM_Results_v5_5_{datetime.now().strftime(\"%Y%m%d_%H%M%S\")}.zip'
                stamped.write_bytes(r['package'])
                r['autosave_path'] = str(latest)
            except Exception as e:
                r['warnings'].append('Automatic result-package save failed: ' + str(e))
"""
replace_once(old_autosave, new_autosave, 'cloud autosave guard')

old_exports = """            # Prepare only the displayed figure. SVG/PDF are ideal for publication;
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
"""
new_exports = """            export_fmt = st.selectbox(
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
"""
replace_once(old_exports, new_exports, 'single-format on-demand export')

p.write_text(s)
print('CLOUD_EXPORT_IO_PATCH_COMPLETE')
