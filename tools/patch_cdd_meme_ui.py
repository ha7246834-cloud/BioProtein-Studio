from pathlib import Path

p = Path('pages/2_Gene_Structure_Domain_Motif.py')
s = p.read_text()

old = """        prog.progress(42, text='Conserved domains...')
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
"""

new = """        prog.progress(42, text='Conserved domains...')
        if do_cdd:
            if len(proteins) > 200:
                r['warnings'].append(
                    f'CDD large-family mode: all {len(proteins)} proteins are analysed in transparent NCBI batches; '
                    'no sequences are subsampled.'
                )
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
                msg = str(e)
                if msg.startswith('CLOUD_RESOURCE_LIMIT:'):
                    r['warnings'].append('MEME cloud safety: ' + msg.split(':', 1)[1].strip())
                else:
                    r['errors'].append('MEME: ' + msg)
"""

if old not in s:
    if new in s:
        print('already patched')
    else:
        raise SystemExit('CDD/MEME UI block not found')
else:
    s = s.replace(old, new, 1)
    p.write_text(s)
    print('CDD_MEME_UI_PATCH_COMPLETE')
