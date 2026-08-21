import json
import platform
from pathlib import Path
from datetime import datetime

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from modules.gdm_common import parse_fasta, fasta_text, protein_qc, translate_cds, looks_ncbi_accession, newick_order, choose_order, zip_files
from modules.gdm_structure import est2genome_ready, gene_structure_batch, ncbi_structures, parse_gene_structure_annotation
from modules.gdm_cdd_meme import run_cdd, collapse_domains, domain_qc, meme_ready, run_meme, motif_qc, parse_cdd, parse_meme_xml
from modules.gdm_phylogeny import build_phylogeny, external_phylogeny_ready, publication_phylogeny_ready
from modules.gdm_plot import gene_structure, missing_structure_figure, architecture, combined, phylogeny_figure, fig_bytes
from modules.gdm_style import STYLE_PRESETS, style_from_preset, assign_colors
from modules.gdm_reference import auto_resolve_gene_structure, auto_reference_ready, datasets_ready, miniprot_ready

st.set_page_config(page_title='Gene Structure, Domains & Motifs | BioProtein Studio', page_icon='🧬', layout='wide')
st.title('🧬 Gene Structure, Domains, Motifs & Phylogeny')
st.caption('Generic multi-gene workflow: reference-backed exon/intron mapping + NCBI CDD + MEME + phylogeny + publication figures.')


def file_text(f):
    return f.getvalue().decode('utf-8', 'replace') if f else ''


def empty_result():
    return dict(
        sequence_qc=pd.DataFrame(), cds_qc=pd.DataFrame(), gene_structures=pd.DataFrame(), gene_qc=pd.DataFrame(),
        domains_raw=pd.DataFrame(), domains=pd.DataFrame(), domain_qc=pd.DataFrame(), motifs=pd.DataFrame(),
        motif_summary=pd.DataFrame(), motif_qc=pd.DataFrame(), phylogeny_qc=pd.DataFrame(),
        reference_meta=pd.DataFrame(), reference_mapping=pd.DataFrame(), annotation_reconciliation=pd.DataFrame(),
        errors=[], warnings=[], figures={}, raw_est={}, ncbi_bundles={}, auto_reference_bundles={},
        tree_text='', alignment_text='', phylogeny_method='', phylogeny_log='', iqtree_report='', phylogeny_command='',
        miniprot_raw='', order=[], package=b'', autosave_path=''
    )


def validation_summary(r):
    rows=[]
    q=r.get('phylogeny_qc',pd.DataFrame())
    if not q.empty:
        x=q.iloc[0]; ok=str(x.get('status','')).startswith('PASS')
        rows.append(dict(analysis='Phylogeny', criterion=f"{x.get('method','')}; {x.get('support','')}", PASS=int(ok), REVIEW=int(not ok)))
    if not r['gene_qc'].empty:
        p=r['gene_qc'].status.astype(str).str.startswith('PASS')
        rows.append(dict(analysis='Gene structure', criterion='Reference-backed mapping; coverage/identity/splice/frame reviewed', PASS=int(p.sum()), REVIEW=int((~p).sum())))
    if not r['domain_qc'].empty:
        p=r['domain_qc'].scientific_status.astype(str).str.startswith('PASS')
        rows.append(dict(analysis='CDD domains', criterion='Specific hits prioritized; raw CDD retained', PASS=int(p.sum()), REVIEW=int((~p).sum())))
    if not r['motif_qc'].empty:
        p=r['motif_qc'].status.eq('PASS')
        rows.append(dict(analysis='MEME motifs', criterion='MEME E<0.05 + >=50% family prevalence for PASS', PASS=int(p.sum()), REVIEW=int((~p).sum())))
    return pd.DataFrame(rows)


def filter_domains(df, mode):
    if df is None or df.empty: return pd.DataFrame()
    if mode=='all': return df.copy()
    s=df[df.hit_type.str.contains('specific',case=False,na=False) & ~df.hit_type.str.contains('non',case=False,na=False)]
    return s if not s.empty else df.copy()


def filter_motifs(df, qc, mode):
    if df is None or df.empty: return pd.DataFrame()
    if mode=='all' or qc is None or qc.empty: return df.copy()
    keep=set(qc[qc.status.eq('PASS')].motif)
    s=df[df.motif.isin(keep)]
    return s if not s.empty else df.copy()


def make_package(r, protein_txt, cds_txt, genomic_txt, structure_txt, params):
    files={
        'inputs/proteins.fasta':protein_txt,
        'METHODS_AND_QC.txt':(
            'Gene structure is reference-backed only: imported annotation, NCBI CDS feature, EMBOSS est2genome, or Protein+Species automatic reference mapping.\n'
            'Automatic mapping: NCBI Datasets reference genome + miniprot + NCBI GFF reconciliation; conservative splice rescue remains REVIEW.\n'
            'Domains: NCBI Batch CD-Search. Motifs: MEME Suite. Phylogeny: uploaded Newick or MAFFT+IQ-TREE/FastTree with NJ fallback.\n'
            'PASS/REVIEW output must be biologically reviewed before publication.\n'
        )
    }
    if cds_txt: files['inputs/cds.fasta']=cds_txt
    if genomic_txt: files['inputs/genomic.fasta']=genomic_txt
    if structure_txt: files['inputs/gene_structure_annotation.txt']=structure_txt
    if r.get('tree_text'): files['phylogeny/tree.nwk']=r['tree_text']
    if r.get('alignment_text'): files['phylogeny/alignment.fasta']=r['alignment_text']
    if r.get('phylogeny_log'): files['phylogeny/phylogeny.log.txt']=r['phylogeny_log']
    if r.get('iqtree_report'): files['phylogeny/iqtree_report.txt']=r['iqtree_report']
    if r.get('phylogeny_command'): files['phylogeny/command.txt']=r['phylogeny_command']
    tables=[
        ('sequence_qc','protein_sequence_qc.csv'),('cds_qc','cds_translation_qc.csv'),('phylogeny_qc','phylogeny_qc.csv'),
        ('gene_structures','gene_structures.csv'),('gene_qc','gene_structure_qc.csv'),('domains_raw','cdd_hits_full.csv'),
        ('domains','cdd_hits_nonredundant.csv'),('domain_qc','domain_validation.csv'),('motifs','motif_sites.csv'),
        ('motif_summary','motif_summary.csv'),('motif_qc','motif_validation.csv'),('reference_meta','auto_reference_metadata.csv'),
        ('reference_mapping','auto_reference_mapping_qc.csv'),('annotation_reconciliation','annotation_reconciliation.csv')]
    for key,name in tables:
        d=r.get(key)
        if isinstance(d,pd.DataFrame) and not d.empty: files['tables/'+name]=d.to_csv(index=False)
    for g,t in r.get('raw_est',{}).items(): files[f'raw/est2genome/{g}.txt']=t
    if r.get('cdd_raw'): files['raw/cdd_output.txt']=r['cdd_raw']
    if r.get('cdd_rid'): files['raw/cdd_search_id.txt']=r['cdd_rid']
    if r.get('meme_xml'): files['raw/meme.xml']=r['meme_xml']
    if r.get('meme_zip'): files['raw/meme_complete_output.zip']=r['meme_zip']
    if r.get('miniprot_raw'): files['raw/miniprot.gff3']=r['miniprot_raw']
    bundles=r.get('auto_reference_bundles',{})
    if bundles:
        files['reference/auto_cds.fasta']=''.join(f'>{g}\n{b.get("cds","")}\n' for g,b in bundles.items())
        files['reference/auto_genomic.fasta']=''.join(f'>{g}|{b.get("reference_record","")}:{b.get("reference_start","")}-{b.get("reference_end","")}({b.get("strand","")})\n{b.get("genomic","")}\n' for g,b in bundles.items())
    for g,b in r.get('ncbi_bundles',{}).items():
        files[f'reference/{g}.cds.fasta']=f'>{g}\n{b["cds"]}\n'
        files[f'reference/{g}.genomic.fasta']=f'>{g}|{b["reference_record"]}\n{b["genomic"]}\n'
    for n,b in r.get('figures',{}).items(): files['figures/'+n]=b
    files['run_manifest.json']=json.dumps(dict(software='BioProtein Studio',module_version='5.5.0',python=platform.python_version(),parameters=params,warnings=r['errors']+r['warnings']),indent=2)
    return zip_files(files)


m1,m2,m3,m4,m5=st.columns(5)
m1.metric('Phylogeny','IQ-TREE Ready' if publication_phylogeny_ready() else ('FastTree Ready' if external_phylogeny_ready() else 'NJ fallback'))
m2.metric('MEME','Ready' if meme_ready() else 'Missing')
m3.metric('miniprot','Ready' if miniprot_ready() else 'Missing')
m4.metric('NCBI Datasets','Ready' if datasets_ready() else 'Missing')
m5.metric('CDD','NCBI remote')
if not publication_phylogeny_ready(): st.info('Publication phylogeny requires MAFFT + IQ-TREE. Auto mode can fall back to FastTree or NJ.')
if not auto_reference_ready(): st.info('Protein + Species gene-structure automation requires miniprot + NCBI Datasets CLI. Manual reference routes remain available.')
if 'gdm_result' not in st.session_state: st.session_state.gdm_result={}

auto,imp,qctab=st.tabs(['🚀 One-click Auto Analysis','📥 Import / Re-validate','✅ Scientific QC Rules'])

with auto:
    mode=st.radio('Primary input',['Protein FASTA','CDS FASTA (auto-translate)'],horizontal=True)
    a,b=st.columns(2)
    with a:
        up=st.file_uploader('Upload primary FASTA',type=['fa','fasta','faa','fna','txt'],key='gdm_primary')
        paste=st.text_area('Or paste FASTA',height=170,key='gdm_paste')
    with b:
        cds_up=st.file_uploader('Optional matching CDS FASTA',type=['fa','fasta','fna','txt'],key='gdm_cds')
        gen_up=st.file_uploader('Optional matching genomic FASTA',type=['fa','fasta','fna','txt'],key='gdm_gen')
        structure_up=st.file_uploader('Optional exon/CDS annotation (CSV/TSV/GFF3/GTF)',type=['csv','tsv','txt','gff','gff3','gtf'],key='gdm_structure')
        tree_up=st.file_uploader('Optional Newick tree',type=['nwk','newick','tree','txt'],key='gdm_tree')

    with st.expander('🧭 Automatic gene structure from Protein + Species',expanded=True):
        auto_reference=st.checkbox('Resolve reference genome automatically when other gene-structure inputs are absent',True)
        r1,r2=st.columns(2)
        with r1: taxon=st.text_input('Species / NCBI taxon',placeholder='e.g. Carica papaya')
        with r2: assembly=st.text_input('Optional assembly accession',placeholder='GCF_... or GCA_...')
        ref_threads=st.slider('Reference mapping threads',1,16,4)
        st.caption('NCBI Datasets → reference genome/GFF3 → miniprot → annotation reconciliation → translation QC. Ambiguous or computationally rescued structures remain REVIEW.')

    with st.expander('NCBI accession route'):
        auto_ncbi=st.checkbox('Try Entrez for real NCBI protein accessions',True)
        email=st.text_input('Entrez email')
        api=st.text_input('NCBI API key (optional)',type='password')

    with st.expander('🌳 Phylogeny settings',expanded=True):
        auto_tree=st.checkbox('Automatically build tree when Newick is absent',True)
        phylo_mode=st.selectbox('Tree method',['auto','publication','fasttree','nj'],format_func=lambda x:{'auto':'Auto — IQ-TREE → FastTree → NJ','publication':'Publication — MAFFT + IQ-TREE + ModelFinder + SH-aLRT + ultrafast bootstrap','fasttree':'Fast screening — MAFFT + FastTree','nj':'Neighbor-Joining — screening fallback'}[x])
        p1,p2,p3=st.columns(3)
        with p1: bootstrap=st.selectbox('Ultrafast bootstrap',[1000,2000,5000])
        with p2: alrt=st.selectbox('SH-aLRT',[1000,2000,5000])
        with p3: phy_threads=st.selectbox('IQ-TREE threads',['AUTO','2','4','8'])

    s1,s2,s3,s4=st.columns(4)
    with s1:
        do_cdd=st.checkbox('Run NCBI CDD',True); cdd_e=st.selectbox('CDD E-value',[0.01,0.001,1e-5])
    with s2:
        do_meme=st.checkbox('Run MEME',True); nm=st.slider('Max motifs',3,20,10)
    with s3:
        minw=st.number_input('Min motif width',3,50,6); maxw=st.number_input('Max motif width',6,100,50)
    with s4:
        model=st.selectbox('MEME model',['zoops','oops','anr'],format_func=lambda x:{'zoops':'ZOOPS — zero/one','oops':'OOPS — exactly one','anr':'ANR — any number'}[x])

    if st.button('🚀 Run complete analysis',type='primary',width='stretch'):
        r=empty_result(); primary=file_text(up) if up else paste.strip(); uploaded_tree=file_text(tree_up)
        if not primary: st.error('Provide FASTA input.'); st.stop()
        try:
            if mode.startswith('Protein'):
                proteins=parse_fasta(primary,'protein'); cds=parse_fasta(file_text(cds_up),'dna') if cds_up else {}
            else:
                cds=parse_fasta(primary,'dna'); proteins,r['cds_qc']=translate_cds(cds)
            genomic=parse_fasta(file_text(gen_up),'dna') if gen_up else {}
        except Exception as e: st.error(f'Sequence parsing failed: {e}'); st.stop()
        r['proteins']=proteins; r['sequence_qc']=protein_qc(proteins); protein_txt=fasta_text(proteins)
        prog=st.progress(5,text='Phylogeny...')

        if uploaded_tree:
            r['tree_text']=uploaded_tree; r['phylogeny_method']='Uploaded Newick'
            r['phylogeny_qc']=pd.DataFrame([dict(method='Uploaded Newick',sequences=len(proteins),model='User supplied',support='Preserved as supplied',status='REFERENCE')])
        elif auto_tree:
            try:
                tr=build_phylogeny(proteins,mode=phylo_mode,bootstrap=bootstrap,alrt=alrt,threads=phy_threads)
                r['tree_text']=tr.get('tree_text',''); r['alignment_text']=tr.get('alignment_text',''); r['phylogeny_method']=tr.get('method',''); r['phylogeny_qc']=tr.get('qc',pd.DataFrame()); r['phylogeny_log']=tr.get('log_text',''); r['iqtree_report']=tr.get('iqtree_report',''); r['phylogeny_command']=tr.get('command','')
                if tr.get('warning'): r['warnings'].append(tr['warning'])
            except Exception as e: r['errors'].append('Phylogeny: '+str(e))

        prog.progress(25,text='Gene structure...')
        structure_text=file_text(structure_up)
        if structure_text:
            try: r['gene_structures'],r['gene_qc']=parse_gene_structure_annotation(structure_text,proteins.keys())
            except Exception as e: r['errors'].append('Gene structure annotation: '+str(e))
        elif cds and genomic:
            try: r['gene_structures'],r['gene_qc'],r['raw_est']=gene_structure_batch(cds,genomic)
            except Exception as e: r['errors'].append('Gene structure: '+str(e))
        elif auto_ncbi and '@' in email:
            ids=[g for g in proteins if looks_ncbi_accession(g)]
            if ids:
                try:
                    r['gene_structures'],r['gene_qc'],r['ncbi_bundles'],missing=ncbi_structures(ids,email,api or None)
                    if missing: r['warnings'].append('NCBI structure unresolved: '+', '.join(missing[:20]))
                except Exception as e: r['errors'].append('NCBI structure: '+str(e))
        if r['gene_structures'].empty and auto_reference and (taxon.strip() or assembly.strip()):
            try:
                ar=auto_resolve_gene_structure(proteins,taxon=taxon,assembly_accession=assembly,api_key=api or '',threads=ref_threads)
                r['gene_structures']=ar['structures']; r['gene_qc']=ar['qc']; r['auto_reference_bundles']=ar['bundles']; r['reference_mapping']=ar['mapping_qc']; r['reference_meta']=ar['reference_table']; r['annotation_reconciliation']=ar.get('annotation_reconciliation',pd.DataFrame()); r['miniprot_raw']=ar['raw_gff']
                acc=ar['reference'].get('accession',''); r['warnings'].append(f'Automatic gene structure used NCBI reference assembly {acc}. Review mapping QC before publication.')
                nrec=int((r['annotation_reconciliation'].decision=='USED').sum()) if not r['annotation_reconciliation'].empty and 'decision' in r['annotation_reconciliation'] else 0
                nres=int(r['gene_qc'].structure_source.astype(str).str.contains('rescue',case=False,na=False).sum()) if not r['gene_qc'].empty and 'structure_source' in r['gene_qc'] else 0
                if nrec: r['warnings'].append(f'Reference annotation reconciled {nrec} mapped gene structure(s) against NCBI GFF3.')
                if nres: r['warnings'].append(f'{nres} gene structure(s) required computational splice rescue and remain REVIEW until annotation/manual confirmation.')
            except Exception as e: r['errors'].append('Automatic reference mapping: '+str(e))
        if r['gene_structures'].empty: r['warnings'].append('Gene structure unavailable for this run. Protein sequence alone cannot define exon/intron boundaries.')

        prog.progress(50,text='NCBI conserved domains...')
        if do_cdd:
            try: r['domains_raw'],r['cdd_rid'],r['cdd_raw']=run_cdd(proteins,float(cdd_e)); r['domains']=collapse_domains(r['domains_raw']); r['domain_qc']=domain_qc(r['domains'])
            except Exception as e: r['errors'].append('CDD: '+str(e))
        prog.progress(70,text='MEME motifs...')
        if do_meme:
            try:
                mm=run_meme(proteins,int(nm),int(minw),int(maxw),model); r['motifs']=mm['sites']; r['motif_summary']=mm['summary']; r['meme_xml']=mm['xml']; r['meme_zip']=mm['zip']
            except Exception as e: r['errors'].append('MEME: '+str(e))
        r['motif_qc']=motif_qc(r['motifs'],r['motif_summary'],len(proteins),r['domains'])

        preferred=newick_order(r['tree_text']) if r['tree_text'] else []
        genes=set(proteins)|set(r['gene_structures'].gene if not r['gene_structures'].empty else [])|set(r['domains'].gene if not r['domains'].empty else [])|set(r['motifs'].gene if not r['motifs'].empty else [])
        r['order']=choose_order(genes,preferred)
        prog.progress(85,text='Publication figures and package...')
        dom_plot=filter_domains(r['domains'],'specific'); mot_plot=filter_motifs(r['motifs'],r['motif_qc'],'pass'); default_style=style_from_preset('Journal Classic')
        specs={
            'phylogenetic_tree':phylogeny_figure(r['tree_text'],r['order'],default_style),
            'gene_structure':gene_structure(r['gene_structures'],r['order'],default_style) if not r['gene_structures'].empty else missing_structure_figure(r['order'],default_style),
            'domain_architecture':architecture(dom_plot,r['order'],'domain','Conserved Domain Architecture',style=default_style),
            'motif_architecture':architecture(mot_plot,r['order'],'motif','Conserved Motif Architecture',style=default_style),
            'integrated_architecture':combined(r['tree_text'],r['gene_structures'],dom_plot,mot_plot,r['order'],default_style)}
        for stem,fig in specs.items():
            if fig:
                for fmt in ['png','svg','pdf','tiff']: r['figures'][f'{stem}.{fmt}']=fig_bytes(fig,fmt,600)
                plt.close(fig)
        params=dict(cdd_evalue=float(cdd_e),meme_nmotifs=int(nm),meme_min_width=int(minw),meme_max_width=int(maxw),meme_model=model,phylogeny_mode=phylo_mode,bootstrap=int(bootstrap),alrt=int(alrt),reference_taxon=taxon,reference_assembly=assembly)
        cds_txt=file_text(cds_up) if cds_up else (primary if mode.startswith('CDS') else '')
        r['package']=make_package(r,protein_txt,cds_txt,file_text(gen_up),structure_text,params)
        try:
            outdir=Path.cwd()/'results'; outdir.mkdir(exist_ok=True); latest=outdir/'BioProtein_Studio_GDM_Results_v5_5_latest.zip'; latest.write_bytes(r['package']); r['autosave_path']=str(latest)
        except Exception as e: r['warnings'].append('Local auto-save failed: '+str(e))
        st.session_state.gdm_result=r; prog.progress(100,text='Complete')

    r=st.session_state.gdm_result
    if r:
        for x in r.get('errors',[]): st.error(x)
        for x in r.get('warnings',[]): st.warning(x)
        st.subheader('Validation dashboard'); st.dataframe(validation_summary(r),width='stretch',hide_index=True)
        with st.expander('Sequence / reference QC'):
            st.dataframe(r['sequence_qc'],width='stretch',hide_index=True)
            if not r['gene_qc'].empty: st.dataframe(r['gene_qc'],width='stretch',hide_index=True)
            if not r['reference_mapping'].empty: st.dataframe(r['reference_mapping'],width='stretch',hide_index=True)
            if not r['annotation_reconciliation'].empty: st.dataframe(r['annotation_reconciliation'],width='stretch',hide_index=True)

        bundles=r.get('auto_reference_bundles',{})
        if bundles:
            c1,c2=st.columns(2)
            auto_cds=''.join(f'>{g}\n{b.get("cds","")}\n' for g,b in bundles.items())
            auto_gen=''.join(f'>{g}|{b.get("reference_record","")}:{b.get("reference_start","")}-{b.get("reference_end","")}({b.get("strand","")})\n{b.get("genomic","")}\n' for g,b in bundles.items())
            with c1: st.download_button('Download auto-generated CDS FASTA',auto_cds,'auto_cds.fasta','text/plain',width='stretch')
            with c2: st.download_button('Download auto-generated genomic FASTA',auto_gen,'auto_genomic.fasta','text/plain',width='stretch')

        v1,v2=st.columns(2)
        with v1: domain_mode=st.radio('Domain display',['specific','all'],horizontal=True,format_func=lambda x:'Specific hits only' if x=='specific' else 'All retained hits')
        with v2: motif_mode=st.radio('Motif display',['pass','all'],horizontal=True,format_func=lambda x:'PASS motifs only' if x=='pass' else 'All motifs')
        dom_plot=filter_domains(r['domains'],domain_mode); mot_plot=filter_motifs(r['motifs'],r['motif_qc'],motif_mode)

        with st.expander('🎨 Graph Studio — auto styles, colors & layout',expanded=True):
            st.caption('Generate publication-ready alternatives, choose a preset, then optionally fine-tune colors and sizing. Figure tabs update live.')
            if st.button('✨ Generate alternative color versions',key='gdm_generate_styles',width='stretch'):
                st.session_state.gdm_show_styles=True
            if st.session_state.get('gdm_show_styles',False):
                cols=st.columns(2)
                for i,pname in enumerate(list(STYLE_PRESETS)[:4]):
                    with cols[i%2]:
                        pf=combined(r['tree_text'],r['gene_structures'],dom_plot,mot_plot,r['order'],style_from_preset(pname))
                        if pf: st.caption(pname); st.pyplot(pf,width='stretch'); plt.close(pf)
            preset=st.radio('Use figure style',list(STYLE_PRESETS.keys()),horizontal=True,key='gdm_style_preset')
            graph_style=style_from_preset(preset)
            if st.checkbox('Fine-tune colors and figure sizing',False,key='gdm_fine_tune'):
                c1,c2,c3,c4=st.columns(4)
                with c1:
                    graph_style['tree_color']=st.color_picker('Tree branches',graph_style['tree_color']); graph_style['support_color']=st.color_picker('Support labels',graph_style['support_color'])
                with c2:
                    graph_style['exon_color']=st.color_picker('Validated exons',graph_style['exon_color']); graph_style['review_exon_color']=st.color_picker('REVIEW exons',graph_style['review_exon_color'])
                with c3:
                    graph_style['intron_color']=st.color_picker('Introns',graph_style['intron_color']); graph_style['backbone_color']=st.color_picker('Protein backbone',graph_style['backbone_color'])
                with c4:
                    graph_style['edge_color']=st.color_picker('Feature borders',graph_style['edge_color']); graph_style['title_size']=st.slider('Title size',10,20,int(graph_style['title_size']))
                graph_style['label_size']=st.slider('Gene label size',6,14,int(graph_style['label_size'])); graph_style['line_width']=st.slider('Line width',0.5,2.5,float(graph_style['line_width']),0.1)
                ml=list(dict.fromkeys(mot_plot.motif.astype(str))) if not mot_plot.empty else []; dl=list(dict.fromkeys(dom_plot.domain.astype(str))) if not dom_plot.empty else []
                mc=assign_colors(ml,graph_style['motif_palette']); dc=assign_colors(dl,graph_style['domain_palette'])
                if ml:
                    st.markdown('**Motif colors**'); cs=st.columns(min(5,len(ml)))
                    for i,lab in enumerate(ml):
                        with cs[i%len(cs)]: graph_style['motif_colors'][lab]=st.color_picker(lab,mc[lab],key='motif_color_'+lab)
                if dl:
                    st.markdown('**Domain colors**'); cs=st.columns(min(4,len(dl)))
                    for i,lab in enumerate(dl):
                        with cs[i%len(cs)]: graph_style['domain_colors'][lab]=st.color_picker(lab,dc[lab],key='domain_color_'+lab)

        tabs=st.tabs(['Phylogeny','Gene Structure','Domains','Motifs','Integrated Figure'])
        figs=[
            phylogeny_figure(r['tree_text'],r['order'],graph_style),
            gene_structure(r['gene_structures'],r['order'],graph_style) if not r['gene_structures'].empty else missing_structure_figure(r['order'],graph_style),
            architecture(dom_plot,r['order'],'domain','Conserved Domain Architecture',style=graph_style),
            architecture(mot_plot,r['order'],'motif','Conserved Motif Architecture',style=graph_style),
            combined(r['tree_text'],r['gene_structures'],dom_plot,mot_plot,r['order'],graph_style)]
        tables=[r['phylogeny_qc'],r['gene_qc'],r['domain_qc'],r['motif_qc'],pd.DataFrame()]
        stems=['phylogenetic_tree','gene_structure','domain_architecture','motif_architecture','integrated_architecture']
        styled={}
        for i,(tab,fig,table) in enumerate(zip(tabs,figs,tables)):
            with tab:
                if fig:
                    st.pyplot(fig,width='stretch'); bytes_now={fmt:fig_bytes(fig,fmt,600) for fmt in ['png','svg','pdf','tiff']}; plt.close(fig)
                    if not table.empty: st.dataframe(table,width='stretch',hide_index=True)
                    for fmt,b in bytes_now.items(): st.download_button(f'Download {fmt.upper()}',b,f'{stems[i]}.{fmt}',key=f'{stems[i]}_{fmt}_{domain_mode}_{motif_mode}')
                    for fmt,b in bytes_now.items(): styled[f'{stems[i]}.{fmt}']=b
        styled['GRAPH_STYLE.json']=json.dumps(graph_style,indent=2)
        st.download_button('🎨 Download customized figure set',zip_files(styled),'BioProtein_Studio_Custom_Figures.zip','application/zip',width='stretch')
        st.download_button('📦 Download complete analysis package',r['package'],'BioProtein_Studio_GDM_Results_v5_5.zip','application/zip',type='primary',width='stretch')
        if r.get('autosave_path'): st.caption('Auto-saved locally: '+r['autosave_path'])

with imp:
    st.write('Import previous NCBI CDD and MEME results and re-validate them.')
    c=st.file_uploader('CDD hit file',type=['txt','tsv'],key='imp_cdd'); m=st.file_uploader('MEME XML',type=['xml'],key='imp_meme')
    if c:
        try: st.dataframe(domain_qc(collapse_domains(parse_cdd(file_text(c)))),width='stretch',hide_index=True)
        except Exception as e: st.error(str(e))
    if m:
        try:
            sites,summary=parse_meme_xml(file_text(m)); st.dataframe(summary,width='stretch',hide_index=True); st.dataframe(sites,width='stretch',hide_index=True)
        except Exception as e: st.error(str(e))

with qctab:
    st.markdown('''**Gene structure:** never inferred from protein alone. Reference CDS/genomic sequence, annotation, NCBI accession, or automatic Protein+Species reference mapping is required. Computational splice rescue is marked REVIEW.\n\n**Phylogeny:** Publication mode uses MAFFT + IQ-TREE + ModelFinder + SH-aLRT + ultrafast bootstrap. FastTree and NJ are screening fallbacks.\n\n**CDD domains:** NCBI specific hits are prioritized; raw full output and Search-ID are retained.\n\n**MEME motifs:** real MEME Suite is used. Motif E < 0.05 and >=50% family prevalence are required for PASS by default.\n\n**Graph Studio:** visual styling changes presentation only; it never changes the underlying biological coordinates, hits, motifs, or tree.''')

st.caption('BioProtein Studio • Phylogeny / Gene Structure / CDD / MEME module • v5.5 Research Release')
