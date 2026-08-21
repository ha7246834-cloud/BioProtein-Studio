import json, platform
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from modules.gdm_common import parse_fasta,fasta_text,protein_qc,translate_cds,looks_ncbi_accession,newick_order,choose_order,zip_files
from modules.gdm_structure import est2genome_ready,gene_structure_batch,ncbi_structures
from modules.gdm_cdd_meme import run_cdd,collapse_domains,domain_qc,meme_ready,run_meme,motif_qc,parse_cdd,parse_meme_xml
from modules.gdm_plot import gene_structure,architecture,combined,fig_bytes

st.set_page_config(page_title='Gene Structure, Domains & Motifs | BioProtein Studio',page_icon='🧬',layout='wide')
st.title('🧬 Gene Structure, Conserved Domains & Motifs')
st.caption('GSDS-like exon–intron mapping + NCBI CDD + MEME motif discovery + transparent scientific QC.')

def file_text(f):return f.getvalue().decode('utf-8','replace') if f else ''
def empty_result():
    return dict(sequence_qc=pd.DataFrame(),cds_qc=pd.DataFrame(),gene_structures=pd.DataFrame(),gene_qc=pd.DataFrame(),domains_raw=pd.DataFrame(),domains=pd.DataFrame(),domain_qc=pd.DataFrame(),motifs=pd.DataFrame(),motif_summary=pd.DataFrame(),motif_qc=pd.DataFrame(),errors=[],warnings=[],figures={},raw_est={},ncbi_bundles={})
def validation_summary(r):
    rows=[]
    if not r['gene_qc'].empty:rows.append(dict(analysis='Gene structure',criterion='CDS coverage ≥95% and identity ≥95%; splice/frame reviewed',PASS=int((r['gene_qc'].status=='PASS').sum()),REVIEW=int((r['gene_qc'].status!='PASS').sum())))
    if not r['domain_qc'].empty:
        p=r['domain_qc'].scientific_status.str.startswith('PASS',na=False);rows.append(dict(analysis='CDD domains',criterion='Specific hits prioritized; E≤1e-5 strong; raw full output retained',PASS=int(p.sum()),REVIEW=int((~p).sum())))
    if not r['motif_qc'].empty:
        p=r['motif_qc'].status=='PASS';rows.append(dict(analysis='MEME motifs',criterion='MEME E<0.05 + ≥50% family prevalence for PASS',PASS=int(p.sum()),REVIEW=int((~p).sum())))
    return pd.DataFrame(rows)
def package(r,protein_txt,cds_txt,gen_txt,params):
    files={'inputs/proteins.fasta':protein_txt,'METHODS_AND_QC.txt':'Gene structure: NCBI reference CDS feature or EMBOSS est2genome.\nDomains: NCBI Batch CD-Search full mode.\nMotifs: MEME Suite protein mode.\nPASS/REVIEW rules are visible in the app and must be biologically reviewed.\n'}
    if cds_txt:files['inputs/cds.fasta']=cds_txt
    if gen_txt:files['inputs/genomic.fasta']=gen_txt
    for k,n in [('sequence_qc','tables/protein_sequence_qc.csv'),('cds_qc','tables/cds_translation_qc.csv'),('gene_structures','tables/gene_structures.csv'),('gene_qc','tables/gene_structure_qc.csv'),('domains_raw','tables/cdd_hits_full.csv'),('domains','tables/cdd_hits_nonredundant.csv'),('domain_qc','tables/domain_validation.csv'),('motifs','tables/motif_sites.csv'),('motif_summary','tables/motif_summary.csv'),('motif_qc','tables/motif_validation.csv')]:
        d=r.get(k)
        if isinstance(d,pd.DataFrame) and not d.empty:files[n]=d.to_csv(index=False)
    for g,t in r.get('raw_est',{}).items():files[f'raw/est2genome/{g}.txt']=t
    if r.get('cdd_raw'):files['raw/cdd_output.txt']=r['cdd_raw']
    if r.get('cdd_rid'):files['raw/cdd_search_id.txt']=r['cdd_rid']
    if r.get('meme_xml'):files['raw/meme.xml']=r['meme_xml']
    if r.get('meme_zip'):files['raw/meme_complete_output.zip']=r['meme_zip']
    for g,b in r.get('ncbi_bundles',{}).items():
        files[f'reference/{g}.cds.fasta']=f'>{g}\n{b["cds"]}\n';files[f'reference/{g}.genomic.fasta']=f'>{g}|{b["reference_record"]}\n{b["genomic"]}\n'
    for n,b in r.get('figures',{}).items():files['figures/'+n]=b
    files['run_manifest.json']=json.dumps(dict(software='BioProtein Studio',module_version='5.0',python=platform.python_version(),parameters=params,warnings=r['errors']+r['warnings']),indent=2)
    return zip_files(files)

c1,c2,c3=st.columns(3);c1.metric('MEME local','Ready' if meme_ready() else 'Not installed');c2.metric('EMBOSS est2genome','Ready' if est2genome_ready() else 'Not installed');c3.metric('CDD','NCBI remote')
if not meme_ready() or not est2genome_ready():st.info('For full one-click mode use the supplied Linux/WSL Conda environment. The app does not replace MEME or spliced alignment with a homemade predictor.')
if 'gdm_result' not in st.session_state:st.session_state.gdm_result={}
auto,imp,qctab=st.tabs(['🚀 One-click Auto Analysis','📥 Import / Re-validate','✅ Scientific QC Rules'])

with auto:
    mode=st.radio('Primary input',['Protein FASTA','CDS FASTA (auto-translate)'],horizontal=True)
    a,b=st.columns(2)
    with a:
        up=st.file_uploader('Upload primary FASTA',type=['fa','fasta','faa','fna','txt'],key='gdm_primary');paste=st.text_area('Or paste FASTA',height=180,key='gdm_paste')
    with b:
        cds_up=st.file_uploader('Optional matching CDS FASTA',type=['fa','fasta','fna','txt'],key='gdm_cds');gen_up=st.file_uploader('Optional matching genomic FASTA',type=['fa','fasta','fna','txt'],key='gdm_gen');tree_up=st.file_uploader('Optional Newick tree',type=['nwk','newick','tree','txt'],key='gdm_tree')
    with st.expander('NCBI reference auto-fetch for gene structure'):
        auto_ncbi=st.checkbox('Try NCBI when CDS + genomic FASTA are absent',True);email=st.text_input('Entrez email');api=st.text_input('NCBI API key (optional)',type='password')
        st.caption('Only FASTA IDs that look like real NCBI protein accessions are queried. Custom labels are never guessed by BLAST.')
    s1,s2,s3,s4=st.columns(4)
    with s1:do_cdd=st.checkbox('Run NCBI CDD',True);cdd_e=st.selectbox('CDD E-value',[0.01,0.001,1e-5])
    with s2:do_meme=st.checkbox('Run MEME',True);nm=st.slider('Max motifs',3,20,10)
    with s3:minw=st.number_input('Min motif width',3,50,6);maxw=st.number_input('Max motif width',6,100,50)
    with s4:model=st.selectbox('MEME model',['zoops','oops','anr'],format_func=lambda x:{'zoops':'ZOOPS — zero/one','oops':'OOPS — exactly one','anr':'ANR — any number'}[x])
    if st.button('🚀 Run complete analysis',type='primary',use_container_width=True):
        r=empty_result();primary=file_text(up) if up else paste.strip()
        if not primary:st.error('Provide FASTA input.');st.stop()
        try:
            if mode.startswith('Protein'):
                proteins=parse_fasta(primary,'protein');cds=parse_fasta(file_text(cds_up),'dna') if cds_up else {}
            else:
                cds=parse_fasta(primary,'dna');proteins,r['cds_qc']=translate_cds(cds)
            genomic=parse_fasta(file_text(gen_up),'dna') if gen_up else {}
        except Exception as e:st.error(f'Sequence parsing failed: {e}');st.stop()
        r['proteins']=proteins;r['sequence_qc']=protein_qc(proteins);protein_txt=fasta_text(proteins)
        if (r['sequence_qc'].status=='REVIEW').any():r['warnings'].append('Sequence-integrity QC flagged one or more proteins.')
        prog=st.progress(10,text='Gene structure...')
        if cds and genomic:
            try:r['gene_structures'],r['gene_qc'],r['raw_est']=gene_structure_batch(cds,genomic)
            except Exception as e:r['errors'].append('Gene structure: '+str(e))
        elif auto_ncbi and '@' in email:
            ids=[x for x in proteins if looks_ncbi_accession(x)]
            if ids:
                try:r['gene_structures'],r['gene_qc'],r['ncbi_bundles'],missing=ncbi_structures(ids,email,api or None);r['warnings'] += ([f'NCBI structure unresolved: {", ".join(missing[:20])}'] if missing else [])
                except Exception as e:r['errors'].append('NCBI structure: '+str(e))
            else:r['warnings'].append('No NCBI accession-style FASTA IDs; gene structure needs matching CDS + genomic FASTA.')
        prog.progress(40,text='Conserved domains...')
        if do_cdd:
            try:r['domains_raw'],r['cdd_rid'],r['cdd_raw']=run_cdd(proteins,float(cdd_e));r['domains']=collapse_domains(r['domains_raw']);r['domain_qc']=domain_qc(r['domains'])
            except Exception as e:r['errors'].append('CDD: '+str(e))
        prog.progress(65,text='MEME motifs...')
        if do_meme:
            try:
                m=run_meme(proteins,int(nm),int(minw),int(maxw),model);r['motifs']=m['sites'];r['motif_summary']=m['summary'];r['meme_xml']=m['xml'];r['meme_zip']=m['zip'];r['logos']=m['logos']
            except Exception as e:r['errors'].append('MEME: '+str(e))
        r['motif_qc']=motif_qc(r['motifs'],r['motif_summary'],len(proteins),r['domains']);order=choose_order(set(proteins)|set(r['gene_structures'].gene if not r['gene_structures'].empty else [])|set(r['domains'].gene if not r['domains'].empty else [])|set(r['motifs'].gene if not r['motifs'].empty else []),newick_order(file_text(tree_up)) if tree_up else [])
        r['order']=order;prog.progress(85,text='Figures and reproducibility package...')
        figs={'gene_structure':gene_structure(r['gene_structures'],order),'domain_architecture':architecture(r['domains'],order,'domain','Conserved Domain Architecture'),'motif_architecture':architecture(r['motifs'],order,'motif','Conserved Motif Architecture'),'integrated_architecture':combined(r['gene_structures'],r['domains'],r['motifs'],order)}
        for stem,fig in figs.items():
            if fig:
                for fmt in ['png','svg','pdf','tiff']:r['figures'][f'{stem}.{fmt}']=fig_bytes(fig,fmt,600)
                plt.close(fig)
        params=dict(cdd_evalue=float(cdd_e),meme_nmotifs=int(nm),meme_min_width=int(minw),meme_max_width=int(maxw),meme_model=model)
        r['package']=package(r,protein_txt,file_text(cds_up) if cds_up else (primary if mode.startswith('CDS') else ''),file_text(gen_up),params);st.session_state.gdm_result=r;prog.progress(100,text='Complete')

    r=st.session_state.gdm_result
    if r:
        for x in r.get('errors',[]):st.error(x)
        for x in r.get('warnings',[]):st.warning(x)
        st.subheader('Validation dashboard');st.dataframe(validation_summary(r),use_container_width=True,hide_index=True)
        with st.expander('Sequence QC'):
            st.dataframe(r['sequence_qc'],use_container_width=True,hide_index=True)
            if not r['cds_qc'].empty:st.dataframe(r['cds_qc'],use_container_width=True,hide_index=True)
        tabs=st.tabs(['Gene Structure','Domains','Motifs','Integrated Figure'])
        items=[(tabs[0],gene_structure(r['gene_structures'],r['order']),r['gene_qc']),(tabs[1],architecture(r['domains'],r['order'],'domain','Conserved Domain Architecture'),r['domain_qc']),(tabs[2],architecture(r['motifs'],r['order'],'motif','Conserved Motif Architecture'),r['motif_qc']),(tabs[3],combined(r['gene_structures'],r['domains'],r['motifs'],r['order']),pd.DataFrame())]
        stems=['gene_structure','domain_architecture','motif_architecture','integrated_architecture']
        for i,(tab,fig,table) in enumerate(items):
            with tab:
                if fig:
                    st.pyplot(fig);plt.close(fig)
                    if not table.empty:st.dataframe(table,use_container_width=True,hide_index=True)
                    for fmt in ['png','svg','pdf','tiff']:st.download_button(f'Download {fmt.upper()}',r['figures'][f'{stems[i]}.{fmt}'],f'{stems[i]}.{fmt}',key=f'{stems[i]}_{fmt}')
                else:st.info('No result available for this analysis.')
        st.download_button('📦 Download complete analysis package',r['package'],'BioProtein_Studio_GDM_Results.zip','application/zip',type='primary',use_container_width=True)

with imp:
    st.write('Import your previous manual NCBI CDD and MEME results and re-validate them.')
    c=st.file_uploader('CDD hit file',type=['txt','tsv'],key='imp_cdd');m=st.file_uploader('MEME XML',type=['xml'],key='imp_meme')
    if c:
        try:st.dataframe(domain_qc(collapse_domains(parse_cdd(file_text(c)))),use_container_width=True,hide_index=True)
        except Exception as e:st.error(str(e))
    if m:
        try:s,u=parse_meme_xml(file_text(m));st.dataframe(u,use_container_width=True,hide_index=True);st.dataframe(s,use_container_width=True,hide_index=True)
        except Exception as e:st.error(str(e))

with qctab:
    st.markdown('''**Gene structure:** reference CDS feature or EMBOSS `est2genome`; PASS requires ≥95% CDS coverage and ≥95% weighted exon identity in the alignment route. Non-canonical splice evidence and frame problems are flagged.\n\n**CDD domains:** NCBI specific hits are strongest. Other hits are ranked by E-value; raw full results and Search-ID are kept. Redundant highly overlapping footprints are collapsed only for visualization.\n\n**MEME motifs:** real MEME Suite is run in protein mode. ZOOPS is the default for gene families. Motif E < 0.05 is required, and the app adds a ≥50% family-prevalence rule for PASS; lower-prevalence significant motifs remain biologically reviewable.\n\n**No fabricated result:** protein sequence alone does not define introns. With custom gene IDs and no genomic/CDS reference, gene structure is deliberately left unavailable.''')

st.caption('BioProtein Studio • Gene Structure / CDD / MEME module • v5.0')
