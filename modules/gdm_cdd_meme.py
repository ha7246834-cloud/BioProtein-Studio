from __future__ import annotations
import io,re,shutil,subprocess,tempfile,time,zipfile
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np,pandas as pd,requests
from .gdm_common import fasta_text,norm_id,read_text
CDD_URL='https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi'

def _rid(t):
    m=re.search(r'(QM\d?-qcdsearch-[A-Za-z0-9\-]+)',t);return m.group(1) if m else None

def _status(t):
    m=re.search(r'#status\s+(\d+)',t);return int(m.group(1)) if m else 0 if 'search completed successfully' in t.lower() else None

def cdd_confidence(hit,e):
    h=str(hit).lower()
    if 'specific' in h and 'non' not in h:return 'HIGH: CDD specific hit'
    if e<=1e-5:return 'HIGH: E≤1e-5'
    if e<=1e-3:return 'MODERATE: E≤1e-3'
    if e<=1e-2:return 'SUPPORTED: E≤0.01'
    return 'WEAK/REVIEW'

def parse_cdd(obj):
    rows=[]
    for line in read_text(obj).splitlines():
        if not line.strip() or line.startswith('#'):continue
        x=line.rstrip().split('\t')
        if len(x)<9 or x[0].lower()=='query':continue
        x += ['']*(11-len(x)); query=x[0].strip(); gene=query.split(' - ',1)[1].strip() if ' - ' in query else query
        try:s,e=int(float(x[3])),int(float(x[4]));ev=float(x[5]);bits=float(x[6])
        except ValueError:continue
        rows.append(dict(gene=norm_id(re.sub(r'\(Warning:.*$','',gene).strip()),query=query,hit_type=x[1].strip().lower(),pssm_id=x[2].strip(),start=min(s,e),end=max(s,e),evalue=ev,bitscore=bits,accession=x[7].strip(),domain=x[8].strip(),incomplete=x[9].strip(),superfamily=x[10].strip(),confidence=cdd_confidence(x[1],ev)))
    return pd.DataFrame(rows)

def run_cdd(proteins,evalue=0.01,timeout=240):
    if not proteins:raise ValueError('Protein FASTA required for CDD.')
    if len(proteins)>1000:raise ValueError('Split jobs larger than 1000 protein sequences.')
    data=dict(queries=fasta_text(proteins),db='cdd',smode='auto',useid1='true',compbasedadj='1',filter='false',evalue=str(evalue),maxhit='500',tdata='hits',dmode='full',cddefl='true')
    h={'User-Agent':'BioProtein-Studio/5.0 research software'}
    r=requests.post(CDD_URL,data=data,timeout=60,headers=h);r.raise_for_status();rid=_rid(r.text)
    if not rid:raise RuntimeError('NCBI CDD did not return a Search-ID.')
    t0=time.time()
    while time.time()-t0<timeout:
        r=requests.get(CDD_URL,params={'cdsid':rid},timeout=60,headers=h);r.raise_for_status();s=_status(r.text)
        if s==0:break
        if s in {1,2,4,5}:raise RuntimeError(f'CDD job status {s}')
        time.sleep(4)
    else:raise TimeoutError(f'CDD job {rid} still running.')
    r=requests.get(CDD_URL,params={'cdsid':rid,'tdata':'hits','dmode':'full','cddefl':'true'},timeout=60,headers=h);r.raise_for_status()
    return parse_cdd(r.text),rid,r.text

def collapse_domains(df):
    if df.empty:return df.copy()
    w=df.copy();p=w[~w.hit_type.str.contains('superfamily',case=False,na=False)]
    if not p.empty:w=p
    w['_r']=np.where(w.hit_type.str.contains('specific',case=False,na=False)&~w.hit_type.str.contains('non',case=False,na=False),0,1)
    w=w.sort_values(['gene','_r','evalue','bitscore'],ascending=[True,True,True,False]);keep=[]
    for _,sub in w.groupby('gene',sort=False):
        accepted=[]
        for _,r in sub.iterrows():
            redundant=False
            for k in accepted:
                ov=max(0,min(r.end,k.end)-max(r.start,k.start)+1);short=min(r.end-r.start+1,k.end-k.start+1)
                if short and ov/short>=0.8:redundant=True;break
            if not redundant:accepted.append(r)
        keep += accepted
    return pd.DataFrame(keep).drop(columns=['_r'],errors='ignore').reset_index(drop=True)

def domain_qc(df):
    if df.empty:return pd.DataFrame()
    o=df.copy();o['length_aa']=o.end-o.start+1
    o['scientific_status']=np.select([o.hit_type.str.contains('specific',case=False,na=False)&~o.hit_type.str.contains('non',case=False,na=False),o.evalue<=1e-5,o.evalue<=1e-3,o.evalue<=1e-2],['PASS: high-confidence specific hit','PASS: strong E-value','SUPPORTED','REVIEW: near default threshold'],default='REVIEW: weak')
    return o

def meme_ready():return shutil.which('meme') is not None

def parse_meme_xml(obj):
    root=ET.fromstring(read_text(obj));tag=lambda e:e.tag.split('}')[-1];names={}
    for e in root.iter():
        if tag(e)=='sequence' and e.attrib.get('id'):names[e.attrib['id']]=norm_id(e.attrib.get('name') or e.attrib['id'])
    meta={};summ=[];n=0
    for e in root.iter():
        if tag(e)!='motif':continue
        n+=1;mid=e.attrib.get('id') or f'motif_{n}';ev=e.attrib.get('e_value',e.attrib.get('evalue','nan'))
        try:ev=float(ev)
        except:ev=np.nan
        m=dict(motif=f'Motif {n}',motif_id=mid,consensus=e.attrib.get('name') or e.attrib.get('alt') or '',width=int(float(e.attrib.get('width',0) or 0)),sites=int(float(e.attrib.get('sites',0) or 0)),evalue=ev)
        meta[mid]=m;summ.append(m.copy())
    rows=[]
    for me in [e for e in root.iter() if tag(e)=='motif']:
        m=meta.get(me.attrib.get('id'),{});w=int(m.get('width',1) or 1)
        for e in me.iter():
            if tag(e)!='contributing_site':continue
            sid=e.attrib.get('sequence_id');pos=e.attrib.get('position')
            if sid is None or pos is None:continue
            p=float(e.attrib.get('pvalue','nan'));start=int(float(pos))+1
            rows.append(dict(gene=names.get(sid,norm_id(sid)),motif=m.get('motif','Motif'),motif_id=me.attrib.get('id'),consensus=m.get('consensus',''),start=start,end=start+w-1,width=w,site_pvalue=p,motif_evalue=m.get('evalue',np.nan)))
    return pd.DataFrame(rows),pd.DataFrame(summ)

def _zip_dir(path):
    b=io.BytesIO()
    with zipfile.ZipFile(b,'w',zipfile.ZIP_DEFLATED) as z:
        for p in path.rglob('*'):
            if p.is_file():z.write(p,p.relative_to(path))
    return b.getvalue()

def run_meme(proteins,nmotifs=10,minw=6,maxw=50,model='zoops',timeout=600):
    if not meme_ready():raise RuntimeError('MEME executable not found. Install MEME Suite with Bioconda in Linux/WSL.')
    if len(proteins)<3:raise ValueError('Use at least 3 related proteins for MEME; 5+ is recommended.')
    if model not in {'zoops','oops','anr'}:raise ValueError('Invalid MEME occurrence model.')
    with tempfile.TemporaryDirectory(prefix='bps_meme_') as td:
        td=Path(td);fa=td/'proteins.fa';out=td/'meme_out';fa.write_text(fasta_text(proteins))
        cmd=['meme',str(fa),'-protein','-oc',str(out),'-nmotifs',str(nmotifs),'-minw',str(minw),'-maxw',str(maxw),'-mod',model]
        p=subprocess.run(cmd,capture_output=True,text=True,timeout=timeout)
        if p.returncode:raise RuntimeError(p.stderr.strip() or 'MEME failed')
        xml=(out/'meme.xml').read_text(errors='replace');sites,summ=parse_meme_xml(xml)
        return dict(sites=sites,summary=summ,xml=xml,zip=_zip_dir(out),logos={x.name:x.read_bytes() for x in out.glob('logo*.png')},command=' '.join(cmd))

def motif_qc(sites,summary,nseq,domains=None):
    if summary.empty:return pd.DataFrame()
    rows=[]
    for _,m in summary.iterrows():
        sub=sites[sites.motif==m.motif] if not sites.empty else pd.DataFrame();ng=sub.gene.nunique() if not sub.empty else 0;prev=100*ng/max(1,nseq)
        overlap=np.nan
        if domains is not None and not domains.empty and not sub.empty:
            vals=[]
            for _,s in sub.iterrows():
                best=0.;L=max(1,s.end-s.start+1)
                for _,d in domains[domains.gene==s.gene].iterrows():best=max(best,100*max(0,min(s.end,d.end)-max(s.start,d.start)+1)/L)
                vals.append(best)
            if vals:overlap=float(np.mean(vals))
        ev=float(m.evalue) if not pd.isna(m.evalue) else np.nan
        status='PASS' if not pd.isna(ev) and ev<0.05 and prev>=50 else 'REVIEW'
        rows.append(dict(motif=m.motif,consensus=m.consensus,width=m.width,meme_evalue=ev,genes_with_motif=ng,prevalence_pct=round(prev,2),mean_domain_overlap_pct=round(overlap,2) if not pd.isna(overlap) else np.nan,status=status))
    return pd.DataFrame(rows)
