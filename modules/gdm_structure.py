from __future__ import annotations
import re, shutil, subprocess, tempfile, time
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple
import numpy as np, pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from .gdm_common import versionless, norm_id

def est2genome_ready()->bool:return shutil.which('est2genome') is not None

def _exact(gene,cds,genomic):
    c=cds.replace('-','').replace('U','T').upper(); g=genomic.replace('-','').replace('U','T').upper()
    i=g.find(c); strand='+'
    if i<0:
        i=g.find(str(Seq(c).reverse_complement()));strand='-'
    if i<0: raise RuntimeError('EMBOSS est2genome is not installed and CDS is not an exact contiguous genomic match.')
    df=pd.DataFrame([dict(gene=gene,feature='CDS',start=i+1,end=i+len(c),identity_pct=100.0,strand=strand,source='Exact contiguous match')])
    qc=dict(gene=gene,method='Exact contiguous match',exons=1,introns=0,cds_coverage_pct=100.0,weighted_identity_pct=100.0,
            canonical_splice_pct=np.nan,frame_multiple_of_3=len(c)%3==0,status='PASS')
    return df,qc,'Exact contiguous match; no intron inferred.'

def est2genome_pair(gene:str,cds:str,genomic:str,timeout=180):
    if not est2genome_ready():return _exact(gene,cds,genomic)
    with tempfile.TemporaryDirectory(prefix='bps_est2genome_') as td:
        td=Path(td); c=td/'cds.fa';g=td/'genomic.fa';o=td/'out.txt'
        c.write_text(f'>{gene}\n{cds}\n');g.write_text(f'>{gene}_genomic\n{genomic}\n')
        p=subprocess.run(['est2genome','-estsequence',str(c),'-genomesequence',str(g),'-outfile',str(o),'-align','N','-auto'],
                         capture_output=True,text=True,timeout=timeout)
        if p.returncode or not o.exists():raise RuntimeError(p.stderr.strip() or 'est2genome failed')
        raw=o.read_text(errors='replace')
    exons=[];introns=[];strand='+'
    for s in map(str.strip,raw.splitlines()):
        if 'REVERSED GENE' in s.upper():strand='-'
        if re.match(r'^Exon\s+',s):
            x=s.split()
            if len(x)>=8:
                try:score=float(x[1]);iden=float(x[2]);gs,ge=int(x[3]),int(x[4]);cs,ce=int(x[6]),int(x[7])
                except ValueError:continue
                exons.append(dict(gene=gene,feature='CDS',start=min(gs,ge),end=max(gs,ge),identity_pct=iden,
                                  score=score,cds_start=min(cs,ce),cds_end=max(cs,ce),strand=strand,source='EMBOSS est2genome'))
        elif re.match(r'^[+\-?]Intron\s+',s):
            x=s.split(); introns.append(x[0])
    df=pd.DataFrame(exons)
    if df.empty:raise RuntimeError('est2genome returned no exon segments')
    covered=set()
    for _,r in df.iterrows():covered.update(range(int(r.cds_start),int(r.cds_end)+1))
    coverage=100*len(covered)/max(1,len(cds.replace('-','')))
    identity=float(np.average(df.identity_pct,weights=df.cds_end-df.cds_start+1))
    canonical=100*np.mean([not x.startswith('?') for x in introns]) if introns else np.nan
    splice_ok=True if not introns else canonical>=80
    qc=dict(gene=gene,method='EMBOSS est2genome',exons=len(df),introns=len(introns),cds_coverage_pct=round(coverage,2),
            weighted_identity_pct=round(identity,2),canonical_splice_pct=round(canonical,2) if not np.isnan(canonical) else np.nan,
            frame_multiple_of_3=len(cds.replace('-',''))%3==0,status='PASS' if coverage>=95 and identity>=95 and splice_ok else 'REVIEW')
    return df,qc,raw

def gene_structure_batch(cds:Dict[str,str],genomic:Dict[str,str]):
    right={versionless(k):k for k in genomic};dfs=[];qcs=[];raw={}
    for ck,cseq in cds.items():
        gk=ck if ck in genomic else right.get(versionless(ck))
        if gk is None:continue
        df,q,t=est2genome_pair(ck,cseq,genomic[gk]);dfs.append(df);qcs.append(q);raw[ck]=t
    if not dfs:raise ValueError('No matching CDS/genomic FASTA IDs.')
    return pd.concat(dfs,ignore_index=True),pd.DataFrame(qcs),raw

def _pid(a,b):
    if not a or not b:return 0.0
    n=min(len(a),len(b));return 100*sum(x==y for x,y in zip(a[:n],b[:n]))/max(len(a),len(b))

def ncbi_structure(accession:str,email:str,api_key:Optional[str]=None,max_links=20):
    if '@' not in email:raise ValueError('Valid email required for NCBI Entrez.')
    Entrez.email=email;Entrez.tool='BioProtein-Studio';Entrez.api_key=api_key or None
    acc=norm_id(accession); ids=Entrez.read(Entrez.esearch(db='protein',term=f'{acc}[Accession]',retmax=3)).get('IdList',[])
    if not ids:return None,None,None
    uid=ids[0]
    with Entrez.efetch(db='protein',id=uid,rettype='fasta',retmode='text') as h:prot=str(SeqIO.read(h,'fasta').seq)
    links=Entrez.read(Entrez.elink(dbfrom='protein',db='nuccore',id=uid));nids=[]
    for b in links:
        for ls in b.get('LinkSetDb',[]):nids += [x.get('Id') for x in ls.get('Link',[]) if x.get('Id')]
    for nid in nids[:max_links]:
        try:
            with Entrez.efetch(db='nuccore',id=nid,rettype='gb',retmode='text') as h: rec=SeqIO.read(h,'genbank')
        except Exception:continue
        for f in rec.features:
            if f.type!='CDS':continue
            pids=f.qualifiers.get('protein_id',[]);trans=f.qualifiers.get('translation',[])
            if not (any(versionless(x)==versionless(acc) for x in pids) or any(_pid(prot,t)>=99 for t in trans)):continue
            parts=list(getattr(f.location,'parts',[f.location]));strand=f.location.strand or 1
            gs=min(int(p.start) for p in parts);ge=max(int(p.end) for p in parts);rows=[]
            for p in parts:
                ps,pe=int(p.start),int(p.end)
                a,b=(ge-pe+1,ge-ps) if strand==-1 else (ps-gs+1,pe-gs)
                rows.append(dict(gene=acc,feature='CDS',start=a,end=b,identity_pct=100.0,strand='-' if strand==-1 else '+',source=f'NCBI {rec.id}'))
            cds=str(f.extract(rec.seq)).upper();tr=trans[0] if trans else str(Seq(cds).translate(to_stop=True));ident=_pid(prot.rstrip('*'),tr.rstrip('*'))
            qc=dict(gene=acc,method='NCBI reference CDS feature',exons=len(rows),introns=max(0,len(rows)-1),cds_coverage_pct=100.0,
                    weighted_identity_pct=100.0,translation_identity_pct=round(ident,2),canonical_splice_pct=np.nan,
                    frame_multiple_of_3=len(cds)%3==0,reference_record=rec.id,status='PASS' if ident>=95 else 'REVIEW')
            bundle=dict(protein=prot,cds=cds,genomic=str(rec.seq[gs:ge]).upper(),reference_record=rec.id)
            return pd.DataFrame(rows).sort_values('start'),qc,bundle
    return None,None,None

def ncbi_structures(accessions:Sequence[str],email:str,api_key:Optional[str]=None):
    dfs=[];q=[];bundles={};missing=[]
    for a in accessions:
        try:df,qc,b=ncbi_structure(a,email,api_key)
        except Exception:df=qc=b=None
        if df is None:missing.append(a)
        else:dfs.append(df);q.append(qc);bundles[a]=b
        time.sleep(0.12)
    return (pd.concat(dfs,ignore_index=True) if dfs else pd.DataFrame(),pd.DataFrame(q),bundles,missing)


def parse_gene_structure_annotation(obj, valid_genes=None):
    """Parse a simple gene-structure CSV/TSV or GFF3/GTF annotation.

    Accepted tabular columns: gene, start, end, optional feature/strand.
    GFF3/GTF exon/CDS rows are accepted when an exact gene/Parent/ID matches.
    Imported annotations are not sequence-realigned, so QC is marked REFERENCE.
    """
    from .gdm_common import read_text
    import io
    text=read_text(obj).lstrip('\ufeff')
    lines=[x for x in text.splitlines() if x.strip() and not x.lstrip().startswith('#')]
    if not lines:
        raise ValueError('Gene-structure annotation is empty.')
    rows=[]
    # GFF3/GTF route
    if any(len(x.split('\t'))>=9 for x in lines[:20]):
        for line in lines:
            x=line.split('\t')
            if len(x)<9: continue
            seqid,source,ftype,start,end,score,strand,phase,attrs=x[:9]
            if ftype.lower() not in {'exon','cds'}: continue
            try:a,b=int(start),int(end)
            except ValueError: continue
            amap={}
            for token in re.split(r';\s*',attrs.strip()):
                if '=' in token:
                    k,v=token.split('=',1);amap[k.strip()]=v.strip().strip('"')
                elif ' ' in token:
                    k,v=token.split(' ',1);amap[k.strip()]=v.strip().strip('"')
            candidates=[]
            for k in ['gene','gene_id','Name','Parent','ID','transcript_id']:
                if amap.get(k): candidates += [z for z in re.split(r',',amap[k]) if z]
            candidates.append(seqid)
            gene=None
            if valid_genes:
                vg={norm_id(g):g for g in valid_genes}; vv={versionless(g):g for g in valid_genes}
                for c in candidates:
                    c=norm_id(c)
                    if c in vg: gene=vg[c];break
                    if versionless(c) in vv: gene=vv[versionless(c)];break
            else:
                gene=norm_id(candidates[0]) if candidates else norm_id(seqid)
            if gene:
                rows.append(dict(gene=gene,feature='CDS' if ftype.lower()=='cds' else 'exon',start=min(a,b),end=max(a,b),strand=strand,source='Imported GFF/GTF annotation'))
    else:
        try:
            df=pd.read_csv(io.StringIO(text),sep=None,engine='python')
        except Exception as e:
            raise ValueError(f'Could not parse gene-structure table: {e}')
        cols={str(c).strip().lower():c for c in df.columns}
        gc=next((cols[k] for k in ['gene','gene_id','id','sequence','seqid'] if k in cols),None)
        sc=next((cols[k] for k in ['start','exon_start','cds_start'] if k in cols),None)
        ec=next((cols[k] for k in ['end','exon_end','cds_end'] if k in cols),None)
        fc=next((cols[k] for k in ['feature','type'] if k in cols),None)
        stc=next((cols[k] for k in ['strand'] if k in cols),None)
        if gc is None or sc is None or ec is None:
            raise ValueError('Structure table needs gene, start and end columns.')
        for _,r in df.iterrows():
            try:a,b=int(float(r[sc])),int(float(r[ec]))
            except Exception: continue
            gene=norm_id(r[gc])
            if valid_genes:
                vg={norm_id(g):g for g in valid_genes}; vv={versionless(g):g for g in valid_genes}
                gene=vg.get(gene,vv.get(versionless(gene)))
            if not gene: continue
            feature=str(r[fc]) if fc is not None and not pd.isna(r[fc]) else 'CDS'
            strand=str(r[stc]) if stc is not None and not pd.isna(r[stc]) else '.'
            rows.append(dict(gene=gene,feature=feature,start=min(a,b),end=max(a,b),strand=strand,source='Imported structure table'))
    out=pd.DataFrame(rows)
    if out.empty:
        raise ValueError('No exon/CDS rows matched the supplied gene IDs.')
    out=out.sort_values(['gene','start','end']).reset_index(drop=True)
    q=[]
    for g,sub in out.groupby('gene',sort=False):
        n=len(sub); q.append(dict(gene=g,method='Imported reference annotation',exons=n,introns=max(0,n-1),cds_coverage_pct=np.nan,weighted_identity_pct=np.nan,canonical_splice_pct=np.nan,frame_multiple_of_3=np.nan,status='REFERENCE'))
    return out,pd.DataFrame(q)
