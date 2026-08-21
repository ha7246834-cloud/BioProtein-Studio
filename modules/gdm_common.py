from __future__ import annotations
import io, re, zipfile
from typing import Dict, Iterable, Sequence, Tuple
import pandas as pd
from Bio import Phylo, SeqIO
from Bio.Seq import Seq

AA_ALLOWED=set('ACDEFGHIKLMNPQRSTVWYBXZJUO*')
DNA_ALLOWED=set('ACGTUNRYKMSWBDHV-')

def norm_id(x:str)->str:
    x=str(x).strip().lstrip('>')
    return x.split()[0] if x else ''

def versionless(x:str)->str:
    return re.sub(r'\.\d+$','',norm_id(x))

def read_text(obj)->str:
    if obj is None:return ''
    if isinstance(obj,str):return obj
    if isinstance(obj,bytes):return obj.decode('utf-8','replace')
    if hasattr(obj,'getvalue'):
        v=obj.getvalue(); return v.decode('utf-8','replace') if isinstance(v,bytes) else str(v)
    return str(obj)

def parse_fasta(obj, alphabet='auto')->Dict[str,str]:
    text=read_text(obj).strip()
    if not text: return {}
    out={}
    for r in SeqIO.parse(io.StringIO(text),'fasta'):
        rid=norm_id(r.id); seq=str(r.seq).replace(' ','').upper()
        if not rid: continue
        if rid in out: raise ValueError(f'Duplicate FASTA identifier: {rid}')
        allowed=AA_ALLOWED if alphabet=='protein' else DNA_ALLOWED if alphabet=='dna' else None
        if allowed is not None:
            bad=sorted(set(seq)-allowed)
            if bad: raise ValueError(f"{rid}: invalid characters: {','.join(bad)}")
        out[rid]=seq
    if not out: raise ValueError('No FASTA records detected.')
    return out

def fasta_text(records:Dict[str,str], width=70)->str:
    lines=[]
    for rid,seq in records.items():
        lines.append(f'>{rid}')
        lines.extend(seq[i:i+width] for i in range(0,len(seq),width))
    return '\n'.join(lines)+'\n'

def protein_qc(records:Dict[str,str])->pd.DataFrame:
    rows=[]; standard=set('ACDEFGHIKLMNPQRSTVWY'); amb=set('BXZJUO')
    for rid,seq in records.items():
        seq=seq.replace('-','').upper(); terminal=seq.endswith('*'); body=seq[:-1] if terminal else seq
        internal=body.count('*'); namb=sum(x in amb for x in body); pct=100*namb/max(1,len(body))
        status='REVIEW' if internal or pct>5 else 'PASS WITH AMBIGUOUS RESIDUES' if namb else 'PASS'
        rows.append(dict(gene=rid,length_aa=len(body),terminal_stop=terminal,internal_stops=internal,
                         ambiguous_residues=namb,ambiguous_pct=round(pct,3),status=status))
    return pd.DataFrame(rows)

def translate_cds(records:Dict[str,str])->Tuple[Dict[str,str],pd.DataFrame]:
    proteins={}; rows=[]
    for rid,seq in records.items():
        s=seq.replace('-','').replace('U','T'); frame=len(s)%3==0; trim=s[:len(s)-len(s)%3]
        aa=str(Seq(trim).translate(to_stop=False)) if trim else ''
        internal=aa[:-1].count('*') if aa.endswith('*') else aa.count('*')
        if aa.endswith('*'): aa=aa[:-1]
        proteins[rid]=aa
        rows.append(dict(gene=rid,cds_length_nt=len(s),frame_multiple_of_3=frame,start_ATG=s.startswith('ATG'),
                         terminal_stop=s[-3:] in {'TAA','TAG','TGA'} if len(s)>=3 else False,
                         internal_stop_count=internal,translated_length_aa=len(aa),
                         status='PASS' if frame and internal==0 else 'REVIEW'))
    return proteins,pd.DataFrame(rows)

def looks_ncbi_accession(x:str)->bool:
    return bool(re.fullmatch(r'(?:[A-Z]{1,4}_?\d{5,}|[A-Z]{3}\d{5,})(?:\.\d+)?',norm_id(x)))

def newick_order(obj)->list[str]:
    t=read_text(obj).strip()
    if not t:return []
    tree=Phylo.read(io.StringIO(t),'newick')
    return [norm_id(x.name) for x in tree.get_terminals() if x.name]

def choose_order(genes:Iterable[str], preferred:Sequence[str]|None=None)->list[str]:
    genes=list(dict.fromkeys(str(x) for x in genes if str(x)))
    if not preferred:return sorted(genes)
    out=[]; used=set()
    for p in preferred:
        for g in genes:
            if g not in used and (g==p or versionless(g)==versionless(p)):
                out.append(g);used.add(g);break
    out.extend(g for g in genes if g not in used)
    return out

def zip_files(files:dict[str,bytes|str])->bytes:
    b=io.BytesIO()
    with zipfile.ZipFile(b,'w',zipfile.ZIP_DEFLATED) as z:
        for n,c in files.items(): z.writestr(n,c.encode() if isinstance(c,str) else c)
    return b.getvalue()
