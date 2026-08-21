import io
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd

def fig_bytes(fig,fmt='png',dpi=600):
    b=io.BytesIO();fig.savefig(b,format=fmt,dpi=dpi,bbox_inches='tight');return b.getvalue()

def _colors():return plt.rcParams['axes.prop_cycle'].by_key().get('color',[None])

def gene_structure(df:pd.DataFrame,order):
    if df.empty:return None
    fig,ax=plt.subplots(figsize=(12,max(3.5,.55*len(order)+1.8)))
    for y,g in enumerate(order):
        s=df[df.gene==g].sort_values('start')
        if s.empty:continue
        ax.plot([s.start.min(),s.end.max()],[y,y],lw=1)
        for _,r in s.iterrows():ax.add_patch(Rectangle((r.start,y-.2),max(1,r.end-r.start+1),.4,edgecolor='black',lw=.7))
    ax.set_yticks(range(len(order)));ax.set_yticklabels(order);ax.invert_yaxis();ax.set_xlabel('Genomic position (bp)');ax.set_title('Gene Structure')
    ax.spines[['top','right','left']].set_visible(False);fig.tight_layout();return fig

def architecture(df,order,col,title):
    if df.empty:return None
    pal=_colors();labs=list(dict.fromkeys(df[col].astype(str)));style={x:pal[i%len(pal)] for i,x in enumerate(labs)};mx=max(1,int(df.end.max()))
    fig,ax=plt.subplots(figsize=(12,max(3.5,.55*len(order)+1.8)))
    for y,g in enumerate(order):
        ax.plot([1,mx],[y,y],lw=.8)
        for _,r in df[df.gene==g].iterrows():ax.add_patch(Rectangle((r.start,y-.2),max(1,r.end-r.start+1),.4,edgecolor='black',lw=.7,facecolor=style[str(r[col])],label=str(r[col])))
    h,l=ax.get_legend_handles_labels();u={}
    for a,b in zip(h,l):u.setdefault(b,a)
    if u:ax.legend(u.values(),u.keys(),bbox_to_anchor=(1.01,1),loc='upper left',frameon=False)
    ax.set_yticks(range(len(order)));ax.set_yticklabels(order);ax.invert_yaxis();ax.set_xlim(0,mx*1.03);ax.set_xlabel('Protein position (aa)');ax.set_title(title)
    ax.spines[['top','right','left']].set_visible(False);fig.tight_layout();return fig

def combined(structures,domains,motifs,order):
    panels=[]
    if not motifs.empty:panels.append(('motif',motifs,'motif'))
    if not domains.empty:panels.append(('domain',domains,'domain'))
    if not structures.empty:panels.append(('structure',structures,None))
    if not panels:return None
    fig,axes=plt.subplots(1,len(panels),figsize=(6*len(panels),max(4,.52*len(order)+2)),squeeze=False);axes=axes[0];pal=_colors()
    for i,(kind,df,col) in enumerate(panels):
        ax=axes[i];mx=max(1,int(df.end.max()))
        if kind=='structure':
            for y,g in enumerate(order):
                s=df[df.gene==g].sort_values('start')
                if s.empty:continue
                ax.plot([s.start.min(),s.end.max()],[y,y],lw=.8)
                for _,r in s.iterrows():ax.add_patch(Rectangle((r.start,y-.19),max(1,r.end-r.start+1),.38,edgecolor='black',lw=.6))
            ax.set_xlabel('Gene position (bp)');ax.set_title('Gene Structure')
        else:
            labs=list(dict.fromkeys(df[col].astype(str)));style={x:pal[j%len(pal)] for j,x in enumerate(labs)}
            for y,g in enumerate(order):
                ax.plot([1,mx],[y,y],lw=.8)
                for _,r in df[df.gene==g].iterrows():ax.add_patch(Rectangle((r.start,y-.19),max(1,r.end-r.start+1),.38,edgecolor='black',lw=.6,facecolor=style[str(r[col])]))
            ax.set_xlabel('Protein position (aa)');ax.set_title('Conserved Motifs' if kind=='motif' else 'Conserved Domains')
        ax.set_xlim(0,mx*1.03);ax.set_yticks(range(len(order)));ax.set_yticklabels(order if i==0 else []);ax.invert_yaxis();ax.spines[['top','right','left']].set_visible(False)
    fig.suptitle('Integrated Gene-Family Architecture',y=1.01);fig.tight_layout();return fig
