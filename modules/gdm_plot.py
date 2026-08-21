import io
from io import StringIO
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import pandas as pd
from Bio import Phylo
from .gdm_style import style_from_preset, assign_colors


def fig_bytes(fig, fmt='png', dpi=600):
    b = io.BytesIO()
    fig.savefig(b, format=fmt, dpi=dpi, bbox_inches='tight')
    return b.getvalue()


def _S(style=None):
    s = style_from_preset('Journal Classic')
    if style:
        s.update(style)
    return s


def _ordered_subset(df, order):
    if df is None or df.empty:
        return pd.DataFrame()
    d = df.copy()
    d['__order__'] = d['gene'].map({g: i for i, g in enumerate(order)})
    return d[d['__order__'].notna()].sort_values(['__order__', 'start', 'end']).drop(columns='__order__')


def _structure_style(row, style):
    status = str(row.get('validation_status', row.get('status', 'PASS'))).upper()
    source = str(row.get('source', ''))
    if status.startswith('REVIEW') or 'rescue' in source.lower():
        return dict(facecolor=style['review_exon_color'], hatch='///')
    return dict(facecolor=style['exon_color'], hatch=None)


def gene_structure(df: pd.DataFrame, order, style=None):
    if df is None or df.empty:
        return None
    s = _S(style)
    fig, ax = plt.subplots(figsize=(12, max(3.8, .62 * len(order) + 1.8)))
    d = _ordered_subset(df, order)
    xmin, xmax = d['start'].min(), d['end'].max()
    for y, g in enumerate(order):
        sub = d[d.gene == g]
        if sub.empty:
            continue
        ax.plot([sub.start.min(), sub.end.max()], [y, y], lw=s['line_width'], color=s['intron_color'], zorder=1)
        for _, r in sub.iterrows():
            sty = _structure_style(r, s)
            ax.add_patch(Rectangle((r.start, y - .20), max(1, r.end-r.start+1), .40,
                                   edgecolor=s['edge_color'], lw=.7, facecolor=sty['facecolor'], hatch=sty['hatch'], zorder=2))
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=s['label_size']); ax.invert_yaxis()
    pad = max(5, (xmax-xmin)*0.02 if xmax > xmin else 5)
    ax.set_xlim(max(0, xmin-pad), xmax+pad)
    ax.set_xlabel('Genomic position (bp)'); ax.set_title('Gene Structure (Exon–Intron)', fontsize=s['title_size'])
    ax.spines[['top','right','left']].set_visible(False)
    leg = [
        Rectangle((0,0),1,1,facecolor=s['exon_color'],edgecolor=s['edge_color'],label='Validated exon'),
        Rectangle((0,0),1,1,facecolor=s['review_exon_color'],edgecolor=s['edge_color'],hatch='///',label='REVIEW exon'),
        Line2D([0],[0],color=s['intron_color'],lw=s['line_width'],label='Intron')
    ]
    ax.legend(handles=leg, loc='upper left', bbox_to_anchor=(0,-0.08), frameon=False, ncol=3)
    fig.tight_layout(rect=[0,0.08,1,1])
    return fig


def missing_structure_figure(order, style=None):
    s = _S(style)
    fig, ax = plt.subplots(figsize=(12,max(3.8,.62*len(order)+1.8)))
    ax.text(.5,.55,'Gene structure unavailable',ha='center',va='center',fontsize=s['title_size'],weight='bold',transform=ax.transAxes)
    ax.text(.5,.45,'Provide matching CDS + genomic FASTA, NCBI-resolvable accessions, or an exon/CDS annotation table.',ha='center',va='center',fontsize=s['label_size'],transform=ax.transAxes,wrap=True)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order,fontsize=s['label_size']); ax.invert_yaxis(); ax.set_xticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title('Gene Structure',fontsize=s['title_size']); fig.tight_layout(); return fig


def _feature_color_map(df, col, s):
    labs = list(dict.fromkeys(df[col].astype(str)))
    custom = s.get('motif_colors',{}) if col == 'motif' else s.get('domain_colors',{})
    palette = s['motif_palette'] if col == 'motif' else s['domain_palette']
    m = assign_colors(labs, palette)
    m.update({str(k):v for k,v in custom.items()})
    return m


def architecture(df, order, col, title, xlabel='Protein position (aa)', style=None):
    if df is None or df.empty:
        return None
    s = _S(style); d = _ordered_subset(df, order); colors = _feature_color_map(d,col,s); mx=max(1,int(d.end.max()))
    fig, ax = plt.subplots(figsize=(12,max(3.8,.62*len(order)+1.8)))
    for y,g in enumerate(order):
        ax.plot([1,mx],[y,y],lw=.8,color=s['backbone_color'],zorder=1)
        for _,r in d[d.gene==g].iterrows():
            lab=str(r[col])
            ax.add_patch(Rectangle((r.start,y-.20),max(1,r.end-r.start+1),.40,edgecolor=s['edge_color'],lw=.6,facecolor=colors[lab],zorder=2,label=lab))
    h,l=ax.get_legend_handles_labels(); u={}
    for a,b in zip(h,l): u.setdefault(b,a)
    if u:
        ax.legend(u.values(),u.keys(),bbox_to_anchor=(1.01,1),loc='upper left',frameon=False,ncol=1 if len(u)<=8 else 2,title='Legend')
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order,fontsize=s['label_size']); ax.invert_yaxis(); ax.set_xlim(0,mx*1.03)
    ax.set_xlabel(xlabel); ax.set_title(title,fontsize=s['title_size']); ax.spines[['top','right','left']].set_visible(False); fig.tight_layout(); return fig


def _tree_layout(tree, order):
    depths=tree.depths()
    if not max(depths.values(),default=0): depths=tree.depths(unit_branch_lengths=True)
    ymap={name:i for i,name in enumerate(order)}; ypos={}
    def assign(clade):
        if clade.is_terminal():
            y=ymap.get(clade.name)
            if y is not None:ypos[clade]=y
            return y
        ys=[assign(c) for c in clade.clades]; ys=[v for v in ys if v is not None]
        if not ys:return None
        y=sum(ys)/len(ys); ypos[clade]=y; return y
    assign(tree.root); return depths,ypos


def _support_label(clade):
    conf=getattr(clade,'confidence',None)
    if conf is not None:
        try:
            v=float(conf); return f'{v:.2f}' if v<=1 else f'{v:.0f}'
        except Exception: pass
    name=str(getattr(clade,'name','') or '').strip()
    if name and any(ch.isdigit() for ch in name) and len(name)<=15:
        return name
    return ''


def _draw_tree_axis(ax, tree_text, order, style=None):
    s=_S(style)
    try: tree=Phylo.read(StringIO(tree_text),'newick')
    except Exception:
        ax.text(.5,.5,'Tree could not be parsed',ha='center',va='center',transform=ax.transAxes); ax.axis('off'); return False
    depths,ypos=_tree_layout(tree,order)
    if not ypos:
        ax.text(.5,.5,'Tree contains no matching terminals',ha='center',va='center',transform=ax.transAxes); ax.axis('off'); return False
    xmax=max(depths.values()) if depths else 1.0
    for clade in tree.find_clades(order='preorder'):
        if clade not in ypos: continue
        x=depths.get(clade,0.0); y=ypos[clade]
        if clade.clades:
            ys=[ypos[c] for c in clade.clades if c in ypos]; xs=[depths.get(c,x) for c in clade.clades if c in ypos]
            if ys:
                ax.plot([x,x],[min(ys),max(ys)],color=s['tree_color'],lw=s['line_width'])
                for cx,cy in zip(xs,ys): ax.plot([x,cx],[cy,cy],color=s['tree_color'],lw=s['line_width'])
        else:
            ax.text(xmax*1.03 if xmax else 1,y,clade.name,va='center',ha='left',fontsize=s['label_size'])
        if not clade.is_terminal():
            label=_support_label(clade)
            if label: ax.text(x+xmax*.012,y-.10,label,fontsize=s['support_size'],color=s['support_color'],ha='left',va='bottom')
    ax.set_xlim(0,xmax*1.38 if xmax else 1); ax.set_ylim(-.5,len(order)-.5); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title('Phylogenetic Tree',fontsize=s['title_size']); [sp.set_visible(False) for sp in ax.spines.values()]; return True


def phylogeny_figure(tree_text, order, style=None):
    if not tree_text:return None
    fig,ax=plt.subplots(figsize=(8.5,max(4,.62*len(order)+1.6)))
    if not _draw_tree_axis(ax,tree_text,order,style): plt.close(fig); return None
    fig.tight_layout(); return fig


def _draw_structure_panel(ax, df, order, show_labels=False, style=None):
    s=_S(style); d=_ordered_subset(df,order); xmin,xmax=d.start.min(),d.end.max(); pad=max(5,(xmax-xmin)*.02 if xmax>xmin else 5)
    for y,g in enumerate(order):
        sub=d[d.gene==g]
        if sub.empty:continue
        ax.plot([sub.start.min(),sub.end.max()],[y,y],lw=s['line_width'],color=s['intron_color'],zorder=1)
        for _,r in sub.iterrows():
            sty=_structure_style(r,s)
            ax.add_patch(Rectangle((r.start,y-.18),max(1,r.end-r.start+1),s['feature_height'],edgecolor=s['edge_color'],lw=.55,facecolor=sty['facecolor'],hatch=sty['hatch'],zorder=2))
    ax.set_xlim(max(0,xmin-pad),xmax+pad); ax.set_yticks(range(len(order))); ax.set_yticklabels(order if show_labels else [],fontsize=s['label_size']); ax.invert_yaxis(); ax.set_xlabel('bp'); ax.set_title('Gene Structure',fontsize=s['title_size']); ax.spines[['top','right','left']].set_visible(False)


def _draw_feature_panel(ax, df, order, col, title, show_labels=False, style=None):
    s=_S(style); d=_ordered_subset(df,order); mx=max(1,int(d.end.max())); colors=_feature_color_map(d,col,s)
    for y,g in enumerate(order):
        ax.plot([1,mx],[y,y],lw=.8,color=s['backbone_color'],zorder=1)
        for _,r in d[d.gene==g].iterrows():
            lab=str(r[col]); ax.add_patch(Rectangle((r.start,y-.18),max(1,r.end-r.start+1),s['feature_height'],edgecolor=s['edge_color'],lw=.55,facecolor=colors[lab],zorder=2,label=lab))
    ax.set_xlim(0,mx*1.03); ax.set_yticks(range(len(order))); ax.set_yticklabels(order if show_labels else [],fontsize=s['label_size']); ax.invert_yaxis(); ax.set_xlabel('aa'); ax.set_title(title,fontsize=s['title_size']); ax.spines[['top','right','left']].set_visible(False)
    h,l=ax.get_legend_handles_labels();u={}
    for a,b in zip(h,l):u.setdefault(b,a)
    return u


def combined(tree_text, structures, domains, motifs, order, style=None):
    s=_S(style); panels=[]
    if tree_text:panels.append('tree')
    panels.append('structure')
    if motifs is not None and not motifs.empty:panels.append('motif')
    if domains is not None and not domains.empty:panels.append('domain')
    widths=[{'tree':1.8,'structure':2.2,'motif':2.0,'domain':2.0}[p] for p in panels]
    fig,axes=plt.subplots(1,len(panels),figsize=(5.2*len(panels),max(4.4,.70*len(order)+1.8)),gridspec_kw={'width_ratios':widths},squeeze=False);axes=axes[0]
    legends={};first_non_tree=True
    for i,p in enumerate(panels):
        ax=axes[i]
        if p=='tree':_draw_tree_axis(ax,tree_text,order,s);continue
        show_labels=first_non_tree and 'tree' not in panels
        if p=='structure':
            if structures is None or structures.empty:
                ax.text(.5,.54,'Gene structure unavailable',ha='center',va='center',fontsize=s['title_size'],weight='bold',transform=ax.transAxes)
                ax.text(.5,.44,'CDS + genomic FASTA, reference mapping, or annotation required',ha='center',va='center',fontsize=s['label_size'],transform=ax.transAxes,wrap=True)
                ax.set_yticks(range(len(order)));ax.set_yticklabels(order if show_labels else [],fontsize=s['label_size']);ax.set_ylim(-.5,len(order)-.5);ax.invert_yaxis();ax.set_xticks([]);ax.set_title('Gene Structure',fontsize=s['title_size']);[sp.set_visible(False) for sp in ax.spines.values()]
            else:_draw_structure_panel(ax,structures,order,show_labels,s)
        elif p=='motif':legends['motif']=_draw_feature_panel(ax,motifs,order,'motif','Conserved Motifs',show_labels,s)
        elif p=='domain':legends['domain']=_draw_feature_panel(ax,domains,order,'domain','Conserved Domains',show_labels,s)
        first_non_tree=False
    if 'structure' in panels:
        exon_leg=[Rectangle((0,0),1,1,facecolor=s['exon_color'],edgecolor=s['edge_color'],label='Validated exon'),Rectangle((0,0),1,1,facecolor=s['review_exon_color'],edgecolor=s['edge_color'],hatch='///',label='REVIEW exon'),Line2D([0],[0],color=s['intron_color'],lw=s['line_width'],label='Intron')]
        fig.legend(handles=exon_leg,loc='lower left',bbox_to_anchor=(.04,-.01),frameon=False,ncol=3,title='Gene structure')
    if legends.get('motif'):fig.legend(legends['motif'].values(),legends['motif'].keys(),loc='lower center',bbox_to_anchor=(.56,-.01),frameon=False,ncol=min(5,max(1,len(legends['motif']))),title='Motifs')
    if legends.get('domain'):fig.legend(legends['domain'].values(),legends['domain'].keys(),loc='lower right',bbox_to_anchor=(.98,-.01),frameon=False,ncol=min(3,max(1,len(legends['domain']))),title='Domains')
    fig.suptitle('Integrated Gene-Family Architecture',y=1.02,fontsize=s['title_size']+1);fig.tight_layout(rect=[0,.08,1,.98]);return fig
