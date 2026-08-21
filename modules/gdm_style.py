STYLE_PRESETS = {
    'Journal Classic': {
        'tree_color': '#222222', 'support_color': '#666666', 'exon_color': '#F2C94C', 'review_exon_color': '#F6B26B',
        'intron_color': '#666666', 'backbone_color': '#D9D9D9', 'edge_color': '#222222',
        'motif_palette': ['#4C78A8','#F58518','#54A24B','#E45756','#B279A2','#FF9DA6','#9D755D','#BAB0AC','#72B7B2','#EECA3B'],
        'domain_palette': ['#2F80ED','#27AE60','#9B51E0','#EB5757','#F2994A','#56CCF2'],
    },
    'Colorblind Safe': {
        'tree_color': '#000000', 'support_color': '#666666', 'exon_color': '#E69F00', 'review_exon_color': '#F0E442',
        'intron_color': '#666666', 'backbone_color': '#CFCFCF', 'edge_color': '#111111',
        'motif_palette': ['#0072B2','#E69F00','#009E73','#D55E00','#CC79A7','#56B4E9','#F0E442','#000000','#999999','#332288'],
        'domain_palette': ['#0072B2','#009E73','#D55E00','#CC79A7','#E69F00','#56B4E9'],
    },
    'Plant & Earth': {
        'tree_color': '#263238', 'support_color': '#607D8B', 'exon_color': '#8BC34A', 'review_exon_color': '#FFB74D',
        'intron_color': '#795548', 'backbone_color': '#D7CCC8', 'edge_color': '#37474F',
        'motif_palette': ['#2E7D32','#66BB6A','#F9A825','#8D6E63','#00897B','#7E57C2','#C62828','#039BE5','#AFB42B','#5D4037'],
        'domain_palette': ['#388E3C','#1565C0','#6A1B9A','#EF6C00','#00838F','#AD1457'],
    },
    'Soft Pastel': {
        'tree_color': '#424242', 'support_color': '#757575', 'exon_color': '#FFD166', 'review_exon_color': '#FFAB91',
        'intron_color': '#8D8D8D', 'backbone_color': '#E5E5E5', 'edge_color': '#4A4A4A',
        'motif_palette': ['#8ECAE6','#FFB703','#90BE6D','#F28482','#B8A1D9','#F6BD60','#84A59D','#F5CAC3','#A8DADC','#CDB4DB'],
        'domain_palette': ['#6FA8DC','#93C47D','#B4A7D6','#E6B8AF','#FFD966','#76A5AF'],
    },
    'High Contrast': {
        'tree_color': '#000000', 'support_color': '#444444', 'exon_color': '#FFD400', 'review_exon_color': '#FF6B00',
        'intron_color': '#000000', 'backbone_color': '#BDBDBD', 'edge_color': '#000000',
        'motif_palette': ['#0057B8','#FF7A00','#00843D','#D22630','#6F2DA8','#00A6D2','#F2C500','#8A1538','#00A499','#7A7A7A'],
        'domain_palette': ['#0057B8','#00843D','#6F2DA8','#D22630','#FF7A00','#00A6D2'],
    },
}


def style_from_preset(name='Journal Classic'):
    base = STYLE_PRESETS.get(name, STYLE_PRESETS['Journal Classic']).copy()
    base.update({
        'name': name,
        'title_size': 13,
        'label_size': 9,
        'support_size': 7,
        'line_width': 1.0,
        'feature_height': 0.36,
        'motif_colors': {},
        'domain_colors': {},
    })
    return base


def assign_colors(labels, palette):
    labels = list(dict.fromkeys(str(x) for x in labels))
    return {lab: palette[i % len(palette)] for i, lab in enumerate(labels)}
