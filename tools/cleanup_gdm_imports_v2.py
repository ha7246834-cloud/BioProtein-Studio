from pathlib import Path

p = Path('pages/2_Gene_Structure_Domain_Motif.py')
s = p.read_text()
lines = s.splitlines()
out = []
for line in lines:
    if line.startswith('from modules.gdm_cdd_meme import '):
        line = ('from modules.gdm_cdd_meme import run_cdd, collapse_domains, domain_qc, meme_ready, '
                'run_meme, motif_qc, parse_cdd, parse_meme_xml, CDD_BATCH_SIZE, CDD_MAX_SEQUENCES, '
                'CLOUD_MEME_MAX_SEQUENCES, CLOUD_MEME_MAX_RESIDUES')
    elif line.startswith('from modules.gdm_phylogeny import '):
        line = ('from modules.gdm_phylogeny import build_phylogeny, external_phylogeny_ready, '
                'publication_phylogeny_ready, phylogeny_tool_status, AUTO_IQTREE_MAX_SEQUENCES, '
                'AUTO_IQTREE_MAX_RESIDUES, CLOUD_AUTO_IQTREE_MAX_SEQUENCES, '
                'CLOUD_AUTO_IQTREE_MAX_RESIDUES, CLOUD_PUBLICATION_MAX_SEQUENCES, '
                'CLOUD_PUBLICATION_MAX_RESIDUES')
    out.append(line)
p.write_text('\n'.join(out) + '\n')
print('GDM_IMPORT_CLEANUP_COMPLETE')
