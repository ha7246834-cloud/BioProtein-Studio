from pathlib import Path

# 1) Engine policy: keep workstation thresholds, add much safer shared-cloud thresholds.
p = Path('modules/gdm_phylogeny.py')
s = p.read_text()
old = """AUTO_IQTREE_MAX_SEQUENCES = 60
AUTO_IQTREE_MAX_RESIDUES = 40000
CLOUD_PUBLICATION_MAX_SEQUENCES = 80
CLOUD_PUBLICATION_MAX_RESIDUES = 60000
"""
new = """# Workstation/local Auto policy.
AUTO_IQTREE_MAX_SEQUENCES = 60
AUTO_IQTREE_MAX_RESIDUES = 40000

# Streamlit Community Cloud is substantially tighter than a workstation.
# A real 40-protein / 8,411-aa family exhausted the shared app while
# ModelFinder + UFBoot + SH-aLRT were running, so Cloud Auto is deliberately
# conservative and falls back to MAFFT + FastTree above this envelope.
CLOUD_AUTO_IQTREE_MAX_SEQUENCES = 20
CLOUD_AUTO_IQTREE_MAX_RESIDUES = 5000
CLOUD_PUBLICATION_MAX_SEQUENCES = 20
CLOUD_PUBLICATION_MAX_RESIDUES = 5000
"""
if old not in s:
    raise SystemExit('threshold block not found')
s = s.replace(old, new, 1)

old = """    if mode == 'auto':
        small_enough_for_iqtree = (
            n <= AUTO_IQTREE_MAX_SEQUENCES
            and total_residues <= AUTO_IQTREE_MAX_RESIDUES
        )
"""
new = """    if mode == 'auto':
        if _is_shared_cloud():
            auto_max_sequences = CLOUD_AUTO_IQTREE_MAX_SEQUENCES
            auto_max_residues = CLOUD_AUTO_IQTREE_MAX_RESIDUES
        else:
            auto_max_sequences = AUTO_IQTREE_MAX_SEQUENCES
            auto_max_residues = AUTO_IQTREE_MAX_RESIDUES

        small_enough_for_iqtree = (
            n <= auto_max_sequences
            and total_residues <= auto_max_residues
        )
"""
if old not in s:
    raise SystemExit('auto threshold block not found')
s = s.replace(old, new, 1)

old = """        if publication_phylogeny_ready():
            return run_mafft_iqtree(
                proteins,
                bootstrap=bootstrap,
                alrt=alrt,
                threads=threads,
                timeout=1800 if _is_shared_cloud() else 7200,
            )
"""
new = """        # Do not silently fall back to a resource-heavy IQ-TREE run on
        # shared Cloud when FastTree is unavailable. A controlled error is
        # safer than killing the whole Streamlit process.
        if publication_phylogeny_ready() and not _is_shared_cloud():
            return run_mafft_iqtree(
                proteins,
                bootstrap=bootstrap,
                alrt=alrt,
                threads=threads,
                timeout=7200,
            )
        if _is_shared_cloud() and publication_phylogeny_ready():
            raise RuntimeError(
                'Cloud Auto selected the safe FastTree route, but MAFFT/FastTree is unavailable. '
                'Do not fall back to IQ-TREE on shared Cloud; restore FastTree or run Publication locally/HPC.'
            )
"""
if old not in s:
    raise SystemExit('auto IQ-TREE fallback block not found')
s = s.replace(old, new, 1)
p.write_text(s)

# 2) UI preflight/progress must use the exact same cloud thresholds.
p = Path('pages/2_Gene_Structure_Domain_Motif.py')
s = p.read_text()
old = "from modules.gdm_phylogeny import build_phylogeny, external_phylogeny_ready, publication_phylogeny_ready, phylogeny_tool_status, AUTO_IQTREE_MAX_SEQUENCES, AUTO_IQTREE_MAX_RESIDUES, CLOUD_PUBLICATION_MAX_SEQUENCES, CLOUD_PUBLICATION_MAX_RESIDUES"
new = "from modules.gdm_phylogeny import build_phylogeny, external_phylogeny_ready, publication_phylogeny_ready, phylogeny_tool_status, AUTO_IQTREE_MAX_SEQUENCES, AUTO_IQTREE_MAX_RESIDUES, CLOUD_AUTO_IQTREE_MAX_SEQUENCES, CLOUD_AUTO_IQTREE_MAX_RESIDUES, CLOUD_PUBLICATION_MAX_SEQUENCES, CLOUD_PUBLICATION_MAX_RESIDUES"
if old not in s:
    raise SystemExit('phylogeny import line not found')
s = s.replace(old, new, 1)

old = """                elif phylo_mode == 'auto':
                    if preview_n <= AUTO_IQTREE_MAX_SEQUENCES and preview_residues <= AUTO_IQTREE_MAX_RESIDUES:
                        tree_plan = 'Auto → MAFFT + IQ-TREE publication-oriented inference'
                    else:
                        tree_plan = 'Auto → MAFFT + FastTree cloud-safe screening'
                    tree_level = 'success'
"""
new = """                elif phylo_mode == 'auto':
                    auto_seq_limit = CLOUD_AUTO_IQTREE_MAX_SEQUENCES if shared_cloud else AUTO_IQTREE_MAX_SEQUENCES
                    auto_residue_limit = CLOUD_AUTO_IQTREE_MAX_RESIDUES if shared_cloud else AUTO_IQTREE_MAX_RESIDUES
                    if preview_n <= auto_seq_limit and preview_residues <= auto_residue_limit:
                        tree_plan = 'Auto → MAFFT + IQ-TREE publication-oriented inference'
                    else:
                        tree_plan = 'Auto → MAFFT + FastTree cloud-safe screening'
                    tree_level = 'success'
"""
if old not in s:
    raise SystemExit('preflight auto route block not found')
s = s.replace(old, new, 1)

old = """        elif phylo_mode == 'auto' and (family_nseq > 60 or family_residues > 40000):
            phylo_progress_text = (
                f'Phylogeny: large family ({family_nseq} proteins, {family_residues:,} aa) → '
                'MAFFT + FastTree cloud-safe screening...'
            )
"""
new = """        elif phylo_mode == 'auto' and (
            family_nseq > (CLOUD_AUTO_IQTREE_MAX_SEQUENCES if Path('/mount/src').exists() else AUTO_IQTREE_MAX_SEQUENCES)
            or family_residues > (CLOUD_AUTO_IQTREE_MAX_RESIDUES if Path('/mount/src').exists() else AUTO_IQTREE_MAX_RESIDUES)
        ):
            phylo_progress_text = (
                f'Phylogeny: cloud-safe route ({family_nseq} proteins, {family_residues:,} aa) → '
                'MAFFT + FastTree screening...'
            )
"""
if old not in s:
    raise SystemExit('progress auto route block not found')
s = s.replace(old, new, 1)
p.write_text(s)

# 3) Regression test from the real failing envelope: 40 proteins on shared cloud must never start IQ-TREE.
p = Path('tests/test_large_family_stability.py')
s = p.read_text()
anchor = """    def test_auto_uses_iqtree_for_small_family(self):
"""
new_test = """    def test_cloud_auto_uses_fasttree_for_40_proteins(self):
        sentinel = {
            'tree_text': '(A:1,B:1,C:1);',
            'alignment_text': '>A\\nAAA\\n',
            'method': 'MAFFT + FastTree',
            'qc': None,
            'warning': '',
            'log_text': '',
        }
        with patch.object(gp, '_is_shared_cloud', return_value=True), \\
             patch.object(gp, 'publication_phylogeny_ready', return_value=True), \\
             patch.object(gp, 'external_phylogeny_ready', return_value=True), \\
             patch.object(gp, 'run_mafft_iqtree', side_effect=AssertionError('IQ-TREE must not run for 40 proteins on shared Cloud')), \\
             patch.object(gp, 'run_mafft_fasttree', return_value=sentinel.copy()) as ft:
            result = gp.build_phylogeny(proteins(40, 211), mode='auto')
        self.assertEqual(result['method'], 'MAFFT + FastTree')
        ft.assert_called_once()

    def test_cloud_publication_guard_covers_real_40_protein_failure(self):
        with patch.object(gp, '_is_shared_cloud', return_value=True), \\
             patch.object(gp, 'publication_phylogeny_ready', return_value=True), \\
             patch.object(gp, 'run_mafft_iqtree') as iq:
            with self.assertRaisesRegex(RuntimeError, 'too resource-intensive'):
                gp.build_phylogeny(proteins(40, 211), mode='publication')
        iq.assert_not_called()

""" + anchor
if anchor not in s:
    raise SystemExit('test insertion anchor not found')
s = s.replace(anchor, new_test, 1)
p.write_text(s)

print('CLOUD_IQTREE_THRESHOLDS_V2_PATCH_COMPLETE')
