from pathlib import Path


def replace_once(text, old, new, label):
    if old not in text:
        raise SystemExit(f'Expected block not found: {label}')
    return text.replace(old, new, 1)

# -----------------------------------------------------------------------------
# Cloud phylogeny policy: IQ-TREE is intentionally local/HPC only.
# -----------------------------------------------------------------------------
p = Path('modules/gdm_phylogeny.py')
s = p.read_text()

s = replace_once(
    s,
    "CLOUD_AUTO_IQTREE_MAX_SEQUENCES = 20\nCLOUD_AUTO_IQTREE_MAX_RESIDUES = 5000\nCLOUD_PUBLICATION_MAX_SEQUENCES = 20\nCLOUD_PUBLICATION_MAX_RESIDUES = 5000\n",
    "CLOUD_AUTO_IQTREE_MAX_SEQUENCES = 0\nCLOUD_AUTO_IQTREE_MAX_RESIDUES = 0\nCLOUD_PUBLICATION_MAX_SEQUENCES = 0\nCLOUD_PUBLICATION_MAX_RESIDUES = 0\n",
    'cloud IQ-TREE constants',
)

old_publication = """    if mode == 'publication':
        if not publication_phylogeny_ready():
            raise RuntimeError('Publication mode requires MAFFT + IQ-TREE.')
        if _is_shared_cloud() and (
            n > CLOUD_PUBLICATION_MAX_SEQUENCES
            or total_residues > CLOUD_PUBLICATION_MAX_RESIDUES
        ):
            raise RuntimeError(
                'Publication IQ-TREE mode is too resource-intensive for this family on shared Streamlit Cloud '
                f'({n} sequences; {total_residues:,} aa). Use Auto/FastTree here, or run Publication mode '
                'locally/HPC and upload the validated Newick tree.'
            )
        return run_mafft_iqtree(
"""
new_publication = """    if mode == 'publication':
        if _is_shared_cloud():
            raise RuntimeError(
                'Publication IQ-TREE is intentionally disabled on shared Streamlit Cloud because '
                'ModelFinder + resampling can terminate the shared app process. Use Auto/FastTree '
                'here, or run Publication mode locally/HPC and upload the validated Newick tree.'
            )
        if not publication_phylogeny_ready():
            raise RuntimeError('Publication mode requires MAFFT + IQ-TREE.')
        return run_mafft_iqtree(
"""
s = replace_once(s, old_publication, new_publication, 'publication cloud guard')

old_auto = """    if mode == 'auto':
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

        if small_enough_for_iqtree and publication_phylogeny_ready():
            result = run_mafft_iqtree(
                proteins,
                bootstrap=bootstrap,
                alrt=alrt,
                threads=threads,
                timeout=1800 if _is_shared_cloud() else 7200,
            )
            return result

        if external_phylogeny_ready():
            result = run_mafft_fasttree(proteins)
            result['warning'] = (
                f'Large-family Auto mode selected MAFFT + FastTree for {n} sequences '
                f'({total_residues:,} aa) to avoid CPU/memory stalls. The full alignment and Newick tree are preserved. '
                'For final publication inference, run IQ-TREE Publication mode locally/HPC or upload a validated Newick tree.'
            )
            return result

        # Do not silently fall back to a resource-heavy IQ-TREE run on
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

    return run_nj_fallback(proteins)
"""
new_auto = """    if mode == 'auto':
        if _is_shared_cloud():
            if external_phylogeny_ready():
                result = run_mafft_fasttree(proteins)
                result['warning'] = (
                    f'Shared-Cloud Auto mode used MAFFT + FastTree for {n} sequences '
                    f'({total_residues:,} aa). IQ-TREE publication inference is intentionally '
                    'disabled on shared Cloud to protect app stability. The full alignment and '
                    'Newick tree are preserved; run IQ-TREE locally/HPC for final publication inference.'
                )
                return result
            return run_nj_fallback(proteins)

        small_enough_for_iqtree = (
            n <= AUTO_IQTREE_MAX_SEQUENCES
            and total_residues <= AUTO_IQTREE_MAX_RESIDUES
        )
        if small_enough_for_iqtree and publication_phylogeny_ready():
            return run_mafft_iqtree(
                proteins,
                bootstrap=bootstrap,
                alrt=alrt,
                threads=threads,
                timeout=7200,
            )
        if external_phylogeny_ready():
            result = run_mafft_fasttree(proteins)
            result['warning'] = (
                f'Large-family Auto mode selected MAFFT + FastTree for {n} sequences '
                f'({total_residues:,} aa). The full alignment and Newick tree are preserved. '
                'For final publication inference, run IQ-TREE Publication mode locally/HPC.'
            )
            return result
        if publication_phylogeny_ready():
            return run_mafft_iqtree(
                proteins,
                bootstrap=bootstrap,
                alrt=alrt,
                threads=threads,
                timeout=7200,
            )

    return run_nj_fallback(proteins)
"""
s = replace_once(s, old_auto, new_auto, 'auto cloud route')
p.write_text(s)

# -----------------------------------------------------------------------------
# Streamlit page: make Cloud mode explicit and remove unsafe controls.
# -----------------------------------------------------------------------------
p = Path('pages/2_Gene_Structure_Domain_Motif.py')
s = p.read_text()

old_metrics = """c1, c2, c3, c4, c5 = st.columns(5)
c1.metric('Phylogeny', 'IQ-TREE publication mode' if publication_phylogeny_ready() else ('MAFFT + FastTree' if external_phylogeny_ready() else 'NJ fallback available'))
"""
new_metrics = """shared_cloud_runtime = bool(phylogeny_tool_status().get('shared_cloud', False))

c1, c2, c3, c4, c5 = st.columns(5)
if shared_cloud_runtime:
    c1.metric('Phylogeny', 'MAFFT + FastTree cloud mode' if external_phylogeny_ready() else 'NJ fallback available')
else:
    c1.metric('Phylogeny', 'IQ-TREE publication mode' if publication_phylogeny_ready() else ('MAFFT + FastTree' if external_phylogeny_ready() else 'NJ fallback available'))
"""
s = replace_once(s, old_metrics, new_metrics, 'cloud metric')

old_info = """if not publication_phylogeny_ready():
    st.info('Publication phylogeny needs MAFFT + IQ-TREE. If IQ-TREE is unavailable, Auto mode falls back to MAFFT + FastTree and then NJ screening.')
"""
new_info = """if shared_cloud_runtime:
    st.info('Shared Streamlit Cloud uses MAFFT + FastTree for automatic phylogeny. IQ-TREE publication inference is intentionally local/HPC-only to prevent shared-process crashes; a validated Newick tree can be uploaded here.')
elif not publication_phylogeny_ready():
    st.info('Publication phylogeny needs MAFFT + IQ-TREE. If IQ-TREE is unavailable, Auto mode falls back to MAFFT + FastTree and then NJ screening.')
"""
s = replace_once(s, old_info, new_info, 'cloud publication info')

s = replace_once(
    s,
    "reference_threads = st.slider('Reference mapping threads', 1, 16, 4)",
    "reference_threads = st.slider('Reference mapping threads', 1, 16, 2 if shared_cloud_runtime else 4)",
    'reference thread default',
)

old_phylo_ui = """        phylo_mode = st.selectbox(
            'Automatic tree method', ['auto', 'publication', 'fasttree', 'nj'],
            format_func=lambda x: {
                'auto': 'Auto — adaptive: IQ-TREE for small families, FastTree for larger families',
                'publication': 'Publication — MAFFT + IQ-TREE + ModelFinder + SH-aLRT + ultrafast bootstrap',
                'fasttree': 'Fast screening — MAFFT + FastTree',
                'nj': 'Internal Neighbor-Joining — quick screening only'
            }[x]
        )
        pc1, pc2, pc3 = st.columns(3)
        with pc1:
            phylo_bootstrap = st.selectbox('Ultrafast bootstrap replicates', [1000, 2000, 5000], index=0)
        with pc2:
            phylo_alrt = st.selectbox('SH-aLRT replicates', [1000, 2000, 5000], index=0)
        with pc3:
            phylo_threads = st.selectbox('IQ-TREE threads', ['2', '4', '8', 'AUTO'], index=0)
        st.caption('Publication mode uses MAFFT alignment followed by IQ-TREE maximum-likelihood inference, automatic ModelFinder selection, SH-aLRT support and ultrafast bootstrap. FastTree remains a screening option.')
"""
new_phylo_ui = """        phylo_options = ['auto', 'fasttree', 'nj'] if shared_cloud_runtime else ['auto', 'publication', 'fasttree', 'nj']
        phylo_mode = st.selectbox(
            'Automatic tree method', phylo_options,
            format_func=lambda x: {
                'auto': 'Auto — MAFFT + FastTree cloud-safe' if shared_cloud_runtime else 'Auto — adaptive: IQ-TREE for small families, FastTree for larger families',
                'publication': 'Publication — MAFFT + IQ-TREE + ModelFinder + SH-aLRT + ultrafast bootstrap',
                'fasttree': 'Fast screening — MAFFT + FastTree',
                'nj': 'Internal Neighbor-Joining — quick screening only'
            }[x]
        )
        if shared_cloud_runtime:
            phylo_bootstrap, phylo_alrt, phylo_threads = 1000, 1000, '2'
            st.caption('Cloud stability mode: IQ-TREE controls are hidden because publication inference is local/HPC-only. Auto uses MAFFT + FastTree without subsampling.')
        else:
            pc1, pc2, pc3 = st.columns(3)
            with pc1:
                phylo_bootstrap = st.selectbox('Ultrafast bootstrap replicates', [1000, 2000, 5000], index=0)
            with pc2:
                phylo_alrt = st.selectbox('SH-aLRT replicates', [1000, 2000, 5000], index=0)
            with pc3:
                phylo_threads = st.selectbox('IQ-TREE threads', ['2', '4', '8', 'AUTO'], index=0)
            st.caption('Publication mode uses MAFFT alignment followed by IQ-TREE maximum-likelihood inference, automatic ModelFinder selection, SH-aLRT support and ultrafast bootstrap. FastTree remains a screening option.')
"""
s = replace_once(s, old_phylo_ui, new_phylo_ui, 'phylogeny UI')

s = s.replace("shared_cloud = Path('/mount/src').exists()", "shared_cloud = shared_cloud_runtime")

old_plan = """                elif phylo_mode == 'auto':
                    auto_seq_limit = CLOUD_AUTO_IQTREE_MAX_SEQUENCES if shared_cloud else AUTO_IQTREE_MAX_SEQUENCES
                    auto_residue_limit = CLOUD_AUTO_IQTREE_MAX_RESIDUES if shared_cloud else AUTO_IQTREE_MAX_RESIDUES
                    if preview_n <= auto_seq_limit and preview_residues <= auto_residue_limit:
                        tree_plan = 'Auto → MAFFT + IQ-TREE publication-oriented inference'
                    else:
                        tree_plan = 'Auto → MAFFT + FastTree cloud-safe screening'
                    tree_level = 'success'
"""
new_plan = """                elif phylo_mode == 'auto':
                    if shared_cloud:
                        tree_plan = 'Auto → MAFFT + FastTree cloud-safe screening'
                    elif preview_n <= AUTO_IQTREE_MAX_SEQUENCES and preview_residues <= AUTO_IQTREE_MAX_RESIDUES:
                        tree_plan = 'Auto → MAFFT + IQ-TREE publication-oriented inference'
                    else:
                        tree_plan = 'Auto → MAFFT + FastTree screening'
                    tree_level = 'success'
"""
s = replace_once(s, old_plan, new_plan, 'preflight route')

old_progress = """        elif phylo_mode == 'auto' and (
            family_nseq > (CLOUD_AUTO_IQTREE_MAX_SEQUENCES if Path('/mount/src').exists() else AUTO_IQTREE_MAX_SEQUENCES)
            or family_residues > (CLOUD_AUTO_IQTREE_MAX_RESIDUES if Path('/mount/src').exists() else AUTO_IQTREE_MAX_RESIDUES)
        ):
            phylo_progress_text = (
                f'Phylogeny: cloud-safe route ({family_nseq} proteins, {family_residues:,} aa) → '
                'MAFFT + FastTree screening...'
            )
"""
new_progress = """        elif phylo_mode == 'auto' and shared_cloud_runtime:
            phylo_progress_text = (
                f'Phylogeny: shared-Cloud safe route ({family_nseq} proteins, {family_residues:,} aa) → '
                'MAFFT + FastTree screening...'
            )
        elif phylo_mode == 'auto' and (
            family_nseq > AUTO_IQTREE_MAX_SEQUENCES or family_residues > AUTO_IQTREE_MAX_RESIDUES
        ):
            phylo_progress_text = (
                f'Phylogeny: large-family route ({family_nseq} proteins, {family_residues:,} aa) → '
                'MAFFT + FastTree screening...'
            )
"""
s = replace_once(s, old_progress, new_progress, 'progress route')

p.write_text(s)

# -----------------------------------------------------------------------------
# Regression tests for cloud policy.
# -----------------------------------------------------------------------------
p = Path('tests/test_large_family_stability.py')
s = p.read_text()
anchor = """    def test_cloud_auto_uses_fasttree_for_40_proteins(self):
"""
if anchor not in s:
    raise SystemExit('Expected test anchor not found')
insert = """    def test_cloud_auto_uses_fasttree_even_for_tiny_family(self):
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
             patch.object(gp, 'run_mafft_iqtree', side_effect=AssertionError('IQ-TREE must never run in shared-Cloud Auto mode')), \\
             patch.object(gp, 'run_mafft_fasttree', return_value=sentinel.copy()) as ft:
            result = gp.build_phylogeny(proteins(5, 120), mode='auto')
        self.assertEqual(result['method'], 'MAFFT + FastTree')
        ft.assert_called_once()

"""
s = s.replace(anchor, insert + anchor, 1)
p.write_text(s)

print('cloud stable runtime v3 patch applied')
