import unittest
from unittest.mock import patch

from modules import gdm_phylogeny as gp


def proteins(n, length=120):
    base = ('ACDEFGHIKLMNPQRSTVWY' * ((length // 20) + 1))[:length]
    return {f'Gene{i:04d}': base for i in range(n)}


class CloudRuntimePolicyTests(unittest.TestCase):
    def test_cloud_auto_never_calls_iqtree(self):
        sentinel = {
            'tree_text': '(A:1,B:1,C:1);',
            'alignment_text': '>A\nAAA\n',
            'method': 'MAFFT + FastTree',
            'qc': None,
            'warning': '',
            'log_text': '',
        }
        with patch.object(gp, '_is_shared_cloud', return_value=True), \
             patch.object(gp, 'external_phylogeny_ready', return_value=True), \
             patch.object(gp, 'publication_phylogeny_ready', return_value=True), \
             patch.object(gp, 'run_mafft_iqtree', side_effect=AssertionError('IQ-TREE must not run on shared Cloud')), \
             patch.object(gp, 'run_mafft_fasttree', return_value=sentinel.copy()) as ft:
            result = gp.build_phylogeny(proteins(3), mode='auto')
        self.assertEqual(result['method'], 'MAFFT + FastTree')
        ft.assert_called_once()

    def test_cloud_publication_is_controlled_error(self):
        with patch.object(gp, '_is_shared_cloud', return_value=True), \
             patch.object(gp, 'run_mafft_iqtree') as iq:
            with self.assertRaisesRegex(RuntimeError, 'intentionally disabled'):
                gp.build_phylogeny(proteins(3), mode='publication')
        iq.assert_not_called()


if __name__ == '__main__':
    unittest.main()
