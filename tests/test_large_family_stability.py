import io
import unittest
from unittest.mock import patch

from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from modules import gdm_phylogeny as gp
from modules.gdm_plot import fig_bytes, MAX_RASTER_PIXELS


def proteins(n, length=120):
    base = ('ACDEFGHIKLMNPQRSTVWY' * ((length // 20) + 1))[:length]
    return {f'Gene{i:04d}': base for i in range(n)}


class LargeFamilyPhylogenyTests(unittest.TestCase):
    def test_auto_uses_fasttree_for_large_family(self):
        sentinel = {
            'tree_text': '(A:1,B:1,C:1);',
            'alignment_text': '>A\nAAA\n',
            'method': 'MAFFT + FastTree',
            'qc': None,
            'warning': '',
            'log_text': '',
        }
        with patch.object(gp, 'publication_phylogeny_ready', return_value=True), \
             patch.object(gp, 'external_phylogeny_ready', return_value=True), \
             patch.object(gp, 'run_mafft_iqtree', side_effect=AssertionError('IQ-TREE should not run')), \
             patch.object(gp, 'run_mafft_fasttree', return_value=sentinel.copy()) as ft:
            result = gp.build_phylogeny(proteins(120), mode='auto')
        self.assertEqual(result['method'], 'MAFFT + FastTree')
        self.assertIn('Large-family Auto mode', result['warning'])
        ft.assert_called_once()

    def test_auto_uses_iqtree_for_small_family(self):
        sentinel = {
            'tree_text': '(A:1,B:1,C:1);',
            'alignment_text': '>A\nAAA\n',
            'method': 'MAFFT + IQ-TREE',
            'qc': None,
            'warning': '',
            'log_text': '',
        }
        with patch.object(gp, 'publication_phylogeny_ready', return_value=True), \
             patch.object(gp, 'run_mafft_iqtree', return_value=sentinel.copy()) as iq:
            result = gp.build_phylogeny(proteins(12), mode='auto')
        self.assertEqual(result['method'], 'MAFFT + IQ-TREE')
        iq.assert_called_once()

    def test_cloud_blocks_oversized_explicit_publication_run(self):
        with patch.object(gp, '_is_shared_cloud', return_value=True), \
             patch.object(gp, 'publication_phylogeny_ready', return_value=True):
            with self.assertRaisesRegex(RuntimeError, 'too resource-intensive'):
                gp.build_phylogeny(proteins(120), mode='publication')

    def test_cloud_auto_threads_are_capped(self):
        with patch.object(gp, '_is_shared_cloud', return_value=True):
            self.assertEqual(gp._safe_iqtree_threads('AUTO'), '2')
            self.assertEqual(gp._safe_iqtree_threads('8'), '4')


class PlotMemoryTests(unittest.TestCase):
    def test_png_raster_is_pixel_capped(self):
        fig, ax = plt.subplots(figsize=(20, 40))
        ax.plot([0, 1], [0, 1])
        data = fig_bytes(fig, 'png', 600)
        plt.close(fig)
        image = Image.open(io.BytesIO(data))
        self.assertLessEqual(image.width * image.height, int(MAX_RASTER_PIXELS * 1.15))

    def test_tiff_serializes(self):
        fig, ax = plt.subplots(figsize=(12, 20))
        ax.plot([0, 1], [1, 0])
        data = fig_bytes(fig, 'tiff', 600)
        plt.close(fig)
        self.assertGreater(len(data), 100)


if __name__ == '__main__':
    unittest.main()
