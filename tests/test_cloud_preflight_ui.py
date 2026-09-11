from pathlib import Path
import unittest


class CloudPreflightUITests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.page = Path('pages/2_Gene_Structure_Domain_Motif.py').read_text()

    def test_preflight_is_visible_before_run_button(self):
        preflight = self.page.index('🧪 Run preflight & resource plan')
        run_button = self.page.index("🚀 Run complete analysis")
        self.assertLess(preflight, run_button)

    def test_preflight_uses_engine_threshold_constants(self):
        for token in [
            'AUTO_IQTREE_MAX_SEQUENCES',
            'AUTO_IQTREE_MAX_RESIDUES',
            'CLOUD_PUBLICATION_MAX_SEQUENCES',
            'CLOUD_PUBLICATION_MAX_RESIDUES',
            'CDD_BATCH_SIZE',
            'CDD_MAX_SEQUENCES',
            'CLOUD_MEME_MAX_SEQUENCES',
            'CLOUD_MEME_MAX_RESIDUES',
        ]:
            self.assertIn(token, self.page)

    def test_preflight_never_claims_subsampling(self):
        self.assertIn('Preflight does not alter or subsample your sequences', self.page)
        self.assertIn('no subsampling', self.page)
        self.assertIn('Nothing will be subsampled', self.page)

    def test_large_auto_plan_is_fasttree(self):
        self.assertIn('Auto → MAFFT + FastTree cloud-safe screening', self.page)
        self.assertIn('Auto → MAFFT + IQ-TREE publication-oriented inference', self.page)


if __name__ == '__main__':
    unittest.main()
