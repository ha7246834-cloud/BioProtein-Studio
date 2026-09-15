import unittest
from unittest.mock import patch

from modules import gdm_reference as gr
from modules import gdm_structure as gs


class ReferenceThreadPolicyTests(unittest.TestCase):
    def test_cloud_miniprot_threads_are_capped(self):
        with patch.object(gr, '_is_shared_cloud', return_value=True):
            self.assertEqual(gr._safe_miniprot_threads(1), 1)
            self.assertEqual(gr._safe_miniprot_threads(4), 2)
            self.assertEqual(gr._safe_miniprot_threads(16), 2)

    def test_local_miniprot_threads_are_preserved(self):
        with patch.object(gr, '_is_shared_cloud', return_value=False):
            self.assertEqual(gr._safe_miniprot_threads(8), 8)


class Est2GenomeFastPathTests(unittest.TestCase):
    def test_exact_contiguous_cds_skips_external_emboss(self):
        cds = 'ATGGCCGCCGCC'
        genomic = 'TTTT' + cds + 'AAAA'
        with patch.object(gs, 'est2genome_ready', return_value=True), \
             patch.object(gs.runtime, 'run_guarded', side_effect=AssertionError('est2genome should not launch')):
            df, qc, raw = gs.est2genome_pair('Gene1', cds, genomic)
        self.assertEqual(qc['method'], 'Exact contiguous match')
        self.assertEqual(qc['status'], 'PASS')
        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]['source'], 'Exact contiguous match')
        self.assertIn('no intron inferred', raw)

    def test_reverse_complement_exact_route_is_preserved(self):
        cds = 'ATGGCCGCCGCC'
        from Bio.Seq import Seq
        genomic = 'TTTT' + str(Seq(cds).reverse_complement()) + 'AAAA'
        with patch.object(gs.runtime, 'run_guarded', side_effect=AssertionError('est2genome should not launch')):
            df, qc, _ = gs.est2genome_pair('Gene2', cds, genomic)
        self.assertEqual(qc['status'], 'PASS')
        self.assertEqual(df.iloc[0]['strand'], '-')


if __name__ == '__main__':
    unittest.main()
