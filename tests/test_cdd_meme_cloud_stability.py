import unittest
from unittest.mock import patch

import pandas as pd

from modules import gdm_cdd_meme as cm


def proteins(n, length=300):
    seq = ('ACDEFGHIKLMNPQRSTVWY' * ((length // 20) + 1))[:length]
    return {f'Gene{i:04d}': seq for i in range(n)}


class CDDBatchingTests(unittest.TestCase):
    def test_large_family_is_batched_without_subsampling(self):
        calls = []
        progress_events = []

        def fake_batch(batch, evalue, deadline):
            calls.append(list(batch))
            gene = next(iter(batch))
            frame = pd.DataFrame([{
                'gene': gene,
                'query': gene,
                'hit_type': 'specific',
                'pssm_id': '1',
                'start': 1,
                'end': 10,
                'evalue': 1e-8,
                'bitscore': 50.0,
                'accession': 'cd00001',
                'domain': 'TEST',
                'incomplete': '',
                'superfamily': '',
                'confidence': 'HIGH',
            }])
            return frame, f'RID{len(calls)}', f'raw-{len(calls)}'

        family = proteins(450, 120)
        with patch.object(cm, '_run_cdd_batch', side_effect=fake_batch), \
             patch.object(cm.time, 'sleep', return_value=None):
            df, rid_text, raw = cm.run_cdd(
                family,
                batch_size=200,
                timeout=60,
                progress=lambda done, total, message: progress_events.append((done, total, message)),
            )

        self.assertEqual([len(x) for x in calls], [200, 200, 50])
        flattened = [gene for batch in calls for gene in batch]
        self.assertEqual(flattened, list(family))
        self.assertEqual(len(df), 3)
        self.assertEqual(rid_text, 'batch1:RID1;batch2:RID2;batch3:RID3')
        self.assertIn('batch 3/3', raw)
        self.assertEqual(progress_events[0][0:2], (0, 3))
        self.assertEqual(progress_events[-1][0:2], (3, 3))
        self.assertIn('Completed NCBI CDD batch 3/3', progress_events[-1][2])

    def test_bad_progress_callback_does_not_break_science(self):
        frame = pd.DataFrame([{
            'gene': 'Gene0000', 'query': 'Gene0000', 'hit_type': 'specific',
            'pssm_id': '1', 'start': 1, 'end': 10, 'evalue': 1e-8,
            'bitscore': 50.0, 'accession': 'cd00001', 'domain': 'TEST',
            'incomplete': '', 'superfamily': '', 'confidence': 'HIGH',
        }])
        with patch.object(cm, '_run_cdd_batch', return_value=(frame, 'RID1', 'raw')):
            df, rid, _ = cm.run_cdd(
                proteins(3, 40), timeout=60,
                progress=lambda *args: (_ for _ in ()).throw(RuntimeError('UI callback failed')),
            )
        self.assertEqual(rid, 'RID1')
        self.assertEqual(len(df), 1)

    def test_more_than_supported_max_is_rejected_clearly(self):
        with self.assertRaisesRegex(ValueError, '1000'):
            cm.run_cdd(proteins(1001, 20))


class MEMECloudPolicyTests(unittest.TestCase):
    def test_small_family_allowed(self):
        allowed, nseq, residues = cm.meme_cloud_policy(proteins(50, 200))
        self.assertTrue(allowed)
        self.assertEqual(nseq, 50)
        self.assertEqual(residues, 10000)

    def test_large_family_rejected_without_subsampling(self):
        family = proteins(200, 400)
        with patch.object(cm, '_is_shared_cloud', return_value=True), \
             patch.object(cm, 'meme_ready', return_value=True), \
             patch.object(cm.subprocess, 'run') as run:
            with self.assertRaisesRegex(RuntimeError, 'CLOUD_RESOURCE_LIMIT'):
                cm.run_meme(family, nmotifs=10)
        run.assert_not_called()

    def test_local_policy_does_not_change_sequences(self):
        family = proteins(200, 400)
        allowed, nseq, residues = cm.meme_cloud_policy(family)
        self.assertFalse(allowed)
        self.assertEqual(nseq, len(family))
        self.assertEqual(residues, 80000)


if __name__ == '__main__':
    unittest.main()
