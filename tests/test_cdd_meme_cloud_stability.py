import os
from pathlib import Path
import subprocess
import tempfile
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
        allowed, nseq, residues = cm.meme_cloud_policy(proteins(20, 200))
        self.assertTrue(allowed)
        self.assertEqual(nseq, 20)
        self.assertEqual(residues, 4000)

    def test_real_40_protein_case_is_allowed_on_cloud(self):
        """Measured MEME memory for 40 proteins / 8,411 aa is ~11 MB, so the
        family is within the evidence-based cloud envelope and must reach the
        guarded runner rather than being blocked."""
        family = proteins(40, 211)
        allowed, nseq, residues = cm.meme_cloud_policy(family)
        self.assertTrue(allowed)
        self.assertEqual(nseq, 40)
        self.assertEqual(residues, 8440)
        with patch.object(cm, '_is_shared_cloud', return_value=True), \
             patch.object(cm, 'meme_ready', return_value=True), \
             patch.object(cm.runtime, 'run_guarded', side_effect=RuntimeError('reached-runner')) as rg:
            with self.assertRaisesRegex(RuntimeError, 'reached-runner'):
                cm.run_meme(family, nmotifs=10)
        rg.assert_called_once()

    def test_oversize_family_rejected_without_subsampling(self):
        """A genuinely oversized family (beyond the wall-time envelope) is still
        refused, and never reaches the runner or subsamples sequences."""
        family = proteins(200, 400)  # 200 seq / 80,000 aa > 120 / 30,000
        with patch.object(cm, '_is_shared_cloud', return_value=True), \
             patch.object(cm, 'meme_ready', return_value=True), \
             patch.object(cm.runtime, 'run_guarded') as rg:
            with self.assertRaisesRegex(RuntimeError, 'CLOUD_RESOURCE_LIMIT'):
                cm.run_meme(family, nmotifs=10)
        rg.assert_not_called()

    def test_policy_counts_sequences_and_residues(self):
        allowed, nseq, residues = cm.meme_cloud_policy(proteins(200, 400))
        self.assertFalse(allowed)
        self.assertEqual(nseq, 200)
        self.assertEqual(residues, 80000)

    def test_vendor_wrapper_allows_40_and_blocks_oversized(self):
        """The shared-cloud MEME wrapper allows the real 40-protein family and
        blocks only genuinely oversized input (BPS_MEME_GATE_ONLY validates the
        gate without launching the slow core)."""
        wrapper = Path('vendor/bin/meme').resolve()
        self.assertTrue(wrapper.exists())
        small = proteins(40, 211)
        big = proteins(200, 400)
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            env = os.environ.copy()
            env['BPS_SHARED_CLOUD'] = '1'
            env['BPS_MEME_GATE_ONLY'] = '1'

            fa = td / 'small.faa'
            fa.write_text(''.join(f'>{n}\n{s}\n' for n, s in small.items()))
            ok = subprocess.run(
                [str(wrapper), str(fa), '-protein', '-oc', str(td / 'o1'), '-nmotifs', '5'],
                capture_output=True, text=True, env=env, timeout=20,
            )
            self.assertEqual(ok.returncode, 0)
            self.assertNotIn('CLOUD_RESOURCE_LIMIT', ok.stderr)

            fb = td / 'big.faa'
            fb.write_text(''.join(f'>{n}\n{s}\n' for n, s in big.items()))
            blocked = subprocess.run(
                [str(wrapper), str(fb), '-protein', '-oc', str(td / 'o2'), '-nmotifs', '5'],
                capture_output=True, text=True, env=env, timeout=20,
            )
        self.assertEqual(blocked.returncode, 75)
        self.assertIn('CLOUD_RESOURCE_LIMIT', blocked.stderr)
        self.assertIn('200 proteins', blocked.stderr)

    def test_vendor_wrapper_stages_embedded_prefix_assets(self):
        """The portable MEME core has a compiled installation prefix. Ensure the
        wrapper stages prior30.plib/template.eps at that exact prefix before the
        core is launched, preventing the production prior-library failure."""
        wrapper = Path('vendor/bin/meme').resolve()
        self.assertTrue(wrapper.exists())
        family = proteins(3, 80)
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            fa = td / 'tiny.faa'
            fa.write_text(''.join(f'>{n}\n{s}\n' for n, s in family.items()))
            env = os.environ.copy()
            env['BPS_MEME_PREFIX_CHECK_ONLY'] = '1'
            proc = subprocess.run(
                [str(wrapper), str(fa), '-protein', '-oc', str(td / 'out'), '-nmotifs', '3'],
                capture_output=True, text=True, env=env, timeout=20,
            )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        prefix = Path('/tmp/bps_meme_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
        self.assertTrue((prefix / 'share/meme-5.5.9/prior30.plib').is_file())
        self.assertTrue((prefix / 'share/meme-5.5.9/template.eps').is_file())


if __name__ == '__main__':
    unittest.main()
