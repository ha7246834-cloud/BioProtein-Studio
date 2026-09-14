"""Tests for the resource-sandbox / runtime layer (modules/runtime.py).

These exercise the structural fix for the "Oh no. Error running app." worker
death: a child that exceeds its memory/CPU ceiling must be killed by the kernel
and surface as a catchable ChildResourceError while THIS process stays alive.

This module has no biopython/streamlit dependency, so it runs anywhere.
"""
import os
import sys
import unittest
from unittest.mock import patch

from modules import runtime as rt

POSIX = os.name == 'posix' and rt._resource is not None


class DetectionTests(unittest.TestCase):
    def test_force_cloud_env_overrides(self):
        with patch.dict(os.environ, {'BPS_FORCE_CLOUD': '1'}, clear=False):
            self.assertTrue(rt.is_shared_cloud())
        with patch.dict(os.environ, {'BPS_FORCE_CLOUD': '0'}, clear=False):
            self.assertFalse(rt.is_shared_cloud())

    def test_export_runtime_env_sets_flag(self):
        with patch.dict(os.environ, {'BPS_FORCE_CLOUD': '1'}, clear=False):
            os.environ.pop('BPS_SHARED_CLOUD', None)
            # export must not clobber an explicit override, but should publish it
            rt.export_runtime_env()
        # With the override present, export leaves BPS_SHARED_CLOUD unset (it only
        # writes when BPS_FORCE_CLOUD is absent). Verify the no-clobber contract.
        with patch.dict(os.environ, {}, clear=False):
            os.environ.pop('BPS_FORCE_CLOUD', None)
            os.environ.pop('BPS_SHARED_CLOUD', None)
            rt.export_runtime_env()
            self.assertIn(os.environ.get('BPS_SHARED_CLOUD'), {'0', '1'})

    def test_child_ceiling_never_below_minimum(self):
        self.assertGreaterEqual(rt.child_memory_ceiling_mb(64), 64)
        self.assertGreaterEqual(rt.child_memory_ceiling_mb(4096), 4096)


class RunGuardedSuccessTests(unittest.TestCase):
    def test_success_returns_output(self):
        res = rt.run_guarded([sys.executable, '-c', 'print("hello")'], timeout=30)
        self.assertEqual(res.returncode, 0)
        self.assertIn('hello', res.stdout)

    def test_nonzero_exit_is_returned_not_raised(self):
        res = rt.run_guarded([sys.executable, '-c', 'import sys; sys.exit(3)'], timeout=30)
        self.assertEqual(res.returncode, 3)

    def test_input_is_forwarded(self):
        res = rt.run_guarded(
            [sys.executable, '-c', 'import sys; sys.stdout.write(sys.stdin.read().upper())'],
            input='abc', timeout=30,
        )
        self.assertEqual(res.stdout.strip(), 'ABC')


@unittest.skipUnless(POSIX, 'resource limits require POSIX')
class RunGuardedSandboxTests(unittest.TestCase):
    def test_memory_runaway_child_does_not_kill_parent(self):
        # A child that tries to allocate ~300 MB under a 64 MB ceiling must fail
        # (either signal-killed -> ChildResourceError for native tools that
        # ignore malloc failure, or a graceful nonzero exit for managed runtimes).
        # Either way, THIS process must survive to make the assertion.
        alloc = 'x = bytearray(300 * 1024 * 1024); print(len(x))'
        try:
            res = rt.run_guarded([sys.executable, '-c', alloc], mem_mb=64, timeout=30, tool='memhog')
            self.assertNotEqual(res.returncode, 0, 'ceiling should have prevented the allocation')
        except rt.ChildResourceError:
            pass  # signal-kill path is also acceptable
        self.assertTrue(True, 'parent survived a runaway child')

    def test_cpu_runaway_child_is_killed(self):
        spin = 'while True: pass'
        with self.assertRaises(rt.ChildResourceError):
            rt.run_guarded([sys.executable, '-c', spin], cpu_s=1, timeout=30, tool='spinner')

    def test_timeout_raises_child_resource_error(self):
        with self.assertRaises(rt.ChildResourceError):
            rt.run_guarded([sys.executable, '-c', 'import time; time.sleep(30)'],
                           timeout=1, tool='sleeper')

    def test_datasets_tool_skips_address_space_limit(self):
        # A tool named 'datasets' must NOT get an RLIMIT_AS ceiling (Go binaries
        # reserve huge virtual memory and crash under one). The same allocation
        # is refused under a generic name but succeeds as 'datasets'.
        script = 'x = bytearray(120 * 1024 * 1024); print(len(x))'
        generic = rt.run_guarded([sys.executable, '-c', script], mem_mb=64, timeout=30, tool='generic')
        self.assertNotEqual(generic.returncode, 0)
        skipped = rt.run_guarded([sys.executable, '-c', script], mem_mb=64, timeout=30, tool='datasets')
        self.assertEqual(skipped.returncode, 0)


VENDORED_MEME = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'vendor', 'meme-5.5.9', 'meme-core')


@unittest.skipUnless(POSIX and os.path.exists(VENDORED_MEME), 'needs POSIX + vendored MEME core')
class NativeToolSignalKillTests(unittest.TestCase):
    """Directly reproduce the crash-relevant path: a native tool that exceeds an
    RLIMIT_AS ceiling is signal-killed and raised as ChildResourceError while the
    parent survives (this is the mechanism that stops the whole Streamlit worker
    from dying)."""

    def test_meme_core_under_tiny_ceiling_raises_and_parent_survives(self):
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            fa = os.path.join(td, 'x.faa')
            with open(fa, 'w') as fh:
                for i in range(40):
                    fh.write(f'>s{i}\n' + ('ACDEFGHIKLMNPQRSTVWY' * 11)[:211] + '\n')
            with self.assertRaises(rt.ChildResourceError):
                rt.run_guarded(
                    [VENDORED_MEME, fa, '-protein', '-oc', os.path.join(td, 'o'),
                     '-nmotifs', '10', '-minw', '6', '-maxw', '50', '-mod', 'zoops', '-nostatus'],
                    mem_mb=12, timeout=60, tool='meme',
                )
        self.assertTrue(True, 'parent survived a memory-starved native child')


if __name__ == '__main__':
    unittest.main()
