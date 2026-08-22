"""BioProtein Studio module package.

Prefer normal system/Conda executables. If a tool is unavailable there,
allow vetted Linux binaries bundled under vendor/bin as a cloud fallback.
"""
from pathlib import Path
import os

_vendor_bin = Path(__file__).resolve().parent.parent / "vendor" / "bin"

if _vendor_bin.is_dir():
    current = os.environ.get("PATH", "")
    entries = current.split(os.pathsep) if current else []
    vendor = str(_vendor_bin)
    if vendor not in entries:
        os.environ["PATH"] = vendor + (os.pathsep + current if current else "")
