from __future__ import annotations

import os
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request
from pathlib import Path


BIN_DIR = Path.home() / '.cache' / 'bioprotein-studio' / 'bin'
SRC_DIR = Path.home() / '.cache' / 'bioprotein-studio' / 'src'


def _prepend_bin():
    BIN_DIR.mkdir(parents=True, exist_ok=True)
    current = os.environ.get('PATH', '')
    b = str(BIN_DIR)
    if b not in current.split(os.pathsep):
        os.environ['PATH'] = b + os.pathsep + current


def _download(url: str, target: Path, timeout: int = 180):
    target.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={'User-Agent': 'BioProtein-Studio/5.5'})
    with urllib.request.urlopen(req, timeout=timeout) as r, open(target, 'wb') as out:
        shutil.copyfileobj(r, out)


def ensure_ncbi_datasets() -> str:
    _prepend_bin()
    found = shutil.which('datasets')
    if found:
        return found
    target = BIN_DIR / 'datasets'
    if not target.exists():
        _download(
            'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets',
            target,
            timeout=240,
        )
        target.chmod(0o755)
    return str(target) if target.exists() else ''


def ensure_miniprot() -> str:
    _prepend_bin()
    found = shutil.which('miniprot')
    if found:
        return found
    target = BIN_DIR / 'miniprot'
    if target.exists():
        target.chmod(0o755)
        return str(target)

    SRC_DIR.mkdir(parents=True, exist_ok=True)
    src = SRC_DIR / 'miniprot-0.18'
    if not src.exists():
        with tempfile.TemporaryDirectory(prefix='bps_miniprot_') as td:
            tgz = Path(td) / 'miniprot-v0.18.tar.gz'
            _download('https://github.com/lh3/miniprot/archive/refs/tags/v0.18.tar.gz', tgz, timeout=240)
            with tarfile.open(tgz, 'r:gz') as tf:
                tf.extractall(SRC_DIR)
    if not src.exists():
        candidates = sorted(SRC_DIR.glob('miniprot-*'))
        if candidates:
            src = candidates[-1]
    if not src.exists():
        return ''

    p = subprocess.run(['make', '-j2'], cwd=src, capture_output=True, text=True, timeout=300)
    built = src / 'miniprot'
    if p.returncode != 0 or not built.exists():
        return ''
    shutil.copy2(built, target)
    target.chmod(0o755)
    return str(target)


def bootstrap_cloud_tools() -> dict:
    """Best-effort bootstrap for binaries absent on Streamlit Community Cloud.

    Existing system/Conda installations always take priority. Failures are intentionally
    non-fatal so manual/reference routes remain usable.
    """
    _prepend_bin()
    status = {}
    try:
        status['datasets'] = ensure_ncbi_datasets()
    except Exception as e:
        status['datasets'] = ''
        status['datasets_error'] = str(e)
    try:
        status['miniprot'] = ensure_miniprot()
    except Exception as e:
        status['miniprot'] = ''
        status['miniprot_error'] = str(e)
    return status
