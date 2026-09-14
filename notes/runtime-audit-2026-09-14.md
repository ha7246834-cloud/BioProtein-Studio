# Runtime audit — 2026-09-14

Observed live failure: shared Streamlit Cloud process terminated during automatic IQ-TREE on a 40-protein / 8,411-aa rice PR-1 family.

Audit conclusions:
- Streamlit Cloud is Python 3.14.x in current deployment logs.
- Previous CI only exercised Python 3.12.
- packages.txt installed IQ-TREE plus a large OpenMPI/ROCm dependency chain on every cold build.
- The live failure occurred during the phylogeny stage, before CDD/MEME/gene-structure stages.
- Shared Cloud should therefore avoid IQ-TREE entirely and use MAFFT + FastTree; local/WSL/HPC can retain IQ-TREE publication inference.
- New CI now runs Streamlit smoke and HTTP health checks under Python 3.12 and 3.14.
