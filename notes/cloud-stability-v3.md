# Streamlit Cloud stability v3

This patch makes shared Streamlit Cloud phylogeny FastTree-only. IQ-TREE publication inference remains available in local/WSL/HPC environments and through uploaded validated Newick trees.

The change was prompted by a real 40-protein / 8,411-aa rice PR-1 run that terminated the shared Streamlit process during IQ-TREE ModelFinder + resampling.
