# BioProtein Studio v5.0

BioProtein Studio is a Streamlit research application for protein characterization and genome-wide gene-family analysis.

## New in v5.0 — Gene Structure, Conserved Domains & Motifs

The new module integrates the workflow often performed manually with GSDS, NCBI CD-Search, MEME Suite and TBtools.

### One workspace

- Protein FASTA or CDS FASTA input
- Protein sequence-integrity QC
- CDS translation QC
- Gene structure from matching CDS + genomic FASTA using **EMBOSS est2genome**, or an NCBI reference CDS feature when a protein accession can be resolved
- Live **NCBI Batch CD-Search** conserved-domain analysis
- Local **MEME Suite** protein motif discovery
- Optional Newick-tree row ordering
- Individual and integrated architecture figures
- 600-DPI PNG/TIFF and vector SVG/PDF output
- Raw tool outputs, cleaned tables, QC tables and a run manifest in one ZIP package

## Scientific safety design

BioProtein Studio does **not** infer exon–intron structure from protein sequence alone. A real genomic/CDS reference is required. If neither matching CDS/genomic sequences nor a resolvable NCBI accession is available, gene structure is reported as unavailable instead of being fabricated.

Validation is rule-based and auditable:

- Gene-structure spliced alignment: CDS coverage >=95% and weighted exon identity >=95% for PASS; splice information and reading frame are retained for review.
- NCBI CDD: specific hits are prioritized; E-values and bit scores are retained; overlapping redundant footprints are reduced while raw full output is preserved.
- MEME: motif E-value is retained; the module additionally reports family prevalence and domain overlap. Default family model is ZOOPS.
- Duplicate FASTA IDs are rejected rather than silently overwritten.
- Raw CDD, MEME and est2genome outputs are included in the result package for reproducibility.

## Full installation — Linux / WSL recommended

Bioconda distributes MEME Suite and EMBOSS for Linux/macOS. On Windows, run the full environment inside WSL2.

```bash
conda env create -f environment-gdm.yml
conda activate bioprotein-studio
streamlit run app.py
```

## Python-only installation

This installs the Streamlit application and remote CDD support, but automatic MEME and spliced est2genome steps remain unavailable until their executables are installed.

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Inputs for the one-click module

### Best case: protein FASTA + matching CDS + genomic FASTA

All three should use the same FASTA identifiers. The application validates sequence data, maps CDS to genomic DNA, runs CDD and MEME, applies QC and creates final figures.

### Protein FASTA with real NCBI protein accession IDs

If CDS/genomic FASTA is not supplied, the module can attempt to retrieve a linked NCBI nucleotide record and only accepts a CDS feature that matches the accession/translation. An Entrez email is required by the interface.

### Protein FASTA with custom IDs only

CDD and MEME can run automatically. Gene structure cannot be reconstructed scientifically from protein sequence alone, so matching CDS + genomic sequence or a valid reference accession is required.

## Output package

The one-click download contains input FASTA files, sequence QC, gene structure coordinates and QC, raw est2genome text, raw and non-redundant CDD hits, CDD Search-ID, MEME XML and original result ZIP, motif coordinates and QC, publication figures, `run_manifest.json`, and methods/QC notes.

## Existing protein module

The original physicochemical protein-analysis workflow remains available in `pages/1_Protein_Analysis.py`.
