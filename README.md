# BioProtein Studio v5.5 Research Release

BioProtein Studio is a Streamlit research application for protein characterization and generic multi-gene / gene-family analysis.

## Main workflows

### Protein Analysis
- FASTA validation and sequence QC
- physicochemical characterization
- statistical summaries and publication outputs

### Gene Structure, Domains, Motifs & Phylogeny
The integrated module combines workflows commonly carried out separately with GSDS, NCBI CD-Search, MEME Suite, phylogeny software and TBtools.

Supported inputs and routes:
- Protein FASTA or CDS FASTA
- matching CDS + genomic FASTA via EMBOSS `est2genome`
- GFF3/GTF or structure-table import
- NCBI protein accession resolution
- **Protein + Species automatic gene structure** using NCBI Datasets + miniprot + NCBI GFF reconciliation
- optional supplied Newick tree
- automatic phylogeny

## Automatic gene structure

Protein sequence alone does not contain intron boundaries, so BioProtein Studio never invents gene structure directly from protein sequence.

When the user supplies Protein FASTA plus a species/taxon, the application can:
1. resolve an annotated NCBI reference assembly;
2. cache the reference genome locally;
3. map proteins to the genome with splice-aware `miniprot`;
4. reconcile mapped loci against the NCBI GFF3 annotation;
5. reconstruct CDS/genomic sequences and exon coordinates;
6. translate the CDS back to protein for validation;
7. flag ambiguous or low-confidence mappings as **REVIEW**.

A conservative canonical one-intron rescue can recover a strongly supported missed splice, but computationally rescued structures remain **REVIEW** until annotation/manual confirmation.

## Conserved domains and motifs

- Domains use the real **NCBI Batch CD-Search** service.
- Specific CDD hits are prioritized for publication display while full raw results are retained.
- Motifs use the real **MEME Suite** executable in protein mode.
- Default publication display shows PASS motifs; all motifs remain available for inspection.

## Phylogeny

Three automatic modes are available:
- **Publication:** MAFFT + IQ-TREE + ModelFinder + SH-aLRT + ultrafast bootstrap.
- **Fast screening:** MAFFT + FastTree.
- **Fallback:** internal Neighbor-Joining for small families.

Uploaded Newick trees are also supported and take priority when supplied.

## Integrated publication figure

The final family architecture can contain four aligned panels:

**Phylogenetic Tree | Gene Structure | Conserved Motifs | Conserved Domains**

Each biological panel keeps its appropriate coordinate system: gene structure in bp and protein motifs/domains in amino-acid coordinates.

## Graph Studio

v5.5 adds interactive figure optimization without changing biological results:
- one-click publication style variants;
- Journal Classic, Colorblind Safe, Plant & Earth, Soft Pastel and High Contrast presets;
- manual colors for tree/support labels, validated exons, REVIEW exons, introns, motif classes and domain classes;
- title size, gene-label size and line-width controls;
- live previews;
- customized PNG, SVG, PDF and TIFF export.

## Reproducibility package

Analysis packages can contain:
- input FASTA files;
- sequence and translation QC;
- phylogeny tree, alignment, command/report/log files;
- gene-structure coordinates and QC;
- automatic reference metadata, mapping QC and annotation reconciliation;
- auto-generated CDS and genomic FASTA;
- raw miniprot GFF3;
- raw/full CDD results and Search-ID;
- MEME XML and complete MEME output;
- motif/domain validation tables;
- 600-DPI PNG/TIFF and SVG/PDF figures;
- run manifest and methods/QC notes.

## Full installation — Linux / WSL recommended

```bash
conda env create -f environment-gdm.yml
conda activate bioprotein-studio
streamlit run app.py
```

The supplied environment includes:
- MEME Suite
- EMBOSS
- MAFFT
- FastTree
- IQ-TREE
- miniprot
- NCBI Datasets CLI

A Python-only `pip install -r requirements.txt` can run the interface and Python components, but external scientific executables are required for the full workflow.

## Scientific principle

BioProtein Studio separates **computed evidence**, **reference evidence**, and **REVIEW** results. Visual customization never changes biological coordinates, domain hits, motif calls, mappings or phylogenetic inference.
