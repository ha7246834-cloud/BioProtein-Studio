import streamlit as st

st.set_page_config(
    page_title="BioProtein Studio",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("🧬 BioProtein Studio")
st.subheader("Integrated Protein & Gene-Family Bioinformatics Analysis Platform")

st.write(
    "Analyze protein sequences, calculate physicochemical properties, run gene-structure, "
    "conserved-domain and MEME motif workflows, generate publication-quality figures, "
    "apply transparent scientific QC, and export reproducible research packages."
)

st.divider()
st.header("Core Capabilities")

col1, col2, col3 = st.columns(3)
with col1:
    st.subheader("🧪 Protein Analysis")
    st.write("Validate raw or multi-FASTA proteins and calculate physicochemical properties with statistical summaries and publication outputs.")
with col2:
    st.subheader("🧬 Gene Structure, Domains & Motifs")
    st.write("GSDS-like CDS-to-genome structure mapping, NCBI Batch CDD conserved-domain search, MEME motif discovery, tree-ordered visualization, and scientific validation.")
with col3:
    st.subheader("📊 Reproducible Research Outputs")
    st.write("Export CSV/Excel data, raw engine outputs, validation tables, run manifests, and 600-DPI PNG/TIFF plus SVG/PDF figures.")

st.divider()
st.header("Recommended Workflow")
st.markdown("""
1. Open **Protein Analysis** for sequence quality and physicochemical characterization.
2. Open **Gene Structure, Domains & Motifs** for the genome-wide family workflow.
3. Provide protein FASTA; optionally provide matching CDS + genomic FASTA for gene structure.
4. If FASTA IDs are real NCBI protein accessions, the module can attempt reference CDS/genomic retrieval.
5. Run **NCBI CDD** and **MEME** from the same workspace.
6. Review every **PASS / REVIEW** QC table before interpreting the biology.
7. Download the complete package containing raw outputs, tables, parameters and publication figures.
""")

st.warning("A protein sequence alone cannot scientifically determine exon–intron structure. BioProtein Studio therefore requires reference genomic/CDS information or a resolvable NCBI accession; it never invents a gene structure from protein sequence alone.")

st.divider()
st.header("Software Information")
info1, info2, info3 = st.columns(3)
info1.metric("Current Version", "5.0")
info2.metric("Core Engines", "Biopython • CDD • MEME • EMBOSS")
info3.metric("Interface", "Streamlit")

st.info("Full automatic motif + spliced gene-structure analysis requires MEME Suite and EMBOSS. Use the supplied Conda environment on Linux/WSL for the complete local workflow.")

st.divider()
developer_col, supervisor_col = st.columns(2)
with developer_col:
    st.subheader("👨‍💻 Developer")
    st.markdown("""
### Muhammad Hammad

Department of Horticulture  
University of the Punjab, Lahore, Pakistan  

**Email:** dr.hammadse@gmail.com
""")
with supervisor_col:
    st.subheader("🎓 Academic Supervisor")
    st.markdown("""
### Dr. Muhammad Shafiq

Associate Professor  
Department of Horticulture  
University of the Punjab, Lahore, Pakistan  

**Email:** shafiq.iags@pu.edu.pk
""")

st.divider()
link1, link2 = st.columns(2)
with link1:
    st.link_button("🌐 Open Live App", "https://bioprotein-studio-hammad.streamlit.app/", use_container_width=True)
with link2:
    st.link_button("💻 View GitHub Repository", "https://github.com/ha7246834-cloud/BioProtein-Studio", use_container_width=True)

st.caption("BioProtein Studio v5.0 — Developed by Muhammad Hammad under the supervision of Dr. Muhammad Shafiq — © 2026")
