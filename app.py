import streamlit as st


st.set_page_config(
    page_title="BioProtein Studio",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("🧬 BioProtein Studio")
st.subheader("Professional Protein Bioinformatics Analysis Platform")

st.write(
    "Analyze protein sequences, calculate physicochemical properties, "
    "generate publication-quality figures, perform statistical analysis, "
    "and export complete research reports."
)

st.divider()

st.header("Core Capabilities")

col1, col2, col3 = st.columns(3)

with col1:
    st.subheader("🧪 Protein Analysis")
    st.write(
        "Analyze raw and multi-FASTA protein sequences with sequence "
        "validation and ambiguous-residue handling."
    )

with col2:
    st.subheader("📊 Publication Figures")
    st.write(
        "Generate bar charts, boxplots, heatmaps, scatterplots, PCA, "
        "radar plots, and vector-quality figures."
    )

with col3:
    st.subheader("📄 Research Reports")
    st.write(
        "Export physicochemical results, statistics, figure data, "
        "automatic interpretations, and Excel reports."
    )

st.divider()

st.header("Analysis Workflow")

st.markdown(
    """
1. Open **Protein Analysis** from the sidebar.
2. Paste or upload one or multiple protein sequences.
3. Select sequence-validation and ambiguous-residue options.
4. Run physicochemical and statistical analysis.
5. Generate publication-quality figures.
6. Download CSV, Excel, PNG, SVG, PDF, or TIFF outputs.
"""
)

st.divider()

st.header("Software Information")

info1, info2, info3 = st.columns(3)

info1.metric("Current Version", "4.2")
info2.metric("Analysis Engine", "Biopython")
info3.metric("Interface", "Streamlit")

st.info(
    "For publication work, use Strict Protein Mode and verify any "
    "ambiguous residues from the original sequence database."
)

st.divider()

developer_col, supervisor_col = st.columns(2)

with developer_col:
    st.subheader("👨‍💻 Developer")
    st.markdown(
        """
### Muhammad Hammad

Department of Horticulture  
University of the Punjab, Lahore, Pakistan  

**Email:** dr.hammadse@gmail.com
"""
    )

with supervisor_col:
    st.subheader("🎓 Academic Supervisor")
    st.markdown(
        """
### Dr. Muhammad Shafiq

Associate Professor  
Department of Horticulture  
University of the Punjab, Lahore, Pakistan  

**Email:** shafiq.iags@pu.edu.pk
"""
    )

st.divider()

link1, link2 = st.columns(2)

with link1:
    st.link_button(
        "🌐 Open Live App",
        "https://bioprotein-studio-hammad.streamlit.app/",
        use_container_width=True,
    )

with link2:
    st.link_button(
        "💻 View GitHub Repository",
        "https://github.com/ha7246834-cloud/BioProtein-Studio",
        use_container_width=True,
    )

st.caption(
    "BioProtein Studio v4.2 — Developed by Muhammad Hammad "
    "under the supervision of Dr. Muhammad Shafiq — © 2026"
)
