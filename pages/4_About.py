import streamlit as st


st.set_page_config(
    page_title="About | BioProtein Studio",
    page_icon="👨‍🔬",
    layout="wide",
)

st.title("👨‍🔬 About BioProtein Studio")

st.write(
    "BioProtein Studio is a web-based research platform for validated "
    "protein sequence analysis, physicochemical characterization, "
    "statistical evaluation, publication-quality visualization, and "
    "research report generation."
)

st.divider()

developer_col, supervisor_col = st.columns(2)

with developer_col:
    st.info("Developer")
    st.markdown(
        """
## Muhammad Hammad

Department of Horticulture  
University of the Punjab, Lahore, Pakistan  

**Email:** dr.hammadse@gmail.com  

**Role:** Developer and Researcher
"""
    )

with supervisor_col:
    st.success("Academic Supervisor")
    st.markdown(
        """
## Dr. Muhammad Shafiq

Associate Professor  
Department of Horticulture  
University of the Punjab, Lahore, Pakistan  

**Email:** shafiq.iags@pu.edu.pk  

**Role:** Academic Supervision and Research Guidance
"""
    )

st.divider()

st.header("Software Features")

st.markdown(
    """
- Raw and multi-FASTA protein input
- Sequence-quality control
- Strict protein and short-peptide validation modes
- Ambiguous-residue handling
- Molecular weight and theoretical pI
- Instability index and stability classification
- Aliphatic index and GRAVY
- Amino-acid composition
- Descriptive statistics and correlation analysis
- IQR outlier analysis
- Publication-quality figures
- PCA, radar plots, heatmaps, and scatterplots
- Automatic results interpretation
- CSV and Excel report export
"""
)

st.divider()

st.header("Software Links")

link_col1, link_col2 = st.columns(2)

with link_col1:
    st.link_button(
        "🌐 Live BioProtein Studio",
        "https://bioprotein-studio-hammad.streamlit.app/",
        use_container_width=True,
    )

with link_col2:
    st.link_button(
        "💻 GitHub Repository",
        "https://github.com/ha7246834-cloud/BioProtein-Studio",
        use_container_width=True,
    )

st.divider()

st.header("Citation")

st.markdown(
    """
**BioProtein Studio, Version 4.2**  
Developed by Muhammad Hammad  
Under the supervision of Dr. Muhammad Shafiq  
Department of Horticulture, University of the Punjab, Lahore, Pakistan  
2026
"""
)

st.caption("© 2026 Muhammad Hammad")
