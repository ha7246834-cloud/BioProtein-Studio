import streamlit as st


st.set_page_config(
    page_title="BioProtein Studio",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


def load_css():
    try:
        with open(
            "assets/style.css",
            "r",
            encoding="utf-8",
        ) as css_file:
            st.markdown(
                f"<style>{css_file.read()}</style>",
                unsafe_allow_html=True,
            )
    except FileNotFoundError:
        pass


load_css()


st.markdown(
    """
    <div class="hero-card">
        <div class="hero-icon">🧬</div>
        <h1>BioProtein Studio</h1>
        <h3>Professional Protein Bioinformatics Analysis Platform</h3>
        <p>
            Analyze protein sequences, calculate physicochemical
            properties, generate publication-quality figures,
            perform statistical analysis, and export complete reports.
        </p>
    </div>
    """,
    unsafe_allow_html=True,
)


st.markdown("## Core Capabilities")

col1, col2, col3 = st.columns(3)

with col1:
    st.markdown(
        """
        <div class="feature-card">
            <h3>🧪 Protein Analysis</h3>
            <p>
                Analyze raw and multi-FASTA protein sequences,
                including sequence validation and ambiguous-residue
                handling.
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

with col2:
    st.markdown(
        """
        <div class="feature-card">
            <h3>📊 Publication Figures</h3>
            <p>
                Generate bar plots, boxplots, heatmaps,
                scatterplots, and vector-quality figures for
                scientific publications.
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

with col3:
    st.markdown(
        """
        <div class="feature-card">
            <h3>📄 Research Reports</h3>
            <p>
                Export CSV files, Excel workbooks, statistical
                summaries, automatic results paragraphs, and
                manuscript-ready captions.
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )


st.markdown("## Analysis Workflow")

st.markdown(
    """
    1. Open **Protein Analysis** from the sidebar.
    2. Paste or upload one or multiple protein sequences.
    3. Select ambiguous-residue handling mode.
    4. Run physicochemical and statistical analysis.
    5. Open **Publication Figures** to inspect figures.
    6. Download reports and publication-ready outputs.
    """
)


st.markdown("## Software Information")

info1, info2, info3 = st.columns(3)

info1.metric(
    "Current Version",
    "4.0 Development",
)

info2.metric(
    "Analysis Engine",
    "Biopython",
)

info3.metric(
    "Interface",
    "Streamlit",
)


st.info(
    "BioProtein Studio is currently under active development. "
    "Results involving removed ambiguous amino-acid residues "
    "should be interpreted as approximate."
)


st.markdown(
    """
    <div class="developer-card">
        <h3>Developed by Muhammad Hammad</h3>
        <p>
            Department of Horticulture<br>
            University of the Punjab, Lahore<br>
            © 2026 Muhammad Hammad
        </p>
    </div>
    """,
    unsafe_allow_html=True,
)