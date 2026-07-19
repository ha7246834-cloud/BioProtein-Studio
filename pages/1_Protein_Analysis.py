import pandas as pd
import streamlit as st

from figures.heatmaps import (
    create_correlation_heatmap,
    create_zscore_heatmap,
)
from figures.publication_plots import (
    create_amino_acid_chart,
    create_bar_chart,
    create_boxplot_panel,
    create_gene_wise_panel,
)
from figures.scatterplots import create_scatterplot
from modules.export_tools import (
    close_figure,
    create_figure_exports,
)
from modules.interpretation import (
    generate_figure_caption,
    generate_results_summary,
)
from modules.physicochemical import analyze_proteins
from modules.sequence_parser import prepare_sequences
from modules.statistics import (
    calculate_correlation_matrix,
    calculate_summary_statistics,
    detect_outliers_iqr,
)
from reports.excel_report import create_excel_report


st.set_page_config(
    page_title="Protein Analysis | BioProtein Studio",
    page_icon="🧬",
    layout="wide",
)

st.title("🧪 Protein Analysis")
st.caption(
    "Validated physicochemical protein analysis, optional statistics, "
    "publication-quality figures, interpretation, and Excel reports."
)


SESSION_KEYS = [
    "properties_df",
    "amino_df",
    "quality_df",
    "statistics_df",
    "correlation_df",
    "outliers_df",
    "interpretation_text",
]

for key in SESSION_KEYS:
    if key not in st.session_state:
        st.session_state[key] = None


def show_figure_with_downloads(
    figure,
    filename,
    caption=None,
):
    exports = create_figure_exports(figure)

    st.image(
        exports["preview_png"],
        width="stretch",
    )

    if caption:
        st.markdown("**Suggested figure caption**")
        st.write(caption)

    cols = st.columns(3)

    with cols[0]:
        st.download_button(
            "PNG 300 DPI",
            data=exports["png_300"],
            file_name=f"{filename}_300dpi.png",
            mime="image/png",
            key=f"{filename}_png300",
        )

    with cols[1]:
        st.download_button(
            "PNG 600 DPI",
            data=exports["png_600"],
            file_name=f"{filename}_600dpi.png",
            mime="image/png",
            key=f"{filename}_png600",
        )

    with cols[2]:
        if "tiff_600" in exports:
            st.download_button(
                "TIFF 600 DPI",
                data=exports["tiff_600"],
                file_name=f"{filename}_600dpi.tiff",
                mime="image/tiff",
                key=f"{filename}_tiff600",
            )

    cols2 = st.columns(2)

    with cols2[0]:
        st.download_button(
            "SVG vector",
            data=exports["svg"],
            file_name=f"{filename}.svg",
            mime="image/svg+xml",
            key=f"{filename}_svg",
        )

    with cols2[1]:
        st.download_button(
            "PDF vector",
            data=exports["pdf"],
            file_name=f"{filename}.pdf",
            mime="application/pdf",
            key=f"{filename}_pdf",
        )

    close_figure(figure)


def reset_results():
    for key in SESSION_KEYS:
        st.session_state[key] = None


st.sidebar.header("Analysis settings")

ambiguous_mode = st.sidebar.radio(
    "Ambiguous residue handling",
    [
        "Strict mode",
        "Remove terminal X only",
        "Remove all ambiguous residues",
    ],
    index=1,
)

st.sidebar.subheader("Sequence validation")

validation_mode = st.sidebar.radio(
    "Sequence type",
    [
        "Strict protein mode",
        "Short peptide mode",
        "Custom minimum length",
    ],
    index=0,
)

if validation_mode == "Strict protein mode":
    minimum_length = 30
    allow_short_peptides = False
    st.sidebar.caption(
        "Sequences shorter than 30 amino acids will be rejected."
    )

elif validation_mode == "Short peptide mode":
    minimum_length = 3
    allow_short_peptides = True
    st.sidebar.caption(
        "Short peptides will be analyzed with a scientific warning."
    )

else:
    minimum_length = st.sidebar.number_input(
        "Minimum sequence length (aa)",
        min_value=3,
        max_value=500,
        value=30,
        step=1,
    )

    allow_short_peptides = st.sidebar.toggle(
        "Allow sequences below this length",
        value=False,
    )

st.sidebar.subheader("Output sections")

show_summary = st.sidebar.toggle(
    "Analysis summary",
    value=True,
)

show_statistics = st.sidebar.toggle(
    "Descriptive statistics",
    value=True,
)

show_outliers = st.sidebar.toggle(
    "IQR outlier analysis",
    value=False,
)

show_excel_report = st.sidebar.toggle(
    "Excel report",
    value=True,
)

show_figures = st.sidebar.toggle(
    "Publication figures",
    value=True,
)

show_amino_acids = st.sidebar.toggle(
    "Amino-acid composition",
    value=True,
)

st.sidebar.info(
    "Strict protein mode is recommended for publication work. "
    "Short-peptide mode is intended for exploratory peptide analysis."
)

if st.sidebar.button("Clear analysis"):
    reset_results()
    st.rerun()


st.header("1. Protein sequence input")

input_method = st.radio(
    "Select input method",
    [
        "Paste sequences",
        "Upload FASTA file",
    ],
    horizontal=True,
)

sequence_text = ""

if input_method == "Paste sequences":
    sequence_text = st.text_area(
        "Paste raw or FASTA protein sequences",
        height=360,
        placeholder=(
            ">Protein_1\n"
            "MALTKLFLVLLAASVSA...\n\n"
            ">Protein_2\n"
            "MKWVTFISLLFLFSSA..."
        ),
    )

else:
    uploaded_file = st.file_uploader(
        "Upload protein FASTA file",
        type=["fasta", "fa", "faa", "txt"],
    )

    if uploaded_file is not None:
        sequence_text = uploaded_file.getvalue().decode(
            "utf-8",
            errors="replace",
        )


if st.button(
    "Analyze proteins",
    type="primary",
):
    if not sequence_text.strip():
        st.error(
            "Please paste or upload at least one protein sequence."
        )

    else:
        try:
            accepted_sequences, quality_records = prepare_sequences(
                input_text=sequence_text,
                ambiguous_mode=ambiguous_mode,
                minimum_length=int(minimum_length),
                allow_short_peptides=allow_short_peptides,
            )

            quality_df = pd.DataFrame(quality_records)
            st.session_state.quality_df = quality_df

            if not accepted_sequences:
                st.session_state.properties_df = None
                st.session_state.amino_df = None
                st.session_state.statistics_df = None
                st.session_state.correlation_df = None
                st.session_state.outliers_df = None
                st.session_state.interpretation_text = None

                st.error(
                    "No sequence could be analyzed. "
                    "Review the quality-control table."
                )

            else:
                properties_df, amino_df = analyze_proteins(
                    accepted_sequences
                )

                statistics_df = calculate_summary_statistics(
                    properties_df
                )

                correlation_df = calculate_correlation_matrix(
                    properties_df
                )

                outliers_df = detect_outliers_iqr(
                    properties_df
                )

                interpretation_text = generate_results_summary(
                    properties_df
                )

                st.session_state.properties_df = properties_df
                st.session_state.amino_df = amino_df
                st.session_state.statistics_df = statistics_df
                st.session_state.correlation_df = correlation_df
                st.session_state.outliers_df = outliers_df
                st.session_state.interpretation_text = (
                    interpretation_text
                )

                st.success(
                    f"{len(properties_df)} protein sequence(s) "
                    "analyzed successfully."
                )

        except Exception as error:
            st.exception(error)


properties_df = st.session_state.properties_df
amino_df = st.session_state.amino_df
quality_df = st.session_state.quality_df
statistics_df = st.session_state.statistics_df
correlation_df = st.session_state.correlation_df
outliers_df = st.session_state.outliers_df
interpretation_text = st.session_state.interpretation_text


if quality_df is not None:
    st.header("2. Sequence quality control")

    st.dataframe(
        quality_df,
        width="stretch",
        hide_index=True,
    )

    st.download_button(
        "Download quality-control CSV",
        data=quality_df.to_csv(index=False).encode("utf-8"),
        file_name="sequence_quality_control.csv",
        mime="text/csv",
    )


if properties_df is not None and not properties_df.empty:
    st.header("3. Physicochemical properties")

    st.dataframe(
        properties_df,
        width="stretch",
        hide_index=True,
    )

    st.download_button(
        "Download physicochemical CSV",
        data=properties_df.to_csv(index=False).encode("utf-8"),
        file_name="physicochemical_properties.csv",
        mime="text/csv",
    )

    if show_summary:
        st.header("4. Analysis summary")

        stable_count = int(
            (properties_df["Stability"] == "Stable").sum()
        )

        unstable_count = int(
            (properties_df["Stability"] == "Unstable").sum()
        )

        hydrophilic_count = int(
            (
                properties_df["Hydropathy_Class"]
                == "Hydrophilic"
            ).sum()
        )

        hydrophobic_count = int(
            (
                properties_df["Hydropathy_Class"]
                == "Hydrophobic"
            ).sum()
        )

        metric1, metric2, metric3, metric4, metric5 = (
            st.columns(5)
        )

        metric1.metric("Proteins", len(properties_df))
        metric2.metric("Stable", stable_count)
        metric3.metric("Unstable", unstable_count)
        metric4.metric("Hydrophilic", hydrophilic_count)
        metric5.metric("Hydrophobic", hydrophobic_count)

        st.subheader("Automatic results paragraph")
        st.write(interpretation_text)

        st.download_button(
            "Download results paragraph",
            data=interpretation_text.encode("utf-8"),
            file_name="physicochemical_results_paragraph.txt",
            mime="text/plain",
        )

    if show_statistics:
        st.header("5. Descriptive statistics")

        if len(properties_df) < 2:
            st.info(
                "Only one sequence was analyzed. Comparative "
                "interpretation of descriptive statistics is limited."
            )

        st.dataframe(
            statistics_df,
            width="stretch",
            hide_index=True,
        )

        st.download_button(
            "Download statistics CSV",
            data=statistics_df.to_csv(
                index=False
            ).encode("utf-8"),
            file_name="summary_statistics.csv",
            mime="text/csv",
        )

    if show_outliers:
        st.header("6. IQR outlier analysis")

        if len(properties_df) < 4:
            st.info(
                "IQR outlier analysis requires at least "
                "four protein sequences."
            )

        elif outliers_df is not None and not outliers_df.empty:
            st.warning(
                f"{len(outliers_df)} potential outlier "
                "observation(s) detected using the 1.5 × IQR method."
            )

            st.dataframe(
                outliers_df,
                width="stretch",
                hide_index=True,
            )

        else:
            st.success(
                "No outliers were detected using the 1.5 × IQR method."
            )

    if show_excel_report:
        st.header("7. Complete Excel report")

        excel_bytes = create_excel_report(
            properties_df=properties_df,
            amino_acid_df=amino_df,
            quality_df=quality_df,
            statistics_df=statistics_df,
            correlation_df=correlation_df,
            outliers_df=outliers_df,
            interpretation_text=interpretation_text,
        )

        st.download_button(
            "Download complete Excel workbook",
            data=excel_bytes,
            file_name="BioProtein_Studio_Report.xlsx",
            mime=(
                "application/vnd.openxmlformats-officedocument."
                "spreadsheetml.sheet"
            ),
        )

    if show_figures:
        st.header("8. Publication figure studio")

        if len(properties_df) == 1:
            st.warning(
                "Only one sequence was analyzed. Comparative "
                "boxplots, correlation heatmaps and scatterplots "
                "are disabled."
            )

            figure_options = [
                "Protein length bar chart",
                "Molecular weight bar chart",
                "Theoretical pI bar chart",
                "Instability index bar chart",
                "Aliphatic index bar chart",
                "GRAVY bar chart",
            ]

        else:
            figure_options = [
                "Paper-style four-panel boxplot",
                "Gene-wise four-panel comparison",
                "Protein length bar chart",
                "Molecular weight bar chart",
                "Theoretical pI bar chart",
                "Instability index bar chart",
                "Aliphatic index bar chart",
                "GRAVY bar chart",
                "Normalized property heatmap",
                "Correlation heatmap",
                "Scatterplot studio",
            ]

        figure_option = st.selectbox(
            "Select figure type",
            figure_options,
        )

        protein_count = len(properties_df)

        if figure_option == "Paper-style four-panel boxplot":
            figure = create_boxplot_panel(properties_df)

            show_figure_with_downloads(
                figure,
                "paper_style_physicochemical_boxplot",
                generate_figure_caption(
                    "boxplot",
                    protein_count,
                ),
            )

        elif figure_option == "Gene-wise four-panel comparison":
            figure = create_gene_wise_panel(properties_df)

            show_figure_with_downloads(
                figure,
                "gene_wise_physicochemical_comparison",
                generate_figure_caption(
                    "gene_wise",
                    protein_count,
                ),
            )

        elif figure_option == "Protein length bar chart":
            figure = create_bar_chart(
                properties_df,
                "Length_aa",
                "Protein length comparison",
                "Protein length (amino acids)",
                decimals=0,
            )

            show_figure_with_downloads(
                figure,
                "protein_length",
            )

        elif figure_option == "Molecular weight bar chart":
            figure = create_bar_chart(
                properties_df,
                "Molecular_Weight_kDa",
                "Molecular weight comparison",
                "Molecular weight (kDa)",
                decimals=2,
            )

            show_figure_with_downloads(
                figure,
                "molecular_weight",
            )

        elif figure_option == "Theoretical pI bar chart":
            figure = create_bar_chart(
                properties_df,
                "Theoretical_pI",
                "Theoretical isoelectric point",
                "Theoretical pI",
                decimals=2,
                reference_line=7,
                reference_label="Neutral pH (7)",
            )

            show_figure_with_downloads(
                figure,
                "theoretical_pi",
            )

        elif figure_option == "Instability index bar chart":
            figure = create_bar_chart(
                properties_df,
                "Instability_Index",
                "Protein instability index",
                "Instability index",
                decimals=2,
                reference_line=40,
                reference_label="Stability threshold (40)",
            )

            show_figure_with_downloads(
                figure,
                "instability_index",
                generate_figure_caption(
                    "instability",
                    protein_count,
                ),
            )

        elif figure_option == "Aliphatic index bar chart":
            figure = create_bar_chart(
                properties_df,
                "Aliphatic_Index",
                "Protein aliphatic index",
                "Aliphatic index",
                decimals=2,
            )

            show_figure_with_downloads(
                figure,
                "aliphatic_index",
            )

        elif figure_option == "GRAVY bar chart":
            figure = create_bar_chart(
                properties_df,
                "GRAVY",
                "Grand average of hydropathicity",
                "GRAVY score",
                decimals=3,
                reference_line=0,
                reference_label="Neutral hydropathy",
            )

            show_figure_with_downloads(
                figure,
                "gravy_score",
                generate_figure_caption(
                    "gravy",
                    protein_count,
                ),
            )

        elif figure_option == "Normalized property heatmap":
            figure = create_zscore_heatmap(properties_df)

            show_figure_with_downloads(
                figure,
                "normalized_physicochemical_heatmap",
                generate_figure_caption(
                    "heatmap",
                    protein_count,
                ),
            )

        elif figure_option == "Correlation heatmap":
            figure = create_correlation_heatmap(
                correlation_df
            )

            show_figure_with_downloads(
                figure,
                "physicochemical_correlation_heatmap",
            )

        elif figure_option == "Scatterplot studio":
            numeric_columns = [
                "Length_aa",
                "Molecular_Weight_kDa",
                "Theoretical_pI",
                "Instability_Index",
                "Aliphatic_Index",
                "GRAVY",
                "Aromaticity",
                "Positive_Residues_KR",
                "Negative_Residues_DE",
                "Charge_Difference",
            ]

            scatter_col1, scatter_col2 = st.columns(2)

            with scatter_col1:
                x_column = st.selectbox(
                    "X-axis property",
                    numeric_columns,
                    index=0,
                )

            with scatter_col2:
                y_column = st.selectbox(
                    "Y-axis property",
                    numeric_columns,
                    index=1,
                )

            if x_column == y_column:
                st.warning(
                    "Select two different properties."
                )

            else:
                figure = create_scatterplot(
                    properties_df,
                    x_column,
                    y_column,
                )

                show_figure_with_downloads(
                    figure,
                    f"scatter_{y_column}_vs_{x_column}",
                )

    if show_amino_acids:
        st.header("9. Amino-acid composition")

        selected_protein = st.selectbox(
            "Select protein",
            properties_df["Protein_ID"].tolist(),
        )

        selected_amino_df = amino_df[
            amino_df["Protein_ID"] == selected_protein
        ]

        st.dataframe(
            selected_amino_df,
            width="stretch",
            hide_index=True,
        )

        amino_figure = create_amino_acid_chart(
            amino_df,
            selected_protein,
        )

        show_figure_with_downloads(
            amino_figure,
            f"{selected_protein}_amino_acid_composition",
        )

        st.download_button(
            "Download amino-acid composition CSV",
            data=amino_df.to_csv(
                index=False
            ).encode("utf-8"),
            file_name="amino_acid_composition.csv",
            mime="text/csv",
        )


elif quality_df is not None:
    st.error(
        "No valid protein result is available. "
        "Review the quality-control section."
    )
