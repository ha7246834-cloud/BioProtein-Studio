import pandas as pd


def get_extreme_record(
    dataframe: pd.DataFrame,
    column: str,
    mode: str,
) -> pd.Series:
    """Return minimum or maximum row for a property."""
    if mode == "minimum":
        index = dataframe[column].idxmin()
    else:
        index = dataframe[column].idxmax()

    return dataframe.loc[index]


def generate_results_summary(
    properties_df: pd.DataFrame,
) -> str:
    """Generate a manuscript-style physicochemical results paragraph."""
    if properties_df is None or properties_df.empty:
        return "No physicochemical results are available."

    number_of_proteins = len(properties_df)

    shortest = get_extreme_record(
        properties_df,
        "Length_aa",
        "minimum",
    )

    longest = get_extreme_record(
        properties_df,
        "Length_aa",
        "maximum",
    )

    lowest_mw = get_extreme_record(
        properties_df,
        "Molecular_Weight_kDa",
        "minimum",
    )

    highest_mw = get_extreme_record(
        properties_df,
        "Molecular_Weight_kDa",
        "maximum",
    )

    lowest_pi = get_extreme_record(
        properties_df,
        "Theoretical_pI",
        "minimum",
    )

    highest_pi = get_extreme_record(
        properties_df,
        "Theoretical_pI",
        "maximum",
    )

    lowest_gravy = get_extreme_record(
        properties_df,
        "GRAVY",
        "minimum",
    )

    highest_gravy = get_extreme_record(
        properties_df,
        "GRAVY",
        "maximum",
    )

    stable_count = int(
        (
            properties_df["Stability"]
            == "Stable"
        ).sum()
    )

    unstable_count = int(
        (
            properties_df["Stability"]
            == "Unstable"
        ).sum()
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

    paragraph = (
        f"A total of {number_of_proteins} protein sequences were "
        f"evaluated for their physicochemical characteristics. "
        f"Protein length ranged from "
        f"{int(shortest['Length_aa'])} amino acids in "
        f"{shortest['Protein_ID']} to "
        f"{int(longest['Length_aa'])} amino acids in "
        f"{longest['Protein_ID']}. Molecular weight varied from "
        f"{lowest_mw['Molecular_Weight_kDa']:.3f} kDa "
        f"({lowest_mw['Protein_ID']}) to "
        f"{highest_mw['Molecular_Weight_kDa']:.3f} kDa "
        f"({highest_mw['Protein_ID']}). The theoretical "
        f"isoelectric point ranged from "
        f"{lowest_pi['Theoretical_pI']:.2f} in "
        f"{lowest_pi['Protein_ID']} to "
        f"{highest_pi['Theoretical_pI']:.2f} in "
        f"{highest_pi['Protein_ID']}. Based on the instability-index "
        f"threshold of 40, {stable_count} proteins were classified "
        f"as stable and {unstable_count} as unstable. "
        f"The GRAVY values ranged from "
        f"{lowest_gravy['GRAVY']:.3f} "
        f"({lowest_gravy['Protein_ID']}) to "
        f"{highest_gravy['GRAVY']:.3f} "
        f"({highest_gravy['Protein_ID']}). "
        f"{hydrophilic_count} proteins displayed negative GRAVY "
        f"values and were therefore classified as hydrophilic, "
        f"whereas {hydrophobic_count} proteins were classified "
        f"as hydrophobic."
    )

    return paragraph


def generate_figure_caption(
    figure_type: str,
    protein_count: int,
) -> str:
    """Generate a journal-style figure caption."""
    captions = {
        "boxplot": (
            f"Distribution of physicochemical properties across "
            f"{protein_count} proteins. The panels show GRAVY, "
            f"aliphatic index, molecular weight, and theoretical "
            f"isoelectric point. Boxes represent the interquartile "
            f"range, central lines indicate medians, and individual "
            f"points represent protein values."
        ),
        "gene_wise": (
            f"Gene-wise comparison of physicochemical properties "
            f"across {protein_count} proteins, including sequence "
            f"length, molecular weight, theoretical pI, and GRAVY."
        ),
        "heatmap": (
            f"Z-score normalized heatmap of physicochemical "
            f"properties across {protein_count} proteins. Positive "
            f"and negative Z-scores indicate values above and below "
            f"the property mean, respectively."
        ),
        "instability": (
            f"Instability-index distribution across {protein_count} "
            f"proteins. The dashed horizontal line at 40 represents "
            f"the conventional threshold separating stable and "
            f"unstable proteins."
        ),
        "gravy": (
            f"GRAVY-score comparison across {protein_count} proteins. "
            f"Negative values indicate hydrophilic tendency, whereas "
            f"positive values indicate hydrophobic tendency."
        ),
    }

    return captions.get(
        figure_type,
        (
            f"Physicochemical comparison across "
            f"{protein_count} proteins."
        ),
    )