import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PCA_COLUMNS = [
    "Length_aa",
    "Molecular_Weight_kDa",
    "Theoretical_pI",
    "Instability_Index",
    "Aliphatic_Index",
    "GRAVY",
    "Aromaticity",
    "Charge_Difference",
]


def _standardize(dataframe: pd.DataFrame) -> pd.DataFrame:
    numeric = dataframe.astype(float)
    standard_deviation = numeric.std(axis=0, ddof=0).replace(0, 1)
    return (numeric - numeric.mean(axis=0)) / standard_deviation


def create_pca_plot(
    properties_df: pd.DataFrame,
):
    """Create a two-dimensional PCA plot using NumPy SVD."""
    available_columns = [
        column
        for column in PCA_COLUMNS
        if column in properties_df.columns
    ]

    matrix = properties_df[available_columns].apply(
        pd.to_numeric,
        errors="coerce",
    )

    matrix = matrix.fillna(matrix.mean())
    standardized = _standardize(matrix)

    u_matrix, singular_values, _ = np.linalg.svd(
        standardized.values,
        full_matrices=False,
    )

    scores = u_matrix * singular_values

    variance = singular_values ** 2
    total_variance = variance.sum()

    if total_variance == 0:
        explained = np.array([0.0, 0.0])
    else:
        explained = variance / total_variance * 100

    if scores.shape[1] < 2:
        scores = np.column_stack(
            [scores[:, 0], np.zeros(len(scores))]
        )
        explained = np.array(
            [explained[0] if len(explained) else 0.0, 0.0]
        )

    figure, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(
        scores[:, 0],
        scores[:, 1],
        s=75,
        edgecolor="black",
        linewidth=0.7,
    )

    for protein_id, x_value, y_value in zip(
        properties_df["Protein_ID"],
        scores[:, 0],
        scores[:, 1],
    ):
        ax.annotate(
            protein_id,
            (x_value, y_value),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )

    ax.axhline(0, linewidth=0.7, linestyle="--")
    ax.axvline(0, linewidth=0.7, linestyle="--")

    ax.set_xlabel(
        f"PC1 ({explained[0]:.1f}%)",
        fontweight="bold",
    )

    ax.set_ylabel(
        f"PC2 ({explained[1]:.1f}%)",
        fontweight="bold",
    )

    ax.set_title(
        "PCA of physicochemical properties",
        fontsize=13,
        fontweight="bold",
        pad=12,
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(linestyle="--", linewidth=0.5, alpha=0.3)

    figure.tight_layout()

    return figure, explained[:2]


def create_radar_plot(
    properties_df: pd.DataFrame,
    selected_proteins: list,
    property_columns: list,
):
    """Create a min-max normalized radar plot."""
    selected = properties_df[
        properties_df["Protein_ID"].isin(selected_proteins)
    ].copy()

    matrix = properties_df[property_columns].apply(
        pd.to_numeric,
        errors="coerce",
    )

    minimum = matrix.min(axis=0)
    maximum = matrix.max(axis=0)
    denominator = (maximum - minimum).replace(0, 1)

    normalized_all = (matrix - minimum) / denominator
    normalized_selected = normalized_all.loc[selected.index]

    labels = [
        column.replace("_", " ")
        for column in property_columns
    ]

    angles = np.linspace(
        0,
        2 * np.pi,
        len(labels),
        endpoint=False,
    ).tolist()

    angles += angles[:1]

    figure, ax = plt.subplots(
        figsize=(8, 7),
        subplot_kw={"polar": True},
    )

    for row_index, protein_id in zip(
        normalized_selected.index,
        selected["Protein_ID"],
    ):
        values = normalized_all.loc[
            row_index,
            property_columns,
        ].tolist()

        values += values[:1]

        ax.plot(
            angles,
            values,
            linewidth=1.5,
            label=protein_id,
        )

        ax.fill(
            angles,
            values,
            alpha=0.08,
        )

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylim(0, 1)

    ax.set_title(
        "Normalized physicochemical radar plot",
        fontsize=13,
        fontweight="bold",
        pad=20,
    )

    ax.legend(
        loc="upper right",
        bbox_to_anchor=(1.3, 1.15),
        frameon=False,
        fontsize=8,
    )

    figure.tight_layout()

    return figure


def create_amino_acid_heatmap(
    amino_acid_df: pd.DataFrame,
):
    """Create a protein × amino-acid percentage heatmap."""
    matrix = amino_acid_df.pivot(
        index="Protein_ID",
        columns="Amino_Acid",
        values="Percentage",
    ).fillna(0)

    matrix = matrix.reindex(
        columns=list("ACDEFGHIKLMNPQRSTVWY"),
        fill_value=0,
    )

    figure_height = max(4.5, len(matrix) * 0.55)

    figure, ax = plt.subplots(
        figsize=(12, figure_height)
    )

    image = ax.imshow(
        matrix.values,
        aspect="auto",
        interpolation="nearest",
    )

    ax.set_xticks(np.arange(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, fontsize=9)

    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_yticklabels(matrix.index, fontsize=9)

    ax.set_xlabel(
        "Amino acid",
        fontweight="bold",
    )

    ax.set_ylabel(
        "Protein",
        fontweight="bold",
    )

    ax.set_title(
        "Amino-acid composition heatmap",
        fontsize=13,
        fontweight="bold",
        pad=12,
    )

    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix.iloc[row, column]
            ax.text(
                column,
                row,
                f"{value:.1f}",
                ha="center",
                va="center",
                fontsize=6,
            )

    colorbar = figure.colorbar(
        image,
        ax=ax,
        shrink=0.85,
    )

    colorbar.set_label(
        "Composition (%)",
        fontweight="bold",
    )

    figure.tight_layout()

    return figure


def create_violin_panel(
    properties_df: pd.DataFrame,
):
    """Create a four-panel violin plot for key properties."""
    specifications = [
        ("GRAVY", "(a) GRAVY"),
        ("Aliphatic_Index", "(b) Aliphatic index"),
        ("Molecular_Weight_kDa", "(c) Molecular weight"),
        ("Theoretical_pI", "(d) Theoretical pI"),
    ]

    figure, axes = plt.subplots(
        1,
        4,
        figsize=(12, 4.8),
    )

    for ax, (column, title) in zip(
        axes,
        specifications,
    ):
        values = pd.to_numeric(
            properties_df[column],
            errors="coerce",
        ).dropna()

        parts = ax.violinplot(
            values,
            showmeans=True,
            showmedians=True,
            showextrema=True,
        )

        jitter = np.linspace(
            0.96,
            1.04,
            len(values),
        )

        ax.scatter(
            jitter,
            values,
            s=24,
            edgecolor="black",
            linewidth=0.5,
            zorder=3,
        )

        ax.set_xticks([1])
        ax.set_xticklabels([column.replace("_", " ")], fontsize=8)
        ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.3)

    figure.suptitle(
        "Distribution of key physicochemical properties",
        fontsize=14,
        fontweight="bold",
    )

    figure.tight_layout(rect=[0, 0, 1, 0.93])

    return figure
