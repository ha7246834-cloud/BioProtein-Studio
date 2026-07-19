import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def create_correlation_heatmap(
    correlation_df: pd.DataFrame,
):
    """Create a publication-quality correlation heatmap."""

    figure_size = max(
        8,
        len(correlation_df.columns) * 0.75,
    )

    figure, ax = plt.subplots(
        figsize=(figure_size, figure_size)
    )

    image = ax.imshow(
        correlation_df.values,
        vmin=-1,
        vmax=1,
        aspect="equal",
        interpolation="nearest",
    )

    labels = [
        column.replace("_", " ")
        for column in correlation_df.columns
    ]

    ax.set_xticks(
        np.arange(len(labels))
    )

    ax.set_xticklabels(
        labels,
        rotation=45,
        ha="right",
        fontsize=8,
    )

    ax.set_yticks(
        np.arange(len(labels))
    )

    ax.set_yticklabels(
        labels,
        fontsize=8,
    )

    ax.set_title(
        "Pearson correlation matrix of physicochemical properties",
        fontsize=14,
        fontweight="bold",
        pad=14,
    )

    for row in range(correlation_df.shape[0]):
        for column in range(correlation_df.shape[1]):

            value = correlation_df.iloc[
                row,
                column,
            ]

            if pd.isna(value):
                label = "NA"
            else:
                label = f"{value:.2f}"

            ax.text(
                column,
                row,
                label,
                ha="center",
                va="center",
                fontsize=7,
            )

    colorbar = figure.colorbar(
        image,
        ax=ax,
        shrink=0.82,
    )

    colorbar.set_label(
        "Pearson correlation coefficient (r)",
        fontsize=10,
        fontweight="bold",
    )

    figure.tight_layout()

    return figure


def create_zscore_heatmap(
    properties_df: pd.DataFrame,
):
    """Create a Z-score normalized protein-property heatmap."""

    columns = [
        "Length_aa",
        "Molecular_Weight_kDa",
        "Theoretical_pI",
        "Instability_Index",
        "Aliphatic_Index",
        "GRAVY",
        "Aromaticity",
    ]

    labels = [
        "Length",
        "MW",
        "pI",
        "Instability",
        "Aliphatic",
        "GRAVY",
        "Aromaticity",
    ]

    available_columns = [
        column
        for column in columns
        if column in properties_df.columns
    ]

    available_labels = [
        labels[columns.index(column)]
        for column in available_columns
    ]

    matrix = properties_df[
        available_columns
    ].astype(float)

    standard_deviation = matrix.std(
        axis=0,
        ddof=0,
    ).replace(0, 1)

    z_scores = (
        matrix - matrix.mean(axis=0)
    ) / standard_deviation

    figure_height = max(
        4.5,
        len(properties_df) * 0.55,
    )

    figure, ax = plt.subplots(
        figsize=(10, figure_height)
    )

    image = ax.imshow(
        z_scores.values,
        aspect="auto",
        interpolation="nearest",
    )

    ax.set_xticks(
        np.arange(len(available_labels))
    )

    ax.set_xticklabels(
        available_labels,
        rotation=45,
        ha="right",
        fontsize=9,
    )

    ax.set_yticks(
        np.arange(len(properties_df))
    )

    ax.set_yticklabels(
        properties_df["Protein_ID"],
        fontsize=9,
    )

    ax.set_xlabel(
        "Physicochemical property",
        fontweight="bold",
    )

    ax.set_ylabel(
        "Protein",
        fontweight="bold",
    )

    ax.set_title(
        "Normalized physicochemical-property heatmap",
        fontsize=14,
        fontweight="bold",
        pad=12,
    )

    for row in range(z_scores.shape[0]):
        for column in range(z_scores.shape[1]):
            ax.text(
                column,
                row,
                f"{z_scores.iloc[row, column]:.2f}",
                ha="center",
                va="center",
                fontsize=7,
            )

    colorbar = figure.colorbar(
        image,
        ax=ax,
        shrink=0.85,
    )

    colorbar.set_label(
        "Z-score",
        fontweight="bold",
    )

    figure.tight_layout()

    return figure