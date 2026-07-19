import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def apply_publication_style(
    ax,
    title,
    x_label="Protein",
    y_label="Value",
):
    """Apply consistent scientific-figure formatting."""
    ax.set_title(
        title,
        fontsize=13,
        fontweight="bold",
        pad=10,
    )

    ax.set_xlabel(
        x_label,
        fontsize=11,
        fontweight="bold",
    )

    ax.set_ylabel(
        y_label,
        fontsize=11,
        fontweight="bold",
    )

    ax.tick_params(
        axis="both",
        labelsize=9,
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.grid(
        axis="y",
        linestyle="--",
        linewidth=0.6,
        alpha=0.3,
    )


def add_value_labels(
    ax,
    bars,
    decimals=2,
):
    """Add numerical values above or below bars."""
    for bar in bars:
        value = bar.get_height()

        if value >= 0:
            offset = 3
            vertical_alignment = "bottom"
        else:
            offset = -3
            vertical_alignment = "top"

        ax.annotate(
            f"{value:.{decimals}f}",
            xy=(
                bar.get_x()
                + bar.get_width() / 2,
                value,
            ),
            xytext=(0, offset),
            textcoords="offset points",
            ha="center",
            va=vertical_alignment,
            fontsize=8,
        )


def set_axis_limits(
    ax,
    values,
):
    """Set safe axis limits for positive and negative values."""
    values = pd.to_numeric(
        values,
        errors="coerce",
    ).dropna()

    minimum = min(
        0,
        float(values.min()),
    )

    maximum = max(
        0,
        float(values.max()),
    )

    data_range = maximum - minimum

    if data_range == 0:
        data_range = 1

    ax.set_ylim(
        minimum - 0.18 * data_range,
        maximum + 0.22 * data_range,
    )


def create_bar_chart(
    properties_df,
    value_column,
    title,
    y_label,
    decimals=2,
    reference_line=None,
    reference_label=None,
):
    """Create one publication-quality bar chart."""
    figure, ax = plt.subplots(
        figsize=(9, 5.5)
    )

    x_positions = np.arange(
        len(properties_df)
    )

    values = pd.to_numeric(
        properties_df[value_column],
        errors="coerce",
    )

    bars = ax.bar(
        x_positions,
        values,
        width=0.67,
        edgecolor="black",
        linewidth=0.8,
    )

    ax.set_xticks(x_positions)

    ax.set_xticklabels(
        properties_df["Protein_ID"],
        rotation=45
        if len(properties_df) > 5
        else 0,
        ha="right"
        if len(properties_df) > 5
        else "center",
    )

    apply_publication_style(
        ax=ax,
        title=title,
        y_label=y_label,
    )

    if reference_line is not None:
        ax.axhline(
            reference_line,
            linestyle="--",
            linewidth=1.2,
            label=reference_label,
        )

        if reference_label:
            ax.legend(
                frameon=False,
                fontsize=9,
            )

    set_axis_limits(
        ax,
        values,
    )

    add_value_labels(
        ax,
        bars,
        decimals,
    )

    figure.tight_layout()

    return figure


def create_boxplot_panel(
    properties_df,
):
    """Create paper-style four-panel boxplots."""
    properties = [
        (
            "GRAVY",
            "GRAVY",
            "(a)",
        ),
        (
            "Aliphatic_Index",
            "Aliphatic index",
            "(b)",
        ),
        (
            "Molecular_Weight_kDa",
            "Molecular weight (kDa)",
            "(c)",
        ),
        (
            "Theoretical_pI",
            "Theoretical pI",
            "(d)",
        ),
    ]

    figure, axes = plt.subplots(
        1,
        4,
        figsize=(12, 4.7),
    )

    for ax, (
        column,
        x_label,
        panel_label,
    ) in zip(
        axes,
        properties,
    ):
        values = pd.to_numeric(
            properties_df[column],
            errors="coerce",
        ).dropna()

        ax.boxplot(
            values,
            widths=0.45,
            patch_artist=True,
            showmeans=True,
            meanprops={
                "marker": "D",
                "markeredgecolor": "black",
                "markerfacecolor": "white",
                "markersize": 5,
            },
            medianprops={
                "linewidth": 1.5,
            },
        )

        jitter = np.linspace(
            0.93,
            1.07,
            len(values),
        )

        ax.scatter(
            jitter,
            values,
            s=28,
            edgecolor="black",
            linewidth=0.5,
            zorder=3,
        )

        ax.set_xticks([1])
        ax.set_xticklabels(
            [x_label],
            fontsize=9,
        )

        ax.set_title(
            panel_label,
            loc="left",
            fontsize=12,
            fontweight="bold",
        )

        ax.tick_params(
            axis="y",
            labelsize=8,
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax.grid(
            axis="y",
            linestyle="--",
            linewidth=0.5,
            alpha=0.3,
        )

    figure.suptitle(
        "Distribution of physicochemical properties",
        fontsize=14,
        fontweight="bold",
    )

    figure.tight_layout(
        rect=[0, 0, 1, 0.93]
    )

    return figure


def create_gene_wise_panel(
    properties_df,
):
    """Create four-panel gene-wise comparison figure."""
    specifications = [
        (
            "Length_aa",
            "(a) Protein length",
            "Length (aa)",
            0,
            None,
        ),
        (
            "Molecular_Weight_kDa",
            "(b) Molecular weight",
            "Molecular weight (kDa)",
            2,
            None,
        ),
        (
            "Theoretical_pI",
            "(c) Theoretical pI",
            "Theoretical pI",
            2,
            7,
        ),
        (
            "GRAVY",
            "(d) GRAVY score",
            "GRAVY",
            3,
            0,
        ),
    ]

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(12, 9),
    )

    axes = axes.flatten()

    x_positions = np.arange(
        len(properties_df)
    )

    for ax, specification in zip(
        axes,
        specifications,
    ):
        (
            column,
            title,
            y_label,
            decimals,
            reference_line,
        ) = specification

        values = pd.to_numeric(
            properties_df[column],
            errors="coerce",
        )

        bars = ax.bar(
            x_positions,
            values,
            width=0.67,
            edgecolor="black",
            linewidth=0.7,
        )

        ax.set_xticks(x_positions)

        ax.set_xticklabels(
            properties_df["Protein_ID"],
            rotation=45,
            ha="right",
            fontsize=8,
        )

        apply_publication_style(
            ax=ax,
            title=title,
            y_label=y_label,
        )

        if reference_line is not None:
            ax.axhline(
                reference_line,
                linestyle="--",
                linewidth=1,
            )

        set_axis_limits(
            ax,
            values,
        )

        add_value_labels(
            ax,
            bars,
            decimals,
        )

    figure.suptitle(
        "Gene-wise physicochemical comparison",
        fontsize=15,
        fontweight="bold",
    )

    figure.tight_layout(
        rect=[0, 0, 1, 0.96]
    )

    return figure


def create_amino_acid_chart(
    amino_acid_df,
    selected_protein,
):
    """Create amino-acid composition bar chart."""
    selected_data = amino_acid_df[
        amino_acid_df["Protein_ID"]
        == selected_protein
    ].copy()

    figure, ax = plt.subplots(
        figsize=(10, 5.5)
    )

    bars = ax.bar(
        selected_data["Amino_Acid"],
        selected_data["Percentage"],
        edgecolor="black",
        linewidth=0.7,
    )

    apply_publication_style(
        ax=ax,
        title=(
            f"Amino-acid composition of "
            f"{selected_protein}"
        ),
        x_label="Amino acid",
        y_label="Composition (%)",
    )

    add_value_labels(
        ax,
        bars,
        decimals=1,
    )

    maximum = pd.to_numeric(
        selected_data["Percentage"],
        errors="coerce",
    ).max()

    ax.set_ylim(
        0,
        maximum * 1.25
        if maximum > 0
        else 1,
    )

    figure.tight_layout()

    return figure