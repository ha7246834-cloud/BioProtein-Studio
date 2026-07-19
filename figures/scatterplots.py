import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def create_scatterplot(
    properties_df: pd.DataFrame,
    x_column: str,
    y_column: str,
):
    """Create a labeled physicochemical scatterplot."""

    x_values = pd.to_numeric(
        properties_df[x_column],
        errors="coerce",
    )

    y_values = pd.to_numeric(
        properties_df[y_column],
        errors="coerce",
    )

    valid_mask = (
        x_values.notna()
        & y_values.notna()
    )

    x_values = x_values[valid_mask]
    y_values = y_values[valid_mask]

    protein_ids = properties_df.loc[
        valid_mask,
        "Protein_ID",
    ]

    figure, ax = plt.subplots(
        figsize=(8, 6)
    )

    ax.scatter(
        x_values,
        y_values,
        s=65,
        edgecolor="black",
        linewidth=0.7,
    )

    for protein_id, x_value, y_value in zip(
        protein_ids,
        x_values,
        y_values,
    ):
        ax.annotate(
            protein_id,
            xy=(x_value, y_value),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )

    if len(x_values) >= 2:
        slope, intercept = np.polyfit(
            x_values,
            y_values,
            1,
        )

        fitted_values = (
            slope * x_values
            + intercept
        )

        order = np.argsort(x_values)

        ax.plot(
            x_values.iloc[order],
            fitted_values.iloc[order],
            linestyle="--",
            linewidth=1.2,
        )

        correlation = np.corrcoef(
            x_values,
            y_values,
        )[0, 1]

        ax.text(
            0.03,
            0.97,
            f"Pearson r = {correlation:.3f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=10,
            bbox={
                "boxstyle": "round",
                "facecolor": "white",
                "alpha": 0.8,
            },
        )

    ax.set_xlabel(
        x_column.replace("_", " "),
        fontsize=11,
        fontweight="bold",
    )

    ax.set_ylabel(
        y_column.replace("_", " "),
        fontsize=11,
        fontweight="bold",
    )

    ax.set_title(
        (
            f"{y_column.replace('_', ' ')} versus "
            f"{x_column.replace('_', ' ')}"
        ),
        fontsize=13,
        fontweight="bold",
        pad=12,
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.grid(
        linestyle="--",
        linewidth=0.5,
        alpha=0.3,
    )

    figure.tight_layout()

    return figure