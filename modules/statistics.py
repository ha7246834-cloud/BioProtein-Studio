import numpy as np
import pandas as pd


NUMERIC_COLUMNS = [
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
    "Alpha_Helix_Fraction",
    "Beta_Turn_Fraction",
    "Beta_Sheet_Fraction",
]


def calculate_summary_statistics(
    properties_df: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate descriptive statistics for physicochemical properties."""
    available_columns = [
        column
        for column in NUMERIC_COLUMNS
        if column in properties_df.columns
    ]

    records = []

    for column in available_columns:
        values = pd.to_numeric(
            properties_df[column],
            errors="coerce",
        ).dropna()

        if values.empty:
            continue

        mean_value = values.mean()
        standard_deviation = values.std(ddof=1)

        coefficient_of_variation = (
            standard_deviation / mean_value * 100
            if mean_value != 0
            else np.nan
        )

        q1 = values.quantile(0.25)
        q3 = values.quantile(0.75)
        iqr = q3 - q1

        records.append(
            {
                "Property": column,
                "Count": int(values.count()),
                "Mean": round(mean_value, 4),
                "Median": round(values.median(), 4),
                "Standard_Deviation": round(
                    standard_deviation,
                    4,
                ),
                "Variance": round(
                    values.var(ddof=1),
                    4,
                ),
                "Minimum": round(values.min(), 4),
                "Q1": round(q1, 4),
                "Q3": round(q3, 4),
                "Maximum": round(values.max(), 4),
                "IQR": round(iqr, 4),
                "CV_Percent": round(
                    coefficient_of_variation,
                    2,
                )
                if not np.isnan(coefficient_of_variation)
                else np.nan,
            }
        )

    return pd.DataFrame(records)


def detect_outliers_iqr(
    properties_df: pd.DataFrame,
) -> pd.DataFrame:
    """Detect outliers using the 1.5 × IQR method."""
    available_columns = [
        column
        for column in NUMERIC_COLUMNS
        if column in properties_df.columns
    ]

    outlier_records = []

    for column in available_columns:
        values = pd.to_numeric(
            properties_df[column],
            errors="coerce",
        )

        q1 = values.quantile(0.25)
        q3 = values.quantile(0.75)
        iqr = q3 - q1

        lower_limit = q1 - 1.5 * iqr
        upper_limit = q3 + 1.5 * iqr

        for index, value in values.items():
            if pd.isna(value):
                continue

            if value < lower_limit or value > upper_limit:
                outlier_records.append(
                    {
                        "Protein_ID": properties_df.loc[
                            index,
                            "Protein_ID",
                        ],
                        "Property": column,
                        "Value": round(float(value), 4),
                        "Lower_Limit": round(
                            float(lower_limit),
                            4,
                        ),
                        "Upper_Limit": round(
                            float(upper_limit),
                            4,
                        ),
                        "Direction": (
                            "Low outlier"
                            if value < lower_limit
                            else "High outlier"
                        ),
                    }
                )

    return pd.DataFrame(outlier_records)


def calculate_correlation_matrix(
    properties_df: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate Pearson correlations among numeric properties."""
    available_columns = [
        column
        for column in NUMERIC_COLUMNS
        if column in properties_df.columns
    ]

    numeric_df = properties_df[
        available_columns
    ].apply(
        pd.to_numeric,
        errors="coerce",
    )

    return numeric_df.corr(method="pearson")