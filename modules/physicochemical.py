import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def calculate_aliphatic_index(sequence: str) -> float:
    """Calculate the aliphatic index of a protein sequence."""
    length = len(sequence)

    if length == 0:
        return 0.0

    percentage_a = sequence.count("A") / length * 100
    percentage_v = sequence.count("V") / length * 100
    percentage_i = sequence.count("I") / length * 100
    percentage_l = sequence.count("L") / length * 100

    return (
        percentage_a
        + 2.9 * percentage_v
        + 3.9 * (percentage_i + percentage_l)
    )


def calculate_extinction_coefficients(
    sequence: str,
) -> tuple[int, int]:
    """Calculate reduced and oxidized extinction coefficients."""
    number_w = sequence.count("W")
    number_y = sequence.count("Y")
    number_c = sequence.count("C")

    reduced = (
        number_w * 5500
        + number_y * 1490
    )

    oxidized = (
        reduced
        + (number_c // 2) * 125
    )

    return reduced, oxidized


def analyze_single_protein(
    protein_id: str,
    sequence: str,
) -> dict:
    """Calculate physicochemical properties for one protein."""
    analysis = ProteinAnalysis(sequence)

    molecular_weight = analysis.molecular_weight()
    theoretical_pi = analysis.isoelectric_point()
    instability_index = analysis.instability_index()
    aromaticity = analysis.aromaticity()
    gravy = analysis.gravy()

    helix, turn, sheet = (
        analysis.secondary_structure_fraction()
    )

    amino_counts = analysis.count_amino_acids()

    positive_residues = (
        amino_counts.get("K", 0)
        + amino_counts.get("R", 0)
    )

    negative_residues = (
        amino_counts.get("D", 0)
        + amino_counts.get("E", 0)
    )

    reduced_extinction, oxidized_extinction = (
        calculate_extinction_coefficients(sequence)
    )

    return {
        "Protein_ID": protein_id,
        "Length_aa": len(sequence),
        "Molecular_Weight_Da": round(
            molecular_weight,
            2,
        ),
        "Molecular_Weight_kDa": round(
            molecular_weight / 1000,
            3,
        ),
        "Theoretical_pI": round(
            theoretical_pi,
            2,
        ),
        "Instability_Index": round(
            instability_index,
            2,
        ),
        "Stability": (
            "Stable"
            if instability_index < 40
            else "Unstable"
        ),
        "Aliphatic_Index": round(
            calculate_aliphatic_index(sequence),
            2,
        ),
        "GRAVY": round(
            gravy,
            3,
        ),
        "Hydropathy_Class": (
            "Hydrophilic"
            if gravy < 0
            else "Hydrophobic"
        ),
        "Aromaticity": round(
            aromaticity,
            3,
        ),
        "Positive_Residues_KR": positive_residues,
        "Negative_Residues_DE": negative_residues,
        "Charge_Difference": (
            positive_residues - negative_residues
        ),
        "Alpha_Helix_Fraction": round(
            helix,
            3,
        ),
        "Beta_Turn_Fraction": round(
            turn,
            3,
        ),
        "Beta_Sheet_Fraction": round(
            sheet,
            3,
        ),
        "Extinction_Reduced": reduced_extinction,
        "Extinction_Oxidized": oxidized_extinction,
    }


def analyze_proteins(
    accepted_sequences: list,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Analyze all accepted sequences.

    Returns:
        physicochemical properties dataframe
        amino-acid composition dataframe
    """
    property_records = []
    amino_acid_records = []

    for record in accepted_sequences:
        protein_id = record["Protein_ID"]
        sequence = record["Sequence"]

        property_records.append(
            analyze_single_protein(
                protein_id,
                sequence,
            )
        )

        analysis = ProteinAnalysis(sequence)
        amino_counts = analysis.count_amino_acids()

        for amino_acid, count in amino_counts.items():
            amino_acid_records.append(
                {
                    "Protein_ID": protein_id,
                    "Amino_Acid": amino_acid,
                    "Count": count,
                    "Percentage": round(
                        count / len(sequence) * 100,
                        3,
                    ),
                }
            )

    properties_df = pd.DataFrame(property_records)
    amino_acid_df = pd.DataFrame(amino_acid_records)

    return properties_df, amino_acid_df