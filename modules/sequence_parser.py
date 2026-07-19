import io
import re
from dataclasses import dataclass

from Bio import SeqIO


STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AMINO_ACIDS = set("XBZJUO")


@dataclass
class ProcessedSequence:
    protein_id: str
    original_sequence: str
    processed_sequence: str
    original_length: int
    analyzed_length: int
    status: str
    message: str


def clean_sequence(sequence: str) -> str:
    """Clean whitespace and terminal stop symbols from a protein sequence."""
    sequence = sequence.upper()
    sequence = re.sub(r"\s+", "", sequence)
    sequence = sequence.rstrip("*")
    return sequence


def parse_sequences(input_text: str) -> list:
    """Parse one raw protein sequence or multiple FASTA sequences."""
    input_text = input_text.strip()

    if not input_text:
        return []

    if input_text.startswith(">"):
        handle = io.StringIO(input_text)
        return list(SeqIO.parse(handle, "fasta"))

    cleaned_sequence = clean_sequence(input_text)

    handle = io.StringIO(
        f">Protein_1\n{cleaned_sequence}"
    )

    return list(SeqIO.parse(handle, "fasta"))


def process_ambiguous_residues(
    sequence: str,
    mode: str,
) -> tuple[str, str, list]:
    """
    Process ambiguous amino-acid residues.

    Returns:
        processed sequence
        warning message
        unsupported residues
    """
    sequence = clean_sequence(sequence)

    unsupported_residues = sorted(
        set(sequence)
        - STANDARD_AMINO_ACIDS
        - AMBIGUOUS_AMINO_ACIDS
    )

    if unsupported_residues:
        return sequence, "", unsupported_residues

    ambiguous_residues = sorted(
        set(sequence) & AMBIGUOUS_AMINO_ACIDS
    )

    if not ambiguous_residues:
        return sequence, "", []

    if mode == "Strict mode":
        return sequence, "", ambiguous_residues

    if mode == "Remove terminal X only":
        original_sequence = sequence
        sequence = sequence.rstrip("X")

        remaining_ambiguous = sorted(
            set(sequence) & AMBIGUOUS_AMINO_ACIDS
        )

        if remaining_ambiguous:
            return sequence, "", remaining_ambiguous

        removed_count = (
            len(original_sequence) - len(sequence)
        )

        warning = (
            f"{removed_count} terminal X residue(s) removed. "
            "Results should be treated as approximate."
        )

        return sequence, warning, []

    if mode == "Remove all ambiguous residues":
        original_length = len(sequence)

        sequence = "".join(
            residue
            for residue in sequence
            if residue in STANDARD_AMINO_ACIDS
        )

        removed_count = original_length - len(sequence)

        warning = (
            f"{removed_count} ambiguous residue(s) removed. "
            "Results should be treated as approximate."
        )

        return sequence, warning, []

    return sequence, "", ambiguous_residues


def prepare_sequences(
    input_text: str,
    ambiguous_mode: str,
) -> tuple[list, list]:
    """
    Parse, clean and validate all supplied protein sequences.

    Returns:
        accepted sequences
        quality-control records
    """
    records = parse_sequences(input_text)

    accepted_sequences = []
    quality_records = []

    for number, record in enumerate(records, start=1):
        protein_id = (
            record.id
            if record.id
            else f"Protein_{number}"
        )

        original_sequence = clean_sequence(
            str(record.seq)
        )

        processed_sequence, warning, invalid = (
            process_ambiguous_residues(
                original_sequence,
                ambiguous_mode,
            )
        )

        if invalid:
            quality_records.append(
                {
                    "Protein_ID": protein_id,
                    "Original_Length": len(original_sequence),
                    "Analyzed_Length": 0,
                    "Status": "Rejected",
                    "Message": (
                        "Unsupported residues: "
                        + ", ".join(invalid)
                    ),
                }
            )
            continue

        if not processed_sequence:
            quality_records.append(
                {
                    "Protein_ID": protein_id,
                    "Original_Length": len(original_sequence),
                    "Analyzed_Length": 0,
                    "Status": "Rejected",
                    "Message": (
                        "No standard amino acids remained "
                        "after sequence processing."
                    ),
                }
            )
            continue

        accepted_sequences.append(
            {
                "Protein_ID": protein_id,
                "Sequence": processed_sequence,
            }
        )

        quality_records.append(
            {
                "Protein_ID": protein_id,
                "Original_Length": len(original_sequence),
                "Analyzed_Length": len(processed_sequence),
                "Status": (
                    "Analyzed with warning"
                    if warning
                    else "Analyzed"
                ),
                "Message": (
                    warning
                    if warning
                    else "No sequence-quality issue detected."
                ),
            }
        )

    return accepted_sequences, quality_records