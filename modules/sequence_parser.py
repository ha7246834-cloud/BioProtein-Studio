import io
import re

from Bio import SeqIO


STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AMINO_ACIDS = set("XBZJUO")


def clean_sequence(sequence: str) -> str:
    sequence = sequence.upper()
    sequence = re.sub(r"\\s+", "", sequence)
    sequence = sequence.rstrip("*")
    return sequence


def parse_sequences(input_text: str) -> list:
    input_text = input_text.strip()

    if not input_text:
        return []

    if input_text.startswith(">"):
        return list(SeqIO.parse(io.StringIO(input_text), "fasta"))

    cleaned = clean_sequence(input_text)
    return list(
        SeqIO.parse(
            io.StringIO(f">Protein_1\\n{cleaned}"),
            "fasta",
        )
    )


def process_ambiguous_residues(
    sequence: str,
    mode: str,
) -> tuple[str, str, list]:
    sequence = clean_sequence(sequence)

    unsupported = sorted(
        set(sequence)
        - STANDARD_AMINO_ACIDS
        - AMBIGUOUS_AMINO_ACIDS
    )

    if unsupported:
        return sequence, "", unsupported

    ambiguous = sorted(
        set(sequence) & AMBIGUOUS_AMINO_ACIDS
    )

    if not ambiguous:
        return sequence, "", []

    if mode == "Strict mode":
        return sequence, "", ambiguous

    if mode == "Remove terminal X only":
        original = sequence
        sequence = sequence.rstrip("X")

        remaining = sorted(
            set(sequence) & AMBIGUOUS_AMINO_ACIDS
        )

        if remaining:
            return sequence, "", remaining

        removed = len(original) - len(sequence)

        return (
            sequence,
            (
                f"{removed} terminal X residue(s) removed. "
                "Results should be treated as approximate."
            ),
            [],
        )

    if mode == "Remove all ambiguous residues":
        original_length = len(sequence)

        sequence = "".join(
            residue
            for residue in sequence
            if residue in STANDARD_AMINO_ACIDS
        )

        removed = original_length - len(sequence)

        return (
            sequence,
            (
                f"{removed} ambiguous residue(s) removed. "
                "Results should be treated as approximate."
            ),
            [],
        )

    return sequence, "", ambiguous


def prepare_sequences(
    input_text: str,
    ambiguous_mode: str,
    minimum_length: int = 30,
    allow_short_peptides: bool = False,
) -> tuple[list, list]:
    records = parse_sequences(input_text)

    accepted = []
    quality = []

    for number, record in enumerate(records, start=1):
        protein_id = record.id or f"Protein_{number}"
        raw_sequence = str(record.seq)

        had_lowercase = any(
            character.islower()
            for character in raw_sequence
        )

        original_sequence = clean_sequence(raw_sequence)

        processed_sequence, warning, invalid = (
            process_ambiguous_residues(
                original_sequence,
                ambiguous_mode,
            )
        )

        messages = []

        if had_lowercase:
            messages.append(
                "Lowercase residues were converted to uppercase."
            )

        if invalid:
            quality.append(
                {
                    "Protein_ID": protein_id,
                    "Original_Length": len(original_sequence),
                    "Analyzed_Length": 0,
                    "Status": "Rejected",
                    "Message": (
                        "Unsupported or ambiguous residues: "
                        + ", ".join(invalid)
                    ),
                }
            )
            continue

        if not processed_sequence:
            quality.append(
                {
                    "Protein_ID": protein_id,
                    "Original_Length": len(original_sequence),
                    "Analyzed_Length": 0,
                    "Status": "Rejected",
                    "Message": (
                        "No standard amino acids remained after processing."
                    ),
                }
            )
            continue

        sequence_length = len(processed_sequence)

        if sequence_length < minimum_length:
            if not allow_short_peptides:
                quality.append(
                    {
                        "Protein_ID": protein_id,
                        "Original_Length": len(original_sequence),
                        "Analyzed_Length": 0,
                        "Status": "Rejected",
                        "Message": (
                            f"Sequence length is {sequence_length} aa. "
                            f"Minimum accepted length is "
                            f"{minimum_length} aa."
                        ),
                    }
                )
                continue

            messages.append(
                f"Short peptide detected ({sequence_length} aa). "
                "Values may not represent a complete protein."
            )

        if warning:
            messages.append(warning)

        accepted.append(
            {
                "Protein_ID": protein_id,
                "Sequence": processed_sequence,
            }
        )

        quality.append(
            {
                "Protein_ID": protein_id,
                "Original_Length": len(original_sequence),
                "Analyzed_Length": sequence_length,
                "Status": (
                    "Analyzed with warning"
                    if messages
                    else "Analyzed"
                ),
                "Message": (
                    " ".join(messages)
                    if messages
                    else "No sequence-quality issue detected."
                ),
            }
        )

    return accepted, quality
