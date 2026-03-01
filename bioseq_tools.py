from modules.dna_rna_tools import (
    is_nucleic_acid,
    transcribe,
    reverse,
    complement,
    reverse_complement
)
from modules.fastq_tools import filter_fastq

def run_dna_rna_tools(*args: str) -> str | list[str] | None:
    """Main function to apply a procedure to one or multiple DNA/RNA sequences.
    Arguments:
        *args: One or more sequences followed by the procedure name as the last argument.
    Returns:
        - Single result if one sequence is provided
        - List of results if multiple sequences are provided
        - None for sequences that are invalid
    Raises:
        ValueError: If no sequences or unknown procedure are provided.
    """
    if not args:
        raise ValueError("No sequences or procedure provided.")

    *sequences, procedure = args
    funcs = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement
    }

    if procedure not in funcs:
        raise ValueError(f"Unknown procedure: {procedure}")

    results: list[str | bool | None] = []
    for sequence in sequences:
        if procedure == "is_nucleic_acid" or is_nucleic_acid(sequence):
            results.append(funcs[procedure](sequence))
        else:
            print(f"Warning: Invalid sequence skipped: {sequence}")
            results.append(None)

    return results[0] if len(results) == 1 else results

