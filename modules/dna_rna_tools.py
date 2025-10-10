def is_nucleic_acid(sequence: str) -> bool:
    "Check if sequence is DNA or RNA."
    dna = {"A", "T", "G", "C"}
    rna = {"A", "U", "G", "C"}
    seq_upper = sequence.upper()
    return set(seq_upper).issubset(dna) or set(seq_upper).issubset(rna)


def transcribe(sequence: str) -> str:
    "Transcribe DNA to RNA."
    return "".join("U" if nuc == "T" else "u" if nuc == "t" else nuc for nuc in sequence)


def reverse(sequence: str) -> str:
    "Return reversed sequence."
    return sequence[::-1]


def complement(sequence: str) -> str:
    "Return complementary sequence (DNA or RNA)."
    comp_map = {}
    if "U" in sequence.upper():
        # RNA
        comp_map = {"A": "U", "a": "u",
                    "U": "A", "u": "a",
                    "G": "C", "g": "c",
                    "C": "G", "c": "g"}
    else:
        # DNA
        comp_map = {"A": "T", "a": "t",
                    "T": "A", "t": "a",
                    "G": "C", "g": "c",
                    "C": "G", "c": "g"}
    return "".join(comp_map.get(nuc, nuc) for nuc in sequence)


def reverse_complement(sequence: str) -> str:
    "Return reverse complementary sequence."
    return reverse(complement(sequence))


def run_dna_rna_tools(*args: str) -> str | list[str] | None:
    "Main function: apply procedure to one or more sequences."
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

    results = []
    for sequence in sequences:
        if procedure == "is_nucleic_acid" or is_nucleic_acid(sequence):
            results.append(funcs[procedure](sequence))
        else:
            print(f"Warning: Invalid sequence skipped: {sequence}")
            results.append(None)

    return results[0] if len(results) == 1 else results

                  
