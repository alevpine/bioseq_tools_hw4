def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str | None = None):
    """
    Convert multiline FASTA sequences into one-line sequences.
    """
    if output_fasta is None:
        output_fasta = input_fasta.replace(".fasta", "_oneline.fasta")

    with open(input_fasta) as f, open(output_fasta, "w") as out:
        seq = ""
        header = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    out.write(header + "\n" + seq + "\n")
                header, seq = line, ""
            else:
                seq += line
        if header:
            out.write(header + "\n" + seq + "\n")

    print(f"FASTA converted and saved to {output_fasta}")


def parse_blast_output(input_file: str, output_file: str):
    """
    Extract best hits from BLAST output and save them alphabetically.
    """
    best_hits = set()

    with open(input_file) as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if "Sequences producing significant alignments:" in line:
            j = i + 2
            if j < len(lines):
                parts = lines[j].split()
                if parts:
                    hit = parts[0]
                    best_hits.add(hit)

    with open(output_file, "w") as out:
        for hit in sorted(best_hits):
            out.write(hit + "\n")

    print(f"Parsed BLAST results saved to {output_file}")

