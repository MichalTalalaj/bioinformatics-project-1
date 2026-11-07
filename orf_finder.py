"""Read a genome and find complete ORFs in all six reading frames."""

import argparse
import csv

from Bio import SeqIO


def load_genome(fasta_path):
    """Return one nonempty FASTA record containing only A, C, G and T."""
    record = SeqIO.read(fasta_path, "fasta")
    if not record.seq:
        raise ValueError("FASTA sequence is empty.")

    record.seq = record.seq.upper()
    if set(str(record.seq)) - set("ACGT"):
        raise ValueError("FASTA sequence must contain only A, C, G and T.")

    return record


def find_orfs(sequence):
    """Find ATG-to-first-STOP ORFs in a validated uppercase Bio.Seq sequence."""
    stop_codons = {"TAA", "TAG", "TGA"}
    sequence_length = len(sequence)
    orfs = []

    for strand, strand_sequence in [(1, sequence), (-1, sequence.reverse_complement())]:
        for offset in range(3):
            for start in range(offset, sequence_length - 2, 3):
                if strand_sequence[start:start + 3] != "ATG":
                    continue

                for stop in range(start + 3, sequence_length - 2, 3):
                    if strand_sequence[stop:stop + 3] not in stop_codons:
                        continue

                    end = stop + 3
                    if strand == 1:
                        reference_start = start + 1
                        reference_end = end
                    else:
                        # Map reverse-complement offsets onto the original reference.
                        reference_start = sequence_length - end + 1
                        reference_end = sequence_length - start

                    orfs.append({
                        "frame": f"{strand * (offset + 1):+d}",
                        "start": reference_start,
                        "end": reference_end,
                        "length": end - start,
                        "protein": str(strand_sequence[start:stop].translate(table=1)),
                    })
                    break

    return orfs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", help="Input FASTA containing one genome.")
    parser.add_argument(
        "--min-len", type=int, default=100,
        help="Minimum ORF length in nucleotides, including STOP (default: 100).",
    )
    parser.add_argument("--out", default="orfs.csv", help="New output CSV file (default: orfs.csv).")
    args = parser.parse_args()
    if args.min_len < 1:
        parser.error("--min-len must be a positive integer.")

    try:
        genome = load_genome(args.fasta)
        orfs = [orf for orf in find_orfs(genome.seq) if orf["length"] >= args.min_len]

        # Exclusive creation protects existing results from accidental overwrites.
        with open(args.out, "x", newline="", encoding="utf-8") as output:
            writer = csv.DictWriter(output, fieldnames=["frame", "start", "end", "length", "protein"])
            writer.writeheader()
            writer.writerows(orfs)
    except (OSError, ValueError) as error:
        parser.exit(2, f"Error: {error}\n")

    print(f"Wrote {len(orfs)} ORFs to {args.out}.")


if __name__ == "__main__":
    main()
