"""Count overlapping k-mers and write counts and relative frequencies to CSV."""

import argparse
import csv
from collections import Counter

from orf_finder import load_genome


def count_kmers(sequence, k):
    """Return (kmer, count, frequency) rows for validated uppercase DNA."""
    sequence = str(sequence)
    total = len(sequence) - k + 1
    if k < 1 or total < 1:
        raise ValueError("k must be between 1 and the sequence length.")

    counts = Counter()
    for start in range(total):
        kmer = sequence[start:start + k]
        counts[kmer] += 1

    def sorting_rule(kmer_count):
        kmer, count = kmer_count
        # Negative counts put frequent k-mers first; the string breaks ties.
        return (-count, kmer)

    sorted_counts = sorted(counts.items(), key=sorting_rule)

    rows = []
    for kmer, count in sorted_counts:
        frequency = count / total
        rows.append((kmer, count, frequency))

    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", help="Input FASTA containing one genome.")
    parser.add_argument("-k", type=int, required=True, help="k-mer length: 1 through genome length.")
    parser.add_argument("--out", required=True, help="New output CSV file.")
    args = parser.parse_args()

    try:
        genome = load_genome(args.fasta)
        rows = count_kmers(genome.seq, args.k)

        # Exclusive creation protects existing results from accidental overwrites.
        with open(args.out, "x", newline="", encoding="utf-8") as output:
            writer = csv.writer(output)
            writer.writerow(["kmer", "count", "frequency"])
            writer.writerows(rows)
    except (OSError, ValueError) as error:
        parser.exit(2, f"Error: {error}\n")

    print(f"Wrote {len(rows)} distinct {args.k}-mers to {args.out}.")


if __name__ == "__main__":
    main()
