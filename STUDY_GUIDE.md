# Private study guide

This is a private learning aid, not a required submission artifact. It covers only approved work already implemented.

## A1.1 — loading a genome

### Purpose and input/output

Read the DNA sequence that the later ORF search will examine. The input is a FASTA file: a header starting with `>` followed by sequence lines. The supplied `genome.fasta` contains one record with ID `NC_001416.1`, description `Enterobacteria phage lambda, complete genome`, and 48,502 nucleotides.

The output of `load_genome(fasta_path)` is a Biopython `SeqRecord`. Its `id` identifies the record, and its `seq` holds the DNA as a `Seq` object. The function returns data in memory; it does not write a report or find ORFs.

### Steps and important functions

1. `SeqIO.read(fasta_path, "fasta")` reads exactly one record. Biopython raises an error for no records or more than one record. A missing file produces a file-not-found error.
2. Reject a record whose sequence is empty.
3. Convert the sequence to uppercase so `acgt` and `ACGT` are treated alike.
4. Compare its symbols with `ACGT`. Reject anything outside that set rather than guessing a base or skipping a position.
5. Return the record for use in the later ORF search.

This keeps the original sequence order and positions. A, C, G and T represent the four DNA bases. An ambiguity code such as `N` means the base is not uniquely specified. Rejecting ambiguity is our approved project convention, not a rule imposed by FASTA or by the assignment.

Biopython is required for FASTA parsing by the assignment. It will not perform ORF discovery for us.

### Environment and verification

`.venv/` isolates project packages from the existing Conda environment. After command-environment network failures, the student installed Biopython manually and supplied the FASTA. Biopython 1.88 was confirmed in `.venv` and pinned in `requirements.txt`. `pip check` reported no broken requirements.

The loader returned a `SeqRecord` containing a `Seq` for the supplied genome. Its ID was `NC_001416.1`, its sequence length was 48,502, and its alphabet contained only `ACGT`. These checks confirm the local header and sequence structure, not an independent checksum match with NCBI.

Small in-memory FASTA examples verified that mixed-case multiline input becomes `ACGTATGC`, and that no records, an empty sequence, multiple records, an ambiguous `N`, an invalid `X` and a gap `-` each raise `ValueError`. A nonexistent path in a temporary directory raised `FileNotFoundError`. The real genome file was not changed by verification.

The loader also passed Python syntax parsing with `ast.parse`; `git diff --check` reported no whitespace errors, and `git check-ignore` confirmed that the virtual environment is ignored. The student accepted A1.1. ORF discovery is covered in A1.2 below.

To inspect the loaded genome again from the project directory in PowerShell:

```powershell
.\.venv\Scripts\python.exe -B -c "from orf_finder import load_genome; record = load_genome('genome.fasta'); print(record.id, len(record.seq))"
```

### Possible instructor questions

1. **Why use Biopython here?** The assignment requires it for FASTA parsing, and `SeqIO.read` already handles the format.
2. **Why `read` rather than `parse`?** `read` requires exactly one record, matching our chosen single-genome input. `parse` iterates over potentially many records.
3. **Why uppercase the DNA?** Case is irrelevant to the base identity, and later codon comparisons must behave consistently.
4. **Why reject `N`?** It represents an unknown base. Our approved convention avoids guessing whether an uncertain codon is a start or STOP.
5. **Does loading discover genes?** No. It reads and validates DNA; the six-frame ORF search belongs to A1.2, and an ORF alone does not prove a biological gene.

## A1.2 — searching all six reading frames

### Purpose and input/output

Find complete open reading frames (ORFs): DNA regions beginning with `ATG` and ending at their first in-frame `TAA`, `TAG` or `TGA`. An ORF is a possible protein-coding region, not evidence by itself that a gene is expressed or functional.

`find_orfs(sequence)` takes the validated uppercase `Seq` from `load_genome(...).seq`. It returns a list of dictionaries with `frame`, `start`, `end`, `length` and `protein`. Length is in nucleotides and includes the STOP codon; protein omits the terminal `*`. This function applies no minimum length and writes no files. Filtering and CSV output are handled by `main()`, described in A1.3.

### Algorithm and important functions

1. Examine the forward sequence and its reverse complement. Complementary bases pair as A-T and C-G; reversing the complemented sequence lets us read the opposite strand in its own 5'-to-3' direction.
2. For each strand, use offsets 0, 1 and 2. A reading frame groups bases into codons (triplets) starting at one of these offsets. Two strands times three offsets gives six frames.
3. Move three bases at a time. At each `ATG`, start a candidate and examine subsequent codons in the same frame.
4. Stop at the first `TAA`, `TAG` or `TGA`. Ignore stops in other frames. A candidate with no complete in-frame STOP is not reported.
5. Record the candidate, then continue checking later start positions. This preserves nested ORFs: `ATGATGTAA` produces `MM` from positions 1-9 and `M` from positions 4-9.
6. Convert offsets into original-reference coordinates and translate the coding slice before STOP using the standard genetic code (`table=1`). Biopython provides `reverse_complement()` and `translate()`; the search itself is our own explicit loops.

The nested loops directly express the assignment's definition and are easy to trace on short DNA strings. They repeat some scanning for nested candidates, which is acceptable for the chosen 48,502-base genome.

### Coordinate conventions

Python uses 0-based, end-exclusive slices: `sequence[start:end]` includes `start` but excludes `end`. Here `end` is set to the position just after the STOP codon, so the nucleotide length is `end - start`.

- Forward-reference positions: `start + 1` through `end`.
- Reverse-strand positions on a reference of length `N`: `N - end + 1` through `N - start`.
- Both reported endpoints are inclusive, and the smaller coordinate is always first. For a reverse ORF, translation begins at the higher reference coordinate.
- Frames are `+1`, `+2`, `+3` for forward offsets 0-2 and `-1`, `-2`, `-3` for reverse-complement offsets 0-2.

Example: the 13-base reference `GGTTATTTCATGG` has reverse complement `CCATGAAATAACC`. Its reverse ORF occupies Python offsets `[2:11]` in the reverse complement, so it is in frame `-3` and maps to reference positions `13 - 11 + 1 = 3` through `13 - 2 = 11`. Its length is 9 nucleotides and its protein is `MK`.

### Verification

Exact hand-checked results passed for each of the six frames, including the reverse coordinate conversion. Additional checks covered all three STOP codons, a later STOP that must be ignored after the first, an out-of-frame STOP that must be skipped, nested starts, separate ORFs, ORFs reaching sequence boundaries, trailing incomplete codons, very short sequences, no start and no complete STOP.

The real lambda input produced 1,800 ORFs before minimum-length filtering: 340 in `+1`, 332 in `+2`, 326 in `+3`, 283 in `-1`, 248 in `-2` and 271 in `-3`. For every ORF, the reported coordinates were used to extract the original-reference segment and reverse-complement it where necessary. Checks confirmed the frame, nucleotide length, initial ATG, terminal STOP, absence of earlier in-frame STOPs, protein sequence and absence of duplicate frame/coordinate records. The count includes short and nested candidates and must not be interpreted as a gene count.

This short Python check can be rerun using the project environment, without creating a test framework:

```python
from Bio.Seq import Seq
from orf_finder import find_orfs

assert find_orfs(Seq("ATGAAATAA")) == [
    {"frame": "+1", "start": 1, "end": 9, "length": 9, "protein": "MK"}
]
assert find_orfs(Seq("GGTTATTTCATGG")) == [
    {"frame": "-3", "start": 3, "end": 11, "length": 9, "protein": "MK"}
]
```

### Possible instructor questions

1. **Why six frames?** Each strand can be divided into triplets at three different offsets.
2. **Why stop at the first in-frame STOP?** That is the ORF definition required by the assignment; an out-of-frame STOP does not terminate this candidate.
3. **Why can ORFs overlap?** We test every ATG separately, including starts inside a previously found candidate, as agreed in S0.
4. **Why do reverse coordinates need conversion?** Reverse-complement offsets refer to a reversed sequence. The assignment requires both strands to use positions on the original forward reference.
5. **Why is protein length shorter than nucleotide length divided by three?** The nucleotide length includes the final STOP codon, which does not encode an amino acid. Protein length is `length / 3 - 1`.

## A1.3 — command line, minimum length and CSV

### Purpose and input/output

Make the ORF search runnable using the assignment's command-line interface and save its results as a CSV table. CSV means comma-separated values: a header names the columns and each following row describes one ORF.

The input is a single-genome FASTA plus a minimum nucleotide length and output filename. The command used to generate the required ORF table was:

```powershell
.\.venv\Scripts\python.exe -B orf_finder.py genome.fasta --min-len 100 --out orfs.csv
```

`--min-len` defaults to 100 and must be a positive integer. `--out` defaults to `orfs.csv`. The minimum includes both ATG and STOP, as agreed in S0. `-B` is a Python option that avoids writing bytecode cache files; it does not change the analysis.

The CSV columns are `frame,start,end,length,protein`, with the same conventions as A1.2. A run with no qualifying ORFs still writes the header. Existing output files are not overwritten; the program reports an error instead. Use a new output path when rerunning if the previous results must be kept.

### Steps and important functions

1. `main()` creates an `argparse.ArgumentParser` to read the FASTA path, `--min-len` and `--out`. Noninteger or nonpositive thresholds are rejected.
2. `load_genome()` reads and validates the sequence, and `find_orfs()` discovers complete ORFs using the unchanged six-frame algorithm.
3. Keep each ORF whose nucleotide `length >= min_len`. Equality is included: a 9-nt ORF passes a threshold of 9.
4. Open a new output file with mode `x`. This prevents accidental replacement of previous results. UTF-8 specifies text encoding; `newline=""` lets the CSV writer manage line endings correctly on Windows.
5. `csv.DictWriter` maps each dictionary to the specified column order, writes the header and then the retained ORFs. Print the number of rows written.
6. Catch likely file and sequence errors and exit with a short message and status 2. A successful run exits with status 0. The `if __name__ == "__main__"` guard runs the CLI only when the file is executed directly, so importing the functions does not run the analysis.

This uses only standard-library modules in addition to the already required Biopython. A length threshold reduces short candidates but is not a statistical test or proof that the remaining ORFs are biological genes. Complete ORFs have lengths divisible by three, so a 100-nt threshold cannot admit a 100-nt or 101-nt ORF; the next possible length is 102 nt.

### Verification and results

Temporary, hand-checkable FASTA inputs tested ORFs of 6, 9 and 12 nucleotides at thresholds 1, 6, 7, 9, 10, 12 and 13. The checks confirmed inclusive filtering and a header-only CSV when nothing qualifies. A separate default-argument check excluded a 99-nt ORF and retained a 102-nt ORF in `orfs.csv`. A lowercase reverse-strand example retained frame `-3`, coordinates 3-11 and protein `MK` after CSV serialization.

Error checks passed for zero, negative and noninteger thresholds; missing FASTA; no records; an empty sequence; multiple records; ambiguous DNA; and a missing output directory. An existing output file remained unchanged. `--help` listed the required options.

The real genome run wrote 708 ORFs, with nucleotide lengths from 102 to 3,399. There were 146 in `+1`, 129 in `+2`, 138 in `+3`, 110 in `-1`, 92 in `-2` and 93 in `-3`. Every row was compared with the filtered search results and reconstructed from the original reference to verify coordinates, frame, nucleotide length, initial ATG, first in-frame STOP and protein translation. The student accepted A1.3.

This Python snippet checks the saved table without replacing it:

```python
import csv

with open("orfs.csv", newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle)
    assert reader.fieldnames == ["frame", "start", "end", "length", "protein"]
    rows = list(reader)

assert len(rows) == 708
assert all(int(row["length"]) >= 100 for row in rows)
assert all(int(row["end"]) - int(row["start"]) + 1 == int(row["length"]) for row in rows)
```

### Possible instructor questions

1. **What does `--min-len 100` measure?** At least 100 nucleotides, including ATG and STOP; it is not a protein-length threshold.
2. **Why is the shortest retained ORF 102 nt?** ORFs consist of complete triplets, so their lengths are divisible by three.
3. **Why use `csv.DictWriter`?** The ORF results are already dictionaries; it writes their values in an explicit, consistent column order.
4. **Why can a valid output contain only a header?** The search may find no complete ORFs meeting the requested threshold. That is a valid empty result, not a missing-column error.
5. **Does importing the script create a CSV?** No. The main guard keeps command-line execution separate from importing functions for checks or reuse.

## A1.4 — counting overlapping k-mers

### Purpose and input/output

A k-mer is a sequence fragment of length `k`. Its count describes how often it occurs, and its relative frequency is its fraction of all counted windows. These are descriptive sequence statistics, not evidence by themselves of a biological function.

The input is the single-genome FASTA, an integer `-k` and an output path `--out`. Both options are required. The command used for lambda was:

```powershell
.\.venv\Scripts\python.exe -B kmer_freq.py genome.fasta -k 3 --out kmers_k3.csv
```

`count_kmers(sequence, k)` returns a list of `(kmer, count, frequency)` tuples. `main()` writes those rows with the CSV header `kmer,count,frequency`. Only observed k-mers are reported. Existing output files are preserved.

### Steps and important functions

1. `main()` reads CLI arguments using `argparse`. It imports `load_genome()` from `orf_finder.py` to reuse the same Biopython parsing, uppercase conversion and input validation. Importing that module does not run the ORF search because of its main guard.
2. `count_kmers()` converts the `Seq` to a string and calculates `total = N - k + 1`, where `N` is the genome length. Reject `k < 1` and `k > N`.
3. Visit start positions 0 through `N - k` and take `sequence[start:start + k]`. Advancing by one base includes overlapping windows. Do not add reverse-complement counts or join the end of the sequence back to its beginning.
4. Start with an empty `collections.Counter`. In an explicit loop, extract each window into `kmer` and increment `counts[kmer]`. Counter gives an unseen fragment an initial count of zero.
5. The named `sorting_rule()` unpacks each `(kmer, count)` pair and returns `(-count, kmer)`. Negative counts put frequent fragments first, and the string breaks equal-count ties alphabetically. `sorted_counts` holds the ordered entries.
6. In a separate loop over `sorted_counts`, calculate `frequency = count / total` and append `(kmer, count, frequency)` to `rows`. The denominator is the number of windows, not the number of distinct k-mers. Return the completed list after the loop.
7. `csv.writer` writes the header and tuple rows to a new UTF-8 file. Mode `x` prevents replacing existing results; `newline=""` lets the CSV writer handle line endings. Likely file/sequence errors produce a short message and exit status 2.

For `AAACAA` with `k=2`, the five windows are `AA`, `AA`, `AC`, `CA`, `AA`. Therefore the sorted rows are `AA: 3, 0.6`, `AC: 1, 0.2` and `CA: 1, 0.2`. Counting nonoverlapping fragments would miss occurrences and give the wrong denominator.

### Verification and results

Hand-counted checks passed for the example above, repeated overlapping `AA` in `AAAA`, and alphabetical ties in `ATATA`. Boundary checks covered `k=1` (individual bases) and `k=N` (one whole-sequence window). An `ATG` input at `k=3` produced only `ATG`, confirming that reverse-complement `CAT` was not added. Lowercase multiline FASTA produced the expected CSV after normalization.

CLI checks passed for missing required options, noninteger `k`, zero, negative and oversized `k`, missing/empty/multiple-record FASTA, an empty sequence, ambiguous DNA, a missing output directory, preservation of an existing output and `--help`.

The lambda CSV contains 64 distinct 3-mers. Counts sum to 48,500, which equals `48,502 - 3 + 1`. Every count was independently verified by checking whether its k-mer begins at each possible sequence position, including overlaps. Frequencies match `count / 48,500`, sum to 1 within floating-point tolerance, and rows follow the required sorting rule. This genome happens to contain all `4**3 = 64` possible 3-mers; reporting all possible strings is not required when some are absent. The student accepted A1.4.

These Python checks can be rerun without changing output files:

```python
import csv
import math
from Bio.Seq import Seq
from kmer_freq import count_kmers

assert count_kmers(Seq("AAACAA"), 2) == [
    ("AA", 3, 0.6), ("AC", 1, 0.2), ("CA", 1, 0.2)
]
assert count_kmers(Seq("AAAA"), 2) == [("AA", 3, 1.0)]

with open("kmers_k3.csv", newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))

assert len(rows) == 64
assert sum(int(row["count"]) for row in rows) == 48500
assert math.isclose(math.fsum(float(row["frequency"]) for row in rows), 1.0)
```

### Possible instructor questions

1. **Why are there `N - k + 1` windows?** The last valid start is `N - k`; including start 0 gives `N - k + 1` positions.
2. **Why count overlapping fragments?** Each position defines a window. For example, `AAAA` contains three overlapping occurrences of `AA`, not two.
3. **What is the frequency denominator?** All counted windows, including repeated fragments, not the number of distinct k-mers.
4. **Why sort by both count and string?** Descending count meets the assignment; alphabetical tie-breaking makes equal-count results deterministic.
5. **Does the k-mer script count only ORFs?** No. It counts the whole supplied forward genome. Only FASTA loading is reused from `orf_finder.py`; ORF discovery is not called.

## A1.5 — Assignment 1 acceptance check

### Purpose and input/output

Check the finished Assignment 1 against the actual PDF, rather than assuming that successful script execution proves every requirement was met. Inputs are page 1 of the PDF, the two scripts, `genome.fasta` and the two saved CSVs. The output of this stage is the recorded audit below; no extra submission artifact or analysis code was needed.

### Requirement mapping

| PDF requirement | Implementation and verification |
| --- | --- |
| Two Python scripts taking FASTA input | `orf_finder.py` and `kmer_freq.py` both accept a positional FASTA path. |
| Biopython parsing and an allowed genome | Both call `load_genome()`, which uses `Bio.SeqIO.read`. The supplied record identifies lambda `NC_001416.1` and contains 48,502 nt. |
| Own ORF function in all six frames | `find_orfs()` uses explicit loops over both orientations and three offsets; exact examples passed for all six frames. Biopython does not discover the ORFs. |
| ATG to the first in-frame TAA, TAG or TGA | The search checks complete codons and stops at the first qualifying STOP. Every saved ORF was checked for no earlier in-frame STOP. |
| Frame, coordinates, length and translated protein | CSV columns are `frame,start,end,length,protein`. Every saved row passed reconstruction and translation checks. |
| 1-based original-reference coordinates on both strands | Forward and reverse examples passed; every saved row's interval and frame were checked against the original sequence. |
| Minimum-length argument and example interface | The CLI accepts `--min-len 100 --out ...`; `orfs.csv` contains 708 qualifying ORFs. Length units, STOP handling and nested starts follow the conventions agreed in S0. |
| CLI-specified k, counts, relative frequencies and descending order | The k-mer CLI accepts `-k 3 --out ...`. All 64 saved counts were independently recounted, with total 48,500; frequencies and sorting passed checks. |
| Both scripts and both requested output tables | `orf_finder.py`, `kmer_freq.py`, `orfs.csv` and `kmers_k3.csv` exist and are nonempty. |

### Audit steps and verification

1. Read the Assignment 1 text directly from the PDF and inspect both scripts. Trace their calls to `load_genome()`, `find_orfs()`, `count_kmers()` and the two `main()` entry points.
2. Run the ORF CLI at `--min-len 100` and the k-mer CLI at `-k 3`, using new output paths in a temporary directory. Both commands exited successfully.
3. Compare each temporary CSV with its saved counterpart. Both pairs were byte-for-byte identical: the same text, values, ordering and line endings.
4. Check all saved rows for their sequence-derived properties, not only file existence. ORFs passed coordinate, frame, length, first-STOP and translation checks. The 3-mers passed independent overlapping occurrence counts, denominator, frequency and sorting checks.
5. Confirm that original file hashes stayed unchanged, both `--help` commands worked, and `pip check` reported no broken requirements. Temporary files were cleaned up.

No implementation gaps or necessary code corrections were identified. This establishes requirement coverage and reproducibility for the chosen input and conventions; the ORF list remains a list of sequence-defined candidates, not a set of experimentally confirmed genes. A1.5 awaits student acceptance.

### Possible instructor questions

1. **Which files does Assignment 1 explicitly request?** `orf_finder.py`, `kmer_freq.py`, `orfs.csv` and `kmers_k3.csv`. The input FASTA supports local reproduction; this study guide is a private learning aid.
2. **Why compare regenerated CSVs with the saved files?** It checks that the saved deliverables actually match what the current scripts produce from the chosen input and arguments.
3. **Why are totals useful checks?** The 3-mer count sum must equal `N - 3 + 1`, and every reported ORF's inclusive interval must agree with its nucleotide length. These catch missing windows and coordinate mistakes that file-existence checks cannot detect.
