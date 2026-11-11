# Private study guide

This is a private learning aid, not a required submission artifact. It covers only approved work already implemented.

For hands-on Assignment 1 verification, open [manual_checks.ipynb](../assignment_1/manual_checks.ipynb), select the project's `.venv` kernel, restart it and run cells top to bottom. This student-requested notebook is also a private study aid. Expected answers appear before each check; an assertion failure means the result needs investigation. The notebook uses the existing scripts and temporary outputs, leaving the saved genome and CSVs unchanged.

Assignment 2 is covered at the end of this guide. Its commands run from `assignment_2/`.

For function-by-function live inspection, open [Assignment 2 manual checks](../assignment_2/manual_checks.ipynb), select `.venv`, restart the kernel and run top to bottom. It displays raw responses and the existing functions' outputs, with online record links. Internet access is required; temporary directories protect your saved report during writing and CLI checks.

Assignment 1 files now live in `assignment_1/`; run the commands and Python examples below from that directory. The virtual environment and `requirements.txt` remain at the project root. Both CLIs overwrite their specified output CSV on a successful run, as requested by the student. Repeated-run checks in temporary copies reproduced the saved tables without changing the originals.

## Using the Assignment 1 manual notebook

The notebook's purpose is to let you inspect small expected answers before checking the real genome and saved results. Inputs are deliberately constructed DNA examples and the existing FASTA/CSVs; outputs are readable values and passed assertions. `load_genome()`, `find_orfs()` and `count_kmers()` are imported from the actual scripts, so this does not maintain another implementation.

Run setup first and check the displayed Python executable. Next check the genome identity and length; compare six-frame ORFs with the explicit expected frame/coordinate/protein rows; check nested starts and first-STOP behavior; hand-count the k-mers; then run both CLIs into temporary files and compare with the saved CSVs. A separate toy CLI check covers minimum-length inclusion/exclusion. The final cell checks totals, basic CSV row properties and original file hashes. Restart the kernel after editing a script so imports are refreshed.

All seven code cells passed when executed directly with the project Python from both supported working directories. The saved outputs were reproduced exactly, and original file hashes were unchanged. No Jupyter installation or graphical session was performed. These checks do not establish correctness for every possible sequence, and an ORF remains a sequence-defined candidate rather than a confirmed gene. A1.5 acceptance remains a separate decision.

Possible instructor questions:

1. **Why use hand-worked examples?** Their expected answers are independent of the implementation and help expose coordinate or counting mistakes.
2. **Does reproducing a saved CSV prove the algorithm is correct?** No. It checks reproducibility; the independently expected toy examples check selected aspects of correctness.
3. **Why use the kernel's Python for CLI checks?** It ensures subprocesses use the same interpreter and packages as the notebook.
4. **Why write temporary outputs?** It exercises the real commands without replacing the student's saved results.

## Using the Assignment 2 manual notebook

Its purpose is to inspect every `genereport.py` function separately and compare extracted annotations with live database records. Inputs are APOE/Homo sapiens and actual API responses. Outputs include raw JSON/GenBank, selected identifiers, parsed fields, Markdown source and temporary files. The notebook imports the real script; it contains no replacement API implementation or fabricated live data.

Follow the numbered cells: inspect NCBI search, ID selection and summary; compare linked RNA IDs and raw GenBank CDS blocks with the two parsers and complete mRNA output; inspect UniProt reviewed/fallback selection, raw comments, text extraction and complete protein fields; inspect Ensembl raw transcript identifiers and parsed fields; then check escaping, field formatting, report assembly, writing and the CLI. Change the selected RNA/protein index to inspect additional records and rerun the dependent parser cells. The CLI cell uses the kernel's Python, performs fresh requests and writes only in a temporary directory.

For online comparison, match the scientific organism, accession/version and assembly before comparing fields. Previously observed transcript/PDB/keyword counts are examples, not assertions that online records can never change. Review any diff against the saved report against current raw data. Missing annotations do not imply biological absence, and extracted disease susceptibility notes must retain their meaning. Raw-to-parsed agreement checks extraction consistency; database content is not independent experimental confirmation.

All notebook cells passed syntax checking, and coverage checking found all 19 helper calls plus the CLI check for `main()`. Ten local parsing/formatting/assembly/writing cells and the final hash check passed with controlled validation inputs; these were not saved as notebook results. Setup worked from both supported directories. This tool session's first live API cell failed with a connection error, so the complete live notebook and graphical Jupyter execution still need the student's verification. No packages were installed, and the script/saved report stayed unchanged.

Possible instructor questions:

1. **Why display raw responses before readable fields?** They expose the exact source paths and values being extracted, making a mapping error visible.
2. **Does a difference from yesterday's saved report necessarily mean a bug?** No. Check current raw data, versions and assembly; online annotations may have changed.
3. **Why are Nucleotide and Ensembl coordinates compared separately?** The CDS locations refer to a transcript, whereas Ensembl gene positions refer to a genomic assembly.
4. **Why call the CLI again after checking individual functions?** It checks argument parsing and the complete retrieval-to-report path together, using fresh requests and temporary output.
5. **What should happen when an API request fails?** Investigate or retry the request; a network failure is not a legitimately missing biological field.

## Assignment 3 — first block: loading GSE40419

Purpose: load raw RPKM expression and the available clinical annotations without starting preprocessing or analysis. Inputs are the local RPKM text file and GSE40419 series matrix under `assignment_3/data/`. Run the first cell in `rnaseq.ipynb` with a pandas-capable kernel. It supports a working directory at the project root or `assignment_3/`.

1. `pd.read_csv()` reads the tab-separated matrix, treating chromosome labels as text. Separate gene/accession/chromosome/start/end/strand into `gene_annotations`; only sample expression columns enter `rpkm`. Repeated gene symbols and original row order are preserved, without aggregation.
2. Standard-library `csv.reader()` reads quoted GEO metadata rows. Collect full sample titles, GEO accessions and source tissue names; extract `gender` and `smoking_status` characteristics into sex and smoking columns. Replace the source's smoking `NA` marker with `pd.NA`, pandas' missing-value marker.
3. Create missing patient and driver-mutation columns. Patient mapping is not verified, so no patient IDs are assigned and LC_C/LC_S prefixes are retained in the sample index. Mutation information is unavailable; no stage or other variable is used as a substitute.
4. Assert unique sample names and identical sample sets, then use `.loc[rpkm.columns]` to reorder metadata. The original files list samples in different orders, so position-only matching would be incorrect. Check numeric, nonmissing, nonnegative expression values and display dimensions, tissue counts, missingness and the first metadata rows.

RPKM is an expression measure adjusted for transcript length and sequencing depth. It is not a raw read count. Gene-annotation coordinates describe genomic features; they must not be treated as expression samples. This block leaves RPKM unchanged: no log transformation, filtering, scaling or imputation of expression values is performed.

The student's approved limitation permits continuing without mutations, using available clinical variables and later testing tumour clusters against smoking status. It is not confirmed teacher approval. The teacher-facing README must document missing mutation annotations and omitted mutation-based analyses at the end. Later analyses still require their own instructions/approval.

Verification with the existing Conda pandas 3.0.6 environment confirmed 36,742 expression rows, 164 samples, six annotation columns, exact name alignment, 87 tumour/77 normal labels and seven missing smoking values. All patient and mutation values remain missing. LC_C1 and LC_S1 retain separate sample identities; repeated gene labels are kept. Source-file hashes were unchanged. Pandas installation in `.venv` has not been performed; environment approval is separate. Patient pairing and cross-dataset compatibility remain unresolved, and no later step has started.

Possible instructor questions:

1. **Why align metadata by name?** The same samples occur in different orders in the two files. Name-based reordering prevents assigning another sample's clinical labels.
2. **Why separate the first six columns?** They describe genes and genomic locations, not RPKM measurements, and would corrupt expression calculations.
3. **Why not strip every sample down to its number?** LC_C1 and LC_S1 would collapse into one identifier without evidence they are the same patient.
4. **Does a missing mutation value mean wild type?** No. It means no mutation annotation is available from the supplied metadata.

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

To inspect the loaded genome again from `assignment_1/` in PowerShell:

```powershell
..\.venv\Scripts\python.exe -B -c "from orf_finder import load_genome; record = load_genome('genome.fasta'); print(record.id, len(record.seq))"
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
..\.venv\Scripts\python.exe -B orf_finder.py genome.fasta --min-len 100 --out orfs.csv
```

`--min-len` defaults to 100 and must be a positive integer. `--out` defaults to `orfs.csv`. The minimum includes both ATG and STOP, as agreed in S0. `-B` is a Python option that avoids writing bytecode cache files; it does not change the analysis.

The CSV columns are `frame,start,end,length,protein`, with the same conventions as A1.2. A run with no qualifying ORFs still writes the header. Existing output files are overwritten on successful reruns, so the assignment's example command can be repeated directly.

### Steps and important functions

1. `main()` creates an `argparse.ArgumentParser` to read the FASTA path, `--min-len` and `--out`. Noninteger or nonpositive thresholds are rejected.
2. `load_genome()` reads and validates the sequence, and `find_orfs()` discovers complete ORFs using the unchanged six-frame algorithm.
3. Keep each ORF whose nucleotide `length >= min_len`. Equality is included: a 9-nt ORF passes a threshold of 9.
4. Open the output file with mode `w`, which creates it or replaces existing contents. UTF-8 specifies text encoding; `newline=""` lets the CSV writer manage line endings correctly on Windows.
5. `csv.DictWriter` maps each dictionary to the specified column order, writes the header and then the retained ORFs. Print the number of rows written.
6. Catch likely file and sequence errors and exit with a short message and status 2. A successful run exits with status 0. The `if __name__ == "__main__"` guard runs the CLI only when the file is executed directly, so importing the functions does not run the analysis.

This uses only standard-library modules in addition to the already required Biopython. A length threshold reduces short candidates but is not a statistical test or proof that the remaining ORFs are biological genes. Complete ORFs have lengths divisible by three, so a 100-nt threshold cannot admit a 100-nt or 101-nt ORF; the next possible length is 102 nt.

### Verification and results

Temporary, hand-checkable FASTA inputs tested ORFs of 6, 9 and 12 nucleotides at thresholds 1, 6, 7, 9, 10, 12 and 13. The checks confirmed inclusive filtering and a header-only CSV when nothing qualifies. A separate default-argument check excluded a 99-nt ORF and retained a 102-nt ORF in `orfs.csv`. A lowercase reverse-strand example retained frame `-3`, coordinates 3-11 and protein `MK` after CSV serialization.

Error checks passed for zero, negative and noninteger thresholds; missing FASTA; no records; an empty sequence; multiple records; ambiguous DNA; and a missing output directory. After the approved overwrite change, repeated runs correctly replaced existing temporary output files, while invalid input was rejected before opening them. `--help` listed the required options.

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
..\.venv\Scripts\python.exe -B kmer_freq.py genome.fasta -k 3 --out kmers_k3.csv
```

`count_kmers(sequence, k)` returns a list of `(kmer, count, frequency)` tuples. `main()` writes those rows with the CSV header `kmer,count,frequency`. Only observed k-mers are reported. Existing output files are overwritten on successful reruns.

### Steps and important functions

1. `main()` reads CLI arguments using `argparse`. It imports `load_genome()` from `orf_finder.py` to reuse the same Biopython parsing, uppercase conversion and input validation. Importing that module does not run the ORF search because of its main guard.
2. `count_kmers()` converts the `Seq` to a string and calculates `total = N - k + 1`, where `N` is the genome length. Reject `k < 1` and `k > N`.
3. Visit start positions 0 through `N - k` and take `sequence[start:start + k]`. Advancing by one base includes overlapping windows. Do not add reverse-complement counts or join the end of the sequence back to its beginning.
4. Start with an empty `collections.Counter`. In an explicit loop, extract each window into `kmer` and increment `counts[kmer]`. Counter gives an unseen fragment an initial count of zero.
5. The named `sorting_rule()` unpacks each `(kmer, count)` pair and returns `(-count, kmer)`. Negative counts put frequent fragments first, and the string breaks equal-count ties alphabetically. `sorted_counts` holds the ordered entries.
6. In a separate loop over `sorted_counts`, calculate `frequency = count / total` and append `(kmer, count, frequency)` to `rows`. The denominator is the number of windows, not the number of distinct k-mers. Return the completed list after the loop.
7. `csv.writer` writes the header and tuple rows to a UTF-8 file. Mode `w` creates the file or replaces existing contents; `newline=""` lets the CSV writer handle line endings. Likely file/sequence errors produce a short message and exit status 2.

For `AAACAA` with `k=2`, the five windows are `AA`, `AA`, `AC`, `CA`, `AA`. Therefore the sorted rows are `AA: 3, 0.6`, `AC: 1, 0.2` and `CA: 1, 0.2`. Counting nonoverlapping fragments would miss occurrences and give the wrong denominator.

### Verification and results

Hand-counted checks passed for the example above, repeated overlapping `AA` in `AAAA`, and alphabetical ties in `ATATA`. Boundary checks covered `k=1` (individual bases) and `k=N` (one whole-sequence window). An `ATG` input at `k=3` produced only `ATG`, confirming that reverse-complement `CAT` was not added. Lowercase multiline FASTA produced the expected CSV after normalization.

CLI checks passed for missing required options, noninteger `k`, zero, negative and oversized `k`, missing/empty/multiple-record FASTA, an empty sequence, ambiguous DNA, a missing output directory and `--help`. After the approved overwrite change, repeated runs correctly replaced existing temporary outputs, while invalid input was rejected before opening them.

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

## A2.1 — NCBI Gene lookup

### Purpose and input/output

Retrieve information about one gene in one organism using NCBI's web API (an interface that returns structured data to a program). From `assignment_2/`, run:

```powershell
..\.venv\Scripts\python.exe genereport.py APOE --organism "Homo sapiens" --out reports/
```

During A2.1–A2.4, database sections were added one at a time and printed to the terminal for verification. A2.5 now assembles those same fields into one Markdown report in `--out`; the terminal prints its path. The earlier stage explanations describe the retrieval logic, which remains unchanged.

### Steps and important functions

1. `main()` uses `argparse` to read the symbol, scientific organism name and future output directory. Empty symbols and organism names are rejected.
2. `get_ncbi_gene()` calls `find_gene_id()`, then `fetch_gene_summary()`, then extracts the report fields. This makes the main sequence of operations visible without mixing it with all the error checks.
3. `find_gene_id()` searches ESearch with the symbol restricted to the Gene Name field and the organism restricted to the Organism field. It requires exactly one matching Gene ID, the database's identifier for that gene. Zero or multiple matches produce a short error.
4. `fetch_gene_summary()` passes that ID to ESummary, then checks the returned Gene ID, official/preferred symbol and scientific organism name. A search can match a synonym or a related organism; an identifier alone does not prove the requested match. Our simple selection convention rejects mismatches rather than guessing. Both search and summary retrieval use `get_ncbi_json()`: `requests.get()` with a 30-second timeout, HTTP status checking and JSON decoding. HTTP is the web request protocol; JSON is a structured text format.
5. Back in `get_ncbi_gene()`, extract `nomenclaturename` (falling back to `description`), `chromosome`, `maplocation`, `otheraliases`, `summary` and `mim`. Convert the list of MIM identifiers to readable text. Optional absent or empty fields become `Unavailable`.
6. Return a dictionary of readable field names and values; `main()` prints each field. Nothing is hardcoded specifically for APOE.

### Biological meaning and field mapping

APOE is a gene symbol; Homo sapiens identifies the species. The student's real NCBI response identifies this combination as Gene ID 348, full name apolipoprotein E, chromosome 19, cytogenetic location 19q13.32 and MIM reference 107741. A cytogenetic location names a chromosome band, not a nucleotide interval. Aliases are alternative names for the same gene. MIM identifiers link to records in Online Mendelian Inheritance in Man (OMIM); a reference alone does not establish a diagnosis or a causal interpretation. The summary is NCBI's description, reproduced without rewriting it as our own biological conclusion.

### Verification and limitations

The real ESummary JSON supplied by the student confirmed the extraction fields and nested `organism.scientificname`. In-memory checks used those supplied values, including the complete summary, and controlled search responses. They verified field extraction, missing-field labels, rejection of ambiguous searches and identity mismatches, query parameters, HTTP status checks, timeout propagation, API errors and CLI argument errors. These checks do not establish that the live search succeeds.

`requests` 2.34.2 imports successfully in `.venv` and is pinned in `requirements.txt`. A live CLI attempt from the assignment directory exited 2 with `Error: NCBI request failed (ConnectionError).` because network connections were blocked in the tool session. No report or output directory was created. The student subsequently ran the full CLI successfully and supplied its output: Gene ID 348, human APOE, with the required fields matching the earlier real response. This verifies the live search-to-summary flow for the example. `pip check` could not start because of a Windows folder-access restriction in this session.

The student then approved separating search, summary retrieval and field extraction for readability. Before-and-after comparisons passed for 21 controlled cases, checking identical returned values and ordering, exception types/messages and request arguments. Successful CLI output and failed-request exit/message also matched. The checks remain explicit `if` statements, now placed with the operation they validate. The student accepted A2.1 including this revision.

### Possible instructor questions

1. **Why use ESearch before ESummary?** ESearch converts the requested symbol and organism into a Gene ID; ESummary retrieves the information for that ID.
2. **Why check the returned symbol and organism?** Searches can match aliases or related organisms. The checks prevent silently reporting another gene or species.
3. **Why use `Unavailable` for absent fields?** Database records differ in annotation coverage; missing information must not be replaced with invented values.
4. **Why use `requests` instead of Biopython here?** Assignment 2 explicitly requires `requests` for database HTTP access and forbids Biopython and ready-made API wrappers.
5. **Does a successful offline check prove the API works?** No. It checks our code against supplied data and controlled cases. A separate live run checks the actual search, service availability and current response structure.

## A2.2 — linked RefSeq mRNA annotations

### Purpose and input/output

Use the Gene ID from A2.1 to obtain all linked RefSeq mRNAs and their coding annotations. An mRNA is a transcript that can be translated into protein; different transcripts of the same gene can have different lengths or coding regions. The student approved including all linked mRNA variants rather than arbitrarily selecting one.

Run the same command from `assignment_2/`:

```powershell
..\.venv\Scripts\python.exe genereport.py APOE --organism "Homo sapiens" --out reports/
```

The Nucleotide section contains each accession/version, nucleotide length, modification date and every CDS with its `protein_id` and `transl_table`. At A2.2 it was printed to the terminal; A2.5 writes it into the Markdown report.

### Steps and important functions

1. `find_refseq_rna_ids()` calls ELink with `dbfrom=gene`, `db=nuccore` and `linkname=gene_nuccore_refseqrna`. Read `linksets`, then `linksetdbs`, then `links`. The student's real response for Gene ID 348 contained five RNA IDs. These numeric database IDs are used to fetch records; they are not transcript accession strings.
2. `get_ncbi_nucleotide()` sends those IDs together to EFetch with `rettype=gb` and `retmode=text`. GenBank is a text format containing record headers, a feature table and a sequence. Split records at their `//` terminators and check that all linked records arrived. Incomplete responses and duplicate mRNA records produce errors.
3. `parse_refseq_mrna()` reads accession/version from VERSION, length and modification date from LOCUS, and organism from ORGANISM. Retain `NM_` or `XM_` accessions whose LOCUS says mRNA; skip noncoding RNAs. `XM_` records are predicted models, not equivalent evidence to curated records. Check the GeneID cross-reference and organism before accepting an mRNA.
4. `parse_cds_features()` reads the fixed feature columns. Each new CDS starts a new annotation; wrapped location lines extend its location. Extract only `protein_id` and `transl_table` qualifiers inside that CDS. Other features and qualifiers are ignored. This is a small parser for the required fields, not a general GenBank parser.
5. Keep locations exactly in GenBank notation. Coordinates are 1-based and inclusive on the transcript, not genomic chromosome coordinates. `join(...)` describes multiple intervals, `complement(...)` specifies reverse orientation, and `<`/`>` mark uncertain boundaries. The code does not reinterpret these annotations.
6. Sort mRNAs by accession/version and display all of them. Missing qualifiers or dates are `Unavailable`; no linked mRNAs produce `RefSeq mRNA: Unavailable`. JSON calls and EFetch pause 0.34 seconds before requesting to stay within NCBI's default request rate. See [NCBI's E-utilities usage guidance](https://www.ncbi.nlm.nih.gov/books/NBK25497/).

### Biological meaning and supplied example

The supplied record is `NM_000041.4`: accession NM_000041, sequence version 4, length 1,166 nt and record modification date `09-SEP-2026`. That date differs from the comment's 2018 date when sequence version 4 replaced version 3; annotations can change without a new sequence version.

Its CDS is `70..1023` on this mRNA, with protein accession/version `NP_000032.1`. CDS means coding sequence, the region annotated for translation. `transl_table` identifies the genetic-code table when explicitly annotated. The supplied record omits that qualifier, so this report displays `Unavailable` rather than presenting an inferred default as a retrieved field. The program retrieves annotations; it does not discover ORFs or translate sequences in Assignment 2.

### Verification and limitations

Checks used the actual five-ID ELink response and the supplied APOE header/CDS fields. A reduced test record preserved those fields while omitting irrelevant references and annotations. Controlled examples checked CRLF line endings, explicit/missing qualifiers, multiple/missing CDSs, wrapped locations, join/complement/partial boundaries, predicted/noncoding RNA filtering, identity mismatches, all-record fetching and sorting, absent links and incomplete/missing/duplicate responses. The controlled five-record case was synthetic and does not establish the contents of the other four real APOE records. Integrated terminal output included every expected field from the supplied example.

The tool session's full CLI run exited 2 with a connection error. The student subsequently supplied successful live output containing all five linked mRNAs: `NM_000041.4`, `NM_001302688.2`, `NM_001302689.2`, `NM_001302690.2` and `NM_001302691.2`. Each showed length, modification date, CDS coordinates and protein ID; all five showed `transl_table: Unavailable`. The first record matched the earlier supplied GenBank annotation. This verifies the live example's retrieval-to-display flow; only the first record was independently compared with its supplied raw annotation. The student explicitly accepted A2.2. No generated report or extra dependency was added, and A2.3 has not started.

### Possible instructor questions

1. **Why link from Gene instead of searching all Nucleotide records for APOE?** ELink uses the database relationship to the verified Gene ID, avoiding unrelated keyword matches.
2. **Why include several mRNAs?** One gene can produce multiple transcript variants. All linked RefSeq mRNAs are included under our approved convention.
3. **What is the difference between record date and sequence version?** The record date reflects updates to the database entry; the accession's version changes when its sequence changes.
4. **Are CDS coordinates chromosome coordinates?** Here they refer to the retrieved transcript and remain in GenBank's 1-based inclusive notation.
5. **Why is `transl_table` unavailable in the example?** There is no explicit qualifier in the supplied record. The report distinguishes retrieved information from an assumption about the genetic code.

## A2.3 — UniProtKB protein annotations

### Purpose and input/output

Retrieve protein information corresponding to the verified gene and organism. The UniProt section follows the NCBI sections, first in the terminal during A2.3 verification and now in the A2.5 report. Each protein entry includes accession, review status, name, function, length in amino acids, ligands, cofactors, diseases, PDB structure IDs and keywords.

### Steps and important functions

1. `get_uniprot_proteins()` first asks `search_uniprot_entries()` for reviewed matches. If there are none, it asks for unreviewed matches. All exact matches in the selected group are retained and sorted by accession. This is our approved selection convention, not a requirement that the PDF specifies.
2. `search_uniprot_entries()` calls the required UniProt search endpoint through `requests`. Its query uses `gene_exact`, the quoted scientific `organism_name` and `reviewed:true` or `reviewed:false`. The `fields` parameter requests only the annotations needed for this stage and identifiers needed to check the match.
3. Verify the returned `genes[].geneName.value` and `organism.scientificName`; ignore records with different primary gene names or species. Follow the response's next-page link until no further page exists, so a page-size limit cannot silently omit matches. See the [UniProt API paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12230682/) for pagination and selected-field retrieval.
4. `parse_uniprot_entry()` reads accession, raw review status, recommended protein name (or a submitted name), sequence length, binding-site ligand names, PDB IDs and keyword names. Repeated binding sites for the same ligand produce one ligand name.
5. `parse_uniprot_comments()` handles FUNCTION, COFACTOR and DISEASE comments. `uniprot_texts()` reads their text-value lists. Preserve function text and disease descriptions/notes, including cautions about susceptibility. For cofactors, preserve names and notes. Absent annotations become `Unavailable` when report fields are assembled.
6. `main()` prints each selected protein's fields. No matching reviewed or unreviewed entry produces `UniProtKB entry: Unavailable`. HTTP/connection failures produce a short database request error, rather than being presented as missing biological annotation.

### Meaning of the fields

The supplied APOE response identifies reviewed protein `P02649`, Apolipoprotein E, with length 317 amino acids. An amino acid is a unit in a protein sequence; this length is not the nucleotide length of an mRNA. A reviewed Swiss-Prot entry has undergone manual curation; an unreviewed TrEMBL entry has not received that level of review. Review status does not make every statement experimental proof.

The supplied response lists heparin at two binding sites. The report lists the ligand name once and does not infer that it is an enzyme cofactor. A cofactor is a non-protein substance required for an enzyme's activity; APOE has no explicit COFACTOR comment in this response. Thus its cofactor field is `Unavailable`, which describes annotation availability rather than proving biological absence.

PDB IDs refer to deposited three-dimensional structures; the report retains all 29 returned IDs without claiming they are all full-length structures. Keywords are UniProt's annotation categories; all 30 supplied names are retained. Disease descriptions and notes remain source annotations. In particular, an association with susceptibility is not rewritten as certain disease causation.

### Verification and limitations

The student supplied the real reviewed human APOE JSON. Checks projected its basic fields, heparin binding sites, 29 PDB IDs and 30 keywords and verified their extraction. Separate checks preserved the supplied microbial function text and complete Alzheimer disease description/note verbatim. Controlled cases covered additional function comments, disease notes without a named disease, missing fields, submitted protein names, cofactor names/notes, pagination, mismatched genes/species, multiple reviewed matches, unreviewed fallback and no matches. Those controlled records are test examples, not database discoveries.

Because APOE lacks a cofactor comment, its structure was verified using UniProt's own [cofactor JSON tests](https://github.com/ebi-uniprot/uniprot-core/blob/main/json-parser/src/test/java/org/uniprot/core/json/parser/uniprot/comment/CofactorCommentTest.java). Submitted-name fields were checked against its [protein-description JSON tests](https://github.com/ebi-uniprot/uniprot-core/blob/main/json-parser/src/test/java/org/uniprot/core/json/parser/uniprot/ProteinDescriptionTest.java). The tool session's live run was blocked, but the student subsequently supplied successful CLI output matching the real APOE JSON: reviewed P02649, 317 aa, both function annotations, heparin, four disease annotations with notes, all 29 PDB IDs, all 30 keywords and no cofactor annotation. The student explicitly accepted A2.3. Ensembl and report assembly have not started.

### Possible instructor questions

1. **Why prefer reviewed entries?** They have been manually curated. Our approved rule uses unreviewed matches when no exact reviewed match exists, with status still shown.
2. **Why check the returned gene and organism after searching?** Query results must identify the requested primary gene name and exact scientific species, avoiding an unintended match.
3. **Why follow pagination?** A search response can contain only one page. Following next-page links retrieves all selected matches.
4. **Is a ligand necessarily a cofactor?** No. Binding a molecule does not establish a requirement for enzyme activity; these annotations are extracted separately.
5. **Why retain disease notes?** They can distinguish causal variants from susceptibility associations. Dropping them could change the biological meaning of the source annotation.

## A2.4 — Ensembl gene lookup

### Purpose and input/output

Retrieve the Ensembl gene annotation for the same symbol and organism. The Ensembl section contains gene ID, assembly, sequence region, start/end, strand, transcript count, canonical transcript and biotype. A2.4 printed it for verification; A2.5 includes it in the report.

### Steps and important functions

1. `fetch_ensembl_gene()` converts the scientific organism name to Ensembl's species form: `Homo sapiens` becomes `homo_sapiens`. URL-encode the species and symbol as path components.
2. Request `/lookup/symbol/{species}/{symbol}` with `expand=1` to include connected transcripts. The `content-type=application/json` query parameter requests JSON, matching the student's successful calls. Connection timeout is 10 seconds; read timeout is 90 seconds because an earlier 30-second read timed out. A timeout reports a request failure, not an unavailable annotation.
3. Check HTTP status, decode JSON and reject an API error or a returned object with the wrong type, gene symbol or species. A valid response must identify the requested gene.
4. `parse_ensembl_gene()` counts the `Transcript` list. A missing list is `Unavailable`, while an explicitly empty list is zero. Count every returned transcript, without filtering by biotype or selecting only the canonical transcript.
5. Preserve `canonical_transcript` exactly, including its version. Transcript objects provide separate stable IDs and version values; compare the canonical ID without its version suffix with the returned transcript IDs. Report an error if it is absent from a supplied list.
6. Extract the remaining fields and print them through `main()`. Keep start/end unchanged on either strand. The assembly and strand accompany the coordinates so their reference and orientation are explicit. Missing optional fields become `Unavailable`.

The endpoint and expanded objects are documented in [Ensembl's symbol lookup documentation](https://rest.ensembl.org/documentation/info/symbol_lookup). Ensembl uses [1-based inclusive coordinates](https://sep2025.archive.ensembl.org/info/docs/api/general_instructions.html); no conversion is needed for this report.

### Meaning of the supplied APOE result

The supplied response reports `ENSG00000130203`, assembly GRCh38, chromosome/sequence region 19, genomic positions 44903787 through 44909396 and strand +1. These are genomic positions, unlike the transcript-relative CDS positions in the Nucleotide section. Strand +1 means forward orientation; -1 would mean reverse orientation while start remains the lower genomic position.

The expanded response contains 39 transcripts and names `ENST00000252486.9` as canonical. A canonical transcript is Ensembl's representative annotation; this does not mean other transcripts do not exist or lack biological importance. `protein_coding` is the gene's biotype, its annotation category. It does not mean all 39 transcript annotations must encode distinct proteins.

The Ensembl transcript count differs from the five linked NCBI RefSeq mRNAs. The report preserves each database's annotation scope; it does not assume their transcript sets or gene boundaries must be identical.

### Verification and limitations

Checked extraction against the supplied real gene fields and all 39 transcript IDs, including the canonical transcript. Controlled examples checked missing versus empty lists, a malformed list, an absent canonical ID, gene/species/type mismatches and preservation of genomic coordinates on a reverse strand. Mocked HTTP checks verified URL construction, expanded JSON parameters, timeouts, API/HTTP failures, integrated terminal output and a short read-timeout error. Syntax checks passed.

An actual call through the new function failed with a connection error in this tool session. The student subsequently supplied successful updated CLI output matching all expected Ensembl fields, including GRCh38 coordinates, 39 transcripts and canonical transcript ENST00000252486.9. This verifies the live example's retrieval-to-display flow. The student explicitly accepted A2.4. No report file or new dependency was added; A2.5 has not started.

### Possible instructor questions

1. **Why request `expand=1`?** It includes the connected transcript objects needed to calculate the transcript count.
2. **Why print the assembly?** Genomic coordinates depend on the reference genome assembly, so numbers without that reference can be ambiguous.
3. **Why does the canonical ID include `.9`?** It is a transcript version suffix. The report preserves it, while membership checking compares the stable ID part.
4. **Why are there 39 Ensembl transcripts but five RefSeq mRNAs?** They are different database annotation sets, and the Ensembl count includes every returned transcript rather than only RefSeq mRNAs.
5. **Does a missing transcript list mean zero transcripts?** No. Missing data is `Unavailable`; zero is used only when the API explicitly returns an empty list.

## A2.5 — assembling the Markdown report

### Purpose and input/output

Combine the four verified database lookups into the assignment's human-readable report. Markdown is the single report format agreed in S0. From `assignment_2/`, run:

```powershell
..\.venv\Scripts\python.exe genereport.py APOE --organism "Homo sapiens" --out reports/
```

On success, the output is `reports/APOE_Homo_sapiens.md`, and the terminal prints the written path. A relative `--out` is resolved from the current working directory. The output directory is created if necessary, and a successful rerun replaces the same report, as explicitly approved. Root `README.md` gives teacher-facing instructions.

### Steps and important functions

1. `main()` performs the NCBI Gene, Nucleotide, UniProtKB and Ensembl lookups. A request or parsing failure stops the program before the report is opened. Retrieval and record-selection rules remain those approved in A2.1–A2.4.
2. `build_report()` assembles four database sections. It includes all extracted Gene fields, every retained mRNA and its CDS annotations, every selected protein's fields and all extracted Ensembl fields. Missing record sets receive explicit `Unavailable` text.
3. `format_markdown_fields()` formats each field as a labelled paragraph. Existing annotation newlines become blank-line paragraph breaks, preserving separate function and disease descriptions without truncating them into table cells.
4. `markdown_text()` escapes punctuation that Markdown could interpret as formatting. For example, the raw file contains `APOE\*4`; a Markdown viewer displays `APOE*4`. This changes presentation syntax, not the biological wording. The same applies to underscores in accession strings.
5. `write_report()` makes the output directory, constructs a filename from symbol and organism and writes the complete report using UTF-8. Characters unsuitable for this simple filename become underscores. The example produces `APOE_Homo_sapiens.md`. No CSV, HTML, JSON report or separate per-database report is generated.
6. `main()` prints a success message only after writing finishes. Database errors, parsing errors and output-path errors produce short messages and exit code 2. A failed lookup preserves any previous report because writing has not started; an interrupted or failed disk write is not guaranteed to preserve an old file.

### Interpretation and requirement coverage

The report presents retrieved annotations rather than a new biological analysis. It preserves missing-field labels, transcript-versus-genomic coordinate conventions, the all-mRNA selection rule and the reviewed-protein preference. Disease susceptibility wording remains intact; a missing annotation is not treated as proof of biological absence. Database annotations can change between runs.

| Required database content | Report location |
| --- | --- |
| Gene ID, full name, chromosome/location, aliases, summary, MIM references | NCBI Gene fields |
| RefSeq mRNA accession/version, length, modification date, CDS location, protein ID, translation table | Each NCBI Nucleotide transcript and CDS block |
| Protein accession, review status, name, function, length, ligands/cofactors, diseases, PDB IDs, keywords | Each UniProtKB protein block |
| Ensembl gene ID, coordinates, transcript count, canonical transcript, biotype | Ensembl fields, with assembly and strand |

### Verification and limitations

Temporary checks verified every expected section, field label and value, all five previously supplied mRNA rows, missing record/CDS output, long and multiple annotation paragraphs and literal Markdown punctuation. Separate checks verified nested directory creation, exact UTF-8 text, replacement/truncation on reruns, a filename remaining inside its output directory, invalid output paths and the complete CLI with controlled lookup results. A forced last-database timeout preserved a pre-existing test report. All test outputs were temporary; controlled descriptions were not saved as real results.

Syntax and CLI help checks passed. The tool session's live run was blocked, but the student subsequently generated the real report locally and accepted it. Inspection of `assignment_2/reports/APOE_Homo_sapiens.md` confirmed all four sections and required fields, all five supplied mRNA/CDS/protein rows, disease context and both function annotations, 29 PDB IDs, 30 keywords, expected missing-field labels and the Ensembl values including 39 transcripts. These checks compared the real saved report with the supplied API/CLI results; it was not generated from controlled test descriptions. A2.5 is accepted. Assignment 3 has not started.

### Possible instructor questions

1. **Why write only after all lookups succeed?** It prevents a failed lookup from replacing an existing report with incomplete data.
2. **Why choose Markdown?** It is the agreed single human-readable format, works as plain text and supports clear sections and paragraphs without another dependency.
3. **Why escape `APOE*4` and accession underscores?** Markdown punctuation can change displayed formatting. Escaping preserves the literal annotation when rendered.
4. **What happens on a successful rerun?** The same gene/organism filename is replaced with the newly retrieved report.
5. **Does a successful controlled CLI check prove the final report is current?** No. It verifies assembly and writing; the real combined API run and generated file still need checking.

## Assignment 3, Part 4 - differential expression and Lab 3 comparison

**Purpose and inputs/outputs.** Compare tumour with normal expression while controlling false discoveries. Inputs are filtered log2(RPKM + 1) GSE40419 records, candidate pairs, log2 RMA GSE10072 values and GPL96 annotation. Outputs are result tables, volcano plots and shared-gene overlap/direction.

**Steps and functions.**
1. Merge GSE40419 tumour and normal metadata on the full LC_C/LC_S stem. Retain 77 candidate pairs in matching column order. Clinical agreement supports consistency, not independently verified identity.
2. scipy.stats.ttest_rel tests each record's mean within-patient difference. Tumour minus normal is calculated on the transformed scale; the +1 pseudocount means it is not an ordinary log2 fold change.
3. multipletests(method="fdr_bh") applies Benjamini-Hochberg correction to valid p-values. Undefined tests stay unavailable. FDR < 0.05 marks significance; Lab 3 Exercise 6 selection also requires absolute difference >= 1.
4. For GSE10072, remove AFFX controls and use all 58/49 samples with ttest_ind(equal_var=False), matching the lab. Do not transform RMA values again. This method ignores dependence between matched tissues.
5. Plot mean difference against -log10(FDR). Red/blue indicate selected higher/lower tumour expression. Grey means not selected, including significant records with small effects.
6. Load/cache GPL96 annotation. First-listed gene symbols follow the lab but are ambiguous for multi-gene probes. Choose the highest-mean-expression record per gene independently of significance. Compare only shared represented genes and count direction agreement. Inherited FDRs remain record-level rather than gene-level.

**Meaning and assumptions.** Paired tests use within-patient differences and assume independent patients and approximately normal differences. Outliers can influence them. Welch's test permits unequal variances but assumes independent observations. FDR controls an expected proportion of false discoveries under its assumptions, not the probability that a particular gene is false. Platform, sample-size and test-design differences limit cross-dataset interpretation. External findings must not guide classifier training or tuning.

**Verification.** Local execution produced 27,190 GSE40419 tests, 20,408 FDR-significant and 2,239 selected records; GSE10072 produced 22,215 tests, 11,047 FDR-significant and 857 selected probes. Counts, FDR ranges and plotting were checked. Official annotation retrieval and complete overlap await local execution because tool network access is blocked; no overlap count is claimed.

**Possible instructor questions.**
1. Why retain pairing for GSE40419? The assignment requests it; valid pairing accounts for between-patient variation.
2. Why Welch for GSE10072? It reproduces Lab 3, with imperfect independence because some tissues are matched.
3. Why distinguish record and gene counts? Several probes or expression records can represent one gene.
4. Why select representatives by mean expression? It avoids choosing records explicitly because they have smaller p-values.
5. Does overlap prove replication? No; annotation ambiguity, different measurements and statistical designs limit that interpretation.

### Part 4 code after simplification

make_de_table(index, difference, p_values) constructs the same result columns for both datasets, corrects only defined p-values with BH, and marks significance and effect-size selection. The paired and Welch tests remain explicit in their own cells. Metadata uses a field-name dictionary to avoid repeated assignments. Gene records are sorted by mean expression, reduced to one record per gene, and combined with an inner join, which retains shared symbols only. These changes preserve the statistical results: exact before/after comparisons passed on the local datasets. Saved outputs and earlier analysis cells were preserved.

The combined Part 4 comparison now reads the submitted data/GPL96_annotation.tsv locally and prints the shared universe, selected counts, overlap, and direction agreement. Verification with that file found 10,305 shared genes and 438 jointly selected representatives, all changing in the same direction. The results remain descriptive cross-platform agreement under the documented mapping and record-selection rules.

## Assignment 3, Part 5 - classification and external ROC diagnosis

**Purpose and inputs/outputs.** Predict tumour versus normal from GSE40419, using patient-grouped validation and external GSE10072 testing. Outputs are fold AUC, held-out predictions, model coefficients and external rankings/class predictions.

**Steps.** Reload unfiltered GSE40419; log-transform once, average records per gene and match measured symbols to annotated GSE10072. These fixed operations do not use DE significance or external labels. Keep full-prefix candidate patients together. Within each outer training fold, StandardScaler and logistic regression are fitted in a pipeline; inner grouped GridSearchCV selects C. Outer patients are excluded from preprocessing and tuning. Tune/refit the final model on GSE40419 only. Absolute coefficients describe importance for training-standardised genes, with correlated and low-variation genes requiring caution.

**Diagnosed artifact.** Probabilities are a sigmoid of decision scores. Large scores can round all probabilities to 1.0, destroying ranking numerically. decision_function retains ordering and is appropriate for ROC AUC. This repairs metric calculation, not calibration or platform shift.

**Verification and limits.** Full execution reproduced mean outer AUC 0.998611 and complete held-out coverage. Both CV levels had zero patient overlap; metadata features and exact duplicate profiles were absent. One shuffled-label diagnostic averaged 0.461955, not a formal significance test. External score AUC is 0.993315, yet all 107 default predictions are tumour. Microarray values can lie thousands of training standard deviations away, so calibrated probabilities and default-threshold performance do not transfer. Candidate identities remain independently unverified. No external result was used to change the model.

**Instructor questions.**
1. Why nested grouped CV? Inner folds tune; outer folds evaluate on excluded candidate patients.
2. Can AUC be high with all predictions tumour? Yes; AUC assesses rankings across thresholds, not one fixed threshold.
3. Why use decision scores? They preserve ordering when probabilities round to 0 or 1.
4. Is platform shift fixed? No; only the numerical loss of ranking is corrected.
5. Are top coefficients causal biomarkers? No; they are model-dependent weights affected by scaling and correlated genes.
