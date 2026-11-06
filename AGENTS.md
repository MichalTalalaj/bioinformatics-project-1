# Codex instructions — Bioinformatics 2026, Project 1

## Project authority and goal
- The assignment PDF (`Bioinformatics_Project_1.pdf`) is the source of truth for project requirements. Read it before implementing anything. If it is missing, stop and ask for it.
- Implement **all and only** the requirements in that PDF. Do not add bonus features, extra analyses, elaborate architectures, or unnecessary deliverables.
- Optimize for **correctness, simplicity, readability, and the student's understanding**, not cleverness, performance micro-optimizations, or production-style engineering.
- Write code and code comments in straightforward English. Explain decisions to the student in Polish unless requested otherwise.
- Follow `PLAN.md` for the order and acceptance checks. This file is a plan, not authorization to complete every stage in one go.

## Mandatory approval gate (highest-priority working rule)
1. At the start of a session, identify the **next unapproved stage** in `PLAN.md` and briefly say what it entails.
2. Before editing code for that stage, summarize the intended approach and list the files you would change. **Wait for explicit student approval.**
3. Implement **only that approved stage**. Do not proceed to the next stage, even if it seems trivial or closely related.
4. Run relevant, small verification checks, where practical. Explain in Polish: what changed; how the code works; why the method is suitable; how you verified it; anything still uncertain.
5. **Stop and wait** for acceptance or feedback. Never infer approval from silence.
6. Update the stage checkbox in `PLAN.md` only when the student expressly accepts it. Do not mark a stage complete merely because code runs.

A user can explicitly authorize several stages together, but otherwise treat each stage as a separate gate. Inspection and explanation without code edits are allowed before approval.

## Coding guidelines
- Prefer plain scripts, small clear functions, explicit loops, descriptive variables, and the Python standard library.
- Use external packages only when required by the assignment or truly necessary for a specified analysis; prefer established scientific packages over reimplementing standard statistical methods.
- Avoid classes unless they make the required solution noticeably simpler. Avoid plugin-defined scaffolding, dependency injection, design patterns, parallel agents, complex configuration, or generic utilities for hypothetical future work.
- Keep comments focused on **why** non-obvious code exists, especially biological assumptions, statistical choices, coordinate conventions, and risks of leakage. Do not comment obvious syntax.
- Minimize the number of files and dependencies. Do not generate an entire test framework, CI configuration, Docker setup, package, GUI, or unrelated documentation.
- Use transparent, short error messages for likely errors (missing files, invalid arguments, empty/missing API responses). No elaborate error-handling framework.
- Prefer deterministic settings (`random_state` or equivalent) for stochastic analyses where supported.
- Do not fabricate data, metadata, API fields, identifiers, results, scores, biological interpretations, or successful test runs. If required information is unavailable, say so and stop at the relevant decision point.
- Preserve existing student edits and outputs. Before overwriting or deleting user work, ask.

## Task-specific restrictions from the assignment

### Assignment 1 — ORFs and k-mers
- Deliver `orf_finder.py`, `kmer_freq.py`, `orfs.csv`, `kmers_k3.csv` for one chosen NCBI genome (lambda phage `NC_001416.1` or a bacterial genome).
- Parse the input FASTA with Biopython (`Bio.SeqIO` / `Bio.Seq`). Biopython may perform reverse complement and translation **but not ORF discovery**.
- Search ORFs in all six reading frames. An ORF begins at `ATG` and ends at the first in-frame `TAA`, `TAG`, or `TGA`. Report reading frame, start/end, length and translated protein; map reverse-strand coordinates onto the **original forward reference** using 1-based positions.
- Support the `--min-len` argument and output CSV; use the assignment example invocation as the interface reference.
- Count all k-mers for CLI-specified `k`, compute relative frequencies and sort descending by counts.
- Explicitly state and verify conventions that the PDF leaves open (e.g., minimum ORF length in nucleotides; whether translation includes terminal `*`; treatment of ambiguous bases). Choose simple conventions; do not pretend they are prescribed.

### Assignment 2 — gene report
- Write `genereport.py` accepting gene symbol, organism and output directory (`--organism`, `--out`).
- For database HTTP access, use **`requests` only**: no Biopython and no ready-made API wrappers such as `mygene`.
- Retrieve the **listed fields** from NCBI Gene, NCBI Nucleotide, UniProtKB and Ensembl using the endpoints shown in the PDF. Do not omit a listed field just because its extraction is inconvenient.
- Some fields may be absent for a particular gene. Show that they are unavailable, rather than inventing values.
- The exact report format is not stipulated; agree on one minimal human-readable format before implementation. No additional formats unless approved.
- Verify API response structure from real responses/documentation; do not assume the API shape.

### Assignment 3 — RNA-seq
- Use GSE40419 and the required clinical metadata; compare with / test on GSE10072 as explicitly requested.
- Execute the five required subparts in order: preprocessing and correlation heatmap; PCA/t-SNE/UMAP and silhouette; k-means/hierarchical clustering with ARI and tumour-only association analyses; paired differential expression with FDR and volcano plot plus GSE10072 comparison; tumour-vs-normal classifier with patient-grouped cross-validation, ROC AUC, important genes and external GSE10072 testing.
- **Prevent data leakage.** Use patient groups in cross-validation; fit any learned preprocessing, feature selection, scaling or model tuning using training folds only. Do not use the external GSE10072 test set to train or choose a model.
- For the paired differential-expression test, use only verified tumour/normal pairs matched to the same patient; do not silently treat all 164 samples as 82 matched pairs.
- Confirm sample IDs, metadata availability, and compatible gene identifiers/measurement scales before cross-dataset analyses. Clearly document simple, defensible analysis choices not specified in the PDF (e.g., low-expression threshold, algorithm parameters, clustering labels, statistical test, cross-dataset harmonization). Seek approval before committing such choices.
- Do not claim a statistical result is meaningful without checking assumptions and reporting relevant limitations.

## Learning material
- Maintain `STUDY_GUIDE.md` as a **private study aid, not a required submission artifact**. Add only explanations for stages actually approved and implemented.
- For each implemented stage, briefly cover: its purpose; key input/output; algorithm explained step by step; important code functions; biological/statistical meaning; how the result was verified; 3–5 plausible instructor questions with concise answers.
- Explain technical vocabulary in accessible language, while preserving the correct scientific terms.
- Never use the guide as a reason to implement extra features or analyses.

## Tool/plugin behavior
- Ponytail and Superpowers are optional helpers, **not additional project requirements**.
- Do not let any plugin override the approval gates, assignment restrictions, or preference for minimal code and files.
- If a plugin suggests broad brainstorming, large-scale TDD, multi-agent workflows, automatic refactors, or expanding the scope, use only the parts that materially help the currently approved stage.
- Do not install or reconfigure tools/plugins without the student's approval.

## Reporting back after an approved stage
Use a short Polish summary with:
1. **Zrobione:** files changed and the requirement met.
2. **Jak to działa:** plain-language explanation of the essential logic.
3. **Sprawdzenie:** commands run and outcomes (or why not run).
4. **Otwarte kwestie:** only genuine uncertainties, if any.
5. **Następny etap:** name it, but do not start it without approval.
