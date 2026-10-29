# External inputs: sources, axes and available preparation

No private cell table is distributed. The full configurations contain expected hashes, seed lists and explicit target filenames. Their aggregate null tables and synthetic fixtures are already included. Sources, conversion steps and remaining gaps follow; availability of original datasets does not by itself recreate historical scoring intermediates.

## Input contracts

- A score TSV has one row per unique `cell` ID and exact nonnegative integer `Neff`, unique integer `RNA_rank` in 1..N, and legal `DNA_rank` in 1..N. The GBML181 tie enumeration additionally needs original textual `RNA_score`; P2 uses its original RNA rank, so its historical score file does not need `RNA_score`. Cohort annotation scores additionally need `truth` (CC true, MC false) and `observed_source_features` (numeric count or missing).
- A threshold CSV has exactly ten rows `N,cutoff,maximum`, one for each cutoff 100,90,...,10. N equals the number of score IDs; maxima are exact integers. `data/bio06/all_thresholds.csv` additionally has `case,null_id`, and `data/cohort_null_maxima.csv` groups by N. The finite maxima are not theoretical cutoffs.
- A nuclear NPZ has string `ids` length N, Boolean `observable` length N in the same axis, and finite `X`, `original`, `alternative` matrices with only `observable=True` rows in that axis order. X columns for P2 are intercept, log1p accepted mitochondrial alignments, log1p nuclear UMI and observed source-feature sites; both outcome matrices have the same `names` string column axis (52 features in the historical P2 input). IDs must match the score IDs exactly. Do not label missing-observation cells as negative.
- A composition NPZ contains `ids,neff,rna_rank,dna_rank,rna_score,truth,observed_source_features,cc_order,anchor,background,untied_endpoints`. The first seven arrays are the parent cell axis; the three ordering arrays are zero-based indices into that axis, and `untied_endpoints` is a Boolean vector on that axis. See `drivers/cohort/prepare_split.py` for the original frozen selection and ID-digest checks.

## Shared conversions

From a privately held saved cell table and an already selected ten-row cutoff table, generate audit-style scores and thresholds. `--dna-column original_DNA_rank` and `--maximum-column random_maxN` are appropriate for the historical P2 global-rank and cutoff-bound tables; use `--include-rna-score` only if that input column is present and needed. Select the correct case's ten rows before invoking the converter:

```sh
python "$REPRO/drivers/prepare_inputs.py" score --cells "$PRIVATE/cells.tsv.gz" \
  --cutoffs "$PRIVATE/case_cutoff_rows.csv" --dna-column original_DNA_rank \
  --maximum-column random_maxN --output "$PRIVATE/audit-input"
```

P2's project-authored prepared NPZ can be converted without changing the nuclear row or feature axes:

```sh
python "$REPRO/drivers/prepare_inputs.py" nuclear --prepared "$PRIVATE/prepared.npz" \
  --score "$PRIVATE/audit-input/score.tsv" --output "$PRIVATE/audit-input/nuclear.npz"
```

The prepared NPZ must contain `barcodes,observable,covariates,original,alternative,names`. These commands check axis compatibility; users must separately check frozen file hashes and scientific provenance. Regeneration can serialize logically identical arrays to different bytes, so do not claim the historical file digest without checking it.

## B1 and B2 experimental benchmarks

The [MERCI v1.0.3 dataset, Zenodo 8065206](https://zenodo.org/records/8065206), provides four processed files per experiment under CC BY 4.0. Replace `bench2` below by `bench1` for the first batch. Preserve the provider's `covearge` spelling.

| Purpose | Provider filename | Axis/meaning |
|---|---|---|
| `--info` | `cell_information(bench2_scRNAseq).txt` | cell IDs, sample, experimental labels and mitochondrial read counts |
| `--expression` | `normalized_exp_matrix(bench2_scRNAseq).txt` | gene rows × cell columns; author-normalized, not raw UMIs; original wrapper maps periods to hyphens in column IDs and checks one-to-one identity |
| `--variants` | `mtSNV_profile(bench2_scRNAseq).txt` | cell, variant ID, position, ALT reads, total reads, allele fraction |
| `--coverage` | `MT_covearge_cells(bench2_scRNAseq).txt` | mitochondrial position rows × cell columns, integer coverage |

### Six-case B1/B2 inputs from the four processed files

Set `REPRO` to this directory and `DATA` to the directory holding the four provider files for the chosen batch. Obtain the upstream [MERCI checkout at commit 79e105aaf13159bee325c3ceceb6ee679e430341](https://github.com/shyhihihi/MERCI/tree/79e105aaf13159bee325c3ceceb6ee679e430341) separately under its academic-use terms. R 4.4.1 with Matrix, data.table, jsonlite, e1071 and proxy in an explicitly supplied user library was used for local checks. The installed Python package and a **new** private output directory are required. For example, use the four bench1 files for B1 frozen split 21:

```sh
python "$REPRO/drivers/benchmark_inputs/run.py" \
  --info "$DATA/cell_information(bench1_scRNAseq).txt" \
  --expression "$DATA/normalized_exp_matrix(bench1_scRNAseq).txt" \
  --variants "$DATA/mtSNV_profile(bench1_scRNAseq).txt" \
  --coverage "$DATA/MT_covearge_cells(bench1_scRNAseq).txt" \
  --upstream "$MERCI_CHECKOUT" --rscript "$RSCRIPT" --r-library "$R_USER_LIBRARY" \
  --experiment bench1 --split 21 --output "$PRIVATE/b1-split21-new"
```

For B2 use the four `bench2_scRNAseq` filenames with `--experiment bench2 --split 1`. Allowed original manifest selections are B1 1–20 (`development_reused`), B1 21–100 (`frozen_evaluation`), and B2 1–20 (`retrospective_second_experiment`); these are technical splits of the two original experiments, not independent patients. The original seeds are `202610260 + split` for B1 and `202610290 + split` for B2, using R Mersenne-Twister/Inversion/Rejection. The published B1 strict-only `drivers/b1/run.py` interface remains available for its original 21–100 range.

The six cases are `original_regular_mixed`, `original_regular_negative`, `original_loo_mixed`, `original_loo_negative`, `strict_regular_mixed` and `strict_regular_negative`. Mixed contains CC plus test MC, negative contains only test MC. Regular uses reference MC for enrichment comparison, background and RNA reference, while original regular visibility also includes candidates. Original LOO uses its own candidate set for feature visibility, comparison and background, and invokes the original leave-one-out RNA function separately for mixed and negative. Strict uses only donor plus reference for source selection and background, excluding test targets. Identical-reference regular RNA weights are computed once for mixed, then subset and reranked for each candidate axis; the prior first B2 negative check called the original regular function independently and found exact weights and ranks. Each case's feature selection and DNA scoring remain separate.

The original 13 mouse mitochondrial gene rows are retained in their frozen order. B2 expression column IDs map periods to hyphens before the one-to-one check. The input coverage must provide integer position-by-cell counts, and variant rows must include `Cell,ID,pos,s_reads,t_reads,AF` with the original variant IDs. Feature selection, depth 5, observability, Yates/BH and background formulas use the project engine; thresholds and RNG are not newly fitted. The original task's saved ten maxima come from `configs/cohort_formal_portable.json`. The ten nested/fixed-N composition conditions instead use the separate BIO-04 matched-N maxima in `data/cohort_null_maxima.csv`; the two tables are not interchanged.

The private `generated/` directory contains six `*_cells.tsv.gz` scores, six feature-selection RDS files, roles and a generation receipt. Score rows use the frozen candidate order and 13 fields: `cell,truth,Neff,DNA_rank,RNA_score,RNA_rank,observed_source_features,source_marker_observable,MT_reads,mean_MT_coverage,observed_MT_positions,nuclear_UMI,raw_top50`. `prepared.npz`, `cohort_manifest.json`, `cohort/` and `run_complete.json` are bound to the freshly generated score and NPZ bytes. The manifest preserves historical hashes separately and checks frozen ordered IDs, split-selection digests and gates. `BYTE_IDENTICAL` means hashes match; `DIFFERENT_BYTES_REQUIRES_COMPARISON` requires explicit value, order and parameter comparison before any equivalence claim. No historical sample is re-signed or reused.

For B1 split 21 and B2 split 1, local full-input executions generated all 12 scores, both composition NPZ files and 32 task/condition member files byte-identically to the historical inputs and outputs. The 13 score columns and ordering, roles, all six feature-selection objects per batch, and common aggregate fields were also compared; [aggregate validation](data/benchmark_inputs/validation.json) records their scope and hashes. Other 118 splits were checked for configuration selection only, not rerun. In these local runs B1 took 82.75 seconds and peaked at 3,629,972 KiB RSS; B2 took 30.92 seconds and peaked at 1,680,184 KiB, each under a 4 GiB address-space limit on one CPU with nice 10. Timings and memory depend on hardware, software and cache state.

The scripts neither download data nor install dependencies. First verify provider file integrity and licensing, and maintain sufficient memory for coverage parsing. Full generated score, role, NPZ and member files contain individual cells and must remain private. Checks of configuration coverage do not claim that every split has been run. The separately available `merci-rank-audit bench2` command remains a complete specific B2 split 1 original-regular negative wrapper, documented in `../docs/bench2.md`.

Small independent engineering checks:

```sh
python "$REPRO/fixtures/benchmark_inputs/check_config.py"
R_LIBS_USER="$R_USER_LIBRARY" "$RSCRIPT" --vanilla "$REPRO/fixtures/b1/check_coverage.R" \
  "$REPRO/drivers/b1/read_coverage.R" "$REPRO/fixtures/b1/coverage.tsv"
python "$REPRO/fixtures/b1/check_prepare.py" --output "$PRIVATE/synthetic-check"
```

For preexisting strict score tables, `drivers/cohort/prepare_split.py` can still regenerate a composition NPZ using the corresponding manifest and frozen seeds; its rebound manifest records the regenerated NPZ hash and retains its historical hash. The 51-case benchmark's 44 synthetic NPZ files are provided; each of seven real cases needs its frozen score and threshold bytes. Create the real NPZ and a locally rebound configuration with `drivers/benchmark/prepare_case.py --config ... --case ... --manifest ... --output-npz ... --output-config ...` after preparing the one-case score manifest. It rejects changed original score/threshold bytes or canonical arrays. The former post-hoc endpoint-versus-old-results/witness comparison is not supplied as a standalone public input chain.

## P2 HCC patient input

Use [GSE242889/GSM7774400](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7774400), provider file [GSM7774400_2T_C24.tar.gz](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7774nnn/GSM7774400/suppl/GSM7774400_2T_C24.tar.gz), and the corresponding raw run [SRR26035078](https://www.ebi.ac.uk/ena/browser/view/SRR26035078), not a different sample. The expression archive supplies barcode, gene and matrix axes. Historical IDs use `GSM7774400:` followed by the exact 24-base GEXSCOPE barcode; preserve exact IDs rather than trimming or guessing a join. The original work retained 3,265 hepatic-compatible scored targets and 2,741 observable nuclear rows. The 52 fixed readouts, two representations, original source-feature observability, covariates and original maxima must be prepared on their original axes before the conversion above. The full frozen input identity is in `configs/nuclear_formal_portable.json`; the P2 benchmark and BIO-06 configs bind the corresponding score and threshold digests.

**Gap:** the public archive and SRA accessions do not themselves supply the frozen role assignments, accepted mitochondrial alignments/UMI consensus, historical source-feature scoring, fixed 52-readout prepared NPZ or serialized score/threshold bytes. The public export here supplies the input converter and analysis drivers, not a standalone complete raw-to-P2 reconstruction. Obtain or reconstruct those project-specific intermediates lawfully before claiming exact reproduction. Nuclear coefficients in this submission use the original threshold; changing the threshold requires a separate conditional calculation.

## GBML181 single-nucleus input

Use [GSE230389/GSM7221332](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7221332) provider files `GSM7221332_GBML181_barcodes.tsv.gz`, `GSM7221332_GBML181_features.tsv.gz`, `GSM7221332_GBML181_matrix.mtx.gz`, and native SRA runs [SRR24283915](https://www.ebi.ac.uk/ena/browser/view/SRR24283915) and [SRR24283916](https://www.ebi.ac.uk/ena/browser/view/SRR24283916). Matrix rows correspond to the features file and columns to barcodes; retain Gene Expression rows (36,601) separately from ATAC peaks (57,101). The two runs are one sample. Preserve each RNA `RNA_score` as original text when reconstructing ties. The prior analysis used 1,343 scored targets and ten matched-N thresholds.

**Gap:** RNA-QC/working role assignments, competitive alignment, original cell/UMI consensus, source-feature selection and score TSV are project-derived intermediates and are not supplied as a portable complete raw-to-score workflow in this bundle. The original study's individual-cell malignant annotations were not obtained. The driver starts with the frozen score and thresholds; it does not imply exact recovery from the public accessions alone. Do not redistribute those private cell files or the upstream academic-only code in this project.
