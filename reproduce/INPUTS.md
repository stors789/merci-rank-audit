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

A complete existing preparation entry for the **B2 split 1 original-regular negative case** is the installed `merci-rank-audit bench2` command, documented in `../docs/bench2.md`. It requires the four files and an independently obtained [MERCI checkout at commit 79e105aaf13159bee325c3ceceb6ee679e430341](https://github.com/shyhihihi/MERCI/tree/79e105aaf13159bee325c3ceceb6ee679e430341) under its separate academic-use license, plus R and the listed packages. The command writes `score.tsv` and `thresholds.csv`; it is not a general B1/B2 generator.

For the 720 annotation cases, bind each historical `*_cells.tsv.gz` file at the path specified in `configs/cohort_formal_portable.json`; that manifest carries its original bytes and ID digest plus the ten maxima. For each of the 120 strict-chain composition splits, use its original frozen mixed and negative score tables and regenerate the selection using its saved seeds and ordered-ID digests:

```sh
python "$REPRO/drivers/cohort/prepare_split.py" --manifest "$REPRO/configs/cohort_formal_portable.json" \
  --experiment bench1 --split 1 --mixed "$PRIVATE/strict_regular_mixed_cells.tsv.gz" \
  --negative "$PRIVATE/strict_regular_negative_cells.tsv.gz" --output "$PRIVATE/prepared.npz" \
  --output-manifest "$PRIVATE/cohort-rebound.json"
```

The optional rebound manifest records the regenerated NPZ byte hash while retaining its historical hash, after verifying original score bytes, ID axes and frozen selection digests. It contains local absolute paths and must remain private. To prepare further splits, pass the latest rebound manifest to the next invocation and write a new one. The 51-case benchmark's 44 synthetic NPZ files are provided; each of seven real cases needs its frozen score and threshold bytes. Create the real NPZ and a locally rebound configuration with `drivers/benchmark/prepare_case.py --config ... --case ... --manifest ... --output-npz ... --output-config ...` after preparing the one-case score manifest. It rejects changed original score/threshold bytes or canonical arrays.

**Gap:** this package does not include a portable reproduction of every B1/B2 split, every original/strict/LOO score chain, or the former post-hoc endpoint-versus-old-results/witness comparison. The single verified B2 wrapper does not establish that all 720 historical cell-score files can be regenerated from the four public files with this release alone. Exact original frozen gzip SHA-256 values may depend on serialization; do not silently substitute different bytes.

## P2 HCC patient input

Use [GSE242889/GSM7774400](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7774400), provider file [GSM7774400_2T_C24.tar.gz](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7774nnn/GSM7774400/suppl/GSM7774400_2T_C24.tar.gz), and the corresponding raw run [SRR26035078](https://www.ebi.ac.uk/ena/browser/view/SRR26035078), not a different sample. The expression archive supplies barcode, gene and matrix axes. Historical IDs use `GSM7774400:` followed by the exact 24-base GEXSCOPE barcode; preserve exact IDs rather than trimming or guessing a join. The original work retained 3,265 hepatic-compatible scored targets and 2,741 observable nuclear rows. The 52 fixed readouts, two representations, original source-feature observability, covariates and original maxima must be prepared on their original axes before the conversion above. The full frozen input identity is in `configs/nuclear_formal_portable.json`; the P2 benchmark and BIO-06 configs bind the corresponding score and threshold digests.

**Gap:** the public archive and SRA accessions do not themselves supply the frozen role assignments, accepted mitochondrial alignments/UMI consensus, historical source-feature scoring, fixed 52-readout prepared NPZ or serialized score/threshold bytes. The public export here supplies the input converter and analysis drivers, not a standalone complete raw-to-P2 reconstruction. Obtain or reconstruct those project-specific intermediates lawfully before claiming exact reproduction. Nuclear coefficients in this submission use the original threshold; changing the threshold requires a separate conditional calculation.

## GBML181 single-nucleus input

Use [GSE230389/GSM7221332](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7221332) provider files `GSM7221332_GBML181_barcodes.tsv.gz`, `GSM7221332_GBML181_features.tsv.gz`, `GSM7221332_GBML181_matrix.mtx.gz`, and native SRA runs [SRR24283915](https://www.ebi.ac.uk/ena/browser/view/SRR24283915) and [SRR24283916](https://www.ebi.ac.uk/ena/browser/view/SRR24283916). Matrix rows correspond to the features file and columns to barcodes; retain Gene Expression rows (36,601) separately from ATAC peaks (57,101). The two runs are one sample. Preserve each RNA `RNA_score` as original text when reconstructing ties. The prior analysis used 1,343 scored targets and ten matched-N thresholds.

**Gap:** RNA-QC/working role assignments, competitive alignment, original cell/UMI consensus, source-feature selection and score TSV are project-derived intermediates and are not supplied as a portable complete raw-to-score workflow in this bundle. The original study's individual-cell malignant annotations were not obtained. The driver starts with the frozen score and thresholds; it does not imply exact recovery from the public accessions alone. Do not redistribute those private cell files or the upstream academic-only code in this project.
