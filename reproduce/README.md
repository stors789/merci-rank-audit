# Reproduce accepted MERCI Rank Audit results

This directory contains project-authored code, aggregate data and synthetic fixtures. The source distribution contains this directory; the wheel supplies the installed numerical library and its quickstart. Real individual-cell scores, private nuclear matrices, source sequencing reads and the separate upstream MERCI code are not distributed. See [INPUTS.md](INPUTS.md) for exact sources, schemas, conversion commands and preparation gaps.

Three reproducibility levels are distinct:

1. The aggregate tables included here redraw all five main figures and five supplementary figures; no original research tree is required.
2. The four analysis drivers accept explicit inputs and original parameters; a small synthetic or frozen-seed prefix runs end to end with independently checked expectations.
3. The original complete 51-case benchmark, 720 annotation tasks, 1,200 composition conditions, P2 reference and BIO-06 scientific runs were performed earlier. The small checks do not repeat those complete raw-data analyses. Exact real-input digests and seeds are retained; producing every historical real input from raw is not fully provided here.

## Install and verify

Exact versions from the accepted B1 reconstruction are listed in
[ENVIRONMENT.md](ENVIRONMENT.md) with a Python 3.12 constraints file.

Use Python 3.10 or newer. From an extracted source distribution or checkout:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install '.[test,reproduce]'
python -m pytest
python reproduce/verify.py --output /tmp/merci-reproduce-check
```

The output directory must not exist. The `reproduce` extra installs Matplotlib ≥3.8,<4, pypdf ≥4,<7 and ReportLab ≥4,<5. NumPy, pandas and SciPy versions are constrained by the installed library. For reproducible resource use, set `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1` before execution. On Linux the previous analysis used one core, nice 10 and an address-space ceiling of 4 GiB; these execution policies are not hidden inside the scripts.

`verify.py` also checks B1 synthetic score binding, preparation and the original float64 observation representation and missing-value handling. It executes the original rare synthetic benchmark's first seed for 100 proposals, a 20-cell synthetic composition task and ten conditions, 100 conditional synthetic nuclear proposals, and four synthetic tie/threshold scenarios. It checks frozen benchmark RMSE, expected confusion counts and transition totals, an independent least-squares coefficient, the 500-valid-beta reporting threshold, and exact hypergeometric gate probabilities. It separately renders and parses all ten figure PDFs and PNGs. Fixture inputs are synthetic; `private_members` outputs can contain IDs when real inputs are supplied and should remain private.

## Figures and data mapping

The following commands can also be run separately from any working directory. Replace `REPRO` with the absolute path to this directory and `OUT` with a destination:

```sh
python "$REPRO/figures/main_figures.py" --data "$REPRO/data/figures" --output "$OUT"
python "$REPRO/figures/supplement_figures.py" --data "$REPRO/data/figures" --output "$OUT"
python "$REPRO/figures/bio06_figure.py" --data "$REPRO/data/bio06" --output "$OUT"
```

| Figure | Exact supplied data and selection |
|---|---|
| 1 | Conceptual diagram; no tabular observations. |
| 2 | `figure2_equal_time.csv`, P2, GBML181 and B2 split 1 negative, two-worker timing batch; uniform and two shared conditional estimates. |
| 3 | `figure3_error_partitions.csv`, `figure3_paired_changes.csv`, B1 frozen and B2 retrospective original-regular mixed results. |
| 4 | `figure4_core_curves.csv`, `figure4_endpoint_transitions.csv`, B1 frozen and B2 retrospective; fixed-N and nested references retain their own denominators. `figure4_channels.csv` is retained as source context. |
| 5 | `figure5_reference_quantiles.csv`, `figure5_reference_signs.csv`, `figure5_selected_witness.csv`, `figure5_fixed_original.csv`, `figure5_fixed_alternative.csv`; two prespecified readouts in both representations. |
| S1 | `FigureS1_gate.csv`, seven real cases, two-worker batch, equal-time gate errors. |
| S2 | `FigureS2_cohorts.csv`, three stages, two designs and three specified state/call series. |
| S3 | `FigureS3_precision.csv`, `FigureS3_reference_quantiles.csv`, four original/alternative outcome combinations; `FigureS3_old_vs_reference.csv` is retained context. |
| S4 | `FigureS4_union.csv`, four fixed cases, one seed per method, four nested proposal budgets. |
| S5 | `data/bio06/rarity_comparison.csv`, both cases at original RNA ranks, eleven threshold tables each; other BIO-06 aggregate source and threshold tables are in `data/bio06/`. |

The selected CSVs preserve source numerical precision. `data/supplementary_data.zip` is the original 59-table aggregate supplement with its own manifest; it contains broader checkpoint, family, composition and nuclear tables. Figures S1–S4 were redrawn using portable code, so their presentation can differ slightly from the original snapshots without changing values. See [SOURCE_MAP.md](SOURCE_MAP.md) for provenance and rights.

## B1 strict input preparation

The source distribution provides a B1 strict entry point from the four public processed input files and the separately obtained pinned MERCI checkout. It supports the frozen B1 split numbers 21–100; split 21 was actually regenerated and compared in the development environment. It produces new private scores and a new cohort manifest before running the two tasks and ten conditions. See [INPUTS.md](INPUTS.md#b1-and-b2-experimental-benchmarks) for the explicit command, historical identity checks, resource and rights boundaries.

## Four full analysis entry points

Paths in a manifest or configuration resolve relative to that file. Output directories must be new. The formal configurations preserve the frozen scientific parameters/seeds and original input digests while replacing private paths; they are portable adaptations, not the original frozen file bytes. To use real cases, provide the exact prepared inputs described in [INPUTS.md](INPUTS.md). Do not treat a changed file as the historical byte-identical input. The complete runs can be computationally expensive and are not part of `verify.py`. The benchmark wrapper runs cases serially; its measured times cannot be substituted for the original two-worker timing batch. The portable benchmark's calendar stop is a 2099 placeholder, while the original per-flow and reference time budgets remain in the config; users should set an appropriate local stop.

```sh
# 51-case formal benchmark (or select one named case).
python "$REPRO/drivers/benchmark/run.py" --config "$REPRO/configs/benchmark_formal_portable.json" \
  --case all --stage both --output "$OUT/benchmark"

# 720 annotation cases and 120 composition splits when their frozen inputs are supplied.
python "$REPRO/drivers/cohort/run.py" --manifest "$REPRO/configs/cohort_formal_portable.json" \
  --output "$OUT/cohort"

# P2 precision reference or one of the 20 frozen repeats.
python "$REPRO/drivers/nuclear/run.py" --manifest "$REPRO/configs/nuclear_manifest_template.json" \
  --config "$REPRO/configs/nuclear_formal_portable.json" --stream reference --output "$OUT/nuclear-reference"
python "$REPRO/drivers/nuclear/run.py" --manifest "$REPRO/configs/nuclear_manifest_template.json" \
  --config "$REPRO/configs/nuclear_formal_portable.json" --stream repeat --repeat 1 --output "$OUT/nuclear-repeat-01"

# BIO-06 exact family scenarios and original-RNA conditional probability streams.
python "$REPRO/drivers/sensitivity/run.py" --config "$REPRO/configs/sensitivity_formal_portable.json" \
  --stage both --output "$OUT/bio06"
```

The P2 precision config names the four historical sequences; `--all-outcomes` runs both representations for every feature in the supplied nuclear matrix. Output raw accepted samples and member files remain private. The composition driver recomputes the core annotation, conditions and transition tables; the original separate retrospective endpoint-versus-old-output checks and witness inventory require the original archived analyses and are not asserted as reproduced by this portable wrapper. BIO-06 reads the supplied frozen aggregate maxima: the original finite R calculation is included at `drivers/sensitivity/nulls.R` for users who explicitly need to regenerate a particular N/seed (the optional fourth argument is a fixture draw count; the original default is 10,000); producing fresh maxima is not needed for the included figures or smoke checks. Original aggregate nulls for BIO-04 are supplied as `data/cohort_null_maxima.csv`.

## Software status

Fixed releases and their distribution assets are available on the [release page](https://github.com/stors789/merci-rank-audit/releases). The CI workflow runs the core tests, small driver/preparation checks, all ten figure redraws and the synthetic R coverage check. Project-authored code is MIT-licensed. The separate upstream MERCI checkout and underlying public data retain their own rights. No biological individual-cell data or restricted upstream source is included.
