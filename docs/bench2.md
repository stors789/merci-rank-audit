# Reproduce the benchmark-2 audit

The optional `bench2` command reads four public processed inputs from
[MERCI v1.0.3, Zenodo 8065206](https://zenodo.org/records/8065206).
Obtain these original files under their CC BY 4.0 terms:

| CLI input | File | Bytes |
|---|---|---:|
| `--info` | `cell_information(bench2_scRNAseq).txt` | 311642 |
| `--expression` | `normalized_exp_matrix(bench2_scRNAseq).txt` | 744245952 |
| `--variants` | `mtSNV_profile(bench2_scRNAseq).txt` | 29309783 |
| `--coverage` | `MT_covearge_cells(bench2_scRNAseq).txt` | 202932422 |

The provider's `covearge` spelling is retained. Obtain the external MERCI
checkout at commit `79e105aaf13159bee325c3ceceb6ee679e430341` from
https://github.com/shyhihihi/MERCI under its own academic-use license.
The R step requires Matrix, data.table, jsonlite, e1071 and proxy. Install
those dependencies into an explicit user library before running the command.
This chain was verified using R 4.4.1 with an independent user library;
R base/recommended libraries and the native runtime were shared.

```sh
merci-rank-audit bench2 --info /path/to/info.txt \
  --expression /path/to/expression.txt --variants /path/to/variants.txt \
  --coverage /path/to/coverage.txt --upstream /path/to/MERCI \
  --rscript /path/to/Rscript --r-library /path/to/R-user-library \
  --output /private/new-output --sample both
```

The output directory must be new. The command generates scores and ten matched-N
thresholds before auditing them. The expected case has 1272 targets and 144
selected features; raw core/union is 306/390, failed-as-empty is 0/390 and
passing-only is 306/390. The supplied order passes, while passage is not certain
across legal Neff orders. At the documented seeds, 5000 uniform orders have
4245 gate hits and 5000 conditional proposals have 4554 accepted orders.
These checks reproduce an existing input from the same original study.
