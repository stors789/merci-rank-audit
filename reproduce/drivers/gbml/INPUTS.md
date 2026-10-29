# GBML181: public inputs, separate dependencies and local reconstruction

This workflow is specific to [GSE230389 / GSM7221332](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7221332), the GBML181 sample, and its two technical runs [SRR24283915](https://www.ncbi.nlm.nih.gov/sra/SRR24283915) and [SRR24283916](https://www.ncbi.nlm.nih.gov/sra/SRR24283916). It reconstructs the previously fixed single-sample RNA working roles and a MERCI-compatible input; it does not recreate the original author's malignant-cell labels, Cell Ranger duplicate flags, physical-transfer truth or a second patient.

Obtain the three GEO attachments `GSM7221332_GBML181_barcodes.tsv.gz`, `GSM7221332_GBML181_features.tsv.gz`, and `GSM7221332_GBML181_matrix.mtx.gz` from that GEO record. The original matrix has 4,379 barcode columns, 36,601 Gene Expression rows and 57,101 ATAC peak rows. Obtain the two **native SRA archives** using [NCBI SRA Toolkit](https://github.com/ncbi/sra-tools) or the archive provider; the ENA single-read FASTQ exports do not preserve the same-spot barcode/UMI technical read needed here. The two native archives used in the local run were 4,754,274,923 and 4,919,265,509 bytes. Archive checksums and provider access conditions should be checked when acquiring new copies.

Obtain the [complete UCSC hg38 FASTA](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) with nuclear sequences and chrM. The original 983,659,424-byte file matched the provider's `hg38.fa.gz` MD5 `1c9dcaddfa41027f17cd8f7a82c7293b`; verify against the provider's [checksum list](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt) when obtaining a copy. The workflow derives chrM and its 21-mer set from this FASTA, while aligning candidates against **all** of hg38. Do not substitute a chrM-only index.

Install [minimap2 2.30-r1287](https://github.com/lh3/minimap2/releases/tag/v2.30). An independently built compatible split index can be supplied with `--index`; it is built from that exact FASTA as follows. Building an index is a separate, resource-intensive step, not part of the raw-data workflow run:

```sh
minimap2 -x splice:sr -t 4 -I 256M -d "$HG38_INDEX" "$HG38_FASTA"
```

The historical index was 15,420,201,719 bytes and used `k=15,w=5`; the workflow checks that size and the minimap version. This is distinct from the 21-mer candidate filter. Plan disk and memory for index construction independently.

The native extractor is the included project-authored `extract_sra.cpp`. Compile it against NCBI VDB headers/libraries and zlib. The locally tested packages were Ubuntu `libncbi-vdb-dev` and `libncbi-vdb3` version `3.0.2+dfsg-2build2` (and SRA Toolkit `3.0.3+dfsg-6ubuntu5` for archive inspection). They may be installed or extracted locally without changing system configuration:

```sh
g++ -std=c++17 -O2 -I"$VDB_INCLUDE" "$REPRO/drivers/gbml/extract_sra.cpp" \
  -L"$VDB_LIB" -Wl,-rpath-link,"$VDB_LIB" -l:libncbi-vdb.so.3 -lz -o "$EXTRACTOR"
```

`VDB_INCLUDE` is the directory containing `vdb/manager.h`; `VDB_LIB` contains `libncbi-vdb.so.3` and its runtime dependencies. When using locally extracted VDB packages, supply their configuration directory through `NCBI_VDB_CONFIG` if it is not installed system-wide. The tested Ubuntu [VDB development package](https://archive.ubuntu.com/ubuntu/pool/universe/n/ncbi-vdb/libncbi-vdb-dev_3.0.2+dfsg-2build2_amd64.deb) and [matching runtime](https://archive.ubuntu.com/ubuntu/pool/universe/n/ncbi-vdb/libncbi-vdb3_3.0.2+dfsg-2build2_amd64.deb) are separate dependencies. Obtain a separate academic-use [MERCI checkout at commit `79e105aaf13159bee325c3ceceb6ee679e430341`](https://github.com/shyhihihi/MERCI/tree/79e105aaf13159bee325c3ceceb6ee679e430341). The workflow checks that commit and the actual files it calls; the upstream source is not included or relicensed. R 4.4.1 with jsonlite 2.0.0, e1071 1.7.17, proxy 0.4.29 and [zipfR](https://cran.r-project.org/package=zipfR) 0.6.70 in an explicit R user library was used. The locally tested Python 3.13.11 dependency closure is in `requirements-py313.txt`; install it from PyPI, and install the `merci-rank-audit` package in that environment or pass a separate `--audit-python` having it installed.

From an extracted source distribution, set the paths above and use a **new private output** directory:

```sh
(
  ulimit -S -v 4194304
  ulimit -H -v 8388608
  nice -n 10 python "$REPRO/drivers/gbml/run.py" \
    --barcodes "$BARCODES" --features "$FEATURES" --matrix "$MATRIX" \
    --sra1 "$SRR24283915" --sra2 "$SRR24283916" --reference "$HG38_FASTA" \
    --index "$HG38_INDEX" --minimap "$MINIMAP" --extractor "$EXTRACTOR" --vdb-lib "$VDB_LIB" \
    --upstream "$MERCI_CHECKOUT" --rscript "$RSCRIPT" --r-library "$R_USER_LIBRARY" \
    --output "$PRIVATE/gbml181-new"
)
```

Allow four CPUs in the parent affinity for the RNA stage; that stage records the original four-thread numerical setup. All other stages are confined to one allowed CPU and numerical thread. Each ordinary child has an address-space limit of 4 GiB; the alignment controller has a 4 GiB soft/8 GiB hard limit solely so its minimap child can use an 8 GiB address-space limit. Limits are not estimates of RSS. Confirm current free memory and disk before running. The process does not automatically download or install anything, write a full FASTQ or a full-genome BAM, delete original inputs, or generate fresh null draws.

Small synthetic checks are separate from the real-data execution:

```sh
python "$REPRO/fixtures/gbml/check_thresholds.py"
python "$REPRO/fixtures/gbml/check_boundaries.py" --extractor "$EXTRACTOR" \
  --vdb-lib "$VDB_LIB" --output "$PRIVATE/gbml-synthetic-check"
"$RSCRIPT" --vanilla "$REPRO/drivers/gbml/qualify.R" "$MERCI_CHECKOUT" /tmp /tmp "$R_USER_LIBRARY" --self-test
"$RSCRIPT" --vanilla "$REPRO/fixtures/gbml/check_score.R" "$MERCI_CHECKOUT" "$R_USER_LIBRARY"
```

The result contains newly generated private roles, candidate reads, accepted chrM BAM, UMI counts, variants, score, `thresholds.csv`, exact-audit manifest and member files. These contain individual-cell data and must remain private. The workflow checks the frozen `reproduce/data/bio06/all_thresholds.csv` byte hash before starting, then selects the ten GBML181 original maxima only when their N matches the newly scored sample; changed maxima are rejected. The saved maxima support the gate and exact-family audit but do not reconstruct the historical full null draws, gate p-values or BH-adjusted gate values; source-feature p/q statistics retain the original calculation. `--audit-python` runs exact families without sampling, and records the actual score/threshold hashes; it does not reuse or re-sign a historic sample. For initial endpoint and component checks, see the aggregate validation supplied with this source distribution.

RNA roles are the original 500–5000-feature and <15% MT QC followed by all non-MT Gene Expression genes, log normalization, randomized SVD 50, 15 neighbors, Leiden resolution 1 and seed 0. Source cluster 13 and target clusters 0/1/2/3/4/6/8 remain fixed. The two archives are processed in their original order, with barcode16+UMI12. Barcode Q20 is counted but not an additional selection gate; UMIs must be ACGT with every base Q≥20. Candidate 21-mers are checked on both orientations; the full hg38 primary chrM competition requires MAPQ ≥30, at least 60 matched bases and no splice. Both runs are merged before the same-cell/UMI/position consensus; base Q≥30 and five-base end trimming apply, and conflicts are discarded. These consensus UMI counts are not identical in definition to the original Cell Ranger read counts. Source selection keeps strict `Cell_reads > 10`, `min_d=2`, 5% observability, original Nmut/OR/p/q defaults and documented shape-only fixes. RNA is normalized by all Gene Expression rows, excluding peaks; original textual `RNA_score` is kept in the score file. Missing source observation is not converted into a covered zero.
