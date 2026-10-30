# P2: public inputs and reconstruction

This workflow reconstructs the previously fixed P2 analysis from the original
count archives and the complete paired FASTQs. The joint RNA identity step uses
all five pairs from [GSE242889](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242889),
while scoring remains within [GSM7774400 / 2T / C24](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7774400).
The resulting roles are working RNA identities, not independently established
malignant labels or direct evidence of physical mitochondrial transfer.

Supply these inputs separately:

* The ten GEO count archives listed by name and URL in [`sample_manifest.json`](sample_manifest.json),
  all in one directory. Their `genes.tsv`, `barcodes.tsv` and `matrix.mtx` are
  interpreted on the original sample axes. Despite their `.tar.gz` names, the
  copies used for the local analysis were uncompressed POSIX tar files.
* The original ENA exports [SRR26035078 mate 1](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR260/078/SRR26035078/SRR26035078_1.fastq.gz)
  and [mate 2](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR260/078/SRR26035078/SRR26035078_2.fastq.gz).
  The previously verified copies were 18,108,503,586 and 21,207,937,884 bytes,
  with provider MD5s `6eaa67ecd2c1f1c8797915a8a7b46913` and
  `a9b893f2a9857e78babec6bc46b8f65c`, respectively. Check the provider's
  current metadata for any newly acquired copy. The original complete pair
  count was 417,883,844; a different run is not a substitute.
* The complete [UCSC hg38 FASTA](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
  and its full-genome minimap2 index. The previously verified FASTA was
  983,659,424 bytes, provider MD5 `1c9dcaddfa41027f17cd8f7a82c7293b`;
  check the [provider checksum list](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt)
  when acquiring it. Use minimap2 2.30-r1287 and the original split index:

  ```sh
  minimap2 -x splice:sr -t 4 -I 256M -d "$HG38_INDEX" "$HG38_FASTA"
  ```

  Building the index is a separate, resource-intensive step. The previously
  used index was 15,420,201,719 bytes; the driver checks that size and the
  minimap version. Size alone does not verify index contents. A chrM-only index
  cannot provide the required whole-genome competition.
* The separately obtained MSigDB file `h.all.v2024.1.Hs.symbols.gmt`, with all
  50 human Hallmark gene sets. Obtain it from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/)
  under its own terms. The original local file has SHA-256
  `ee2463540042078bfa3f67828e1e223bb354446d9fbb4d22845866835ba5c772`;
  this is its local identity, not a publisher-supplied digest. The wrapper
  requires those exact bytes, including gene-set order. It is not bundled or
  relicensed here.
* A separately obtained [MERCI checkout at commit `79e105aaf13159bee325c3ceceb6ee679e430341`](https://github.com/shyhihihi/MERCI/tree/79e105aaf13159bee325c3ceceb6ee679e430341)
  under its academic-use license. The driver checks the revision and the
  required source files. R 4.4.1 with jsonlite 2.0.0, e1071 1.7.17,
  proxy 0.4.29 and zipfR 0.6.70 in an explicit R user library was used.
  The tested Linux/Python 3.13.11 dependency closure is in
  [`../gbml/requirements-py313.txt`](../gbml/requirements-py313.txt). Install
  `merci-rank-audit` in that environment or provide a separate audit Python.

The included source distribution also supplies the project-authored fixed C5
role decisions and the published fixed threshold table. The workflow verifies
the table's SHA-256 and requires the original matched N of 3,265; it does not
generate fresh thresholds or reuse historical samples. The fixed maxima allow
the gate and exact member-family audit, but do not reconstruct the historical
full null draws, gate p-values or gate BH-adjusted values. Source-feature p/q
statistics retain the original calculation.

The real inputs were regenerated locally through the stages in this guide, and
the aggregate comparisons are in [P2 validation](../../data/p2/validation.json).
The combined `run.py` wrapper has not independently repeated the full FASTQs.
The original complete extraction ran with a stricter plus-line check; after that
process exited, the code was corrected to accept the original plus-prefix rule
and checked synthetically. The first nuclear models used row-major arrays and
showed up to 5.55×10⁻¹⁶ differences in standard deviations; restoring the
historical column-major layout and rerunning only the nuclear stage made its
input NPZ and both model tables byte-identical to the historical outputs.

From an extracted source distribution, supply explicit paths and a new private
output directory:

```sh
(
  ulimit -S -v 4194304
  ulimit -H -v 8388608
  nice -n 10 python "$REPRO/drivers/p2/run.py" \
    --archives "$GEO_ARCHIVE_DIR" --r1 "$R1_FASTQ" --r2 "$R2_FASTQ" \
    --reference "$HG38_FASTA" --index "$HG38_INDEX" --minimap "$MINIMAP" \
    --gmt "$HALLMARK_GMT" --upstream "$MERCI_CHECKOUT" --rscript "$RSCRIPT" \
    --r-library "$R_USER_LIBRARY" --audit-python "$AUDIT_PYTHON" \
    --output "$PRIVATE/p2-new"
)
```

Give the parent access to four CPUs. Ordinary children have a 4 GiB address-space
limit; the alignment controller has a 4 GiB soft and 8 GiB hard limit so its
minimap child can use an 8 GiB ceiling. These limits are not estimates of resident
memory. Check memory and disk before execution, including space for new private
candidate chunks, BAMs and consensus files. Each stage logs to the output folder;
on a failure its partial output is retained and no success receipt is written.
The wrapper does not resume failed stages automatically or download inputs.
If alignment alone was interrupted, its completed chunks can be resumed with
the same arguments and an explicit `--resume`:

```sh
python "$REPRO/drivers/p2/align_candidates.py" \
  --extraction "$PRIVATE/p2-new/extraction" --roles "$PRIVATE/p2-new/roles/P2_roles.tsv.gz" \
  --index "$HG38_INDEX" --minimap "$MINIMAP" --output "$PRIVATE/p2-new/alignment" --resume
```

Use the same soft/hard address-space limits and four-CPU allocation. The
alignment driver checks its own source, extraction/chunk manifest, role digest, candidate
file size/mtime, index size/mtime and minimap executable/version before reusing
a completed chunk. A partial BAM or log is renamed and kept before its chunk is
retried; a completed alignment is never overwritten. This stage-specific resume
does not continue the remaining wrapper stages automatically.

After that recovered alignment completes, the remaining stages can be run in
order from the same private directory. The output paths below must be new:

```sh
(
  set -e
  ulimit -S -v 4194304
  ulimit -H -v 8388608
  export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4 NUMBA_NUM_THREADS=4
  P2_RUN="$PRIVATE/p2-new"
  nice -n 10 python "$REPRO/drivers/p2/consensus.py" --extraction "$P2_RUN/extraction" \
    --alignment "$P2_RUN/alignment" --output "$P2_RUN/consensus"
  nice -n 10 python "$REPRO/drivers/p2/qualify.py" --roles "$P2_RUN/roles/P2_roles.tsv.gz" \
    --consensus "$P2_RUN/consensus" --alignment "$P2_RUN/alignment" --reference "$P2_RUN/reference/chrM.fa" \
    --upstream "$MERCI_CHECKOUT" --r-library "$R_USER_LIBRARY" --rscript "$RSCRIPT" --output "$P2_RUN/qualification"
  nice -n 10 python "$REPRO/drivers/p2/prepare_merci.py" --qualification "$P2_RUN/qualification" \
    --counts "$P2_RUN/counts/raw_HCC_counts.h5ad" --output "$P2_RUN/merci_inputs"
  nice -n 10 "$RSCRIPT" --vanilla "$REPRO/drivers/p2/score_merci.R" "$P2_RUN/merci_inputs" "$P2_RUN/scores" \
    "$MERCI_CHECKOUT" "$R_USER_LIBRARY" "$REPRO/data/bio06/all_thresholds.csv"
  nice -n 10 python "$REPRO/drivers/p2/nuclear.py" --counts "$P2_RUN/counts/raw_HCC_counts.h5ad" \
    --metadata "$P2_RUN/qualification/full_metadata.tsv.gz" --cells "$P2_RUN/scores/cells.tsv" \
    --gmt "$HALLMARK_GMT" --output "$P2_RUN/nuclear"
  nice -n 10 python "$REPRO/drivers/p2/prepare_audit.py" --cells "$P2_RUN/scores/cells.tsv" \
    --threshold-table "$REPRO/data/bio06/all_thresholds.csv" --nuclear "$P2_RUN/nuclear/nuclear.npz" \
    --output "$P2_RUN/audit_input"
  export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1
  nice -n 10 "$AUDIT_PYTHON" -m merci_rank_audit audit --manifest "$P2_RUN/audit_input/manifest.json" \
    --output "$P2_RUN/exact_audit" --sample none
  nice -n 10 python "$REPRO/drivers/p2/short_prefix.py" --manifest "$P2_RUN/audit_input/manifest.json" \
    --audit-python "$AUDIT_PYTHON" --output "$P2_RUN/reference_prefix"
)
```

Inspect the individual stage receipts after a manual continuation. The wrapper's
`run_complete.json` is produced only by a successful uninterrupted wrapper run.

The separate synthetic checks exercise boundaries and failure evidence:

```sh
python "$REPRO/fixtures/p2/check_gexscope.py"
python "$REPRO/fixtures/p2/check_half_counts.py"
python "$REPRO/fixtures/p2/check_joint_concat.py"
python "$REPRO/fixtures/p2/check_preflight.py" --output "$PRIVATE/p2-preflight-check"
python "$REPRO/fixtures/p2/check_incomplete_extraction.py" --output "$PRIVATE/p2-incomplete-check"
python "$REPRO/fixtures/p2/check_alignment_resume.py" --output "$PRIVATE/p2-alignment-resume-check"
python "$REPRO/fixtures/p2/check_alignment_consensus.py" --output "$PRIVATE/p2-alignment-boundary-check"
python "$REPRO/fixtures/p2/check_nuclear.py" --output "$PRIVATE/p2-nuclear-check"
"$RSCRIPT" --vanilla "$REPRO/drivers/p2/source_features.R" /tmp /tmp \
  "$MERCI_CHECKOUT" "$R_USER_LIBRARY" --self-test
python "$REPRO/fixtures/p2/check_qualification.py" --upstream "$MERCI_CHECKOUT" \
  --r-library "$R_USER_LIBRARY" --rscript "$RSCRIPT" --output "$PRIVATE/p2-qualification-check"
python "$REPRO/fixtures/p2/check_score_input.py" --upstream "$MERCI_CHECKOUT" \
  --r-library "$R_USER_LIBRARY" --rscript "$RSCRIPT" --output "$PRIVATE/p2-score-check"
```

The last check runs the scorer on artificial inputs and expects rejection of a
nonmatching N; add `--full-size` with a different output directory to check its
complete 3,265-target artificial path. The incomplete-extraction check confirms
that readable partial candidate chunks remain without a completion receipt;
it is not an automatic FASTQ restart facility. The alignment check uses a fake
minimap executable and a sparse placeholder index; it tests recovery, not mapping.

The fixed GEXSCOPE barcode is R1 positions 1–8, 25–32 and 49–56 concatenated,
with the exact `GSM7774400:` namespace in the RNA matrix. It is not the GBML181
16-base layout. UMI bases are R1 positions 58–69 with ACGT and Q20, and the
following poly-T window requires at least ten T bases. Both-orientation chrM
21-mer candidates are aligned against all of hg38. The original primary chrM
filter uses MAPQ ≥30, at least 60 matched bases and no splice; joint cell/UMI
consensus trims five read-end bases, requires base Q30 and discards conflicting
UMI positions. Discovery and evaluation halves are fixed by the original
cell/UMI hash. Strict source and target eligibility requires more than 1,000
accepted mitochondrial reads, with the original discovery, observability and
source-enrichment criteria preserved.

RNA score inputs are read from the newly generated original-GEO count matrix,
not a different run's expression estimate. Nuclear output contains the fixed
50 Hallmark sets plus the original six-contact and Human10 readouts, in both
representations, on the original 3,428-cell percentile domain. The driver writes
new scores, the original-label HC3 models, an exact three-family audit and the
first 100 proposals of the original conditional reference stream, binding each
to the new score, threshold and nuclear bytes.

Covariates follow the original `std > 0` selection rule. In a synthetic
constant-UMI example, floating-point rounding produced a tiny nonzero standard
deviation and retained that constant column; the model's SVD projection then
determines its numerical rank. The fixture records this behavior rather than
changing the historical rule. Compare the actual generated P2 covariate matrix
with the historical input before claiming equivalence.

All real-cell outputs, candidate reads, BAMs and derived matrices are private
artifacts. The source, fixtures and aggregate validation do not include the
original reads, third-party GMT, academic-use upstream code or individual-cell
results.
