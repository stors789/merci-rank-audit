# MERCI Rank Audit

Audit fixed integer Neff scores, unique RNA ranks and ten matched-N gate maxima.
The package gives exact raw, failed-as-empty, and passing-only member families;
it also supports requested conditional samples and optional nuclear outcomes.
These are properties of the specified ranking mechanism, not probabilities of
physical mitochondrial transfer.

[Source code](https://github.com/stors789/merci-rank-audit) · [Versioned releases](https://github.com/stors789/merci-rank-audit/releases)

Install the release wheel in an isolated Python environment:

```sh
python -m pip install merci_rank_audit-0.1.0-py3-none-any.whl
merci-rank-audit init-example --output ./example
merci-rank-audit audit --manifest ./example/manifest.json --output ./exact
merci-rank-audit audit --manifest ./example/manifest.json --output ./sampled --sample both
```

The example has synthetic nuclear outcomes, a certain gate and an impossible
gate. The exact-only run does not sample. For a specific member's certificate:

```sh
merci-rank-audit audit --manifest ./example/manifest.json --case synthetic_certain \
  --member synthetic_19 --decision selected --output ./one-member
```

When sampling is requested, accepted conditional orders are drawn with the
1/M overlap correction. Sampled frequencies are not exact family membership.
Nuclear beta quantiles require at least 500 valid beta values. Reuse with
`--reuse-samples OLD_OUTPUT` requires identical input identities, budget,
methods, seeds, package version, sampling-code and runtime identities, and
every expected output file and digest; unbound old samples are rejected.

## Input manifest

The manifest has a positive integer `budget` and a `cases` list. Each case
has a unique name, `score`, `thresholds`, `seed_uniform` and
`seed_conditional`, plus optional `nuclear`. Paths are absolute or relative
to the manifest. Score TSV requires unique nonempty `cell`, integer
`Neff`, and unique `RNA_rank` 1..N; `DNA_rank` and `observable` are
optional. Missing DNA rank leaves the original gate unevaluated.

Threshold CSV has exactly ten distinct rows `N,cutoff,maximum` at cutoffs
100, 90, ..., 10 and the same actual N. Neff, ranks and thresholds must be
exact integers in 0..2^53−1; direct APIs reject float values.
Optional nuclear NPZ contains `ids`, Boolean `observable`, `X`,
`original`, `alternative`, and `names`; ids must exactly match the score
set, X/Y rows follow ids[observable], and the two Y matrices match names.

## Outputs

Each case has `summary.json`, `members.tsv`, and optional
`witness.json`/`witness.tsv`. Probability source is exact for certain or
single-random-gate cases; a general multi-gate case requires requested MC
sampling. Impossible G has an undefined passing-only family and probability.
Exact rational numerators and denominators use the explicitly reported integer
base (10 normally, 16 with a `0x` prefix for very large integers). The decimal
field is a 40-digit approximation and the log field remains finite for a tiny
positive probability. MC numerators and denominators are observed counts, not
exact probabilities. Exact member rationals are deduplicated in
`probability_values.json` and referenced by `probability_exact_id`; unusually
long rational components are omitted from individual TSV rows to avoid repeating
large integers for thousands of members.
Sampling outputs include actual accepted counts and prefixes; the nuclear
status distinguishes not supplied, not requested, and no accepted order.
Existing output directories are never overwritten.

## Optional bench2 reproduction

```sh
merci-rank-audit bench2 --info /path/info.txt --expression /path/expression.txt \
  --variants /path/variants.txt --coverage /path/coverage.txt \
  --upstream /path/to/MERCI-checkout --rscript /path/to/Rscript \
  --r-library /path/to/isolated/R/library --output /private/new-output --sample both
```

This runs the external R scorer using the four user-provided processed inputs
and the pinned upstream checkout, then audits the newly generated score. The
upstream code and real inputs are not in the wheel. Their licenses apply separately; see [NOTICE.md](NOTICE.md). This package's own code is [MIT licensed](LICENSE). The release includes synthetic test
inputs; real-cell inputs and the external MERCI checkout are obtained separately.

## Testing and reproducibility

The test suite checks exact family membership against independent enumeration,
conditional sampling, integer and identity boundaries, witness orders, optional
nuclear outcomes and sample-reuse provenance. After installing the wheel:

```sh
python -m pip install pytest
python -m pytest -q tests
```

The source archive contains the tests and synthetic inputs under
`src/merci_rank_audit/examples/`. The release also supplies a separate
synthetic-test-data ZIP and SHA256 checksums. Bench2 input acquisition and R
requirements are described in [the reproduction guide](docs/bench2.md).
GitHub Actions runs the same tests on a standard Linux runner.

The implementation was tested with Python 3.13. Results describe fixed-score
ranking families and the specified tie mechanism. Nuclear Monte Carlo intervals
and experimental annotations retain their stated scopes.
