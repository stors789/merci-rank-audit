# Changelog

## 0.1.4 (candidate)

- Add an explicit GBML181 reconstruction workflow from its public matrix and
  native SRA inputs, subject to the documented reference and upstream licenses.

## 0.1.3

- Add B1/B2 six-case input generation from explicitly supplied processed files,
  preserving historical roles, scoring cases, thresholds and private input binding.
- Add static checks of all 720 historical case selectors and retain the existing
  B1 strict-only and B2 one-case entry points.

## 0.1.2

- Add an explicit B1 frozen strict split preparation from the four public processed inputs and pinned upstream MERCI checkout.
- Preserve the historical float64 observed-source-feature representation in cohort preparation.
- Record newly generated and historical identities separately; do not reuse old samples.

## 0.1.1

- Include portable aggregate plotting and analysis-reproduction materials in the source distribution.
- Keep the 0.1.0 numerical implementation unchanged.
- Add four analysis-driver checks and redraw all ten main and supplementary figures in CI.

## 0.1.0

Initial public software release based on the validated 0.1.0rc6 candidate.

- Exact raw, failed-as-empty and passing-only member families and requested witnesses.
- Direct analytic member marginals for certain and single-random-gate cases.
- Conditional order sampling, optional nuclear propagation and content-bound reuse.
- Synthetic examples, independent-enumeration tests and explicit bench2 reproduction.
- MIT license for project-authored code; separate upstream source and data terms.
