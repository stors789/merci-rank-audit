# Provenance and rights

The included plotting selections and supplement archive are project-authored aggregates from the fixed author revision01. Their original project source mapping is as follows; paths below refer to the private historical research repository and are provenance labels, not runtime dependencies.

| Data | Historical aggregate source |
|---|---|
| Figure 2 / S1 | `reports/merci_bio03_formal_benchmark_20261028/figures/real_equal_time_data.csv` and formal benchmark summaries. |
| Figure 3 / S2 | `reports/merci_bio04_cohort_composition_20261029/summary/task_error_summary.csv`, `task_paired_summary.csv`, `condition_summary.csv`. |
| Figure 4 | `reports/merci_bio04_cohort_composition_20261029/summary/condition_summary.csv`, `transition_summary.csv`, `channel_summary.csv`. |
| Figure 5 / S3 | `reports/merci_bio05_nuclear_precision_20261029/summary/reference_quantiles.csv`, `reference_observed_signs.csv`, precision summaries, and fixed P2 nuclear aggregates. |
| S4 | `reports/merci_installable_release_20261029/representative_results_rc5/union_prefixes.csv`. |
| S5 | `reports/merci_bio06_fixed_input_sensitivity_20261029/rarity/rarity_comparison.csv`; exact thresholds and scenario summaries in the same BIO-06 report. |
| 59-table ZIP | `manuscripts/merci_rank_audit_bioinformatics_20261029_revision01/supplementary_data.zip`, including `MANIFEST.json` with source paths and original byte hashes. |
| 44 synthetic NPZ fixtures | Frozen BIO-03 synthetic prepared arrays; identifiers are synthetic, not biological cell IDs. |
| Drivers | Project-authored BIO-03 benchmark workers, BIO-04 exact family/cohort functions, BIO-05 nuclear worker/FWL and BIO-06 scenario/rarity functions, adapted to explicit input and output paths. |
| B1 strict preparation | Project-authored `scripts/merci_methodology_20261026/{prepare.R,engine.R}` and `scripts/merci_frozen_evaluation_20261026/benchmark.R` adapted to four explicit public inputs and the frozen split selection; pinned upstream MERCI functions are loaded from a separately obtained checkout. |
| B1/B2 six-case preparation | The same project-authored scoring engine and benchmark loops, `scripts/merci_frozen_evaluation_20261026/prepare.R`, and original role/seed registries, adapted to explicit inputs; original regular and two candidate-specific LOO calls use the separately obtained pinned MERCI functions. |
| GBML181 input preparation | Project-authored `scripts/merci_GBML181_public_recompute_20261028/role_assessment.py`, `scripts/merci_GBML181_full_recompute_20261028/` and the prior competitive-alignment/UMI consensus functions adapted to explicit public inputs; upstream MERCI functions remain in the separately obtained checkout. |
| P2 input preparation | Project-authored `scripts/hcc_state_structure_20261019/`, `scripts/hcc_p2_complete_rna_library_20261022/` and `scripts/merci_nuclear_robustness_20261026/` adapted to explicit public archives, FASTQs, reference, GMT and output paths; pinned MERCI and MSigDB remain separate dependencies. |

Project-authored software is MIT-licensed as stated in the package LICENSE. The MERCI benchmark dataset is CC BY 4.0 from [Zenodo 8065206](https://doi.org/10.5281/zenodo.8065206); cite Zhang et al., Cancer Cell 2023, DOI 10.1016/j.ccell.2023.09.003. The underlying P2 source is Li et al., Hepatology 2024, DOI 10.1097/HEP.0000000000000673; the GBML181 source is Ravn-Boess et al., Cell Reports 2023, DOI 10.1016/j.celrep.2023.113374 and GEO GSE230389. The figures and aggregates do not relicense underlying biological datasets; users must check and honor the original source terms. Upstream MERCI is separately academic-use licensed and is not included. Project-authored code follows the existing NOTICE for the general rare-event mixture method. No original reads, real-cell score/nuclear tables or third-party full texts are included.
