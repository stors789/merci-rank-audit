# Provenance and upstream rights

MERCI Rank Audit contains project-authored ranking bounds, constrained member
queries, finite hypergeometric calculations, overlap-corrected sampling, input
validation and Frisch–Waugh–Lovell code. Its own code is licensed under MIT;
see LICENSE. The small R helpers perform explicit-input orchestration and
scoring preparation.

The general mixture and inverse-overlap method for sampling a union of rare
events is due to Owen, Maximov and Chertkov (2019),
https://doi.org/10.1214/18-EJS1527. This software adapts it to a discrete
ranking and nested-gate setting.

The original MERCI source is not included. The optional bench2 command uses
a separately obtained checkout at commit
79e105aaf13159bee325c3ceceb6ee679e430341 from
https://github.com/shyhihihi/MERCI. Its academic-research-only license,
including commercial-use restrictions, applies separately. The associated
benchmark data at https://doi.org/10.5281/zenodo.8065206 are licensed CC BY 4.0.
Neither original MERCI source nor real-cell input files are distributed here.
The MIT license for this package does not relicense those external resources.
