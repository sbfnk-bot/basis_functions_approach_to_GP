# Deriving the m/c/lengthscale relationship for Matern nu = 1/2

Riutort-Mayol et al. (2023, Stat Comput, Sec 4.3.1) give empirical constants
linking the number of HSGP basis functions `m`, the boundary factor `c` and
the lengthscale `l` (normalised by the input half-range `S`) for the squared
exponential and Matern 3/2 and 5/2 kernels, but not for Matern 1/2 (the
Ornstein-Uhlenbeck covariance). `EpiNow2` (PR epiforecasts/EpiNow2#1562,
issue epiforecasts/EpiNow2#1557) falls back to the Matern 3/2 constants for
nu = 1/2 and the `"ou"` kernel, describing them as "the most conservative
published values".

This reimplements the paper's empirical procedure (Eqs 15/18 for the
accuracy criterion, Eqs 6-7 for the HSGP eigenfunctions, the Sec 2.1 spectral
densities) from scratch, to validate it against the published constants and
then apply it to nu = 1/2.

## Method

`derive_m_c_relationship.R`: for a grid of boundary factors `c` and
lengthscales `l` (with `S = 1`), finds the minimum number of basis functions
`m` for which the relative total-variation distance between the exact and
HSGP-approximate covariance functions falls below 1%, integrated over
`tau in [0, min(L, 2S)]` (two points can only be `2S` apart within the actual
data domain `[-S, S]`).

`fit_constants.R`: fits `m ~ m_factor * c / l` (no intercept) to the grid, in
the region where `c / l` is well clear of the feasibility boundary, mirroring
the paper's own near-linear relationship (Sec 4.2).

`c_min_root.R`: estimates `c_factor` by root-finding the smallest `c`, for
each `l`, at which a large (400) `m` still reaches the 1% threshold at all.

## Validation (squared exponential, Matern 3/2, Matern 5/2)

| kernel | m_factor (this sim) | m_factor (published) | m_factor (draft notebook*) |
|---|---|---|---|
| squared exponential | 1.95 | 1.75 | 2.0 |
| Matern 3/2 | 3.83 | 3.42 | 3.32 |
| Matern 5/2 | 2.99 | 2.65 | not in this draft |

\* `Numeric_diagnosis_function/Diagnosis_HSGP_v2.Rmd` in this repo, an
earlier review-round draft of the paper's own derivation, reports these
intermediate constants before further refinement for publication.

This simulation lands within about 12% of the published constants, and
within about 2% of the paper's own earlier draft for the squared exponential
case — closer to their own intermediate numbers than to the final rounded
ones, which is a reasonable sanity check given the published constants are
themselves described as a practical simplification of a more complex fitted
curve (see the draft notebook's commented-out quadratic-in-log(c) models).

`c_factor` estimates from the root-finding approach ran about 30-40% below
published (e.g. 2.17 vs 3.20 for squared exponential). This is a weaker
replication than `m_factor`: `c_factor` in the paper marks where the
*simple linear formula* becomes a good approximation to the true curve (the
draft notebook uses an affine condition, e.g. `c >= 0.97 + 1.45*l` for
squared exponential, not a pure `c >= c_factor * l`), not the literal
feasibility boundary this script estimates. Treat the `c_factor` numbers
below as indicative only.

## Result for Matern 1/2

| | m_factor | c_factor |
|---|---|---|
| raw simulation | 9.28 | 3.23 (indicative only, see above) |
| bias-corrected (divide by the mean 1.12 ratio seen for the three validated kernels) | 8.28 | - |

Matern 1/2 (the Ornstein-Uhlenbeck covariance) is the roughest of the four
kernels the paper and `EpiNow2` consider (non-differentiable at `tau = 0`,
vs. once/twice/infinitely differentiable for 3/2, 5/2 and squared
exponential), and its spectral density decays the slowest at high frequency.
Both the published constants' own ordering (SE 1.75 < Matern 5/2 2.65 <
Matern 3/2 3.42, increasing as smoothness decreases) and this simulation
predict it needs substantially more basis functions than Matern 3/2 for the
same accuracy, not fewer or similar.

**This means `EpiNow2`'s current fallback of reusing the Matern 3/2
`m_factor` (3.42) for nu = 1/2 and the `"ou"` kernel is not actually
conservative for this kernel** — the simulated `m_factor` for 1/2 is roughly
2.4-2.7x larger (8.3-9.3 vs 3.42). A user fitting `gp_opts(kernel = "ou")`
or `matern_order = 1/2` with automatic basis-function selection is likely
getting noticeably too few basis functions for the stated accuracy target.
