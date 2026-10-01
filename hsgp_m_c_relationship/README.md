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

The draft notebook (`Numeric_diagnosis_function/Diagnosis_HSGP_v2.Rmd` in
this repo, an earlier review-round version of the paper's own derivation)
reports two different numbers per kernel: the actual `lm(log(m) ~ log(l) +
log(c))` regression coefficient fitted to their simulated grid, and a
further hand-rounded "approximation" of that fit chosen for a cleaner
formula (e.g. SE: fitted `1.73 * c / l^1.05`, rounded afterwards to a
cruder `2 * c / l`). The first round of this validation compared against
the rounded numbers and looked like a near-exact match; re-checking against
their actual regression coefficients (the fairer comparison, and close to
the final published constants) shows my simulation runs consistently
8-13% high instead:

| kernel | m_factor (this sim) | m_factor (their precise regression) | m_factor (published) |
|---|---|---|---|
| squared exponential | 1.89 | 1.73 | 1.75 |
| Matern 3/2 | 3.67 | 3.32 | 3.42 |
| Matern 5/2 | 2.86 | not in this draft | 2.65 |

I traced part of that gap to the integration domain for the accuracy
criterion (Eq 18 gives no explicit bounds): integrating the covariance
error over `tau in [0, 2S]` (any two points at most `2S` apart within the
actual data domain `[-S, S]`) gave m_factor 12-15% high; narrowing to
`tau in [0, S]` (used above) brought it down to 8-10% high, consistently
across all three kernels. Loosening or tightening the "near-linear regime"
cutoff used for the final `m = m_factor * c / l` fit moved the estimate by
less than 2% either way, so that isn't the remaining source. The residual
8-10% gap is unresolved: the paper's actual simulation that produced
`area_gp.rData`/`area_diff.rData` (the raw accuracy computation the draft
notebook's regression is fitted to) isn't published in this repo, only the
downstream curve-fit is, so some detail of their exact criterion (possibly
a different normalisation, grid resolution, or numerical quadrature) isn't
independently checkable from here.

`c_factor` estimates (root-finding the smallest `c`, per `l`, at which a
large `m` still reaches the 1% threshold) ran 30-45% below published
regardless of integration domain. This is expected to be a weaker
replication than `m_factor`: `c_factor` in the paper marks where the
*simple proportional formula* becomes a good approximation to a more
complex true curve (the draft notebook uses an affine condition, e.g.
`c >= 0.97 + 1.45*l` for squared exponential, not a pure `c >= c_factor *
l`), not the literal feasibility boundary this script estimates. Treat the
`c_factor` numbers below as indicative only.

## Result for Matern 1/2

| | m_factor | c_factor |
|---|---|---|
| raw simulation | 8.62 | 2.87 (indicative only, see above) |
| bias-corrected (divide by the mean ~1.08 ratio seen for the three validated kernels, tau in [0,S]) | ~8.0 | - |

The bias-corrected estimate is stable across both integration-domain
choices tested (8.3 with `tau in [0,2S]`, 8.0 with `tau in [0,S]`), which is
reassuring given the domain choice is the main unresolved methodological
uncertainty here.

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
