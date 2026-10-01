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
`tau in [0, min(L, S)]` in the numerator (the error), against the exact
kernel's full `[0, Inf)` area in the denominator (see "Validation" below for
why `[0, S]` rather than the wider `[0, 2S]` tried first).

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
error (the numerator only; the denominator stays the full `[0, Inf)` kernel
area throughout, an untested axis of this same uncertainty) over
`tau in [0, 2S]` (any two points at most `2S` apart within the actual data
domain `[-S, S]`) gave m_factor 12-15% high; narrowing to `tau in [0, S]`
(used throughout this script, and in the table above) brought it down to
8-10% high, consistently across all three kernels. The residual 8-10% gap
is unresolved: the paper's actual simulation that produced
`area_gp.rData`/`area_diff.rData` (the raw accuracy computation the draft
notebook's regression is fitted to) isn't published in this repo, only the
downstream curve-fit is, so some detail of their exact criterion (possibly
a different normalisation, grid resolution, or numerical quadrature) isn't
independently checkable from here.

The "near-linear regime" cutoff used for the final `m = m_factor * c / l`
fit (`c / l >= 5`) is **not** a negligible source of uncertainty on its own,
contrary to an earlier version of this document: refitting at cutoffs from
3 to 20 moves the squared-exponential/Matern-3/2/5-2 estimates by 1.5-2.8%,
and the Matern 1/2 estimate (the kernel this document is about) by up to
4.6% (8.54 at cutoff 3 vs 9.01 at cutoff 20; cutoff 30 leaves too few grid
points - under 10 - for a reliable fit for Matern 1/2, so wasn't used). This
doesn't change the qualitative conclusion below, but it's a real, checkable
contributor to the overall uncertainty, not a negligible one.

`c_factor` estimates (root-finding the smallest `c`, per `l`, at which a
large `m` still reaches the 1% threshold, via `c_min_root.R`) ran about
45-52% below published with the current `tau in [0, S]` domain (se: 1.73 vs
3.20 published; Matern 3/2: 2.17 vs 4.50; Matern 5/2: 2.03 vs 4.10). This is
expected to be a weaker replication than `m_factor`: `c_factor` in the paper
marks where the *simple proportional formula* becomes a good approximation
to a more complex true curve (the draft notebook uses an affine condition,
e.g. `c >= 0.97 + 1.45*l` for squared exponential, not a pure
`c >= c_factor * l`), not the literal feasibility boundary this script
estimates. Treat the `c_factor` numbers below as indicative only.

## Result for Matern 1/2

| | m_factor | c_factor |
|---|---|---|
| raw simulation | 8.62 | 2.68 (root-find, `c_min_root.R`; indicative only, see above) |
| bias-corrected (divide by the mean ~1.08 ratio of this sim to published, seen for the three validated kernels) | ~8.0 | - |

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
2.5x larger (raw: 8.6; bias-corrected for the 8% high bias seen on the three
validated kernels: 8.0 -> about 2.3-2.5x across that range). A user fitting
`gp_opts(kernel = "ou")` or `matern_order = 1/2` with automatic
basis-function selection is likely getting noticeably too few basis
functions for the stated accuracy target.
