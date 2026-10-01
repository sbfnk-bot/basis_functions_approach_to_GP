# Reproduces the empirical m/c/(l/S) relationship of Riutort-Mayol et al.
# (2023, Stat Comput, Sec 4.3), for the squared exponential and Matern
# kernels, as a validation step before deriving the equivalent relationship
# for the (uncovered) Matern nu = 1/2 kernel.
#
# S (half-range of the input domain) is fixed to 1 throughout, so l/S == l
# and L == c.

# --- Exact covariance functions (paper, Sec 2.1 / page 4) -------------------

matern_cov <- function(tau, ell, nu, alpha = 1) {
  if (is.infinite(nu)) {
    return(alpha * exp(-0.5 * tau^2 / ell^2))
  }
  if (nu == 1 / 2) {
    return(alpha * exp(-abs(tau) / ell))
  }
  if (nu == 3 / 2) {
    r <- sqrt(3) * abs(tau) / ell
    return(alpha * (1 + r) * exp(-r))
  }
  if (nu == 5 / 2) {
    r <- sqrt(5) * abs(tau) / ell
    return(alpha * (1 + r + r^2 / 3) * exp(-r))
  }
  stop("unsupported nu")
}

# --- Spectral density, D = 1 (paper, Sec 2.1, general Matern formula; the
# published Eqs 1-3 absorb the "ordinary frequency" 4*pi^2 factor into the
# omega convention used throughout the rest of the paper, i.e. omega is
# evaluated at sqrt(lambda_j) directly, as confirmed against the SE case
# (Eq 1) and the companion repo's own spd() helper for SE) -------------------

spectral_density <- function(omega, ell, nu, alpha = 1) {
  if (is.infinite(nu)) {
    return(alpha * sqrt(2 * pi) * ell * exp(-0.5 * ell^2 * omega^2))
  }
  num <- alpha * 2 * sqrt(pi) * gamma(nu + 1 / 2) * (2 * nu)^nu
  den <- gamma(nu) * ell^(2 * nu)
  (num / den) * (2 * nu / ell^2 + omega^2)^(-(nu + 1 / 2))
}

# Sanity check: nu = 1/2 should reduce to the standard OU/exponential
# spectral density S(omega) = 2 * ell / (1 + ell^2 * omega^2).
stopifnot(all.equal(
  spectral_density(seq(0, 5, 0.5), ell = 1.3, nu = 1 / 2),
  2 * 1.3 / (1 + 1.3^2 * seq(0, 5, 0.5)^2)
))

# --- HSGP eigenfunctions/eigenvalues (paper, Eqs 6-7) ------------------------

eigenvalue <- function(j, L) (j * pi / (2 * L))^2
eigenfunction <- function(j, L, x) sqrt(1 / L) * sin(sqrt(eigenvalue(j, L)) * (x + L))

# Approximate covariance with m basis functions (paper, unnumbered eq before
# Eq 15, page 7): k_tilde_m(tau) = sum_j s(sqrt(lambda_j)) phi_j(tau) phi_j(0)
approx_cov <- function(tau, m, L, ell, nu, alpha = 1) {
  j <- seq_len(m)
  lambda_j <- eigenvalue(j, L)
  s_j <- spectral_density(sqrt(lambda_j), ell = ell, nu = nu, alpha = alpha)
  phi0 <- eigenfunction(j, L, 0)
  # outer(tau, j) matrix of phi_j(tau)
  phi_tau <- outer(tau, j, function(x, jj) eigenfunction(jj, L, x))
  as.numeric(phi_tau %*% (s_j * phi0))
}

# --- Accuracy criterion (paper, Eqs 15/18): total variation distance between
# exact and approximate covariance, as a proportion of the exact kernel's
# area. Integrated over tau in [0, min(L, S)] (narrower than the [0, 2S]
# first tried, see README): this reduced an initial 12-15% overshoot versus
# the paper's own fitted constants to 8-10%, so it is the domain actually
# used here, not [0, 2S].

# Closed-form area under the exact kernel over tau in [0, Inf) (alpha = 1),
# used as the reference area in Eq 18 rather than a truncated one.
kernel_area <- function(ell, nu) {
  if (is.infinite(nu)) {
    return(ell * sqrt(pi / 2))
  }
  if (nu == 1 / 2) {
    return(ell)
  }
  if (nu == 3 / 2) {
    return(2 * ell / sqrt(3))
  }
  if (nu == 5 / 2) {
    return(8 * ell / (3 * sqrt(5)))
  }
  stop("unsupported nu")
}

relative_error <- function(m, L, ell, nu, S, n_grid = 2000) {
  tau_max <- min(L, S)
  tau <- seq(0, tau_max, length.out = n_grid)
  k <- matern_cov(tau, ell, nu)
  k_tilde <- approx_cov(tau, m, L, ell, nu)
  # trapezoidal integration
  trapz <- function(y) {
    dx <- diff(tau)
    sum((head(y, -1) + tail(y, -1)) / 2 * dx)
  }
  trapz(abs(k - k_tilde)) / kernel_area(ell, nu)
}

# Minimum m (<= m_max) such that relative_error(m) < 0.01, via doubling then
# bisection on the integer grid. Returns NA if even m_max is not accurate
# enough (c too small for this ell, per the paper's "minimum valid c" finding).
min_m <- function(L, ell, nu, S, m_max = 400, threshold = 0.01) {
  if (relative_error(m_max, L, ell, nu, S) >= threshold) {
    return(NA_integer_)
  }
  lo <- 1L
  hi <- m_max
  # exponential search for an upper bound that satisfies the threshold
  m <- 4L
  while (m < hi && relative_error(m, L, ell, nu, S) >= threshold) {
    lo <- m
    m <- min(m * 2L, hi)
  }
  hi <- m
  while (hi - lo > 1L) {
    mid <- (lo + hi) %/% 2L
    if (relative_error(mid, L, ell, nu, S) < threshold) {
      hi <- mid
    } else {
      lo <- mid
    }
  }
  hi
}

# --- Grid (paper, Sec 4.3, Fig 5 analogue) -----------------------------------

S <- 1
c_grid <- c(1.2, 1.5, 2, 2.5, 3, 4, 5, 6)
ell_grid <- seq(0.05, 2, by = 0.05)

derive_relationship <- function(nu, label) {
  grid <- expand.grid(c = c_grid, ell = ell_grid)
  grid$m <- mapply(
    function(c, ell) min_m(L = c * S, ell = ell, nu = nu, S = S),
    grid$c, grid$ell
  )
  grid$nu <- label
  grid
}

results <- list(
  se = derive_relationship(Inf, "se"),
  m32 = derive_relationship(3 / 2, "matern32"),
  m52 = derive_relationship(5 / 2, "matern52"),
  m12 = derive_relationship(1 / 2, "matern12")
)

saveRDS(results, "hsgp_m_c_relationship/min_m_results.rds")

for (nm in names(results)) {
  n_na <- sum(is.na(results[[nm]]$m))
  cat(nm, ": ", nrow(results[[nm]]), " grid points, ", n_na,
    " unreachable within m_max\n",
    sep = ""
  )
}
