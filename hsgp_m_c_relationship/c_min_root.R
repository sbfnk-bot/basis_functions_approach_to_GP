# Precise minimum-c-for-a-given-ell boundary, via root-finding on c (rather
# than reading off the coarse c_grid used for min_m()), with m fixed at a
# large m_max standing in for "m large enough that only c limits accuracy".

source("hsgp_m_c_relationship/derive_m_c_relationship.R")

c_min_for_ell <- function(ell, nu, m_max = 400, threshold = 0.01,
                          c_lo = 1.01, c_hi = 10) {
  f <- function(c) relative_error(m_max, L = c * S, ell = ell, nu = nu, S = S) - threshold
  if (f(c_hi) > 0) {
    return(NA_real_) # not reachable even at c_hi
  }
  if (f(c_lo) < 0) {
    return(c_lo) # already accurate at the smallest c considered
  }
  uniroot(f, c(c_lo, c_hi), tol = 1e-4)$root
}

ell_grid_fine <- seq(0.05, 1.5, by = 0.05)

fit_c_factor <- function(nu, label) {
  c_min <- vapply(ell_grid_fine, c_min_for_ell, numeric(1), nu = nu)
  df <- data.frame(ell = ell_grid_fine, c_min = c_min)
  df <- df[!is.na(df$c_min) & df$c_min > 1.21, ] # drop the c >= 1.2 floor
  fit <- lm(c_min ~ 0 + ell, data = df)
  c_factor <- unname(coef(fit)[1])
  cat(sprintf("%-10s  c_factor (root-find) = %.3f\n", label, c_factor))
  invisible(list(df = df, c_factor = c_factor))
}

c_fits <- list(
  se = fit_c_factor(Inf, "se"),
  matern32 = fit_c_factor(3 / 2, "matern32"),
  matern52 = fit_c_factor(5 / 2, "matern52"),
  matern12 = fit_c_factor(1 / 2, "matern12")
)

saveRDS(c_fits, "hsgp_m_c_relationship/c_fits.rds")
