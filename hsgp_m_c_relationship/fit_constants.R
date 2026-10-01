# Fit m ~ m_factor * c / ell (paper Eqs 19,21,23: a direct proportionality,
# not a general log-log regression) to the simulated minimum-m grid, in the
# region where c is comfortably above the minimum valid boundary (paper Sec
# 4.2: the near-linear relationship only holds once c is "large enough"), and
# estimate c_factor as the smallest c (per ell) for which any m is accurate.

results <- readRDS("hsgp_m_c_relationship/min_m_results.rds")

published <- list(
  se = c(m_factor = 1.75, c_factor = 3.2),
  matern32 = c(m_factor = 3.42, c_factor = 4.5),
  matern52 = c(m_factor = 2.65, c_factor = 4.1)
)

fit_one <- function(df, label) {
  df <- df[!is.na(df$m), ]
  df$ratio <- df$c / df$ell

  # m_factor: linear fit through the origin, restricted to the regime where
  # c/ell is generously above the minimum (well clear of the "stairs" bend)
  fit_df <- df[df$ratio >= 5 & df$m >= 5, ]
  m_fit <- lm(m ~ 0 + ratio, data = fit_df)
  m_factor <- unname(coef(m_fit)[1])

  # c_factor: minimum c (per ell, from the discrete c_grid) for which the
  # threshold is reachable at all, regressed against ell through the origin
  c_min_per_ell <- aggregate(c ~ ell, data = df, FUN = min)
  c_fit <- lm(c ~ 0 + ell, data = c_min_per_ell)
  c_factor <- unname(coef(c_fit)[1])

  pub <- published[[label]]
  cat(sprintf(
    "%-10s  m_factor = %.3f  (published %s)   c_factor ~= %.3f  (published %s)\n",
    label, m_factor,
    if (is.null(pub)) "n/a" else sprintf("%.2f", pub["m_factor"]),
    c_factor,
    if (is.null(pub)) "n/a" else sprintf("%.2f", pub["c_factor"])
  ))
  invisible(list(m_factor = m_factor, c_factor = c_factor, fit_df = fit_df))
}

fits <- list(
  se = fit_one(results$se, "se"),
  matern32 = fit_one(results$m32, "matern32"),
  matern52 = fit_one(results$m52, "matern52"),
  matern12 = fit_one(results$m12, "matern12")
)

saveRDS(fits, "hsgp_m_c_relationship/fits.rds")
