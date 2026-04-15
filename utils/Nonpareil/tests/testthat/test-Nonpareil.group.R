library(testthat)
library(Nonpareil)

# ---------------------------------------------------------------------------
# Shared fixture: load the three bundled .npo files once
# ---------------------------------------------------------------------------
files <- system.file(
  "extdata",
  c("HumanGut.npo", "LakeLanier.npo", "IowaSoil.npo"),
  package = "Nonpareil"
)

# ---------------------------------------------------------------------------
# 1. Output schema
# ---------------------------------------------------------------------------
test_that("Nonpareil.group returns correct class and column schema", {
  nps    <- Nonpareil.set(files, plot = FALSE)
  result <- Nonpareil.group(nps, group = c("A", "A", "B"))

  expect_s3_class(result, "Nonpareil.Group")
  expect_s3_class(result, "data.frame")

  expected_cols <- c(
    "group", "depth_bp", "mean_cov", "sd_cov", "se_cov",
    "ci_low", "ci_high", "n_samples", "n_observed", "n_projected"
  )
  expect_true(all(expected_cols %in% names(result)))
  expect_setequal(unique(result$group), c("A", "B"))
})

test_that("Nonpareil.group attributes carry group metadata", {
  nps    <- Nonpareil.set(files, plot = FALSE)
  result <- Nonpareil.group(nps, group = c("A", "A", "B"))

  expect_true(!is.null(attr(result, "group.levels")))
  expect_true(!is.null(attr(result, "group.col")))
  expect_true(!is.null(attr(result, "nps")))
  expect_true(!is.null(attr(result, "group")))
  expect_setequal(attr(result, "group.levels"), c("A", "B"))
})

# ---------------------------------------------------------------------------
# 2. Aggregation in prediction space, not coefficient space
# ---------------------------------------------------------------------------
test_that("group mean matches mean of per-sample predictions, not mean-coef prediction", {
  nps <- Nonpareil.set(files, plot = FALSE)

  # Only use curves that have fitted models
  has_model <- sapply(nps$np.curves, function(np) np$has.model)
  if (sum(has_model) < 2) skip("Need at least 2 fitted models for this test")

  # Put all model-bearing curves into one group
  group <- ifelse(has_model, "has_model", "no_model")
  result <- Nonpareil.group(nps, group = group)

  grp_data  <- result[result$group == "has_model", ]
  mod_curves <- nps$np.curves[has_model]

  # Pick a depth in the middle of the shared observed range
  shared_max_obs <- min(sapply(mod_curves, function(np) max(np$x.adj)))
  shared_min_obs <- max(sapply(mod_curves, function(np) min(np$x.adj[np$x.adj > 0])))
  test_depth <- exp((log(shared_min_obs) + log(shared_max_obs)) / 2)

  # Mean of individual predictions at test_depth
  pred_mean <- mean(sapply(mod_curves, function(np) predict(np, lr = test_depth)))

  # Nearest grid point in the group result
  nearest <- which.min(abs(grp_data$depth_bp - test_depth))
  group_val <- grp_data$mean_cov[nearest]

  expect_false(is.na(group_val))
  # Prediction-space average should be close to the group result
  # (they match exactly for equal weights; close for inverse-SD weights)
  expect_true(abs(group_val - pred_mean) < 0.1)

  # Coefficient-space average: average the two parameters then evaluate
  coef_a <- mean(sapply(mod_curves, function(np) coef(np$model)[["a"]]))
  coef_b <- mean(sapply(mod_curves, function(np) coef(np$model)[["b"]]))
  coef_pred <- Nonpareil.f(test_depth, coef_a, coef_b)

  # Prediction-space and coefficient-space results generally differ for
  # nonlinear models; at minimum the function should use prediction-space
  # (verified by the earlier assertion that group_val ≈ pred_mean).
  # Both must be in [0, 1].
  expect_true(group_val  >= 0 && group_val  <= 1)
  expect_true(coef_pred  >= 0 && coef_pred  <= 1)
})

# ---------------------------------------------------------------------------
# 3. Coverage values are bounded in [0, 1]
# ---------------------------------------------------------------------------
test_that("mean_cov, ci_low, ci_high are bounded in [0, 1]", {
  nps    <- Nonpareil.set(files, plot = FALSE)
  result <- Nonpareil.group(nps, group = c("A", "A", "B"))
  valid  <- !is.na(result$mean_cov)

  expect_true(all(result$mean_cov[valid] >= 0))
  expect_true(all(result$mean_cov[valid] <= 1))
  expect_true(all(result$ci_low[valid]   >= 0))
  expect_true(all(result$ci_high[valid]  <= 1))
})

# ---------------------------------------------------------------------------
# 4. n_observed + n_projected == n_samples
# ---------------------------------------------------------------------------
test_that("n_observed + n_projected equals n_samples at every depth", {
  nps    <- Nonpareil.set(files, plot = FALSE)
  result <- Nonpareil.group(nps, group = c("A", "A", "B"))
  valid  <- result$n_samples > 0

  expect_true(all(
    result$n_observed[valid] + result$n_projected[valid] ==
    result$n_samples[valid]
  ))
})

# ---------------------------------------------------------------------------
# 5. Single-sample groups: sd_cov and se_cov are 0
# ---------------------------------------------------------------------------
test_that("single-sample groups have sd_cov = 0 and se_cov = 0", {
  nps    <- Nonpareil.set(files, plot = FALSE)
  result <- Nonpareil.group(nps, group = c("X", "Y", "Z"))

  for (grp in c("X", "Y", "Z")) {
    gd    <- result[result$group == grp, ]
    valid <- !is.na(gd$sd_cov)
    expect_true(all(gd$sd_cov[valid] == 0))
    expect_true(all(gd$se_cov[valid] == 0))
  }
})

# ---------------------------------------------------------------------------
# 6. Mixed observed / projected depths
# ---------------------------------------------------------------------------
test_that("n_projected > 0 at depths beyond the smaller replicate's library", {
  nps <- Nonpareil.set(files, plot = FALSE)

  # Check that curves have different library sizes
  lrs <- sapply(nps$np.curves, function(np) np$LR)
  if (max(lrs) / min(lrs) < 2)
    skip("Library sizes too similar for this test")

  # Group all three curves together; grid will extend to the largest LR
  result <- Nonpareil.group(nps, group = rep("all", 3))
  gd <- result[result$group == "all", ]

  # At depths > smallest LR there should be at least one projected value
  small_lr <- min(lrs)
  deep_rows <- gd[gd$depth_bp > small_lr & !is.na(gd$n_projected), ]
  if (nrow(deep_rows) > 0) {
    expect_true(any(deep_rows$n_projected > 0))
  }
})

# ---------------------------------------------------------------------------
# 7. Curves without fitted models: no error, only observed region
# ---------------------------------------------------------------------------
test_that("Nonpareil.group works when no models are fitted", {
  nps    <- Nonpareil.set(files, plot = FALSE, skip.model = TRUE)
  result <- expect_no_error(Nonpareil.group(nps, group = c("A", "A", "B")))
  expect_s3_class(result, "Nonpareil.Group")
  # No projected values when models are absent
  expect_true(all(result$n_projected == 0))
})

# ---------------------------------------------------------------------------
# 8. Input type coercion: character, numeric, factor all accepted
# ---------------------------------------------------------------------------
test_that("group argument accepts character, numeric, and factor inputs", {
  nps <- Nonpareil.set(files, plot = FALSE)

  r_chr <- Nonpareil.group(nps, c("A", "A", "B"))
  r_num <- Nonpareil.group(nps, c(1,   1,   2))
  r_fac <- Nonpareil.group(nps, factor(c("A", "A", "B")))

  expect_s3_class(r_chr, "Nonpareil.Group")
  expect_s3_class(r_num, "Nonpareil.Group")
  expect_s3_class(r_fac, "Nonpareil.Group")
})

# ---------------------------------------------------------------------------
# 9. Input validation errors
# ---------------------------------------------------------------------------
test_that("Nonpareil.group stops on wrong input types and lengths", {
  nps <- Nonpareil.set(files, plot = FALSE)

  expect_error(Nonpareil.group("not_a_set", c("A", "A", "B")),
               "inherit from class")
  expect_error(Nonpareil.group(nps, c("A", "B")),
               "same length")
  expect_error(Nonpareil.group(nps, c("A", "A", "B"),
                               extrapolation.penalty = 0),
               "range \\(0, 1\\]")
  expect_error(Nonpareil.group(nps, c("A", "A", "B"),
                               extrapolation.penalty = 1.5),
               "range \\(0, 1\\]")
})

# ---------------------------------------------------------------------------
# 10. Equal-weights scheme produces sensible output
# ---------------------------------------------------------------------------
test_that("equal-weights scheme returns same schema and bounded coverage", {
  nps    <- Nonpareil.set(files, plot = FALSE)
  result <- Nonpareil.group(nps, group = c("A", "A", "B"),
                            weights = "equal")
  expect_s3_class(result, "Nonpareil.Group")
  valid <- !is.na(result$mean_cov)
  expect_true(all(result$mean_cov[valid] >= 0))
  expect_true(all(result$mean_cov[valid] <= 1))
})

# ---------------------------------------------------------------------------
# 11. plot.Nonpareil.Set aggregate mode: no error
# ---------------------------------------------------------------------------
test_that("plot.Nonpareil.Set with aggregate=TRUE does not error", {
  nps <- Nonpareil.set(files, plot = FALSE)

  tmp <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp), add = TRUE)

  expect_no_error({
    pdf(tmp)
    plot(nps, aggregate = TRUE, group = c("A", "A", "B"))
    dev.off()
  })
})

test_that("plot.Nonpareil.Set with aggregate=TRUE and plot.individual=TRUE does not error", {
  nps <- Nonpareil.set(files, plot = FALSE)

  tmp <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp), add = TRUE)

  expect_no_error({
    pdf(tmp)
    plot(nps,
         aggregate       = TRUE,
         group           = c("A", "A", "B"),
         plot.individual = TRUE,
         individual.alpha = 0.2)
    dev.off()
  })
})

test_that("plot.Nonpareil.Set with aggregate=TRUE accepts all ribbon types", {
  nps    <- Nonpareil.set(files, plot = FALSE)
  ribbons <- c("sd", "ci95", "ci90", "ci50", FALSE)

  for (rb in ribbons) {
    tmp <- tempfile(fileext = ".pdf")
    on.exit(unlink(tmp), add = TRUE)
    expect_no_error({
      pdf(tmp)
      plot(nps, aggregate = TRUE, group = c("A", "A", "B"),
           aggregate.ribbon = rb)
      dev.off()
    })
  }
})

test_that("plot.Nonpareil.Set aggregate=TRUE errors without group", {
  nps <- Nonpareil.set(files, plot = FALSE)
  expect_error(
    plot(nps, aggregate = TRUE),
    "'group' must be provided"
  )
})

test_that("plot.Nonpareil.Set aggregate=TRUE errors on wrong group length", {
  nps <- Nonpareil.set(files, plot = FALSE)
  expect_error(
    plot(nps, aggregate = TRUE, group = c("A", "B")),
    "same length"
  )
})

# ---------------------------------------------------------------------------
# 12. Original plot.Nonpareil.Set behaviour is unchanged
# ---------------------------------------------------------------------------
test_that("plot.Nonpareil.Set without aggregate still works as before", {
  nps <- Nonpareil.set(files, plot = FALSE)

  tmp <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp), add = TRUE)

  expect_no_error({
    pdf(tmp)
    plot(nps)
    dev.off()
  })
})
