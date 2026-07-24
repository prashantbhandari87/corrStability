# ---- calculate_metabolite_correlations --------------------------------------

test_that("correlations return one row per metabolite with the expected columns", {
  ms_t  <- transform_metabolites(make_metab_set(), seed = 123)
  corrs <- calculate_metabolite_correlations(ms_t)

  expect_s3_class(corrs, "data.frame")
  expect_equal(nrow(corrs), 3L)                       # glucose, lactate, citrate
  expect_setequal(corrs$metabolite, c("glucose", "lactate", "citrate"))

  expected_cols <- c(
    "metabolite",
    "original_vs_fixed", "original_vs_random",
    "original_vs_log_fixed", "original_vs_log_random",
    "fixed_vs_random", "fixed_vs_log_fixed", "fixed_vs_log_random",
    "random_vs_log_fixed", "random_vs_log_random",
    "log_fixed_vs_log_random"
  )
  expect_setequal(names(corrs), expected_cols)
})

test_that("for a zero-free metabolite, original_vs_fixed is exactly 1", {
  ms_t  <- transform_metabolites(make_metab_set(), seed = 123)
  corrs <- calculate_metabolite_correlations(ms_t)
  citrate <- corrs[corrs$metabolite == "citrate", ]
  # citrate had no zeros, so fixed_ == original -> perfect correlation
  expect_equal(citrate$original_vs_fixed, 1)
})

test_that("correlating an untransformed set warns and yields no rows", {
  ms <- make_metab_set()                              # not transformed
  expect_warning(res <- calculate_metabolite_correlations(ms),
                 "transform_metabolites")
  expect_null(res)
})

# ---- summarize_correlations -------------------------------------------------

test_that("summary has one row per correlation type and the right stats", {
  ms_t  <- transform_metabolites(make_metab_set(), seed = 123)
  corrs <- calculate_metabolite_correlations(ms_t)
  summ  <- summarize_correlations(corrs)

  expect_equal(nrow(summ), 10L)                       # 10 pairwise types
  expect_setequal(names(summ),
                  c("correlation_type", "mean", "median", "min", "max", "sd"))
  # min never exceeds max, for every row
  expect_true(all(summ$min <= summ$max))
})

# ---- run_metabolite_analysis ------------------------------------------------

test_that("the workflow returns transformed set, correlations, and summary", {
  res <- run_metabolite_analysis(make_metab_set(), seed = 123)
  expect_named(res, c("transformed", "correlations", "summary"))
  expect_s3_class(res$transformed, "metabolite_set")
  expect_s3_class(res$correlations, "data.frame")
  expect_s3_class(res$summary, "data.frame")
})

test_that("the workflow equals running the steps by hand", {
  ms <- make_metab_set()
  res <- run_metabolite_analysis(ms, seed = 123)

  manual_t <- transform_metabolites(ms, seed = 123)
  manual_c <- calculate_metabolite_correlations(manual_t)
  manual_s <- summarize_correlations(manual_c)

  expect_equal(res$correlations, manual_c)
  expect_equal(res$summary, manual_s)
})

test_that("the whole workflow is reproducible under a fixed seed", {
  a <- run_metabolite_analysis(make_metab_set(), seed = 123)
  b <- run_metabolite_analysis(make_metab_set(), seed = 123)
  expect_equal(a$correlations, b$correlations)
  expect_equal(a$summary, b$summary)
})
