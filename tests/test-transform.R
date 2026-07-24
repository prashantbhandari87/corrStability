# ---- structure & type stability ---------------------------------------------

test_that("transform returns a metabolite_set and leaves metabolite_cols intact", {
  ms   <- make_metab_set()
  ms_t <- transform_metabolites(ms, seed = 123)
  expect_s3_class(ms_t, "metabolite_set")
  expect_identical(ms_t$metabolite_cols, ms$metabolite_cols)
})

test_that("transform adds the four derived columns per metabolite", {
  ms_t <- transform_metabolites(make_metab_set(), seed = 123)
  for (col in c("glucose", "lactate", "citrate")) {
    expect_true(all(paste0(c("fixed_", "random_", "log_fixed_", "log_random_"), col)
                    %in% names(ms_t$data)))
  }
})

# ---- fixed_ imputation is deterministic and correct -------------------------

test_that("fixed_ replaces zeros with min-nonzero / 3 and leaves the rest alone", {
  ms_t <- transform_metabolites(make_metab_set(), seed = 123)
  # glucose min nonzero = 3, so the leading zero -> 1
  expect_equal(ms_t$data$fixed_glucose, c(1, 3, 6, 9))
  # lactate min nonzero = 2, so the row-2 zero -> 2/3
  expect_equal(ms_t$data$fixed_lactate, c(2, 2/3, 4, 8))
  # citrate has no zeros, so fixed_ equals the original
  expect_equal(ms_t$data$fixed_citrate, ms_t$data$citrate)
})

test_that("log_fixed_ is the log of fixed_", {
  ms_t <- transform_metabolites(make_metab_set(), seed = 123)
  expect_equal(ms_t$data$log_fixed_glucose, log(ms_t$data$fixed_glucose))
})

# ---- random_ imputation: reproducible, and design intent --------------------

test_that("the same seed gives identical random imputation", {
  a <- transform_metabolites(make_metab_set(), seed = 123)
  b <- transform_metabolites(make_metab_set(), seed = 123)
  expect_equal(a$data, b$data)
})

test_that("fixed_ is seed-independent but random_ is seed-dependent", {
  a <- transform_metabolites(make_metab_set(), seed = 1)
  b <- transform_metabolites(make_metab_set(), seed = 2)
  # deterministic imputation is unaffected by the seed
  expect_equal(a$data$fixed_glucose, b$data$fixed_glucose)
  # random imputation is not
  expect_false(isTRUE(all.equal(a$data$random_glucose, b$data$random_glucose)))
})

test_that("random_ imputes only the zero rows, leaving nonzeros equal to original", {
  ms_t <- transform_metabolites(make_metab_set(), seed = 123)
  # glucose zero is row 1; rows 2:4 must be untouched
  expect_equal(ms_t$data$random_glucose[2:4], c(3, 6, 9))
})

# ---- RNG hygiene: the global stream must be left untouched -------------------

test_that("a seeded transform does not disturb the caller's RNG state", {
  set.seed(999)
  before <- get(".Random.seed", envir = .GlobalEnv)
  invisible(transform_metabolites(make_metab_set(), seed = 123))
  after <- get(".Random.seed", envir = .GlobalEnv)
  expect_identical(before, after)
})

test_that("seed = NULL uses the ambient RNG (results follow the outer seed)", {
  ms <- make_metab_set()
  set.seed(1); a <- transform_metabolites(ms, seed = NULL)
  set.seed(1); b <- transform_metabolites(ms, seed = NULL)
  expect_equal(a$data$random_glucose, b$data$random_glucose)
})

# ---- edge case: all-zero columns are skipped --------------------------------

test_that("an all-zero metabolite column is skipped, not imputed", {
  df <- make_metab_df()
  df$empty <- c(0, 0, 0, 0)
  ms   <- metabolite_set(df, id_cols = c("sample_id", "group"))
  ms_t <- transform_metabolites(ms, seed = 123)
  expect_true("empty" %in% ms_t$metabolite_cols)          # still a metabolite
  expect_false("fixed_empty" %in% names(ms_t$data))       # but no derived cols
})
