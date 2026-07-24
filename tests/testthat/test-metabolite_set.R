# ---- construction & the three resolution strategies -------------------------

test_that("id_cols keeps the numeric remainder as metabolites", {
  ms <- metabolite_set(make_metab_df(), id_cols = c("sample_id", "group"))
  expect_s3_class(ms, "metabolite_set")
  expect_setequal(ms$metabolite_cols, c("glucose", "lactate", "citrate"))
})

test_that("metabolite columns are found regardless of naming (no ^X reliance)", {
  df <- data.frame(
    sample_id   = "s1",
    glucose     = 1,
    lactate     = 2,
    HMDB0000122 = 3,
    `100.5`     = 4,       # a name the old ^X code would never have matched
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  ms <- metabolite_set(df, id_cols = "sample_id")
  expect_setequal(ms$metabolite_cols,
                  c("glucose", "lactate", "HMDB0000122", "100.5"))
})

test_that("pattern strategy matches a regex", {
  df <- make_metab_df()
  names(df)[names(df) == "glucose"] <- "HMDB0001"
  names(df)[names(df) == "lactate"] <- "HMDB0002"
  ms <- metabolite_set(df, pattern = "^HMDB")
  expect_setequal(ms$metabolite_cols, c("HMDB0001", "HMDB0002"))
})

test_that("explicit metabolites vector is used verbatim", {
  ms <- metabolite_set(make_metab_df(), metabolites = c("glucose", "citrate"))
  expect_identical(ms$metabolite_cols, c("glucose", "citrate"))
})

test_that("id_cols strategy silently drops non-numeric non-id columns", {
  df <- make_metab_df()
  df$note <- c("a", "b", "c", "d")          # character, not an id column
  ms <- metabolite_set(df, id_cols = c("sample_id", "group"))
  expect_false("note" %in% ms$metabolite_cols)
})

# ---- error paths ------------------------------------------------------------

test_that("supplying no strategy errors", {
  expect_error(metabolite_set(make_metab_df()), "exactly one")
})

test_that("supplying more than one strategy errors", {
  expect_error(
    metabolite_set(make_metab_df(), pattern = "^gl", id_cols = "sample_id"),
    "only one"
  )
})

test_that("explicit metabolites that don't exist error", {
  expect_error(
    metabolite_set(make_metab_df(), metabolites = c("glucose", "nope")),
    "not found"
  )
})

test_that("id_cols that don't exist error", {
  expect_error(
    metabolite_set(make_metab_df(), id_cols = "nope"),
    "not found"
  )
})

test_that("a pattern matching nothing errors", {
  expect_error(
    metabolite_set(make_metab_df(), pattern = "ZZZ_no_match"),
    "No metabolite columns"
  )
})

# ---- validation -------------------------------------------------------------

test_that("non-numeric metabolite columns are rejected at construction", {
  df <- make_metab_df()
  df$note <- c("a", "b", "c", "d")
  expect_error(
    metabolite_set(df, metabolites = c("glucose", "note")),
    "numeric"
  )
})

# ---- print / format ---------------------------------------------------------

test_that("print returns the object invisibly", {
  ms <- make_metab_set()
  expect_output(out <- print(ms), "metabolite_set")
  expect_identical(out, ms)
})

test_that("format reports the sample and metabolite counts", {
  ms <- make_metab_set()
  txt <- format(ms)
  expect_true(any(grepl("samples:\\s+4", txt)))
  expect_true(any(grepl("metabolites:\\s+3", txt)))
})
