# Auto-sourced by testthat before the test files run.
#
# A tiny, fully hand-computable fixture. Metabolite columns deliberately use
# non-"X" names (glucose/lactate/citrate) so the tests double as proof that
# column identification no longer depends on the old "^X" prefix. Zeros are
# placed so imputation values are known exactly:
#
#   glucose: min nonzero = 3  -> zero imputed to 3/3 = 1
#   lactate: min nonzero = 2  -> zero imputed to 2/3
#   citrate: no zeros         -> fixed_ == original, so cor(original, fixed) = 1
make_metab_df <- function() {
  data.frame(
    sample_id = c("s1", "s2", "s3", "s4"),
    group     = c("ctrl", "case", "ctrl", "case"),
    glucose   = c(0, 3, 6, 9),
    lactate   = c(2, 0, 4, 8),
    citrate   = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )
}

# A ready-made metabolite_set for tests that don't care about construction.
make_metab_set <- function() {
  metabolite_set(make_metab_df(), id_cols = c("sample_id", "group"))
}
