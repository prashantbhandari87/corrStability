# =============================================================================
# Analysis methods for metabolite_set.
#
# Each generic dispatches on the metabolite_set object and reads
# x$metabolite_cols -- no column-name guessing anywhere below.
# =============================================================================

# --- helper: derived column names for one metabolite -------------------------
# Kept in one place so transform and correlate agree on the naming scheme.
derived_names <- function(col) {
  list(
    fixed      = paste0("fixed_",      col),
    random     = paste0("random_",     col),
    log_fixed  = paste0("log_fixed_",  col),
    log_random = paste0("log_random_", col)
  )
}

# =============================================================================
# transform_metabolites
# =============================================================================

#' Impute zeros and add log-transformed metabolite variants
#'
#' For each metabolite column, appends four derived columns:
#' \code{fixed_} (zeros replaced by min-nonzero / 3), \code{random_} (zeros
#' replaced by that value plus uniform jitter scaled by the column SD), and the
#' natural logs \code{log_fixed_} and \code{log_random_}.
#'
#' @param x A \code{\link{metabolite_set}}.
#' @param seed Integer seed for the random imputation, or \code{NULL} to leave
#'   the RNG state untouched. Defaults to \code{123}.
#' @param ... Unused; for method consistency.
#' @return The \code{metabolite_set} with derived columns added to its
#'   \code{data}. \code{metabolite_cols} is unchanged (it still lists the
#'   original columns).
#' @export
transform_metabolites <- function(x, ...) {
  UseMethod("transform_metabolites")
}

#' @rdname transform_metabolites
#' @importFrom stats sd runif
#' @export
transform_metabolites.metabolite_set <- function(x, seed = 123, ...) {
  # Restore the caller's RNG stream on exit rather than clobbering it globally.
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  
  data <- x$data
  
  for (col in x$metabolite_cols) {
    non_zero <- data[[col]][data[[col]] > 0]
    if (length(non_zero) == 0) next
    min_nz <- min(non_zero)
    nm <- derived_names(col)
    
    # fixed_: zeros -> min_nonzero / 3
    data[[nm$fixed]] <- data[[col]]
    data[[nm$fixed]][data[[nm$fixed]] == 0] <- min_nz / 3
    
    # random_: zeros -> min_nonzero / 3 + jitter (SD-scaled, vectorised)
    data[[nm$random]] <- data[[col]]
    zero_idx <- which(data[[nm$random]] == 0)
    if (length(zero_idx) > 0) {
      base_value  <- min_nz / 3
      sd_original <- sd(data[[col]], na.rm = TRUE)
      jitter <- runif(length(zero_idx), 0, 0.2)
      data[[nm$random]][zero_idx] <- base_value + sd_original * jitter
    }
    
    # log transforms
    data[[nm$log_fixed]]  <- log(data[[nm$fixed]])
    data[[nm$log_random]] <- log(data[[nm$random]])
  }
  
  x$data <- data
  x
}

# =============================================================================
# calculate_metabolite_correlations
# =============================================================================

#' Pairwise correlations among metabolite variants
#'
#' @param x A \code{metabolite_set} that has been passed through
#'   \code{\link{transform_metabolites}}.
#' @param ... Unused; for method consistency.
#' @return A data.frame: one row per metabolite, one column per variant pair.
#' @export
calculate_metabolite_correlations <- function(x, ...) {
  UseMethod("calculate_metabolite_correlations")
}

#' @rdname calculate_metabolite_correlations
#' @importFrom stats cor
#' @export
calculate_metabolite_correlations.metabolite_set <- function(x, ...) {
  data <- x$data
  
  rows <- lapply(x$metabolite_cols, function(col) {
    nm <- derived_names(col)
    if (!all(unlist(nm) %in% names(data))) {
      warning("Skipping '", col,
              "': derived columns missing. Run transform_metabolites() first.",
              call. = FALSE)
      return(NULL)
    }
    
    v <- c(list(original = col), nm)   # named list of the five column names
    cc <- function(a, b) cor(data[[v[[a]]]], data[[v[[b]]]], use = "complete.obs")
    
    data.frame(
      metabolite              = col,
      original_vs_fixed       = cc("original",  "fixed"),
      original_vs_random      = cc("original",  "random"),
      original_vs_log_fixed   = cc("original",  "log_fixed"),
      original_vs_log_random  = cc("original",  "log_random"),
      fixed_vs_random         = cc("fixed",     "random"),
      fixed_vs_log_fixed      = cc("fixed",     "log_fixed"),
      fixed_vs_log_random     = cc("fixed",     "log_random"),
      random_vs_log_fixed     = cc("random",    "log_fixed"),
      random_vs_log_random    = cc("random",    "log_random"),
      log_fixed_vs_log_random = cc("log_fixed", "log_random"),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, rows)
}

# =============================================================================
# summarize_correlations  (plain function -- input is a data.frame, not the set)
# =============================================================================

#' Summarise correlation results across metabolites
#'
#' @param correlation_results Output of
#'   \code{\link{calculate_metabolite_correlations}}.
#' @return A data.frame of mean / median / min / max / sd per correlation type.
#' @importFrom stats median sd
#' @export
summarize_correlations <- function(correlation_results) {
  correlation_cols <- setdiff(names(correlation_results), "metabolite")
  agg <- function(f) vapply(correlation_cols,
                            function(cn) f(correlation_results[[cn]], na.rm = TRUE),
                            numeric(1))
  data.frame(
    correlation_type = correlation_cols,
    mean   = agg(mean),
    median = agg(median),
    min    = agg(min),
    max    = agg(max),
    sd     = agg(sd),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# run_metabolite_analysis  (full workflow)
# =============================================================================

#' Run the full transform-and-correlate workflow
#'
#' @param x A \code{\link{metabolite_set}}.
#' @param seed Passed to \code{\link{transform_metabolites}}.
#' @param ... Unused; for method consistency.
#' @return A list with elements \code{transformed} (the transformed
#'   \code{metabolite_set}), \code{correlations}, and \code{summary}.
#' @export
#' @examples
#' \dontrun{
#'   ms  <- metabolite_set(my_data, id_cols = c("sample_id", "group"))
#'   res <- run_metabolite_analysis(ms)
#'   res$summary
#' }
run_metabolite_analysis <- function(x, ...) {
  UseMethod("run_metabolite_analysis")
}

#' @rdname run_metabolite_analysis
#' @export
run_metabolite_analysis.metabolite_set <- function(x, seed = 123, ...) {
  transformed  <- transform_metabolites(x, seed = seed)
  correlations <- calculate_metabolite_correlations(transformed)
  list(
    transformed  = transformed,
    correlations = correlations,
    summary      = summarize_correlations(correlations)
  )
}