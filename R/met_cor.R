#' Impute zeros and create log-transformed metabolite variants
#'
#' For every column whose name begins with \code{"X"}, appends four derived
#' columns: \code{fixed_} (zeros replaced by min-nonzero / 3), \code{random_}
#' (zeros replaced by that value plus random jitter), and the log transforms
#' \code{log_fixed_} and \code{log_random_}.
#'
#' @param data A data.frame with metabolite columns prefixed by \code{"X"}.
#' @param seed Integer seed for the random imputation, or \code{NULL} to leave
#'   the RNG state untouched. Defaults to \code{123} to match the original.
#' @return \code{data} with the derived columns appended.
#' @importFrom stats sd runif
#' @export
transform_metabolites <- function(data, seed = 123) {
  # Restore the caller's RNG stream on exit instead of clobbering it globally
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  
  metabolite_cols <- grep("^X", names(data), value = TRUE)
  
  for (col in metabolite_cols) {
    non_zero_values <- data[[col]][data[[col]] > 0]
    if (length(non_zero_values) == 0) next
    min_nonzero <- min(non_zero_values)
    
    fixed_col  <- paste0("fixed_", col)
    random_col <- paste0("random_", col)
    
    data[[fixed_col]] <- data[[col]]
    data[[fixed_col]][data[[fixed_col]] == 0] <- min_nonzero / 3
    
    data[[random_col]] <- data[[col]]
    zero_idx <- which(data[[random_col]] == 0)
    if (length(zero_idx) > 0) {
      base_value  <- min_nonzero / 3
      sd_original <- sd(data[[col]], na.rm = TRUE)
      jitter <- runif(length(zero_idx), 0, 0.2)
      data[[random_col]][zero_idx] <- base_value + sd_original * jitter
    }
    
    data[[paste0("log_fixed_",  col)]] <- log(data[[fixed_col]])
    data[[paste0("log_random_", col)]] <- log(data[[random_col]])
  }
  
  data
}

#' Pairwise correlations among metabolite variants
#'
#' @param data A data.frame produced by \code{\link{transform_metabolites}}.
#' @return A data.frame with one row per metabolite and one column per
#'   variant pair.
#' @importFrom stats cor
#' @export
calculate_metabolite_correlations <- function(data) {
  metabolite_cols <- grep("^X", names(data), value = TRUE)
  
  rows <- lapply(metabolite_cols, function(col) {
    cols <- c(
      original   = col,
      fixed      = paste0("fixed_", col),
      random     = paste0("random_", col),
      log_fixed  = paste0("log_fixed_", col),
      log_random = paste0("log_random_", col)
    )
    if (!all(cols[-1] %in% names(data))) return(NULL)
    
    cc <- function(a, b) cor(data[[cols[[a]]]], data[[cols[[b]]]],
                             use = "complete.obs")
    
    data.frame(
      metabolite              = col,
      original_vs_fixed       = cc("original", "fixed"),
      original_vs_random      = cc("original", "random"),
      original_vs_log_fixed   = cc("original", "log_fixed"),
      original_vs_log_random  = cc("original", "log_random"),
      fixed_vs_random         = cc("fixed", "random"),
      fixed_vs_log_fixed      = cc("fixed", "log_fixed"),
      fixed_vs_log_random     = cc("fixed", "log_random"),
      random_vs_log_fixed     = cc("random", "log_fixed"),
      random_vs_log_random    = cc("random", "log_random"),
      log_fixed_vs_log_random = cc("log_fixed", "log_random"),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, rows)
}

#' Summarise correlation results across metabolites
#'
#' @param correlation_results Output of
#'   \code{\link{calculate_metabolite_correlations}}.
#' @return A data.frame of mean/median/min/max/sd per correlation type.
#' @importFrom stats median sd
#' @export
summarize_correlations <- function(correlation_results) {
  correlation_cols <- setdiff(names(correlation_results), "metabolite")
  
  data.frame(
    correlation_type = correlation_cols,
    mean   = sapply(correlation_cols, function(x) mean(correlation_results[[x]],   na.rm = TRUE)),
    median = sapply(correlation_cols, function(x) median(correlation_results[[x]], na.rm = TRUE)),
    min    = sapply(correlation_cols, function(x) min(correlation_results[[x]],    na.rm = TRUE)),
    max    = sapply(correlation_cols, function(x) max(correlation_results[[x]],    na.rm = TRUE)),
    sd     = sapply(correlation_cols, function(x) sd(correlation_results[[x]],     na.rm = TRUE)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Run the full metabolite transform-and-correlate workflow
#'
#' @param data A data.frame with metabolite columns prefixed by \code{"X"}.
#' @param seed Passed to \code{\link{transform_metabolites}}.
#' @return A list with elements \code{transformed_data}, \code{correlations},
#'   and \code{summary}.
#' @export
#' @examples
#' \dontrun{
#'   res <- run_metabolite_analysis(my_data)
#'   res$summary
#' }
run_metabolite_analysis <- function(data, seed = 123) {
  transformed <- transform_metabolites(data, seed = seed)
  correlations <- calculate_metabolite_correlations(transformed)
  list(
    transformed_data = transformed,
    correlations     = correlations,
    summary          = summarize_correlations(correlations)
  )
}