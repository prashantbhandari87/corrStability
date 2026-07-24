# =============================================================================
# metabolite_set : an S3 class holding a wide metabolomics table plus the
# explicitly-resolved set of metabolite (feature) columns.
#
# Identification happens ONCE, at construction, and is validated there.
# Every downstream method reads x$metabolite_cols instead of re-guessing.
# =============================================================================

#' Resolve which columns are metabolites
#'
#' Three mutually exclusive strategies, in priority order:
#' \describe{
#'   \item{\code{metabolites}}{An explicit character vector of column names.}
#'   \item{\code{pattern}}{A regex matched against \code{names(data)}
#'     (e.g. \code{"^HMDB"}).}
#'   \item{\code{id_cols}}{Names of metadata columns; every \emph{other} numeric
#'     column is treated as a metabolite. Usually the most robust choice.}
#' }
#'
#' @param data A data.frame.
#' @param metabolites,pattern,id_cols See Details; supply exactly one.
#' @return A character vector of column names.
#' @keywords internal
#' @noRd
resolve_metabolite_cols <- function(data, metabolites = NULL,
                                    pattern = NULL, id_cols = NULL) {
  supplied <- c(!is.null(metabolites), !is.null(pattern), !is.null(id_cols))
  if (sum(supplied) == 0)
    stop("Specify exactly one of `metabolites`, `pattern`, or `id_cols`.",
         call. = FALSE)
  if (sum(supplied) > 1)
    stop("Specify only one of `metabolites`, `pattern`, or `id_cols`.",
         call. = FALSE)
  
  if (!is.null(metabolites)) {
    missing <- setdiff(metabolites, names(data))
    if (length(missing))
      stop("Metabolite columns not found in data: ",
           paste(missing, collapse = ", "), call. = FALSE)
    cols <- metabolites
  } else if (!is.null(pattern)) {
    cols <- grep(pattern, names(data), value = TRUE)
  } else {
    missing <- setdiff(id_cols, names(data))
    if (length(missing))
      stop("`id_cols` not found in data: ",
           paste(missing, collapse = ", "), call. = FALSE)
    candidate <- setdiff(names(data), id_cols)
    cols <- candidate[vapply(data[candidate], is.numeric, logical(1))]
  }
  
  if (length(cols) == 0)
    stop("No metabolite columns matched the given specification.", call. = FALSE)
  cols
}

#' Bare constructor (internal, no validation)
#'
#' @param data A data.frame.
#' @param metabolite_cols Character vector of resolved column names.
#' @return An object of class \code{"metabolite_set"}.
#' @keywords internal
#' @noRd
new_metabolite_set <- function(data, metabolite_cols) {
  stopifnot(is.data.frame(data), is.character(metabolite_cols))
  structure(
    list(data = data, metabolite_cols = metabolite_cols),
    class = "metabolite_set"
  )
}

#' Validate a metabolite_set's invariants
#'
#' @param x A candidate \code{metabolite_set}.
#' @return \code{x}, invisibly, if valid; otherwise an error.
#' @keywords internal
#' @noRd
validate_metabolite_set <- function(x) {
  cols <- x$metabolite_cols
  
  missing <- setdiff(cols, names(x$data))
  if (length(missing))
    stop("metabolite_cols refer to absent columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  
  non_numeric <- cols[!vapply(x$data[cols], is.numeric, logical(1))]
  if (length(non_numeric))
    stop("Metabolite columns must be numeric; these are not: ",
         paste(non_numeric, collapse = ", "), call. = FALSE)
  
  x
}

#' Create a metabolite_set
#'
#' User-facing constructor. Resolves the metabolite columns from one of three
#' strategies, then validates the result.
#'
#' @param data A data.frame in wide format (one row per sample, one column per
#'   metabolite, plus any metadata columns).
#' @param metabolites Optional explicit character vector of metabolite columns.
#' @param pattern Optional regex matched against column names.
#' @param id_cols Optional names of metadata columns; the remaining numeric
#'   columns are taken to be metabolites.
#' @return An object of class \code{"metabolite_set"}.
#' @export
#' @examples
#' df <- data.frame(
#'   sample_id = c("a", "b", "c"),
#'   group     = c("x", "y", "x"),
#'   glucose   = c(0, 2.1, 3.4),
#'   lactate   = c(1.2, 0, 0.9)
#' )
#' ms <- metabolite_set(df, id_cols = c("sample_id", "group"))
#' ms
metabolite_set <- function(data, metabolites = NULL,
                           pattern = NULL, id_cols = NULL) {
  cols <- resolve_metabolite_cols(data, metabolites, pattern, id_cols)
  validate_metabolite_set(new_metabolite_set(data, cols))
}

#' @export
format.metabolite_set <- function(x, ...) {
  eg <- utils::head(x$metabolite_cols, 3L)
  more <- length(x$metabolite_cols) - length(eg)
  eg_str <- paste(eg, collapse = ", ")
  if (more > 0) eg_str <- paste0(eg_str, ", ... (", more, " more)")
  c(
    "<metabolite_set>",
    paste0("  samples:     ", nrow(x$data)),
    paste0("  metabolites: ", length(x$metabolite_cols)),
    paste0("  columns:     ", eg_str)
  )
}

#' @export
print.metabolite_set <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}