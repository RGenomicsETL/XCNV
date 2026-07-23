#' Write X-CNV-compatible output
#'
#' Write a prediction data frame as a comma-separated X-CNV output table. A
#' file is written only when an explicit output path is supplied.
#'
#' @param x A data frame returned by `predict_cnv()`.
#' @param output Explicit output path.
#' @return The normalized output path, invisibly.
#' @export
write_xcnv <- function(x, output) {
  if (!is.data.frame(x) || !"MVP_score" %in% names(x)) {
    stop("'x' must be a prediction data frame with MVP_score", call. = FALSE)
  }
  if (!is.character(output) || length(output) != 1L || is.na(output) || !nzchar(output)) {
    stop("'output' must be one explicit file path", call. = FALSE)
  }
  parent <- dirname(output)
  if (!dir.exists(parent)) {
    stop("output directory does not exist: ", parent, call. = FALSE)
  }
  output <- normalizePath(output, winslash = "/", mustWork = FALSE)
  temporary <- tempfile("xcnv-output-")
  on.exit(unlink(temporary), add = TRUE)
  # data.table::fwrite in the upstream script emits two columns named Type;
  # readers that make names unique expose the numeric one as Type.1. Keep that
  # file-level header while the R data frame uses the unambiguous Type.1 name.
  header <- names(x)
  header[header == "Type.1"] <- "Type"
  utils::write.table(
    x, file = temporary, sep = ",", quote = FALSE, row.names = FALSE,
    col.names = header, na = ""
  )
  if (file.exists(output) && !unlink(output)) {
    stop("cannot replace existing output: ", output, call. = FALSE)
  }
  if (!file.rename(temporary, output)) {
    stop("could not move output into place: ", output, call. = FALSE)
  }
  invisible(output)
}

#' Run X-CNV with a file-compatible interface
#'
#' Read CNVs, compute annotations and MVP scores, and optionally write the
#' output table. Unlike the historical executable, this function never writes
#' beside the input by default and never installs software or downloads data.
#'
#' @param input A CNV path, data frame, or matrix.
#' @param resources An `xcnv_resources` object or explicit resource directory.
#' @param model A model object/path, or `NULL` to use a resource-bundle
#'   override when present and otherwise the bundled published model.
#' @param output `NULL` to return results without writing, or one explicit CSV
#'   path.
#' @param overlap_backend Passed to `annotate_cnv()`.
#' @return A prediction data frame; the same object is returned invisibly when
#'   `output` is non-`NULL`.
#' @export
run_xcnv <- function(
  input, resources, model = NULL, output = NULL,
  overlap_backend = c("duckdb", "reference")
) {
  overlap_backend <- match.arg(overlap_backend)
  result <- predict_cnv(
    input, resources = resources, model = model,
    overlap_backend = overlap_backend
  )
  if (!is.null(output)) write_xcnv(result, output)
  result
}

#' @rdname run_xcnv
#' @export
xcnv <- run_xcnv
