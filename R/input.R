# Internal coercion is deliberately positional: the published executable uses
# the first four tab-separated fields and ignores later fields.
.coerce_cnv_table <- function(x) {
  if (is.character(x) && length(x) == 1L && !is.na(x) && file.exists(x)) {
    if (isTRUE(file.info(x)$isdir)) {
      stop("'x' must be a CNV file, not a directory", call. = FALSE)
    }
    if (identical(file.size(x), 0)) {
      raw <- data.frame(matrix(nrow = 0L, ncol = 4L))
    } else {
      raw <- data.table::fread(
        x, header = FALSE, sep = "\t", quote = "", fill = TRUE,
        data.table = FALSE, showProgress = FALSE
      )
    }
  } else if (is.data.frame(x) || is.matrix(x)) {
    raw <- as.data.frame(x, stringsAsFactors = FALSE)
  } else {
    stop("'x' must be a path to a CNV table, data.frame, or matrix", call. = FALSE)
  }

  if (ncol(raw) < 4L) {
    stop("CNV input must contain at least four columns", call. = FALSE)
  }
  raw <- raw[, seq_len(4L), drop = FALSE]
  names(raw) <- c("chr", "start", "end", "type")
  raw$chr <- as.character(raw$chr)
  raw$type <- as.character(raw$type)

  # This intentionally follows the executable: chromosome normalization is
  # triggered by the first record and then removes the literal string 'chr'.
  if (nrow(raw) > 0L && !is.na(raw$chr[[1L]]) &&
      grepl("chr", raw$chr[[1L]])) {
    raw$chr <- gsub("chr", "", raw$chr, fixed = TRUE)
  }

  start <- suppressWarnings(as.integer(as.character(raw$start)))
  end <- suppressWarnings(as.integer(as.character(raw$end)))
  raw$start <- start
  raw$end <- end
  raw
}

.validate_cnv_table <- function(x) {
  allowed_chr <- c(as.character(seq_len(22L)), "X", "Y")
  bad_chr <- is.na(x$chr) | !(x$chr %in% allowed_chr)
  bad_pos <- is.na(x$start) | is.na(x$end)
  bad_type <- is.na(x$type) | !(x$type %in% c("gain", "loss"))
  bad_width <- !bad_pos & x$start > x$end
  if (any(bad_chr) || any(bad_pos) || any(bad_type) || any(bad_width)) {
    stop(
      "invalid CNV input: chromosomes must be 1-22, X, or Y; positions must be present; type must be 'gain' or 'loss'; start must not exceed end",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate CNV records
#'
#' Validate the first four columns of a CNV table using the input rules of the
#' pinned X-CNV executable. Chromosome names may use the executable's leading
#' `chr` convention, and coordinates are retained as integer values.
#'
#' @param x A path, data frame, or matrix containing chromosome, start, end,
#'   and CNV type in its first four columns.
#' @return `TRUE` invisibly when validation succeeds; otherwise an ordinary R
#'   error is raised.
#' @export
validate_cnv <- function(x) {
  .validate_cnv_table(.coerce_cnv_table(x))
}

#' Read and normalize a CNV table
#'
#' Read a tab-separated BED-like CNV table or normalize an in-memory table.
#' The first four columns are used, matching the published command-line tool;
#' additional columns are ignored.
#'
#' @param x A path, data frame, or matrix containing chromosome, start, end,
#'   and CNV type in its first four columns.
#' @return A four-column data frame with columns `chr`, `start`, `end`, and
#'   `type`, classed as `xcnv_cnv`.
#' @export
read_cnv <- function(x) {
  out <- .coerce_cnv_table(x)
  .validate_cnv_table(out)
  class(out) <- c("xcnv_cnv", "data.frame")
  out
}

.prepare_cnv <- function(x) {
  cnvs <- if (inherits(x, "xcnv_cnv")) x else read_cnv(x)
  n_input <- nrow(cnvs)
  drop_first <- n_input == 1L

  # The executable duplicates a single record, changes the first duplicate to
  # a 50-bp record, sorts, and later removes the first sorted record. Keeping
  # this quirk is required for compatibility with short single-row inputs.
  if (drop_first) {
    altered <- cnvs
    altered$end[[1L]] <- altered$start[[1L]] + 49L
    cnvs <- rbind(altered, cnvs)
  }

  if (nrow(cnvs) > 0L) {
    chromosome_order <- match(cnvs$chr, c(as.character(seq_len(22L)), "X", "Y"))
    cnvs <- cnvs[order(chromosome_order, cnvs$start, cnvs$end, seq_len(nrow(cnvs))), , drop = FALSE]
    rownames(cnvs) <- NULL
    cnvs$row_id <- seq_len(nrow(cnvs))
  } else {
    cnvs$row_id <- integer()
  }
  list(data = cnvs, drop_first = drop_first, n_input = n_input)
}
