.validate_overlap_tables <- function(a, b) {
  if (!all(c("chr", "start", "end") %in% names(a)) ||
      !all(c("chr", "start", "end") %in% names(b))) {
    stop("BED interval tables need chr, start, and end columns", call. = FALSE)
  }
  if (anyNA(a[, c("chr", "start", "end")]) ||
      anyNA(b[, c("chr", "start", "end")])) {
    stop("BED interval coordinates must not be missing", call. = FALSE)
  }
  invisible(TRUE)
}

.bed_overlap_pairs_reference <- function(a, b, min_a = 0, min_b = 0) {
  .validate_overlap_tables(a, b)
  out <- vector("list", nrow(a))
  for (i in seq_len(nrow(a))) {
    hit <- which(
      b$chr == a$chr[[i]] & b$start < a$end[[i]] & b$end > a$start[[i]]
    )
    if (length(hit)) {
      overlap <- pmax(
        0L,
        pmin(a$end[[i]], b$end[hit]) - pmax(a$start[[i]], b$start[hit])
      )
      keep <- overlap > 0L
      if (min_a > 0 || min_b > 0) {
        len_a <- a$end[[i]] - a$start[[i]]
        len_b <- b$end[hit] - b$start[hit]
        keep <- keep & overlap / len_a >= min_a & overlap / len_b >= min_b
      }
      hit <- hit[keep]
      overlap <- overlap[keep]
    } else {
      hit <- integer()
      overlap <- numeric()
    }
    if (!length(hit)) {
      out[[i]] <- data.frame(a_id = i, b_id = NA_integer_, overlap = 0,
                             stringsAsFactors = FALSE)
    } else {
      out[[i]] <- data.frame(a_id = rep.int(i, length(hit)), b_id = hit,
                             overlap = overlap, stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) {
    return(data.frame(a_id = integer(), b_id = integer(), overlap = numeric()))
  }
  do.call(rbind, out)
}

.new_overlap_context <- function(backend = c("duckdb", "reference")) {
  backend <- match.arg(backend)
  context <- new.env(parent = emptyenv())
  context$backend <- backend
  context$counter <- 0L
  context$con <- NULL
  if (identical(backend, "duckdb")) {
    context$con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  }
  context
}

.close_overlap_context <- function(context) {
  if (!is.null(context$con)) {
    DBI::dbDisconnect(context$con, shutdown = TRUE)
    context$con <- NULL
  }
  invisible(NULL)
}

.bed_overlap_pairs_duckdb <- function(a, b, min_a = 0, min_b = 0, context) {
  .validate_overlap_tables(a, b)
  if (!nrow(a)) {
    return(data.frame(a_id = integer(), b_id = integer(), overlap = numeric()))
  }
  if (!nrow(b)) {
    return(data.frame(
      a_id = seq_len(nrow(a)), b_id = rep.int(NA_integer_, nrow(a)),
      overlap = numeric(nrow(a))
    ))
  }
  context$counter <- context$counter + 1L
  a_name <- paste0("xcnv_a_", context$counter)
  b_name <- paste0("xcnv_b_", context$counter)
  a_sql <- DBI::dbQuoteIdentifier(context$con, a_name)
  b_sql <- DBI::dbQuoteIdentifier(context$con, b_name)
  a_input <- data.frame(
    a_id = seq_len(nrow(a)), chr = as.character(a$chr),
    start = as.numeric(a$start), end = as.numeric(a$end)
  )
  b_input <- data.frame(
    b_id = seq_len(nrow(b)), chr = as.character(b$chr),
    start = as.numeric(b$start), end = as.numeric(b$end)
  )
  DBI::dbWriteTable(context$con, a_name, a_input, temporary = TRUE)
  DBI::dbWriteTable(context$con, b_name, b_input, temporary = TRUE)
  on.exit({
    DBI::dbRemoveTable(context$con, a_name)
    DBI::dbRemoveTable(context$con, b_name)
  }, add = TRUE)
  overlap <- 'least(a."end", b."end") - greatest(a."start", b."start")'
  conditions <- c(
    'a.chr = b.chr', 'a."start" < b."end"', 'a."end" > b."start"'
  )
  if (min_a > 0) {
    conditions <- c(
      conditions,
      sprintf('(%s) / nullif(a."end" - a."start", 0) >= %.17g', overlap, min_a)
    )
  }
  if (min_b > 0) {
    conditions <- c(
      conditions,
      sprintf('(%s) / nullif(b."end" - b."start", 0) >= %.17g', overlap, min_b)
    )
  }
  query <- paste0(
    "SELECT a.a_id, b.b_id, CASE WHEN b.b_id IS NULL THEN 0 ELSE ", overlap,
    " END AS overlap FROM ", a_sql, " AS a LEFT JOIN ", b_sql, " AS b ON ",
    paste(conditions, collapse = " AND "),
    " ORDER BY a.a_id, b.b_id NULLS LAST"
  )
  out <- DBI::dbGetQuery(context$con, query)
  out$a_id <- as.integer(out$a_id)
  out$b_id <- as.integer(out$b_id)
  out$overlap <- as.numeric(out$overlap)
  out
}

.bed_overlap_pairs <- function(a, b, min_a = 0, min_b = 0, context = NULL) {
  own_context <- is.null(context)
  if (own_context) context <- .new_overlap_context("duckdb")
  if (own_context) on.exit(.close_overlap_context(context), add = TRUE)
  if (identical(context$backend, "reference")) {
    .bed_overlap_pairs_reference(a, b, min_a = min_a, min_b = min_b)
  } else {
    .bed_overlap_pairs_duckdb(a, b, min_a = min_a, min_b = min_b, context = context)
  }
}

.read_table <- function(path, header = FALSE, comment.char = "") {
  data.table::fread(
    path, header = header, sep = "\t", quote = "", fill = TRUE,
    comment.char = comment.char, data.table = FALSE, showProgress = FALSE
  )
}

.read_bed <- function(path, ncol, names) {
  x <- .read_table(path, header = FALSE, comment.char = "#")
  if (ncol(x) < ncol) {
    stop("resource has too few columns: ", basename(path), call. = FALSE)
  }
  x <- x[, seq_len(ncol), drop = FALSE]
  names(x) <- names
  x$chr <- as.character(x$chr)
  x$start <- suppressWarnings(as.integer(as.character(x$start)))
  x$end <- suppressWarnings(as.integer(as.character(x$end)))
  if (anyNA(x$start) || anyNA(x$end)) {
    stop("resource has non-numeric BED coordinates: ", basename(path), call. = FALSE)
  }
  x
}

.read_genome <- function(path) {
  x <- .read_table(path, header = FALSE, comment.char = "#")
  if (ncol(x) < 2L) stop("genome.txt must have two columns", call. = FALSE)
  x <- x[, 1:2, drop = FALSE]
  names(x) <- c("chr", "length")
  x$chr <- as.character(x$chr)
  x$length <- suppressWarnings(as.numeric(as.character(x$length)))
  x
}

.read_ljb_sites <- function(path) {
  x <- .read_table(path, header = FALSE, comment.char = "#")
  if (ncol(x) < 4L) {
    stop("the LJB26 sites file must expose its numeric site index in column 4", call. = FALSE)
  }
  x <- x[, 1:4, drop = FALSE]
  names(x) <- c("chr", "pos", "id", "site_id")
  x$chr <- as.character(x$chr)
  pos <- suppressWarnings(as.integer(as.character(x$pos)))
  if (anyNA(pos)) stop("the LJB26 sites file has invalid positions", call. = FALSE)
  # The upstream awk command selects field 9 of the combined -wao row: after
  # the five input fields this is field 4 of the converted sites file. The
  # published file is point-like; accepting BED-like third coordinates also
  # makes the boundary explicit without changing the point case.
  third <- suppressWarnings(as.integer(as.character(x$id)))
  point_like <- is.na(third) | third == pos
  if (any(!point_like & third < pos)) {
    stop("the LJB26 sites file has an interval end before its start", call. = FALSE)
  }
  starts <- ifelse(point_like, pos - 1L, pos)
  ends <- ifelse(point_like, pos, third)
  site_id <- suppressWarnings(as.integer(as.character(x$site_id)))
  if (anyNA(site_id) || any(site_id < 1L)) {
    stop("the LJB26 sites file has invalid numeric site indices", call. = FALSE)
  }
  data.frame(
    chr = x$chr, start = starts, end = ends,
    site_id = site_id,
    stringsAsFactors = FALSE
  )
}

.read_ljb_scores <- function(path) {
  x <- .read_table(path, header = TRUE, comment.char = "#")
  if (ncol(x) < 15L) {
    stop("the LJB26 score table must contain at least 15 feature columns", call. = FALSE)
  }
  x <- x[, seq_len(15L), drop = FALSE]
  names(x) <- names(x)
  x
}

.read_cdts <- function(path) {
  .read_bed(path, 4L, c("chr", "start", "end", "value"))
}

.read_encode <- function(path) {
  .read_bed(path, 4L, c("chr", "start", "end", "annotation"))
}

.read_hi <- function(path) {
  .read_bed(path, 7L, c("chr", "start", "end", "gene", "pLI", "Episcore", "GHIS"))
}

.read_merged_sites <- function(path) {
  .read_bed(path, 5L, c("chr", "start", "end", "type", "site_id"))
}

.read_sample_info <- function(path) {
  x <- .read_table(path, header = TRUE)
  if (!"ethnicity_abbr" %in% names(x)) {
    stop("sample.info.txt must contain an ethnicity_abbr column", call. = FALSE)
  }
  x$ethnicity_abbr <- as.character(x$ethnicity_abbr)
  x
}

.read_merged_sample_info <- function(path) {
  x <- .read_table(path, header = FALSE)
  if (ncol(x) < 3L) {
    stop("merged.cnv.sample.info.txt must have site, sample, and ethnicity columns", call. = FALSE)
  }
  x <- x[, 1:3, drop = FALSE]
  names(x) <- c("site_id", "sample_id", "ethnicity")
  if (nrow(x) && identical(
    tolower(as.character(unlist(x[1L, ], use.names = FALSE))),
    c("site_id", "sample_id", "ethnicity")
  )) {
    x <- x[-1L, , drop = FALSE]
  }
  x$site_id <- as.character(x$site_id)
  x$sample_id <- as.character(x$sample_id)
  x$ethnicity <- as.character(x$ethnicity)
  x
}

.load_annotation_resources <- function(resources) {
  files <- resources$files
  list(
    genome = .read_genome(files$genome),
    ljb_sites = .read_ljb_sites(files$ljb_sites),
    ljb_scores = .read_ljb_scores(files$ljb_scores),
    cdts = .read_cdts(files$cdts),
    merged_sites = .read_merged_sites(files$merged_sites),
    sample_info = .read_sample_info(files$sample_info),
    merged_sample_info = .read_merged_sample_info(files$merged_sample_info),
    encode = .read_encode(files$encode),
    hi = .read_hi(files$hi)
  )
}

.compute_ljb <- function(a, r, context) {
  pairs <- .bed_overlap_pairs(a, r$ljb_sites, context = context)
  score_names <- names(r$ljb_scores)
  out <- matrix(0, nrow(a), 15L, dimnames = list(NULL, score_names))
  missing_values <- c(0, 0, -12.30, -11.958, -20.000, 0.0003)
  names(missing_values) <- score_names[10:15]
  for (i in seq_len(nrow(a))) {
    idx <- pairs$b_id[pairs$a_id == i]
    score_idx <- r$ljb_sites$site_id[idx]
    if (any(score_idx > nrow(r$ljb_scores), na.rm = TRUE)) {
      stop("the LJB26 sites file references a missing score row", call. = FALSE)
    }
    for (j in seq_len(9L)) {
      values <- suppressWarnings(as.numeric(r$ljb_scores[[j]][score_idx]))
      # The executable counts positive records and divides by every -wao row.
      out[i, j] <- sum(!is.na(values) & values > 0) / length(values)
    }
    for (j in 10:15) {
      values <- suppressWarnings(as.numeric(r$ljb_scores[[j]][score_idx]))
      values[is.na(values)] <- missing_values[[j - 9L]]
      out[i, j] <- mean(values)
    }
  }
  as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
}

.compute_cdts <- function(a, r, context) {
  pairs <- .bed_overlap_pairs(a, r$cdts, context = context)
  out <- matrix(0, nrow(a), 2L, dimnames = list(NULL, c("CDTS_1st", "CDTS_5th")))
  lengths <- a$end - a$start + 1L
  for (i in seq_len(nrow(a))) {
    values <- as.character(r$cdts$value[pairs$b_id[pairs$a_id == i]])
    values[is.na(values) | values == "."] <- "100"
    values <- suppressWarnings(as.numeric(values))
    out[i, 1L] <- min(1, sum(values < 1) * 10 / lengths[[i]])
    out[i, 2L] <- min(1, sum(values < 5) * 10 / lengths[[i]])
  }
  as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
}

.compute_encode <- function(a, r, context) {
  cols <- c("pELS", "CTCF-bound", "PLS", "dELS", "CTCF-only", "DNase-H3K4me3")
  pairs <- .bed_overlap_pairs(a, r$encode, context = context)
  out <- matrix(0, nrow(a), length(cols), dimnames = list(NULL, cols))
  lengths <- a$end - a$start + 1L
  for (i in seq_len(nrow(a))) {
    row_pairs <- pairs[pairs$a_id == i, , drop = FALSE]
    for (k in seq_len(nrow(row_pairs))) {
      b <- row_pairs$b_id[[k]]
      if (is.na(b)) next
      labels <- strsplit(as.character(r$encode$annotation[[b]]), ",", fixed = TRUE)[[1L]]
      labels <- labels[labels != "."]
      unknown <- setdiff(labels, cols)
      if (length(unknown)) {
        stop("unknown ENCODE annotation label(s): ", paste(unknown, collapse = ", "), call. = FALSE)
      }
      out[i, match(labels, cols)] <- out[i, match(labels, cols)] + row_pairs$overlap[[k]]
    }
    out[i, ] <- out[i, ] / lengths[[i]]
  }
  as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
}

.compute_hi <- function(a, r, context) {
  cols <- c("pLI", "Episcore", "GHIS")
  pairs <- .bed_overlap_pairs(a, r$hi, context = context)
  out <- matrix(0, nrow(a), 3L, dimnames = list(NULL, cols))
  for (i in seq_len(nrow(a))) {
    idx <- pairs$b_id[pairs$a_id == i]
    for (j in seq_along(cols)) {
      values <- as.character(r$hi[[cols[[j]]]][idx])
      values[is.na(values) | values == "."] <- "0"
      values <- suppressWarnings(as.numeric(values))
      out[i, j] <- max(values)
    }
  }
  as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
}

.compute_frequency <- function(a, r, context) {
  pairs <- .bed_overlap_pairs(
    a, r$merged_sites, min_a = 0.5, min_b = 0.5, context = context
  )
  populations <- sort(unique(r$sample_info$ethnicity_abbr))
  populations <- populations[!is.na(populations)]
  n <- nrow(a)
  gain <- matrix(0, n, length(populations), dimnames = list(NULL, paste0("gain_freq_", populations)))
  loss <- matrix(0, n, length(populations), dimnames = list(NULL, paste0("loss_freq_", populations)))
  gain_overall <- numeric(n)
  loss_overall <- numeric(n)
  denom <- table(factor(r$sample_info$ethnicity_abbr, levels = populations))

  for (k in seq_len(nrow(pairs))) {
    b <- pairs$b_id[[k]]
    if (is.na(b)) next
    content <- r$merged_sample_info[r$merged_sample_info$site_id == as.character(r$merged_sites$site_id[[b]]), , drop = FALSE]
    if (nrow(content)) content <- content[!duplicated(content$sample_id), , drop = FALSE]
    if (!nrow(content)) next
    i <- pairs$a_id[[k]]
    type <- as.character(r$merged_sites$type[[b]])
    for (ethnicity in content$ethnicity) {
      if (is.na(ethnicity) || !ethnicity %in% populations) next
      if (identical(type, "gain")) {
        gain[i, ethnicity == populations] <- gain[i, ethnicity == populations] + 1 / denom[[ethnicity]]
        gain_overall[[i]] <- gain_overall[[i]] + 1 / nrow(r$sample_info)
      } else if (identical(type, "loss")) {
        loss[i, ethnicity == populations] <- loss[i, ethnicity == populations] + 1 / denom[[ethnicity]]
        loss_overall[[i]] <- loss_overall[[i]] + 1 / nrow(r$sample_info)
      }
    }
  }
  out <- cbind(as.data.frame(gain, check.names = FALSE),
               as.data.frame(loss, check.names = FALSE),
               gain.freq = gain_overall, loss.freq = loss_overall)
  out
}

#' Annotate CNV records with X-CNV features
#'
#' Compute the coding, genome-wide, regulatory, haploinsufficiency, and
#' population-frequency features used by X-CNV. Overlap is computed in-process
#' using BED half-open interval arithmetic; the reciprocal 0.5 thresholds used
#' by the upstream frequency step are retained. CNV lengths use the upstream
#' inclusive `end - start + 1` convention.
#'
#' @param x A path, data frame, or matrix containing CNV records.
#' @param resources An `xcnv_resources` object or an explicit resource directory.
#' @param overlap_backend `"duckdb"` for the production inequality-join path,
#'   or `"reference"` for the small brute-force correctness oracle.
#' @return A data frame in legacy feature order, with the input columns first,
#'   the model and auxiliary annotation columns next, and `Type.1` and
#'   `Length` at the end.
#' @export
annotate_cnv <- function(x, resources, overlap_backend = c("duckdb", "reference")) {
  overlap_backend <- match.arg(overlap_backend)
  resources <- if (inherits(resources, "xcnv_resources")) {
    resources
  } else {
    xcnv_resources(resources, require = "annotations")
  }
  prepared <- .prepare_cnv(x)
  a <- prepared$data
  annotation <- .load_annotation_resources(resources)
  context <- .new_overlap_context(overlap_backend)
  on.exit(.close_overlap_context(context), add = TRUE)
  ljb <- .compute_ljb(a, annotation, context)
  cdts <- .compute_cdts(a, annotation, context)
  frequency <- .compute_frequency(a, annotation, context)
  encode <- .compute_encode(a, annotation, context)
  hi <- .compute_hi(a, annotation, context)
  out <- cbind(
    data.frame(Chr = a$chr, Start = a$start, End = a$end, Type = a$type,
               stringsAsFactors = FALSE, check.names = FALSE),
    ljb, cdts, frequency, encode, hi,
    Type.1 = as.integer(a$type == "gain"),
    Length = a$end - a$start + 1L
  )
  if (prepared$drop_first) out <- out[-1L, , drop = FALSE]
  rownames(out) <- NULL
  out
}
