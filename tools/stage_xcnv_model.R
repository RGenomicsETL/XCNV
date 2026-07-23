#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: stage_xcnv_model.R INPUT_RDATA OUTPUT_TSV", call. = FALSE)
}
input <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
if (!requireNamespace("xgboost", quietly = TRUE)) {
  stop("model staging requires xgboost; the installed XCNV package does not", call. = FALSE)
}
if (!requireNamespace("openssl", quietly = TRUE)) {
  stop("model staging requires openssl for the source receipt", call. = FALSE)
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("model staging requires data.table to validate the portable artifact", call. = FALSE)
}
if (!identical(as.character(utils::packageVersion("xgboost")), "1.2.0.1")) {
  stop(
    "the pinned serialized X-CNV model must be opened with xgboost 1.2.0.1 before conversion; do not deserialize it with a newer XGBoost",
    call. = FALSE
  )
}

environment <- new.env(parent = emptyenv())
loaded <- load(input, envir = environment)
if (!"xcnv.model" %in% loaded || !inherits(environment$xcnv.model, "xgb.Booster")) {
  stop("INPUT_RDATA must contain an xgb.Booster named xcnv.model", call. = FALSE)
}
model <- environment$xcnv.model
feature_names <- model$feature_names
if (!length(feature_names) || anyDuplicated(feature_names)) {
  stop("source model feature names are missing or duplicated", call. = FALSE)
}

dump_path <- tempfile("xcnv-dump-")
on.exit(unlink(dump_path), add = TRUE)
xgboost::xgb.dump(model, fname = dump_path, with_stats = FALSE)
lines <- readLines(dump_path, warn = FALSE)
tree_id <- -1L
rows <- list()
for (line in lines) {
  if (grepl("^booster\\[[0-9]+\\]:$", line)) {
    tree_id <- as.integer(sub("^booster\\[([0-9]+)\\]:$", "\\1", line))
    next
  }
  if (tree_id < 0L) tree_id <- 0L
  text <- trimws(line)
  branch <- regexec(
    "^([0-9]+):\\[([^<]+)<([^]]+)\\] yes=([0-9]+),no=([0-9]+),missing=([0-9]+)",
    text
  )
  match <- regmatches(text, branch)[[1L]]
  if (length(match)) {
    feature <- match[[3L]]
    if (grepl("^f[0-9]+$", feature)) {
      index <- as.integer(sub("^f", "", feature)) + 1L
      if (index < 1L || index > length(feature_names)) stop("split feature index is out of range", call. = FALSE)
      feature <- feature_names[[index]]
    }
    rows[[length(rows) + 1L]] <- data.frame(
      tree_id = tree_id, node_id = as.integer(match[[2L]]), feature = feature,
      threshold = as.numeric(match[[4L]]), yes_id = as.integer(match[[5L]]),
      no_id = as.integer(match[[6L]]), missing_id = as.integer(match[[7L]]),
      leaf_value = NA_real_, stringsAsFactors = FALSE
    )
    next
  }
  leaf <- regexec("^([0-9]+):leaf=([^,]+)", text)
  match <- regmatches(text, leaf)[[1L]]
  if (length(match)) {
    rows[[length(rows) + 1L]] <- data.frame(
      tree_id = tree_id, node_id = as.integer(match[[2L]]), feature = NA_character_,
      threshold = NA_real_, yes_id = NA_integer_, no_id = NA_integer_,
      missing_id = NA_integer_, leaf_value = as.numeric(match[[3L]]),
      stringsAsFactors = FALSE
    )
    next
  }
  if (nzchar(text)) stop("unrecognized XGBoost dump line: ", text, call. = FALSE)
}
nodes <- do.call(rbind, rows)
nodes <- nodes[order(nodes$tree_id, nodes$node_id), , drop = FALSE]

config <- xgboost::xgb.config(model)
base_match <- regexec('"base_score":"([^"]+)"', config)
base_parts <- regmatches(config, base_match)[[1L]]
if (length(base_parts) != 2L) stop("could not obtain base_score from XGBoost config", call. = FALSE)
base_score <- as.numeric(base_parts[[2L]])
objective <- if (!is.null(model$params$objective)) model$params$objective else "reg:logistic"
if (!identical(objective, "reg:logistic")) stop("only reg:logistic is supported", call. = FALSE)

parent <- dirname(output)
if (!dir.exists(parent)) dir.create(parent, recursive = TRUE)
temporary <- tempfile("xcnv-portable-", tmpdir = parent)
on.exit(unlink(temporary), add = TRUE)
source_connection <- file(input, open = "rb")
source_hash <- as.character(openssl::sha256(source_connection))
close(source_connection)
header <- c(
  "# format=xcnv-tree-v1",
  paste0("# objective=", objective),
  paste0("# base_score=", format(base_score, digits = 17, scientific = FALSE)),
  paste0("# feature_names=", paste(feature_names, collapse = ",")),
  paste0("# source_sha256=", source_hash),
  paste0("# source_basename=", basename(input)),
  paste0("# source_xgboost_version=", as.character(utils::packageVersion("xgboost")))
)
writeLines(header, temporary, useBytes = TRUE)
utils::write.table(
  nodes, temporary, sep = "\t", quote = FALSE, row.names = FALSE,
  col.names = TRUE, na = "NA", append = TRUE
)

# Validate the stable artifact against the source booster before publication.
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("cannot locate staging script", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_argument), mustWork = TRUE)
repository_root <- dirname(dirname(script_path))
source(file.path(repository_root, "R", "model.R"))
portable <- .read_portable_model(temporary)
set.seed(20260721)
probe <- matrix(stats::runif(1024L * length(feature_names)), ncol = length(feature_names))
colnames(probe) <- feature_names
probe[1L, ] <- 0
probe[2L, ] <- 1
probe[3L, seq.int(1L, length(feature_names), by = 3L)] <- NA_real_
source_prediction <- as.numeric(xgboost::predict(model, probe))
portable_prediction <- .predict_portable_tree(portable, probe)
max_abs_difference <- max(abs(source_prediction - portable_prediction))
if (!is.finite(max_abs_difference) || max_abs_difference > 1e-6) {
  stop("portable-model differential failed: max absolute difference ", max_abs_difference, call. = FALSE)
}

if (file.exists(output) && !unlink(output)) stop("cannot replace output: ", output, call. = FALSE)
if (!file.rename(temporary, output)) stop("could not publish portable model: ", output, call. = FALSE)
receipt <- data.frame(
  format = "xcnv-tree-v1", source_basename = basename(input),
  source_sha256 = source_hash, xgboost_version = as.character(utils::packageVersion("xgboost")),
  trees = length(unique(nodes$tree_id)), nodes = nrow(nodes), probe_rows = nrow(probe),
  max_abs_difference = max_abs_difference, stringsAsFactors = FALSE
)
utils::write.table(
  receipt, paste0(output, ".receipt.tsv"), sep = "\t", quote = FALSE,
  row.names = FALSE, na = ""
)
