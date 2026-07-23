#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 5L) {
  stop(
    "usage: build_real_example_resources.R SOURCE_DIR OUTPUT_DIR CHROM START END",
    call. = FALSE
  )
}

source_dir <- normalizePath(arguments[[1L]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(arguments[[2L]], winslash = "/", mustWork = FALSE)
chrom <- sub("^chr", "", arguments[[3L]])
start <- suppressWarnings(as.integer(arguments[[4L]]))
end <- suppressWarnings(as.integer(arguments[[5L]]))
if (is.na(start) || is.na(end) || start >= end) {
  stop("START and END must define one non-empty half-open interval", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source_path <- function(name) {
  path <- file.path(source_dir, name)
  if (!file.exists(path)) stop("missing source resource: ", path, call. = FALSE)
  path
}

copy_resource <- function(name) {
  if (!file.copy(source_path(name), file.path(output_dir, name), overwrite = TRUE)) {
    stop("could not copy resource: ", name, call. = FALSE)
  }
}

stream_lines <- function(path, callback, compressed = FALSE) {
  input <- if (compressed) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(input))
  repeat {
    lines <- readLines(input, n = 100000L, warn = FALSE)
    if (!length(lines)) break
    callback(lines)
  }
}

subset_bed <- function(input_name, output_name, columns = 3L) {
  output <- file(file.path(output_dir, output_name), open = "wt")
  on.exit(close(output))
  stream_lines(source_path(input_name), function(lines) {
    fields <- strsplit(lines[!startsWith(lines, "#")], "\t", fixed = TRUE)
    keep <- vapply(fields, function(row) {
      if (length(row) < columns || row[[1L]] != chrom) return(FALSE)
      row_start <- suppressWarnings(as.integer(row[[2L]]))
      row_end <- suppressWarnings(as.integer(row[[3L]]))
      !is.na(row_start) && !is.na(row_end) && row_start < end && row_end > start
    }, logical(1))
    if (any(keep)) writeLines(lines[!startsWith(lines, "#")][keep], output)
  })
}

combined_name <- "hg19_ljb26_all_converted.vcf.gz"
sites_output <- file(
  file.path(output_dir, "hg19_ljb26_all_converted_sites.vcf"),
  open = "wt"
)
scores_output <- file(
  file.path(output_dir, "hg19_ljb26_all_converted_scores.txt"),
  open = "wt"
)
site_id <- 0L
header_seen <- FALSE
stream_lines(source_path(combined_name), function(lines) {
  if (!header_seen) {
    header <- lines[startsWith(lines, "#CHROM\t")]
    if (length(header)) {
      fields <- strsplit(sub("^#", "", header[[1L]]), "\t", fixed = TRUE)[[1L]]
      if (length(fields) < 18L) stop("combined LJB26 header has too few columns", call. = FALSE)
      writeLines(paste(fields[4:18], collapse = "\t"), scores_output)
      header_seen <<- TRUE
    }
  }
  records <- lines[!startsWith(lines, "#")]
  fields <- strsplit(records, "\t", fixed = TRUE)
  for (row in fields) {
    if (length(row) < 18L || row[[1L]] != chrom) next
    pos <- suppressWarnings(as.integer(row[[2L]]))
    row_end <- suppressWarnings(as.integer(row[[3L]]))
    if (is.na(pos) || is.na(row_end) || pos <= start || pos > end) next
    site_id <<- site_id + 1L
    writeLines(paste(chrom, pos, row_end, site_id, sep = "\t"), sites_output)
    writeLines(paste(row[4:18], collapse = "\t"), scores_output)
  }
}, compressed = TRUE)
close(sites_output)
close(scores_output)
if (!header_seen) stop("combined LJB26 header was not found", call. = FALSE)

subset_bed("CDTS_percentile.txt", "CDTS_percentile.txt", columns = 4L)
subset_bed("merged.cnv.sites.bed", "merged.cnv.sites.bed", columns = 5L)
subset_bed("hg19-ccREs.bed", "hg19-ccREs.bed", columns = 4L)
subset_bed("gencode_v19_features.bed", "gencode_v19_features.bed", columns = 7L)

site_lines <- readLines(file.path(output_dir, "merged.cnv.sites.bed"), warn = FALSE)
site_fields <- strsplit(site_lines, "\t", fixed = TRUE)
site_ids <- vapply(site_fields, `[[`, character(1), 5L)
sample_output <- file(
  file.path(output_dir, "merged.cnv.sample.info.txt"),
  open = "wt"
)
stream_lines(source_path("merged.cnv.sample.info.txt"), function(lines) {
  ids <- sub("\t.*$", "", lines)
  keep <- ids %in% site_ids
  if (any(keep)) writeLines(lines[keep], sample_output)
})
close(sample_output)

copy_resource("genome.txt")
copy_resource("sample.info.txt")
