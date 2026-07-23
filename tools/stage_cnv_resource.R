#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    "usage: stage_cnv_resource.R SOURCE_ID URL RELEASE ASSEMBLY EXPECTED_SHA256 OUTPUT",
    call. = FALSE
  )
}
names(args) <- c("source_id", "url", "release", "assembly", "expected_sha256", "output")
if (!requireNamespace("openssl", quietly = TRUE)) stop("staging requires openssl", call. = FALSE)
if (!grepl("^[0-9a-fA-F]{64}$", args[["expected_sha256"]])) {
  stop("EXPECTED_SHA256 must contain 64 hexadecimal characters", call. = FALSE)
}
output <- normalizePath(args[["output"]], winslash = "/", mustWork = FALSE)
if (!dir.exists(dirname(output))) dir.create(dirname(output), recursive = TRUE)
temporary <- tempfile("xcnv-resource-", tmpdir = dirname(output))
on.exit(unlink(temporary), add = TRUE)
utils::download.file(args[["url"]], temporary, mode = "wb", quiet = FALSE)
connection <- file(temporary, open = "rb")
actual_sha256 <- as.character(openssl::sha256(connection))
close(connection)
if (!identical(tolower(actual_sha256), tolower(args[["expected_sha256"]]))) {
  stop("SHA-256 mismatch: expected ", args[["expected_sha256"]], ", observed ", actual_sha256, call. = FALSE)
}
if (file.exists(output) && !unlink(output)) stop("cannot replace output: ", output, call. = FALSE)
if (!file.rename(temporary, output)) stop("could not publish staged resource", call. = FALSE)
receipt <- data.frame(
  source_id = args[["source_id"]], release = args[["release"]],
  assembly = args[["assembly"]], source_url = args[["url"]],
  sha256 = tolower(actual_sha256), bytes = file.info(output)$size,
  retrieved_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)
utils::write.table(
  receipt, paste0(output, ".receipt.tsv"), sep = "\t", quote = FALSE,
  row.names = FALSE, na = ""
)
