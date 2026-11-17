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
# on.exit() does not run at the top level of an Rscript; remove the partial
# download explicitly when the transfer or the checksum fails. R's timeout
# bounds the whole transfer (60 s by default), which a gigabyte source exceeds.
options(timeout = max(3600, getOption("timeout")))
tryCatch(
  utils::download.file(args[["url"]], temporary, mode = "wb", quiet = FALSE),
  error = function(condition) {
    unlink(temporary)
    stop("download failed: ", conditionMessage(condition), call. = FALSE)
  }
)
connection <- file(temporary, open = "rb")
actual_sha256 <- tolower(unclass(as.character(openssl::sha256(connection))))
close(connection)
if (!identical(actual_sha256, tolower(args[["expected_sha256"]]))) {
  unlink(temporary)
  stop("SHA-256 mismatch: expected ", args[["expected_sha256"]], ", observed ", actual_sha256, call. = FALSE)
}
# The receipt is written before the resource is published, and the previous
# receipt is removed first, so that a published resource never sits beside
# the receipt of another download. A failure after the removal leaves a
# resource without a receipt, which its consumers reject.
receipt_file <- paste0(output, ".receipt.tsv")
receipt_temporary <- paste0(temporary, ".receipt.tsv")
receipt <- data.frame(
  source_id = args[["source_id"]], release = args[["release"]],
  assembly = args[["assembly"]], source_url = args[["url"]],
  sha256 = actual_sha256, bytes = file.info(temporary)$size,
  retrieved_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)
tryCatch(
  utils::write.table(
    receipt, receipt_temporary, sep = "\t", quote = FALSE,
    row.names = FALSE, na = ""
  ),
  error = function(condition) {
    unlink(c(temporary, receipt_temporary))
    stop("could not write receipt: ", conditionMessage(condition), call. = FALSE)
  }
)
if (file.exists(receipt_file) && !file.remove(receipt_file)) {
  unlink(c(temporary, receipt_temporary))
  stop("could not remove the previous receipt", call. = FALSE)
}
# file.rename() replaces an existing output.
if (!file.rename(temporary, output)) {
  unlink(c(temporary, receipt_temporary))
  stop("could not publish staged resource", call. = FALSE)
}
if (!file.rename(receipt_temporary, receipt_file)) {
  unlink(receipt_temporary)
  stop("could not publish the receipt of ", output, call. = FALSE)
}
