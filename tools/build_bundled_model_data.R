#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) > 1L) {
  stop("usage: build_bundled_model_data.R [PORTABLE_MODEL_TSV]", call. = FALSE)
}

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) {
  stop("cannot locate build_bundled_model_data.R", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", script_argument), mustWork = TRUE)
repository_root <- dirname(dirname(script_path))
source(file.path(repository_root, "R", "model.R"))

model_path <- if (length(arguments)) {
  normalizePath(arguments[[1L]], winslash = "/", mustWork = TRUE)
} else {
  file.path(repository_root, "tools", "xcnv.portable.model.tsv")
}
xcnv_portable_model <- .read_portable_model(model_path)

output <- file.path(repository_root, "data", "xcnv_portable_model.rda")
save(xcnv_portable_model, file = output, compress = "xz", version = 3)
