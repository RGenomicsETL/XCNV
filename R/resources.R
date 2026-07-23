.xcnv_resource_spec <- data.frame(
  key = c(
    "genome", "ljb_sites", "ljb_scores", "cdts", "merged_sites",
    "sample_info", "merged_sample_info", "encode", "hi", "model"
  ),
  file = c(
    "genome.txt", "hg19_ljb26_all_converted_sites.vcf",
    "hg19_ljb26_all_converted_scores.txt", "CDTS_percentile.txt",
    "merged.cnv.sites.bed", "sample.info.txt",
    "merged.cnv.sample.info.txt", "hg19-ccREs.bed",
    "gencode_v19_features.bed", "xcnv.portable.model.tsv"
  ),
  group = c(
    rep("annotations", 9L), "model"
  ),
  stringsAsFactors = FALSE
)
.xcnv_model_candidates <- c("xcnv.portable.model.tsv", "xcnv.fixture.model.tsv")

#' Show the X-CNV resource contract
#'
#' Return the file names and resource groups expected by the package. The
#' trained model and all large annotations are external resources; they are not
#' downloaded or installed automatically.
#'
#' @return A data frame describing required resource keys, file names, and
#'   groups.
#' @export
xcnv_resource_manifest <- function() {
  .xcnv_resource_spec
}

.resolve_resource_root <- function(path = NULL) {
  if (is.null(path)) {
    path <- Sys.getenv("XCNV_RESOURCE_DIR", unset = "")
    if (!nzchar(path)) {
      stop(
        "an explicit resource directory is required (or set XCNV_RESOURCE_DIR)",
        call. = FALSE
      )
    }
  }
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("'path' must be one existing resource directory", call. = FALSE)
  }
  if (!dir.exists(path)) {
    stop("resource directory does not exist: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

#' Validate an X-CNV resource bundle
#'
#' Check that a local directory contains the named X-CNV annotation and/or
#' model resources. Validation never downloads files and never changes the
#' working directory.
#'
#' @param path Resource directory. If omitted, `XCNV_RESOURCE_DIR` must be set.
#' @param require Resource groups to require: `"annotations"`, `"model"`, or
#'   `"all"`.
#' @return An object of class `xcnv_resources`, containing the normalized root
#'   and resolved resource paths.
#' @export
validate_xcnv_resources <- function(path = NULL, require = c("annotations", "model")) {
  require <- match.arg(require, c("annotations", "model", "all"), several.ok = TRUE)
  if ("all" %in% require) {
    require <- c("annotations", "model")
  }
  root <- .resolve_resource_root(path)
  needed <- .xcnv_resource_spec$key[.xcnv_resource_spec$group %in% require]
  files <- vector("list", length(.xcnv_resource_spec$key))
  names(files) <- .xcnv_resource_spec$key
  missing <- character()
  for (key in .xcnv_resource_spec$key) {
    if (key == "model") {
      candidates <- file.path(root, .xcnv_model_candidates)
      hit <- candidates[file.exists(candidates)]
      files[[key]] <- if (length(hit)) hit[[1L]] else NA_character_
    } else {
      hit <- file.path(root, .xcnv_resource_spec$file[.xcnv_resource_spec$key == key])
      files[[key]] <- if (file.exists(hit)) hit else NA_character_
    }
    if (key %in% needed && is.na(files[[key]])) {
      missing <- c(missing, .xcnv_resource_spec$file[.xcnv_resource_spec$key == key])
    }
  }
  if (length(missing)) {
    stop(
      "resource bundle is missing required file(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  structure(list(root = root, files = files, required = require), class = "xcnv_resources")
}

#' Discover an X-CNV resource bundle
#'
#' Resolve and validate a resource directory. With no `path`, only the
#' explicit `XCNV_RESOURCE_DIR` environment variable is consulted; the current
#' working directory is never searched implicitly.
#'
#' @param path Resource directory, or `NULL` to use `XCNV_RESOURCE_DIR`.
#' @param require Resource groups to require.
#' @return An object of class `xcnv_resources`.
#' @export
xcnv_resources <- function(path = NULL, require = c("annotations", "model")) {
  validate_xcnv_resources(path = path, require = require)
}

#' Locate the bundled deterministic fixture
#'
#' Resolve the small, non-scientific fixture shipped for tests and examples.
#' It is not a substitute for the Zhang et al. annotations or model.
#'
#' @return An object of class `xcnv_resources`.
#' @export
xcnv_fixture_resources <- function() {
  path <- system.file("extdata", "fixture", package = "XCNV")
  if (!nzchar(path)) {
    stop("the XCNV fixture is not available in this installation", call. = FALSE)
  }
  xcnv_resources(path, require = "all")
}
