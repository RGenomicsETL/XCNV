#' Recommended source authorities for constitutional CNV annotation
#'
#' Return landing pages and release policies for resources used to construct a
#' provenance-labelled CNV annotation pack. This function does not download
#' data. Staging records the selected release, assembly, and source URL.
#'
#' @return A data frame of source authorities and current pinned releases where
#'   the upstream publishes numbered releases.
#' @export
cnv_resource_sources <- function() {
  path <- system.file("extdata", "cnv_resource_sources.tsv", package = "XCNV")
  if (!nzchar(path)) stop("CNV resource source manifest is missing", call. = FALSE)
  data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE, showProgress = FALSE)
}
