#' Load the bundled real X-CNV example
#'
#' Return ClinVar allele 72461, a pathogenic deletion containing `SOX2`, with
#' the case-limited GRCh37 annotation relations needed to run X-CNV and the
#' objective ACMG/ClinGen 2019 content-evidence lane. The annotation files are
#' subsets for this case and must not be used for unrelated CNVs.
#'
#' @return A list containing `cnv`, `clinvar`, `resources`, `genes`, and
#'   `dosage`.
#' @export
xcnv_real_example <- function() {
  root <- system.file("extdata", "real_case", package = "XCNV")
  if (!nzchar(root)) {
    stop("the real X-CNV example is not available in this installation", call. = FALSE)
  }
  case <- utils::read.delim(
    file.path(root, "case.tsv"), check.names = FALSE,
    stringsAsFactors = FALSE
  )
  genes <- utils::read.delim(
    file.path(root, "case.genes.tsv"), check.names = FALSE,
    stringsAsFactors = FALSE
  )
  dosage <- utils::read.delim(
    file.path(root, "case.dosage.tsv"), check.names = FALSE,
    stringsAsFactors = FALSE
  )
  cnv <- data.frame(
    chr = as.character(case$grch37_chrom),
    start = as.integer(case$grch37_start0),
    end = as.integer(case$grch37_end0),
    type = as.character(case$cnv_type),
    stringsAsFactors = FALSE
  )
  list(
    cnv = cnv,
    clinvar = case,
    resources = xcnv_resources(root, require = "annotations"),
    genes = genes,
    dosage = dosage
  )
}
