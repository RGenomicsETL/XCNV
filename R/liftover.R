.xcnv_liftover_reasons <- c(
  "source_contig_absent", "start_unmapped", "end_unmapped",
  "end_before_start", "contig_changed", "strand_reversed"
)

.xcnv_existing_file <- function(path, label) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("'", label, "' must be one existing file path", call. = FALSE)
  }
  if (!file.exists(path) || isTRUE(file.info(path)$isdir)) {
    stop("'", label, "' does not exist: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.xcnv_label <- function(value, label) {
  if (!(is.character(value) || is.factor(value)) || length(value) != 1L ||
      is.na(value) || !nzchar(trimws(as.character(value)))) {
    stop("'", label, "' must be one nonempty label", call. = FALSE)
  }
  as.character(value)
}

#' Declare a coordinate liftover for X-CNV input
#'
#' Declare the chain and reference sequences that lift CNV calls from another
#' assembly into the assembly of the X-CNV resources. XCNV never downloads,
#' bundles, or infers a liftover: the caller stages the files and states their
#' assemblies and release. The declaration is passed to `predict_cnv()`.
#'
#' The lift is computed by the DuckHTS `bcftools_liftover()` function, a port
#' of the bcftools `+liftover` plugin. Each CNV becomes one symbolic VCF record:
#' `<DEL>` for a loss and `<DUP>` for a gain. The record position is the first
#' affected base, `start + 1` in VCF's 1-based convention, and `END` is the last
#' affected base, which equals the BED `end`. A lifted record with position `p`
#' and `END` `e` is returned as the BED interval `[p - 1, e)`.
#'
#' @param chain Chain file, in the UCSC chain format, from `source_assembly`
#'   to `destination_assembly`; plain or gzip-compressed.
#' @param source_fasta FASTA reference of `source_assembly`.
#' @param destination_fasta FASTA reference of `destination_assembly`.
#' @param source_assembly Assembly of the CNV calls, for example `"GRCh38"`.
#' @param destination_assembly Assembly of the X-CNV resources, for example
#'   `"GRCh37"`.
#' @param release Release label of the chain, carried into every output row.
#' @return An object of class `xcnv_liftover`.
#' @export
xcnv_liftover <- function(
  chain, source_fasta, destination_fasta, source_assembly,
  destination_assembly, release
) {
  source_assembly <- .xcnv_canonical_assembly(
    .xcnv_label(source_assembly, "source_assembly")
  )
  destination_assembly <- .xcnv_canonical_assembly(
    .xcnv_label(destination_assembly, "destination_assembly")
  )
  if (identical(source_assembly, destination_assembly)) {
    stop("a liftover needs different source and destination assemblies", call. = FALSE)
  }
  chain <- .xcnv_existing_file(chain, "chain")
  structure(
    list(
      chain = chain,
      chain_id = basename(chain),
      source_fasta = .xcnv_existing_file(source_fasta, "source_fasta"),
      destination_fasta = .xcnv_existing_file(destination_fasta, "destination_fasta"),
      source_assembly = source_assembly,
      destination_assembly = destination_assembly,
      release = .xcnv_label(release, "release")
    ),
    class = "xcnv_liftover"
  )
}

.xcnv_check_liftover <- function(liftover, resources) {
  if (!inherits(liftover, "xcnv_liftover")) {
    stop("'liftover' must be an xcnv_liftover() declaration", call. = FALSE)
  }
  if (!identical(liftover$destination_assembly, resources$assembly)) {
    stop(
      "liftover destination assembly ", liftover$destination_assembly,
      " does not match the X-CNV resource assembly ", resources$assembly,
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Lift validated CNV rows with DuckHTS bcftools_liftover() on a connection
# that has the extension loaded. Returns one row per input row, in input order.
.xcnv_lift_cnv <- function(cnvs, liftover, con) {
  n <- nrow(cnvs)
  out <- data.frame(
    chr = as.character(cnvs$chr), start = as.integer(cnvs$start),
    end = as.integer(cnvs$end), type = as.character(cnvs$type),
    lifted_assembly = rep.int(liftover$destination_assembly, n),
    lifted_chr = rep.int(NA_character_, n),
    lifted_start = rep.int(NA_integer_, n),
    lifted_end = rep.int(NA_integer_, n),
    liftover_chain = rep.int(liftover$chain_id, n),
    liftover_release = rep.int(liftover$release, n),
    liftover_length_ratio = rep.int(NA_real_, n),
    liftover_reason = rep.int(NA_character_, n),
    status = rep.int("supported", n),
    stringsAsFactors = FALSE
  )
  if (!n) return(out)
  records <- data.frame(
    row_id = seq_len(n),
    chrom = out$chr,
    pos = as.numeric(out$start) + 1,
    alt = ifelse(out$type == "loss", "<DEL>", "<DUP>"),
    end_pos = as.numeric(out$end),
    stringsAsFactors = FALSE
  )
  table_name <- "xcnv_liftover_records"
  DBI::dbWriteTable(con, table_name, records, temporary = TRUE, overwrite = TRUE)
  on.exit(DBI::dbRemoveTable(con, table_name), add = TRUE)
  # max_snp_gap and max_indel_inc are NULL: they govern sequence-resolved
  # alleles, not symbolic records, so the DuckHTS defaults apply unchanged.
  query <- paste0(
    "SELECT row_id, unnest(lifted) FROM (SELECT row_id, bcftools_liftover(",
    "chrom, CAST(pos AS BIGINT), 'N', alt, ",
    DBI::dbQuoteString(con, liftover$chain), ", ",
    DBI::dbQuoteString(con, liftover$destination_fasta), ", ",
    DBI::dbQuoteString(con, liftover$source_fasta), ", ",
    "CAST(NULL AS INTEGER), CAST(NULL AS INTEGER), false, ",
    "CAST(end_pos AS BIGINT), true) AS lifted FROM ",
    DBI::dbQuoteIdentifier(con, table_name), ") ORDER BY row_id"
  )
  lifted <- DBI::dbGetQuery(con, query)
  reject <- ifelse(is.na(lifted$reject_reason), "", lifted$reject_reason)
  note <- ifelse(is.na(lifted$note), "", lifted$note)
  mapped <- !is.na(lifted$mapped) & lifted$mapped
  dest_chr <- sub("^chr", "", lifted$dest_chrom)
  reason <- rep.int(NA_character_, n)
  reason[!mapped] <- ifelse(
    grepl("MissingContig", reject[!mapped], fixed = TRUE),
    "source_contig_absent", "start_unmapped"
  )
  start_lifted <- is.na(reason)
  end_failed <- start_lifted & is.na(lifted$dest_end)
  reason[end_failed] <- ifelse(
    grepl("EndLiftOrderViolation", note[end_failed], fixed = TRUE),
    "end_before_start", "end_unmapped"
  )
  reason[is.na(reason) & dest_chr != out$chr] <- "contig_changed"
  reason[is.na(reason) & lifted$reverse_complemented] <- "strand_reversed"
  ok <- is.na(reason)
  out$liftover_reason <- reason
  out$status[!ok] <- "unsupported"
  out$lifted_chr[ok] <- dest_chr[ok]
  out$lifted_start[ok] <- as.integer(lifted$dest_pos[ok] - 1)
  out$lifted_end[ok] <- as.integer(lifted$dest_end[ok])
  out$liftover_length_ratio[ok] <-
    (out$lifted_end[ok] - out$lifted_start[ok]) / (out$end[ok] - out$start[ok])
  out
}

# Annotate the supported lifted rows and place their features at the input
# rows. Unsupported rows keep missing features and receive no model score.
.xcnv_lifted_features <- function(lift, resources, context) {
  supported <- which(lift$status == "supported")
  # Present the rows in X-CNV order so that annotation keeps their order.
  supported <- supported[order(
    match(lift$lifted_chr[supported], c(as.character(seq_len(22L)), "X", "Y")),
    lift$lifted_start[supported], lift$lifted_end[supported], supported
  )]
  annotated <- .annotate_cnv_with(
    data.frame(
      chr = lift$lifted_chr[supported], start = lift$lifted_start[supported],
      end = lift$lifted_end[supported], type = lift$type[supported],
      stringsAsFactors = FALSE
    ),
    resources, context
  )
  features <- annotated[rep.int(NA_integer_, nrow(lift)), , drop = FALSE]
  rownames(features) <- NULL
  features[supported, ] <- annotated
  features$Chr <- lift$chr
  features$Start <- lift$start
  features$End <- lift$end
  features$Type <- lift$type
  features
}
