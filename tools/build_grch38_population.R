#!/usr/bin/env Rscript

# Build the GRCh38 population CNV relation of criterion 4O of the ACMG/ClinGen
# 2019 lane from the gnomAD-SV 4.1 sites BED, pinned by SHA-256 and staged by
# stage_cnv_resource.R. gnomAD releases are numbered, so this relation is
# built apart from the gene and dosage relations, whose ClinGen source
# changes nightly.
#
# One row per PASS deletion (SVTYPE DEL, mobile-element deletions included)
# or duplication (SVTYPE DUP). `frequency` is the highest allele frequency of
# the continental genetic-ancestry groups with more than 2000 alleles tested:
# African/African American (afr), admixed American (amr), East Asian (eas),
# non-Finnish European (nfe) and South Asian (sas). This is the population
# rule of the ClinGen CNV FAQ for 4O, which applies the ClinGen SVI BA1
# recommendation (Ghosh et al. 2018, PMID 30311383): any continental
# population with more than 2000 alleles. The founder and unassigned groups
# (ami, asj, fin, mid, rmi) are left out, as in that recommendation and in
# gnomAD's own group maximum. A record without such a group, or with
# frequency 0, is left out. BED coordinates are kept as published.
#
# Output, under OUTPUT_DIR:
#   sources/gnomad.v4.1.sv.sites.bed.gz and its receipt   staged input
#   grch38/population.parquet                            population relation
#   grch38/population-sources.tsv                        receipt of the input

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("usage: build_grch38_population.R OUTPUT_DIR", call. = FALSE)
}
for (package in c("DBI", "duckdb", "XCNV")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("building requires ", package, call. = FALSE)
  }
}
output_dir <- normalizePath(arguments[[1L]], winslash = "/", mustWork = FALSE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
stage_script <- file.path(dirname(normalizePath(script)), "stage_cnv_resource.R")

gnomad_release <- "4.1"
url <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.bed.gz"
sha256 <- "5bd4b9855f7b4b23238d76a3c58e098ba419e7c030f0e32970394eec26c73100"
source_file <- file.path(output_dir, "sources", basename(url))
status <- system2(file.path(R.home("bin"), "Rscript"), c(
  shQuote(stage_script), "gnomad_sv", shQuote(url), gnomad_release, "GRCh38",
  sha256, shQuote(source_file)
))
if (!identical(status, 0L)) stop("staging failed: ", url, call. = FALSE)

groups <- c("afr", "amr", "eas", "nfe", "sas")
minimum_alleles <- 2000L

relation_dir <- file.path(output_dir, "grch38")
dir.create(relation_dir, recursive = TRUE, showWarnings = FALSE)
# The connection is XCNV's own overlap connection, with its DuckDB settings.
context <- XCNV:::.new_overlap_context("duckdb")
con <- context$con
quote <- function(x) as.character(DBI::dbQuoteString(con, x))

# Every column is read as text and the used ones are cast by name, so that a
# renamed column stops the build.
eligible <- vapply(groups, function(group) {
  sprintf(
    "CASE WHEN TRY_CAST(\"AN_%1$s\" AS INTEGER) > %2$d THEN TRY_CAST(\"AF_%1$s\" AS DOUBLE) END",
    group, minimum_alleles
  )
}, character(1))
alleles <- vapply(groups, function(group) {
  sprintf("TRY_CAST(\"AN_%s\" AS INTEGER)", group)
}, character(1))
DBI::dbExecute(con, sprintf(
  "CREATE TEMP TABLE sites AS
   SELECT \"#chrom\" AS chr, CAST(\"start\" AS BIGINT) AS start, CAST(\"end\" AS BIGINT) AS \"end\",
     name AS record_id, SVTYPE AS svtype, FILTER AS filter,
     list_value(%s) AS eligible, list_value(%s) AS alleles
   FROM read_csv(%s, delim = '\t', header = true, quote = '', escape = '',
     all_varchar = true, nullstr = 'NA', sample_size = -1)",
  paste(eligible, collapse = ", "), paste(alleles, collapse = ", "), quote(source_file)
))
DBI::dbExecute(con, sprintf(
  "COPY (
     SELECT chr, start, \"end\", record_id,
       CASE svtype WHEN 'DEL' THEN 'loss' ELSE 'gain' END AS cnv_type,
       list_max(eligible) AS frequency,
       list_value(%s)[list_position(eligible, list_max(eligible))] AS population,
       alleles[list_position(eligible, list_max(eligible))] AS allele_number,
       'gnomad_sv' AS source_id, %s AS source_release, 'GRCh38' AS assembly
     FROM sites
     WHERE filter = 'PASS' AND svtype IN ('DEL', 'DUP') AND list_max(eligible) > 0
     ORDER BY chr, start, record_id
   ) TO %s (FORMAT parquet)",
  paste(vapply(groups, quote, character(1)), collapse = ", "),
  quote(gnomad_release), quote(file.path(relation_dir, "population.parquet"))
))
counts <- DBI::dbGetQuery(con,
  "SELECT count(*) AS sites,
     count(*) FILTER (WHERE filter = 'PASS' AND svtype IN ('DEL', 'DUP')) AS pass_cnv
   FROM sites"
)
population <- DBI::dbGetQuery(con, sprintf(
  "SELECT * FROM read_parquet(%s)", quote(file.path(relation_dir, "population.parquet"))
))
invisible(XCNV::derive_cnv_2019_evidence(
  data.frame(cnv_id = character(), chr = character(), start = integer(),
             end = integer(), cnv_type = character()),
  data.frame(chr = "chr1", start = 0L, end = 1L, gene_id = "probe", biotype = "protein_coding",
             source_id = "probe", source_release = "probe", assembly = "GRCh38"),
  data.frame(chr = "chr1", start = 0L, end = 1L, record_id = "probe", hi_score = 0,
             ts_score = 0, source_id = "probe", source_release = "probe", assembly = "GRCh38"),
  population
))
XCNV:::.close_overlap_context(context)
if (!file.copy(
  paste0(source_file, ".receipt.tsv"), file.path(relation_dir, "population-sources.tsv"),
  overwrite = TRUE
)) {
  stop("could not copy the receipt to ", relation_dir, call. = FALSE)
}
message(
  "gnomAD-SV sites: ", counts$sites, "; PASS DEL and DUP: ", counts$pass_cnv,
  "; population rows: ", nrow(population), " (", sum(population$cnv_type == "loss"),
  " losses, ", sum(population$cnv_type == "gain"), " gains; ",
  sum(population$frequency >= 0.01), " at 1 % or more)"
)
