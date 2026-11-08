#!/usr/bin/env Rscript

# Build the GRCh38 gene and ClinGen dosage relations of the ACMG/ClinGen 2019
# lane from pinned sources. Each source is staged by stage_cnv_resource.R,
# which checks its SHA-256 and writes a receipt. ClinGen publishes its
# curation lists nightly at fixed URLs: a new download needs new pins and a new
# release date below. The Ensembl GRCh38 to GRCh37 chain is staged for the
# liftover of the X-CNV lane and copied unchanged to liftover/.
#
# Output, under OUTPUT_DIR:
#   sources/<file> and sources/<file>.receipt.tsv   staged inputs
#   grch38/genes.parquet                           gene relation
#   grch38/dosage.parquet                          dosage relation
#   grch38/sources.tsv                             receipts of all inputs
#   liftover/GRCh38_to_GRCh37.chain.gz             Ensembl chain

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("usage: build_grch38_resources.R OUTPUT_DIR", call. = FALSE)
}
for (package in c("DBI", "duckdb", "XCNV")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("building requires ", package, call. = FALSE)
  }
}
output_dir <- normalizePath(arguments[[1L]], winslash = "/", mustWork = FALSE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
stage_script <- file.path(dirname(normalizePath(script)), "stage_cnv_resource.R")

gencode_release <- "50"
clingen_release <- "2026-11-08"
pins <- data.frame(
  source_id = c("gencode", "clingen_dosage", "clingen_dosage", "ensembl_assembly_mapping"),
  release = c(gencode_release, clingen_release, clingen_release, "76"),
  assembly = c("GRCh38", "GRCh38", "GRCh38", "GRCh38,GRCh37"),
  url = c(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz",
    "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv",
    "https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv",
    "https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh38_to_GRCh37.chain.gz"
  ),
  sha256 = c(
    "5b9ad9aa71f6ab5d55123c51f74723644a895bc0cbaa9cee8fbc11bb9dec0bf7",
    "2b3f8318df9114b01095576cae3b90683d5a295353fe01f6e51b8e7d6184590a",
    "06bd0333e2c94585c3f823a69be220e07e0c0e1e82dc47e6a325fa734ab8d454",
    "6bee5bfdd6defc0fb2bd388ef51819c32bbfae1d3e7a3345aa5ebbdc3dd492b1"
  ),
  stringsAsFactors = FALSE
)
pins$file <- file.path(output_dir, "sources", basename(pins$url))
for (i in seq_len(nrow(pins))) {
  status <- system2(file.path(R.home("bin"), "Rscript"), c(
    shQuote(stage_script), pins$source_id[[i]], shQuote(pins$url[[i]]),
    pins$release[[i]], shQuote(pins$assembly[[i]]), pins$sha256[[i]], shQuote(pins$file[[i]])
  ))
  if (!identical(status, 0L)) stop("staging failed: ", pins$url[[i]], call. = FALSE)
}
source_file <- stats::setNames(pins$file, basename(pins$url))

# ClinGen lists carry five comment lines and a commented header. Check the
# columns read by position so that a format change stops the build.
clingen_columns <- function(path, first) {
  header <- strsplit(sub("^#", "", readLines(path, n = 6L)[[6L]]), "\t", fixed = TRUE)[[1L]]
  expected <- c(first, "Genomic Location", "Haploinsufficiency Score", "Triplosensitivity Score")
  if (!identical(header[c(1L, 2L, 4L, 5L, 13L)], expected)) {
    stop("unexpected ClinGen header in ", path, call. = FALSE)
  }
}
clingen_columns(source_file[["ClinGen_gene_curation_list_GRCh38.tsv"]], c("Gene Symbol", "Gene ID"))
clingen_columns(
  source_file[["ClinGen_region_curation_list_GRCh38.tsv"]], c("ISCA ID", "ISCA Region Name")
)

relation_dir <- file.path(output_dir, "grch38")
dir.create(relation_dir, recursive = TRUE, showWarnings = FALSE)
# The connection is XCNV's own overlap connection, with its DuckDB settings.
context <- XCNV:::.new_overlap_context("duckdb")
con <- context$con
quote <- function(x) as.character(DBI::dbQuoteString(con, x))

# GTF is 1-based closed; BED start is GTF start - 1.
DBI::dbExecute(con, sprintf(
  "COPY (
     SELECT chr, gtf_start - 1 AS start, gtf_end AS \"end\",
       regexp_extract(attributes, 'gene_id \"([^\"]+)\"', 1) AS gene_id,
       regexp_extract(attributes, 'gene_name \"([^\"]+)\"', 1) AS gene_name,
       regexp_extract(attributes, 'gene_type \"([^\"]+)\"', 1) AS biotype,
       strand, 'gencode' AS source_id, %s AS source_release, 'GRCh38' AS assembly
     FROM read_csv(%s, delim = '\t', header = false, comment = '#', quote = '',
       escape = '', auto_detect = false, columns = {'chr': 'VARCHAR', 'source': 'VARCHAR',
       'feature': 'VARCHAR', 'gtf_start': 'INTEGER', 'gtf_end': 'INTEGER',
       'score': 'VARCHAR', 'strand': 'VARCHAR', 'frame': 'VARCHAR',
       'attributes': 'VARCHAR'})
     WHERE feature = 'gene'
     ORDER BY chr, start, gene_id
   ) TO %s (FORMAT parquet)",
  quote(gencode_release), quote(source_file[["gencode.v50.basic.annotation.gtf.gz"]]),
  quote(file.path(relation_dir, "genes.parquet"))
))

# ClinGen locations are chrN:start-end, 1-based closed. Records without a
# GRCh38 location ('tbd' or empty) cannot be placed and are left out.
clingen_select <- function(path, record_type, record_id, name) {
  sprintf(
    "SELECT regexp_extract(column03, '^(chr[^:]+):([0-9]+)-([0-9]+)$', 1) AS chr,
       CAST(regexp_extract(column03, '^(chr[^:]+):([0-9]+)-([0-9]+)$', 2) AS INTEGER) - 1 AS start,
       CAST(regexp_extract(column03, '^(chr[^:]+):([0-9]+)-([0-9]+)$', 3) AS INTEGER) AS \"end\",
       %s AS record_id, %s AS record_type, %s AS name,
       TRY_CAST(column04 AS DOUBLE) AS hi_score, TRY_CAST(column12 AS DOUBLE) AS ts_score,
       'clingen_dosage' AS source_id, %s AS source_release, 'GRCh38' AS assembly
     FROM read_csv(%s, delim = '\t', header = false, skip = 6, quote = '',
       escape = '', all_varchar = true, null_padding = true)
     WHERE regexp_matches(column03, '^chr[^:]+:[0-9]+-[0-9]+$')",
    record_id, quote(record_type), name, quote(clingen_release), quote(path)
  )
}
DBI::dbExecute(con, sprintf(
  "COPY (SELECT * FROM (%s UNION ALL %s) ORDER BY chr, start, record_id) TO %s (FORMAT parquet)",
  clingen_select(
    source_file[["ClinGen_gene_curation_list_GRCh38.tsv"]], "gene", "'NCBIGene:' || column01", "column00"
  ),
  clingen_select(source_file[["ClinGen_region_curation_list_GRCh38.tsv"]], "region", "column00", "column01"
  ),
  quote(file.path(relation_dir, "dosage.parquet"))
))
read_relation <- function(name) {
  DBI::dbGetQuery(con, sprintf(
    "SELECT * FROM read_parquet(%s)", quote(file.path(relation_dir, name))
  ))
}
genes <- read_relation("genes.parquet")
dosage <- read_relation("dosage.parquet")
invisible(XCNV::derive_cnv_2019_evidence(
  data.frame(cnv_id = character(), chr = character(), start = integer(),
             end = integer(), cnv_type = character()),
  genes, dosage
))
XCNV:::.close_overlap_context(context)
receipts <- do.call(rbind, lapply(paste0(pins$file, ".receipt.tsv"), utils::read.delim,
  colClasses = "character"
))
utils::write.table(
  receipts, file.path(relation_dir, "sources.tsv"), sep = "\t", quote = FALSE,
  row.names = FALSE, na = ""
)
liftover_dir <- file.path(output_dir, "liftover")
dir.create(liftover_dir, showWarnings = FALSE)
chain <- source_file[["GRCh38_to_GRCh37.chain.gz"]]
if (!file.copy(chain, file.path(liftover_dir, basename(chain)), overwrite = TRUE)) {
  stop("could not copy the chain to ", liftover_dir, call. = FALSE)
}
message(
  "genes: ", nrow(genes), " rows; dosage: ", nrow(dosage), " rows (",
  sum(dosage$record_type == "gene"), " genes, ", sum(dosage$record_type == "region"), " regions)"
)
