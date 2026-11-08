library(XCNV)

fixture <- xcnv_fixture_resources()
duckhts_data <- system.file("extdata", package = "Rduckhts")
duckhts_liftover <- function(set, source_assembly = "GRCh38") {
  xcnv_liftover(
    chain = file.path(duckhts_data, paste0("liftover_", set, ".chain")),
    source_fasta = file.path(duckhts_data, paste0("liftover_", set, "_src.fa")),
    destination_fasta = file.path(duckhts_data, paste0("liftover_", set, "_dst.fa")),
    source_assembly = source_assembly, destination_assembly = "GRCh37",
    release = paste0("duckhts-test-", set)
  )
}

# The declaration is explicit: every file must exist and the two assemblies
# must differ.
expect_error(xcnv_liftover(
  tempfile(), tempfile(), tempfile(), "GRCh38", "GRCh37", "r1"
))
expect_error(duckhts_liftover("chr23_alias", source_assembly = "GRCh37"))
clip_pad <- duckhts_liftover("clip_pad")
expect_true(inherits(clip_pad, "xcnv_liftover"))
expect_identical(clip_pad$chain_id, "liftover_clip_pad.chain")

# Without a declaration GRCh38 input is still refused.
expect_error(predict_cnv(
  data.frame(chr = "1", start = 1100L, end = 1199L, type = "gain"),
  xcnv_resources(fixture$root, require = "all", assembly = "GRCh38")
))
# The declaration must lift into the resource assembly.
expect_error(predict_cnv(
  data.frame(chr = "1", start = 1100L, end = 1199L, type = "gain"),
  fixture,
  liftover = xcnv_liftover(
    clip_pad$chain, clip_pad$source_fasta, clip_pad$destination_fasta,
    "GRCh37", "GRCh38", "r1"
  )
))

lift_rows <- function(cnvs, liftover) {
  context <- XCNV:::.new_overlap_context("reference", duckhts = TRUE)
  on.exit(XCNV:::.close_overlap_context(context), add = TRUE)
  XCNV:::.xcnv_lift_cnv(cnvs, liftover, context$con)
}

# BED [s, e) becomes a symbolic record at POS s + 1 with END e. The DuckHTS
# clip_pad chain maps source bases 1-8 to 1-8, drops source base 9, and maps
# 10-13 to 9-12. Its contigs are renamed to chromosome 1 on both sides here so
# that the same geometry lifts within one contig.
renamed <- tempfile("xcnv-clip-pad-")
dir.create(renamed)
rename_copy <- function(file) {
  out <- file.path(renamed, file)
  writeLines(gsub("chr[SD]", "1", readLines(file.path(duckhts_data, file))), out)
  out
}
clip_pad_1 <- xcnv_liftover(
  chain = rename_copy("liftover_clip_pad.chain"),
  source_fasta = rename_copy("liftover_clip_pad_src.fa"),
  destination_fasta = rename_copy("liftover_clip_pad_dst.fa"),
  source_assembly = "GRCh38", destination_assembly = "GRCh37",
  release = "duckhts-test-clip_pad"
)
clip_rows <- lift_rows(
  data.frame(
    chr = "1", start = c(1L, 1L, 8L, 5L, 1L), end = c(12L, 7L, 12L, 8L, 9L),
    type = c("loss", "gain", "loss", "loss", "loss"),
    stringsAsFactors = FALSE
  ),
  clip_pad_1
)
expect_identical(
  clip_rows$status,
  c("supported", "supported", "unsupported", "supported", "unsupported")
)
expect_identical(clip_rows$lifted_chr, c("1", "1", NA, "1", NA))
expect_identical(clip_rows$lifted_start, c(1L, 1L, NA, 5L, NA))
# [1, 12) loses the dropped base: [1, 11) with a length ratio of 10/11.
expect_identical(clip_rows$lifted_end, c(11L, 7L, NA, 8L, NA))
expect_equal(clip_rows$liftover_length_ratio, c(10 / 11, 1, NA, 1, NA))
# POS 9 and END 9 fall in the dropped source base.
expect_identical(
  clip_rows$liftover_reason,
  c(NA, NA, "start_unmapped", NA, "end_unmapped")
)
expect_identical(clip_rows$liftover_release, rep("duckhts-test-clip_pad", 5L))

# The unrenamed chain lifts chrS to another contig, chrD.
expect_identical(
  lift_rows(data.frame(chr = "chrS", start = 1L, end = 12L, type = "loss"), clip_pad)$liftover_reason,
  "contig_changed"
)

# The DuckHTS chr23_alias chain lifts X to a different contig, chrLiftX, and
# has no chain for chromosome 1.
alias_rows <- lift_rows(
  data.frame(chr = c("X", "1"), start = c(1L, 1L), end = c(8L, 8L),
             type = c("gain", "gain"), stringsAsFactors = FALSE),
  duckhts_liftover("chr23_alias")
)
expect_identical(alias_rows$liftover_reason, c("contig_changed", "source_contig_absent"))
expect_true(all(is.na(alias_rows$lifted_chr)))

# End to end on the synthetic fixture. Source chromosome 1 [1000, 1400) maps
# forward to destination chromosome 1 [0, 400); chromosome 2 maps to the
# reverse strand; chromosome 3 maps to chromosome 4.
staged <- tempfile("xcnv-liftover-")
dir.create(staged)
chain <- file.path(staged, "test.chain")
writeLines(c(
  "chain 1000 1 2000 + 1000 1400 1 2000 + 0 400 1", "400", "",
  "chain 1000 2 1000 + 0 1000 2 1000 - 0 1000 2", "1000", "",
  "chain 1000 3 1000 + 0 1000 4 1000 + 0 1000 3", "1000", ""
), chain)
fasta <- function(path, contigs, sizes) {
  writeLines(unlist(lapply(seq_along(contigs), function(i) {
    c(paste0(">", contigs[[i]]), strrep("A", sizes[[i]]))
  })), path)
  path
}
liftover <- xcnv_liftover(
  chain = chain,
  source_fasta = fasta(file.path(staged, "source.fa"), c("1", "2", "3"), c(2000, 1000, 1000)),
  destination_fasta = fasta(file.path(staged, "destination.fa"), c("1", "2", "4"), c(2000, 1000, 1000)),
  source_assembly = "GRCh38", destination_assembly = "GRCh37",
  release = "test-chain-1"
)
grch38 <- data.frame(
  chr = c("3", "1", "2", "1", "1"),
  start = c(100L, 1300L, 100L, 1100L, 1500L),
  end = c(200L, 1340L, 200L, 1199L, 1600L),
  type = c("gain", "loss", "gain", "gain", "loss"),
  stringsAsFactors = FALSE
)
lifted <- predict_cnv(grch38, fixture, liftover = liftover)
direct <- predict_cnv(
  data.frame(chr = c("1", "1"), start = c(100L, 300L), end = c(199L, 340L),
             type = c("gain", "loss")),
  fixture
)
expect_identical(lifted$Chr, c("1", "1", "1", "2", "3"))
expect_identical(lifted$Start, c(1100L, 1300L, 1500L, 100L, 100L))
expect_identical(lifted$End, c(1199L, 1340L, 1600L, 200L, 200L))
expect_true(all(lifted$assembly == "GRCh38"))
expect_true(all(lifted$lifted_assembly == "GRCh37"))
expect_identical(lifted$status, c("supported", "supported", rep("unsupported", 3L)))
expect_identical(
  lifted$liftover_reason,
  c(NA, NA, "start_unmapped", "strand_reversed", "contig_changed")
)
expect_identical(lifted$lifted_chr, c("1", "1", NA, NA, NA))
expect_identical(lifted$lifted_start, c(100L, 300L, NA, NA, NA))
expect_identical(lifted$lifted_end, c(199L, 340L, NA, NA, NA))
expect_equal(lifted$liftover_length_ratio, c(1, 1, NA, NA, NA))
expect_true(all(lifted$liftover_chain == "test.chain"))
expect_true(all(lifted$liftover_release == "test-chain-1"))
# Supported rows score exactly as the same interval given in GRCh37; the
# unsupported rows have no features and no MVP_score.
expect_equal(lifted$MVP_score[1:2], direct$MVP_score, tolerance = 1e-12)
expect_equal(lifted$pLI[1:2], direct$pLI)
expect_true(all(is.na(lifted$MVP_score[3:5])))
expect_true(all(is.na(lifted$pLI[3:5])))
expect_identical(names(lifted)[length(names(lifted))], "MVP_score")
expect_true(all(lifted$resource_release == direct$resource_release[[1L]]))

# The reference overlap backend gives the same lifted result.
expect_equal(
  predict_cnv(grch38, fixture, overlap_backend = "reference", liftover = liftover),
  lifted
)
# A call set with no liftable row keeps its rows and scores none of them.
none <- predict_cnv(grch38[c(3L, 5L), ], fixture, liftover = liftover)
expect_identical(none$status, c("unsupported", "unsupported"))
expect_true(all(is.na(none$MVP_score)))
