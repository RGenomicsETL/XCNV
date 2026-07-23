library(XCNV)

a <- data.frame(
  chr = c("1", "1", "2", "3"), start = c(10L, 30L, 1L, 8L),
  end = c(20L, 50L, 10L, 8L)
)
b <- data.frame(
  chr = c("1", "1", "2"), start = c(15L, 45L, 5L),
  end = c(25L, 55L, 8L)
)

context <- XCNV:::.new_overlap_context("duckdb")
observed <- XCNV:::.bed_overlap_pairs(a, b, context = context)
expected <- XCNV:::.bed_overlap_pairs_reference(a, b)
expect_equal(observed, expected)

observed_reciprocal <- XCNV:::.bed_overlap_pairs(
  a, b, min_a = 0.5, min_b = 0.5, context = context
)
expected_reciprocal <- XCNV:::.bed_overlap_pairs_reference(
  a, b, min_a = 0.5, min_b = 0.5
)
expect_equal(observed_reciprocal, expected_reciprocal)

fixture <- xcnv_fixture_resources()
input <- data.frame(
  chr = c("1", "1"), start = c(100L, 300L), end = c(199L, 340L),
  type = c("gain", "loss")
)
expect_equal(
  annotate_cnv(input, fixture, overlap_backend = "duckdb"),
  annotate_cnv(input, fixture, overlap_backend = "reference")
)

published_sites <- tempfile()
writeLines(c("1\t120\t120\t2", "1\t140\t145\t1"), published_sites)
sites <- XCNV:::.read_ljb_sites(published_sites)
expect_equal(sites$start, c(119L, 140L))
expect_equal(sites$end, c(120L, 145L))
expect_equal(sites$site_id, c(2L, 1L))

score_names <- c(
  "SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred",
  "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred",
  "RadialSVM_pred", "LR_pred", "VEST3_score", "CADD_phred", "GERP++_RS",
  "phyloP46way_placental", "phyloP100way_vertebrate", "SiPhy_29way_logOdds"
)
scores <- as.data.frame(matrix(c(rep(0, 15L), rep(1, 15L)), nrow = 2L, byrow = TRUE))
names(scores) <- score_names
point_result <- XCNV:::.compute_ljb(
  data.frame(chr = "1", start = 119L, end = 120L),
  list(ljb_sites = sites[1L, ], ljb_scores = scores),
  context
)
expect_equal(point_result$SIFT_pred, 1)

headerless_samples <- tempfile()
writeLines(c("1\tS1\tAFR", "2\tS2\tEAS"), headerless_samples)
sample_rows <- XCNV:::.read_merged_sample_info(headerless_samples)
expect_equal(nrow(sample_rows), 2L)
expect_equal(sample_rows$sample_id, c("S1", "S2"))

XCNV:::.close_overlap_context(context)
