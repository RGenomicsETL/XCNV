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
XCNV:::.close_overlap_context(context)
