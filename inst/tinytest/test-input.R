library(XCNV)

expect_true(validate_cnv(data.frame(chr = "chr1", start = 1, end = 2, type = "gain")))
expect_equal(
  read_cnv(data.frame(chr = "chr1", start = "1", end = "2", type = "gain"))$chr,
  "1"
)
expect_error(validate_cnv(data.frame(chr = "23", start = 1, end = 2, type = "gain")))
expect_error(validate_cnv(data.frame(chr = "1", start = 3, end = 2, type = "gain")))
expect_error(validate_cnv(data.frame(chr = "1", start = 1, end = 2, type = "dup")))

fixture <- xcnv_fixture_resources()
expect_true(inherits(fixture, "xcnv_resources"))
expect_equal(nrow(xcnv_resource_manifest()), 10L)
expect_error(xcnv_resources(tempdir(), require = "all"))
sources <- cnv_resource_sources()
expect_true(all(c("clingen_dosage", "gencode", "gnomad_sv") %in% sources$source_id))
expect_equal(sources$current_release[sources$source_id == "gencode"], "50")

# The old executable's single-row workaround is intentionally observable for a
# record shorter than 50 bases: the altered second row survives the final drop.
short <- annotate_cnv(
  data.frame(chr = "1", start = 100L, end = 110L, type = "gain"), fixture
)
expect_equal(short$Start, 100L)
expect_equal(short$End, 149L)

normal <- read_cnv(data.frame(chr = "1", start = 100L, end = 199L, type = "gain"))
expect_equal(normal$type, "gain")
