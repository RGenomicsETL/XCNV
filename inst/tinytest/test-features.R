library(XCNV)
fixture <- xcnv_fixture_resources()
input <- data.frame(
  chr = c("1", "1"), start = c(100L, 300L), end = c(199L, 340L),
  type = c("gain", "loss")
)
features <- annotate_cnv(input, fixture)

expected_model <- c(
  "SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred",
  "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred",
  "RadialSVM_pred", "LR_pred", "VEST3_score", "CADD_phred", "GERP++_RS",
  "phyloP46way_placental", "phyloP100way_vertebrate", "SiPhy_29way_logOdds",
  "CDTS_1st", "CDTS_5th", "gain.freq", "loss.freq", "pELS", "CTCF-bound",
  "PLS", "dELS", "CTCF-only", "DNase-H3K4me3", "pLI", "Episcore", "GHIS",
  "Type.1", "Length"
)
expect_true(all(expected_model %in% names(features)))
expect_equal(features$Type, c("gain", "loss"))
expect_equal(features$Type.1, c(1L, 0L))
expect_equal(features$Length, c(100L, 41L))
expect_equal(features$SIFT_pred, c(0.5, 1))
expect_equal(features$CDTS_1st, c(0.1, 0))
expect_equal(features$CDTS_5th, c(0.2, 10 / 41), tolerance = 1e-12)
expect_equal(features$pELS, c(0.1, 0))
expect_equal(features[["CTCF-bound"]], c(0.1, 0), tolerance = 1e-12)
expect_equal(features$PLS, c(0.1, 0))
expect_equal(features$dELS, c(0, 10 / 41), tolerance = 1e-12)
expect_equal(features[["CTCF-only"]], c(0, 10 / 41), tolerance = 1e-12)
expect_equal(features$pLI, c(0.8, 0.1))
expect_equal(features$Episcore, c(0.4, 0.9))
expect_equal(features$GHIS, c(0.2, 0.3))
expect_equal(features$gain_freq_AFR, c(0.5, 0))
expect_equal(features$gain_freq_EAS, c(1, 0))
expect_equal(features$loss_freq_EAS, c(0, 1))
expect_equal(features$gain.freq, c(2 / 3, 0), tolerance = 1e-12)
expect_equal(features$loss.freq, c(0, 1 / 3), tolerance = 1e-12)

predictions <- predict_cnv(input, fixture)
expect_equal(predictions$MVP_score, stats::plogis(c(0.5, -0.59)), tolerance = 1e-12)
expect_equal(names(predictions)[length(names(predictions))], "MVP_score")

model <- load_xcnv_model(fixture$files$model)
expect_true(inherits(model, "xcnv_portable_model"))
expect_equal(model$model_type, "xcnv-tree-v1")
bundled_model <- xcnv_bundled_model()
expect_true(inherits(bundled_model, "xcnv_portable_model"))
expect_equal(nrow(bundled_model$nodes), 41L)
expect_identical(load_xcnv_model(), bundled_model)
model_data <- new.env(parent = emptyenv())
utils::data("xcnv_portable_model", package = "XCNV", envir = model_data)
expect_identical(model_data$xcnv_portable_model, bundled_model)

real_example <- xcnv_real_example()
expect_identical(real_example$clinvar$clinvar_allele_id, 72461L)
real_prediction <- predict_cnv(real_example$cnv, real_example$resources)
expect_equal(real_prediction$MVP_score, 0.8537559, tolerance = 1e-7)
expect_equal(real_prediction$pLI, 0.396696528099993, tolerance = 1e-12)
real_cnv <- transform(
  real_example$cnv,
  cnv_id = real_example$clinvar$cnv_id,
  cnv_type = type
)
real_evidence <- derive_cnv_2019_evidence(
  real_cnv, real_example$genes, real_example$dosage
)
expect_identical(sort(real_evidence$criterion), c("1A", "2A", "3A"))
expect_identical(
  score_cnv_2019(real_evidence)$classification,
  "pathogenic"
)

probe <- matrix(c(0, 1, NA_real_), ncol = 1L, dimnames = list(NULL, "Type"))
expect_equal(
  XCNV:::.predict_portable_tree(model, probe),
  stats::plogis(c(-0.59, 0.5, -0.59)), tolerance = 1e-12
)
expect_error(load_xcnv_model(system.file("extdata", "fixture", "genome.txt", package = "XCNV")))

out <- tempfile(fileext = ".csv")
on.exit(unlink(out), add = TRUE)
run_xcnv(input, fixture, output = out)
expect_true(file.exists(out))
written <- read.csv(out, check.names = FALSE)
expect_equal(written$MVP_score, predictions$MVP_score, tolerance = 1e-12)
header <- strsplit(readLines(out, n = 1L), ",", fixed = TRUE)[[1L]]
type_code_column <- which(names(predictions) == "Type.1")
expect_equal(header[[type_code_column]], "Type")
