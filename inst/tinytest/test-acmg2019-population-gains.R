library(XCNV)

# Partial-gene gains require transcript and breakpoint evidence.
gain_cnvs <- data.frame(
  cnv_id = c(
    "internal", "left_partial", "right_partial", "gene_start_partial",
    "gene_end_partial", "whole_gene", "flanking", "left_adjacent", "right_adjacent"
  ),
  chr = "chr1",
  start = c(140L, 90L, 140L, 100L, 140L, 100L, 90L, 80L, 200L),
  end = c(160L, 160L, 210L, 160L, 200L, 200L, 210L, 100L, 220L),
  cnv_type = "gain", assembly = "GRCh38"
)
gain_genes <- data.frame(
  chr = "chr1", start = c(100L, 50L), end = c(200L, 250L),
  gene_id = c("HI_GENE", "RNA_GENE"), biotype = c("protein_coding", "lncRNA"),
  source_id = "GENCODE", source_release = "50", assembly = "GRCh38"
)
gain_dosage <- data.frame(
  chr = "chr1", start = 100L, end = 200L, record_id = "HI_GENE",
  hi_score = 3, ts_score = 0, source_id = "ClinGen",
  source_release = "2026-11-08", assembly = "GRCh38"
)
gain_population <- data.frame(
  chr = "chr1", start = 50L, end = 250L, record_id = "whole_gene_dup",
  cnv_type = "gain", frequency = 0.02, source_id = "population_source",
  source_release = "1", assembly = "GRCh38"
)
gain_evidence <- derive_cnv_2019_evidence(
  gain_cnvs, gain_genes, gain_dosage, population = gain_population
)
expect_equal(
  gain_evidence$cnv_id[gain_evidence$criterion == "4O"],
  c("whole_gene", "flanking", "left_adjacent", "right_adjacent")
)
gain_summary <- score_cnv_2019(gain_evidence)
expect_equal(
  as.character(gain_summary$classification[
    match(gain_cnvs$cnv_id[1:5], gain_summary$cnv_id)
  ]),
  rep("uncertain_significance", 5L)
)
expect_equal(
  as.character(gain_summary$classification[
    match(gain_cnvs$cnv_id[6:9], gain_summary$cnv_id)
  ]),
  rep("benign", 4L)
)

# An interval match alone does not establish a partial gain's transcript effect.
internal_population <- gain_population
internal_population$start <- 140L
internal_population$end <- 160L
internal_evidence <- derive_cnv_2019_evidence(
  gain_cnvs[1, ], gain_genes, gain_dosage, population = internal_population
)
expect_false("4O" %in% internal_evidence$criterion)

# Every overlapping functional gene must be completely contained in a gain.
extra_gene <- gain_genes[1, ]
extra_gene$gene_id <- "SECOND_GENE"
extra_gene$start <- 180L
extra_gene$end <- 230L
multiple_gene_dosage <- gain_dosage
multiple_gene_dosage$ts_score <- 3
multiple_gene_evidence <- derive_cnv_2019_evidence(
  gain_cnvs[gain_cnvs$cnv_id == "whole_gene", ],
  rbind(gain_genes, extra_gene), multiple_gene_dosage, population = gain_population
)
expect_false("4O" %in% multiple_gene_evidence$criterion)
expect_equal(multiple_gene_evidence$score[multiple_gene_evidence$criterion == "2A"], 1)

# A containing population deletion can support a smaller partial-gene loss.
partial_loss <- gain_cnvs[1, ]
partial_loss$cnv_type <- "loss"
loss_dosage <- gain_dosage
loss_dosage$hi_score <- 0
loss_population <- gain_population
loss_population$cnv_type <- "loss"
loss_evidence <- derive_cnv_2019_evidence(
  partial_loss, gain_genes, loss_dosage, population = loss_population
)
expect_equal(loss_evidence$score[loss_evidence$criterion == "4O" & loss_evidence$selected], -1)

# Functional biotypes and the established score accept positional arguments.
positional_cnv <- gain_cnvs[gain_cnvs$cnv_id == "whole_gene", ]
positional_cnv$cnv_type <- "loss"
expect_equal(
  derive_cnv_2019_evidence(positional_cnv, gain_genes, gain_dosage, "protein_coding"),
  derive_cnv_2019_evidence(positional_cnv, gain_genes, gain_dosage)
)
expect_equal(
  derive_cnv_2019_evidence(
    positional_cnv, gain_genes, gain_dosage, c("protein_coding", "lncRNA"), 2
  ),
  derive_cnv_2019_evidence(
    positional_cnv, gain_genes, gain_dosage,
    functional_biotypes = c("protein_coding", "lncRNA"), established_score = 2
  )
)
expect_equal(
  derive_cnv_2019_evidence(
    gain_cnvs, gain_genes, gain_dosage, "protein_coding", 3,
    population = gain_population
  ),
  gain_evidence
)
