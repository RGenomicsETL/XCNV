library(XCNV)

expect_equal(
  as.character(classify_cnv_score(c(-1, -0.95, 0, 0.95, 1))),
  c("benign", "likely_benign", "uncertain_significance", "likely_pathogenic", "pathogenic")
)
expect_true(all(c("loss", "gain") %in% unique(cnv_2019_criteria()$cnv_type)))

evidence <- data.frame(
  cnv_id = c("loss-1", "loss-1", "loss-1", "gain-1"),
  cnv_type = c("loss", "loss", "loss", "gain"),
  criterion = c("2A", "4A", "4A", "1B"),
  score = c(1, 0.45, 0.45, -0.60),
  source_id = c("ClinGen", "PMID", "PMID", "GENCODE"),
  source_release = c("2026-07-21", "31690835", "31690835", "50"),
  evidence_id = c("gene:1", "case:1", "case:2", "gene-count:0")
)
groups <- score_cnv_2019(evidence, detail = "groups")
expect_equal(groups$applied_score[groups$cnv_id == "loss-1" & groups$cap_group == "4ABC"], 0.90)
summary <- score_cnv_2019(evidence)
expect_equal(summary$score[summary$cnv_id == "loss-1"], 1.90)
expect_equal(summary$classification[summary$cnv_id == "loss-1"], "pathogenic")
expect_equal(summary$classification[summary$cnv_id == "gain-1"], "uncertain_significance")

bad <- evidence[1, , drop = FALSE]
bad$score <- 0.5
expect_error(score_cnv_2019(bad))

cnvs <- data.frame(
  cnv_id = c("del", "dup"), chr = c("1", "1"), start = c(90L, 290L),
  end = c(220L, 420L), cnv_type = c("loss", "gain")
)
genes <- data.frame(
  chr = c("1", "1"), start = c(100L, 300L), end = c(200L, 400L),
  gene_id = c("GENE1", "GENE2"), biotype = "protein_coding",
  source_id = "GENCODE", source_release = "50"
)
dosage <- data.frame(
  chr = c("1", "1"), start = c(100L, 300L), end = c(200L, 400L),
  record_id = c("HI1", "TS1"), hi_score = c(3, 0), ts_score = c(0, 3),
  source_id = "ClinGen", source_release = "2026-07-21"
)
derived <- derive_cnv_2019_evidence(cnvs, genes, dosage)
expect_true(all(c("1A", "2A", "3A") %in% derived$criterion[derived$cnv_id == "del"]))
expect_true(all(c("1A", "2A", "3A") %in% derived$criterion[derived$cnv_id == "dup"]))
expect_equal(score_cnv_2019(derived)$score, c(1, 1))
bad$criterion <- "not-a-criterion"
expect_error(score_cnv_2019(bad))
bad <- evidence[1, setdiff(names(evidence), "source_release"), drop = FALSE]
expect_error(score_cnv_2019(bad))
