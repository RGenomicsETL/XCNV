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
  assembly = "GRCh37",
  provider_id = c("provider-a", "provider-a", "provider-a", "provider-b"),
  evidence_id = c("gene:1", "case:1", "case:2", "gene-count:0")
)
groups <- score_cnv_2019(evidence, detail = "groups")
expect_equal(groups$applied_score[groups$cnv_id == "loss-1" & groups$cap_group == "4ABC"], 0.90)
summary <- score_cnv_2019(evidence)
expect_equal(summary$score[summary$cnv_id == "loss-1"], 1.90)
expect_equal(summary$classification[summary$cnv_id == "loss-1"], "pathogenic")
expect_equal(summary$classification[summary$cnv_id == "gain-1"], "uncertain_significance")
expect_true(all(summary$provider_id == "xcnv-acmg-clingen-2019"))
expect_true(!any(summary$provider_id == "conflict"))
expect_equal(summary$resource_release[summary$cnv_id == "gain-1"], "50")
expect_true(is.na(summary$resource_release[summary$cnv_id == "loss-1"]))
expect_equal(summary$evidence_provider_id[summary$cnv_id == "loss-1"], "provider-a")
expect_equal(summary$evidence_provider_id[summary$cnv_id == "gain-1"], "provider-b")
expect_true(is.na(summary$evidence_source_release[summary$cnv_id == "loss-1"]))
expect_equal(summary$evidence_release_count[summary$cnv_id == "loss-1"], 2L)

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
expect_true(all(derived$provider_id[derived$source_id == "GENCODE"] == "GENCODE"))
expect_true(all(derived$provider_id[derived$source_id == "ClinGen"] == "ClinGen"))
expect_equal(score_cnv_2019(derived)$score, c(1, 1))
unspecified_summary <- score_cnv_2019(derived)
expect_true(all(unspecified_summary$status == "unspecified"))
expect_true(all(is.na(unspecified_summary$classification)))
policy_summary <- score_cnv_2019(derived)
expect_true(all(policy_summary$provider_id == "xcnv-acmg-clingen-2019"))
expect_true(all(is.na(policy_summary$evidence_provider_id)))
expect_equal(policy_summary$evidence_provider_count, c(2L, 2L))
expect_true(all(policy_summary$policy_id == "acmg-clingen-2019-constitutional-cnv"))
expect_true(all(is.na(policy_summary$model_id)))
expect_true(!"MVP_score" %in% names(policy_summary))

cnvs38 <- transform(cnvs, assembly = "GRCh38")
genes38 <- transform(genes, assembly = "GRCh38")
dosage38 <- transform(dosage, assembly = "GRCh38")
evidence38 <- derive_cnv_2019_evidence(cnvs38, genes38, dosage38)
expect_true(all(evidence38$assembly == "GRCh38"))
classification38 <- score_cnv_2019(evidence38)
expect_equal(classification38$assembly, c("GRCh38", "GRCh38"))
expect_equal(classification38$status, c("supported", "supported"))

mixed_assembly <- genes38
mixed_assembly$assembly <- "GRCh37"
expect_error(derive_cnv_2019_evidence(cnvs38, mixed_assembly, dosage38))

unsupported <- transform(cnvs, assembly = "T2T-CHM13")
unsupported_genes <- transform(genes, assembly = "T2T-CHM13")
unsupported_dosage <- transform(dosage, assembly = "T2T-CHM13")
unsupported_evidence <- derive_cnv_2019_evidence(
  unsupported, unsupported_genes, unsupported_dosage
)
expect_true(all(unsupported_evidence$status == "unsupported"))
unsupported_result <- score_cnv_2019(unsupported_evidence)
expect_true(all(unsupported_result$status == "unsupported"))
expect_true(all(is.na(unsupported_result$classification)))

conflicting_status <- evidence[1:2, , drop = FALSE]
conflicting_status$status <- c("supported", "conflict")
conflict_result <- score_cnv_2019(conflicting_status)
expect_true(all(conflict_result$status == "conflict"))
expect_true(all(is.na(conflict_result$classification)))

per_cnv_conflict <- rbind(
  transform(evidence[1, , drop = FALSE], cnv_id = "bad-1", policy_id = "policy-a"),
  transform(evidence[2, , drop = FALSE], cnv_id = "bad-1", assembly = "GRCh38", policy_id = "policy-b"),
  transform(evidence[4, , drop = FALSE], cnv_id = "good-1", policy_id = "acmg-clingen-2019-constitutional-cnv")
)
per_cnv_summary <- score_cnv_2019(per_cnv_conflict)
expect_equal(per_cnv_summary$status[per_cnv_summary$cnv_id == "bad-1"], "conflict")
expect_true(is.na(per_cnv_summary$classification[per_cnv_summary$cnv_id == "bad-1"]))
expect_equal(per_cnv_summary$status[per_cnv_summary$cnv_id == "good-1"], "supported")
expect_true(!is.na(per_cnv_summary$classification[per_cnv_summary$cnv_id == "good-1"]))
bad$criterion <- "not-a-criterion"
expect_error(score_cnv_2019(bad))
bad <- evidence[1, setdiff(names(evidence), "source_release"), drop = FALSE]
expect_error(score_cnv_2019(bad))
