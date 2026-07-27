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

# A foreign policy never reaches a summary/classification, including when it is
# mixed with otherwise valid native-policy evidence.
foreign_policy <- evidence[1, , drop = FALSE]
foreign_policy$policy_id <- "foreign-policy"
expect_error(score_cnv_2019(foreign_policy))
mixed_policy <- evidence[1:2, , drop = FALSE]
mixed_policy$policy_id <- c("acmg-clingen-2019-constitutional-cnv", "foreign-policy")
expect_error(score_cnv_2019(mixed_policy))

# Scores are parsed through factor/character labels rather than factor codes.
# In particular, a factor label "0" cannot become the factor code 1 for 2A.
factor_zero_score <- evidence[1, , drop = FALSE]
factor_zero_score$score <- factor("0")
expect_error(score_cnv_2019(factor_zero_score))
factor_one_score <- evidence[1, , drop = FALSE]
factor_one_score$score <- factor("1")
expect_equal(score_cnv_2019(factor_one_score)$score, 1)
nonnumeric_score <- evidence[1, , drop = FALSE]
nonnumeric_score$score <- "nonsense"
expect_error(score_cnv_2019(nonnumeric_score))

# Every required evidence identity/provenance value is present on every row;
# no summary may obtain it from a different row.
for (column in c(
  "cnv_id", "criterion", "source_id", "source_release", "evidence_id", "provider_id"
)) {
  for (value in list(NA_character_, "")) {
    malformed_evidence <- evidence[1, , drop = FALSE]
    malformed_evidence[[column]] <- value
    expect_error(score_cnv_2019(malformed_evidence))
  }
}

# model_id is a supported evidence-model provenance field, not a discarded
# input: one recognized identity reaches both audit surfaces.
modelled_evidence <- evidence[1, , drop = FALSE]
modelled_evidence$model_id <- "evidence-model-v1"
expect_equal(score_cnv_2019(modelled_evidence, detail = "groups")$model_id, "evidence-model-v1")
expect_equal(score_cnv_2019(modelled_evidence)$model_id, "evidence-model-v1")
empty_model_id <- modelled_evidence
empty_model_id$model_id <- ""
expect_error(score_cnv_2019(empty_model_id))
conflicting_model_id <- evidence[1:2, , drop = FALSE]
conflicting_model_id$model_id <- c("evidence-model-a", "evidence-model-b")
expect_error(score_cnv_2019(conflicting_model_id))

# Model attribution is based only on selected evidence that contributes to an
# output grain. An unselected model cannot be broadcast to the selected 2A
# group or the CNV/type summary.
unselected_model <- evidence[1:2, , drop = FALSE]
unselected_model$selected <- c(TRUE, FALSE)
unselected_model$model_id <- c(NA_character_, "m-unselected")
unselected_model_groups <- score_cnv_2019(unselected_model, detail = "groups")
unselected_model_summary <- score_cnv_2019(unselected_model)
expect_true(is.na(unselected_model_groups$model_id[
  unselected_model_groups$cap_group == "2A"
]))
expect_true(is.na(unselected_model_summary$model_id))

# Missing model provenance in one selected contributor prevents an attribution
# for both the group and its summary, rather than ignoring the missing row.
mixed_missing_models <- evidence[2:3, , drop = FALSE]
mixed_missing_models$model_id <- c("m-selected", NA_character_)
mixed_missing_groups <- score_cnv_2019(mixed_missing_models, detail = "groups")
mixed_missing_summary <- score_cnv_2019(mixed_missing_models)
expect_true(is.na(mixed_missing_groups$model_id[
  mixed_missing_groups$cap_group == "4ABC"
]))
expect_true(is.na(mixed_missing_summary$model_id))

# Evidence identifiers, provenance, and supplied contract labels are labels,
# not values to stringify. TRUE and FALSE must not create fake identities.
explicit_contract_evidence <- evidence[1, , drop = FALSE]
explicit_contract_evidence$model_id <- "evidence-model-v1"
explicit_contract_evidence$policy_id <- "acmg-clingen-2019-constitutional-cnv"
explicit_contract_evidence$status <- "supported"
for (column in c(
  "cnv_id", "cnv_type", "criterion", "source_id", "source_release",
  "evidence_id", "provider_id", "assembly", "model_id", "policy_id", "status"
)) {
  for (value in list(TRUE, FALSE)) {
    logical_label <- explicit_contract_evidence
    logical_label[[column]] <- value
    expect_error(score_cnv_2019(logical_label))
  }
}
for (value in list(1, I(list("provider")), I(matrix("provider", nrow = 1L)))) {
  nonlabel_provider <- explicit_contract_evidence
  nonlabel_provider$provider_id <- value
  expect_error(score_cnv_2019(nonlabel_provider))
}

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

# Dosage labels are parsed safely: NA/nonnumeric labels cannot become 2A, while
# a factor label "3" is the established score rather than its factor code.
no_2a <- function(x) {
  derived_x <- derive_cnv_2019_evidence(cnvs[1, , drop = FALSE], genes, x)
  !"2A" %in% derived_x$criterion
}
dosage_na_score <- dosage
dosage_na_score$hi_score[1] <- NA_real_
expect_true(no_2a(dosage_na_score))
dosage_nonnumeric_score <- dosage
dosage_nonnumeric_score$hi_score <- as.character(dosage_nonnumeric_score$hi_score)
dosage_nonnumeric_score$hi_score[1] <- "nonsense"
expect_true(no_2a(dosage_nonnumeric_score))
dosage_factor_score <- dosage
dosage_factor_score$hi_score <- factor(c("3", "0"))
expect_true("2A" %in% derive_cnv_2019_evidence(
  cnvs[1, , drop = FALSE], genes, dosage_factor_score, established_score = factor("3")
)$criterion)
expect_error(derive_cnv_2019_evidence(cnvs, genes, dosage, established_score = NA_real_))

# Required automatic-provider identifiers/provenance are row-level: missing
# values cannot borrow the relation's identity from another row.
missing_cnv_id <- cnvs
missing_cnv_id$cnv_id[1] <- NA_character_
expect_error(derive_cnv_2019_evidence(missing_cnv_id, genes, dosage))
missing_gene_id <- genes
missing_gene_id$gene_id[1] <- ""
expect_error(derive_cnv_2019_evidence(cnvs, missing_gene_id, dosage))
missing_record_id <- dosage
missing_record_id$record_id[1] <- NA_character_
expect_error(derive_cnv_2019_evidence(cnvs, genes, missing_record_id))
missing_gene_source <- genes
missing_gene_source$source_id[2] <- NA_character_
expect_error(derive_cnv_2019_evidence(cnvs, missing_gene_source, dosage))
missing_gene_release <- genes
missing_gene_release$source_release[2] <- ""
expect_error(derive_cnv_2019_evidence(cnvs, missing_gene_release, dosage))
missing_dosage_source <- dosage
missing_dosage_source$source_id[2] <- NA_character_
expect_error(derive_cnv_2019_evidence(cnvs, genes, missing_dosage_source))
missing_dosage_release <- dosage
missing_dosage_release$source_release[2] <- ""
expect_error(derive_cnv_2019_evidence(cnvs, genes, missing_dosage_release))
mixed_gene_provider <- genes
mixed_gene_provider$source_id[2] <- "OTHER-GENCODE"
expect_error(derive_cnv_2019_evidence(cnvs, mixed_gene_provider, dosage))
mixed_dosage_release <- dosage
mixed_dosage_release$source_release[2] <- "other-release"
expect_error(derive_cnv_2019_evidence(cnvs, genes, mixed_dosage_release))

# CNV, gene, and dosage identifiers/provenance are labels, so logical values
# (including record_id) cannot become fabricated TRUE/FALSE identities.
for (value in list(TRUE, FALSE)) {
  for (column in c("cnv_id", "cnv_type")) {
    logical_cnv_label <- cnvs
    logical_cnv_label[[column]] <- value
    expect_error(derive_cnv_2019_evidence(logical_cnv_label, genes, dosage))
  }
  for (column in c("gene_id", "source_id", "source_release")) {
    logical_gene_label <- genes
    logical_gene_label[[column]] <- value
    expect_error(derive_cnv_2019_evidence(cnvs, logical_gene_label, dosage))
  }
  for (column in c("record_id", "source_id", "source_release")) {
    logical_dosage_label <- dosage
    logical_dosage_label[[column]] <- value
    expect_error(derive_cnv_2019_evidence(cnvs, genes, logical_dosage_label))
  }
}

# Functional filtering has an explicit biotype vocabulary. Incomplete,
# padded, or placeholder values fail before filtering, so they cannot fabricate
# 1B/3A absence rows.
for (value in list(
  NA_character_, "", " protein_coding ", TRUE, I(list("protein_coding"))
)) {
  incomplete_biotype <- genes
  incomplete_biotype$biotype <- value
  expect_error(derive_cnv_2019_evidence(cnvs[1, , drop = FALSE], incomplete_biotype, dosage))
}
for (value in list(
  NA_character_, "", " protein_coding ", TRUE, I(list("protein_coding")), character()
)) {
  expect_error(derive_cnv_2019_evidence(
    cnvs[1, , drop = FALSE], genes, dosage, functional_biotypes = value
  ))
}
for (placeholder in c("unknown", ".", "NA", "N-A", "N/A", "null", "unspecified", "not_provided", "UNKNOWN")) {
  placeholder_biotype <- genes
  placeholder_biotype$biotype <- placeholder
  expect_error(derive_cnv_2019_evidence(cnvs[1, , drop = FALSE], placeholder_biotype, dosage))
  expect_error(derive_cnv_2019_evidence(
    cnvs[1, , drop = FALSE], genes, dosage, functional_biotypes = placeholder
  ))
}
expect_error(derive_cnv_2019_evidence(
  cnvs[1, , drop = FALSE], genes, dosage,
  functional_biotypes = c("protein_coding", "protein_coding")
))
nonfunctional_biotype <- genes
nonfunctional_biotype$biotype <- "lncRNA"
nonfunctional_evidence <- derive_cnv_2019_evidence(
  cnvs[1, , drop = FALSE], nonfunctional_biotype, dosage
)
expect_identical(sort(nonfunctional_evidence$criterion), c("1B", "2A", "3A"))
factor_biotype <- genes
factor_biotype$biotype <- factor(factor_biotype$biotype)
expect_true("1A" %in% derive_cnv_2019_evidence(
  cnvs[1, , drop = FALSE], factor_biotype, dosage,
  functional_biotypes = factor("protein_coding")
)$criterion)
factor_record_id <- dosage
factor_record_id$record_id <- factor(factor_record_id$record_id)
expect_true("dosage:HI1" %in% derive_cnv_2019_evidence(
  cnvs[1, , drop = FALSE], genes, factor_record_id
)$evidence_id)

cnvs38 <- transform(cnvs, assembly = "GRCh38")
genes38 <- transform(genes, assembly = "GRCh38")
dosage38 <- transform(dosage, assembly = "GRCh38")
evidence38 <- derive_cnv_2019_evidence(cnvs38, genes38, dosage38)
expect_true(all(evidence38$assembly == "GRCh38"))
classification38 <- score_cnv_2019(evidence38)
expect_equal(classification38$assembly, c("GRCh38", "GRCh38"))
expect_equal(classification38$status, c("supported", "supported"))

# Every native BED relation rejects zero-width, reversed, fractional, missing,
# infinite, negative, and nameless intervals before the overlap query. In
# particular, a zero-width GRCh38 CNV is never scored.
expect_invalid_geometry <- function(relation, invalid) {
  inputs <- list(cnvs = cnvs38, genes = genes38, dosage = dosage38)
  inputs[[relation]] <- invalid
  expect_error(do.call(derive_cnv_2019_evidence, inputs))
}
invalid_interval <- function(x, kind) {
  x <- x[1, , drop = FALSE]
  switch(
    kind,
    zero_width = x$end <- x$start,
    reversed = x$start <- x$end + 1,
    fractional = x$start <- x$start + 0.5,
    missing = x$start <- NA_real_,
    infinite = x$end <- Inf,
    negative = x$start <- -1,
    empty_contig = x$chr <- ""
  )
  x
}
for (relation in c("cnvs", "genes", "dosage")) {
  source <- list(cnvs = cnvs38, genes = genes38, dosage = dosage38)[[relation]]
  for (kind in c(
    "zero_width", "reversed", "fractional", "missing", "infinite",
    "negative", "empty_contig"
  )) {
    expect_invalid_geometry(relation, invalid_interval(source, kind))
  }
}

# Contigs are explicit character/factor labels on every BED relation. They
# cannot be stringified from another type or padded into a different key.
contig_counterexamples <- list(
  TRUE, 1, I(list("1")), I(matrix("1", nrow = 1L)),
  NA_character_, "", " 1", "1 "
)
for (relation in c("cnvs", "genes", "dosage")) {
  source <- list(cnvs = cnvs38, genes = genes38, dosage = dosage38)[[relation]]
  for (value in contig_counterexamples) {
    invalid_contig <- source[1, , drop = FALSE]
    invalid_contig$chr <- value
    expect_invalid_geometry(relation, invalid_contig)
  }
}

# BED strings are checked lexically before double conversion. The accepted
# human-genomic range is [0, 2^53), so these would-be rounding counterexamples
# never reach an overlap query.
safe_cnvs <- cnvs38[1, , drop = FALSE]
safe_genes <- genes38[1, , drop = FALSE]
safe_dosage <- dosage38[1, , drop = FALSE]
safe_cnvs$start <- "9007199254740990"
safe_cnvs$end <- "9007199254740991"
safe_genes$start <- "9007199254740990"
safe_genes$end <- "9007199254740991"
safe_dosage$start <- "9007199254740990"
safe_dosage$end <- "9007199254740991"
expect_true("2A" %in% derive_cnv_2019_evidence(safe_cnvs, safe_genes, safe_dosage)$criterion)
for (outside_safe_range in c("9007199254740992", "9007199254740993", "9007199254740995")) {
  precision_counterexample <- safe_cnvs
  precision_counterexample$end <- outside_safe_range
  expect_error(derive_cnv_2019_evidence(precision_counterexample, safe_genes, safe_dosage))
}

# A zero-row CNV relation (including one declaring assembly) is neutral. Empty
# gene/dosage relations are not treated as absence evidence for nonempty CNVs.
empty_cnvs <- cnvs38[FALSE, , drop = FALSE]
empty_evidence <- derive_cnv_2019_evidence(empty_cnvs, genes38, dosage38)
expect_equal(nrow(empty_evidence), 0L)
expect_identical(names(empty_evidence), c(
  "cnv_id", "cnv_type", "criterion", "score", "source_id", "source_release",
  "evidence_id", "selected", "provider_id", "assembly", "model_id", "policy_id", "status"
))
expect_equal(nrow(derive_cnv_2019_evidence(
  empty_cnvs, genes38[FALSE, , drop = FALSE], dosage38[FALSE, , drop = FALSE]
)), 0L)
expect_error(derive_cnv_2019_evidence(cnvs38[1, , drop = FALSE], genes38[FALSE, , drop = FALSE], dosage38))
expect_error(derive_cnv_2019_evidence(cnvs38[1, , drop = FALSE], genes38, dosage38[FALSE, , drop = FALSE]))
conflicting_empty_cnv_resources <- dosage38
conflicting_empty_cnv_resources$assembly <- "GRCh37"
expect_error(derive_cnv_2019_evidence(empty_cnvs, genes38, conflicting_empty_cnv_resources))

# Native interval matching uses exact caller/provider contig names. A valid
# no-overlap contig is admissible and does not borrow resources from "1".
no_overlap_cnvs <- cnvs38[1, , drop = FALSE]
no_overlap_cnvs$cnv_id <- "no-overlap"
no_overlap_cnvs$chr <- "chr1"
no_overlap_evidence <- derive_cnv_2019_evidence(no_overlap_cnvs, genes38, dosage38)
expect_identical(sort(no_overlap_evidence$criterion), c("1B", "3A"))

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

bad$criterion <- "not-a-criterion"
expect_error(score_cnv_2019(bad))
bad <- evidence[1, setdiff(names(evidence), "source_release"), drop = FALSE]
expect_error(score_cnv_2019(bad))
