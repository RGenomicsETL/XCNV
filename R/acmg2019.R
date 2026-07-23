.cnv_2019_rows <- function(cnv_type) {
  common <- data.frame(
    criterion = c(
      "1A", "1B", "3A", "3B", "3C",
      "4A", "4B", "4C", "4D", "4E", "4F", "4G", "4H",
      "4I", "4J", "4K", "4L", "4M", "4N", "4O",
      "5A", "5B", "5C", "5D", "5E", "5F", "5G", "5H"
    ),
    default_score = c(
      0, -0.60, 0, 0.45, 0.90,
      0.45, 0.30, 0.15, 0, 0.10, 0.15, 0.30, 0.45,
      -0.45, -0.30, -0.15, 0.45, 0.30, -0.90, -1,
      NA, -0.30, -0.15, NA, NA, 0, 0.10, 0.15
    ),
    min_score = c(
      0, -0.60, 0, 0.45, 0.90,
      0.15, 0, 0, -0.30, 0, 0.15, 0.30, 0.45,
      -0.45, -0.30, -0.15, 0, 0, -0.90, -1,
      0, -0.45, -0.30, 0, -0.45, 0, 0, 0
    ),
    max_score = c(
      0, -0.60, 0, 0.45, 0.90,
      0.45, 0.45, 0.30, 0, 0.15, 0.15, 0.30, 0.45,
      0, 0, 0, 0.45, 0.30, 0, 0,
      0.45, 0, 0, 0.45, 0, 0, 0.15, 0.30
    ),
    cap_group = c(
      "1A", "1B", "3A", "3B", "3C",
      rep("4ABC", 3), "4D", "4E", rep("4FGH", 3),
      "4I", "4J", "4K", "4L", "4M", "4N", "4O",
      "5A", "5B", "5C", "5D", "5E", "5F", "5G", "5H"
    ),
    group_min = c(
      0, -0.60, 0, 0.45, 0.90,
      rep(0, 3), -0.30, 0, rep(0, 3),
      -0.90, -0.90, -0.30, 0, 0, -0.90, -1,
      0, -0.45, -0.30, 0, -0.45, 0, 0, 0
    ),
    group_max = c(
      0, -0.60, 0, 0.45, 0.90,
      rep(0.90, 3), 0, 0.30, rep(0.45, 3),
      0, 0, 0, 0.45, 0.45, 0, 0,
      0.45, 0, 0, 0.45, 0, 0, 0.15, 0.30
    ),
    repeatable = c(
      rep(FALSE, 5), rep(TRUE, 15), rep(FALSE, 8)
    ),
    stringsAsFactors = FALSE
  )
  if (identical(cnv_type, "loss")) {
    specific <- data.frame(
      criterion = c(
        "2A", "2B", "2C-1", "2C-2", "2D-1", "2D-2", "2D-3", "2D-4",
        "2E-PVS1", "2E-PVS1_Strong", "2E-PVS1_Moderate",
        "2E-PVS1_Supporting", "2E-NA", "2F", "2G", "2H"
      ),
      default_score = c(1, 0, .90, 0, 0, .90, .30, .90, .90, .45, .30, .15, 0, -1, 0, .15),
      min_score = c(1, 0, .45, 0, 0, .45, 0, .45, .45, .30, .15, 0, 0, -1, 0, .15),
      max_score = c(1, 0, 1, .45, 0, .90, .45, 1, .90, .90, .45, .30, 0, -1, 0, .15),
      stringsAsFactors = FALSE
    )
  } else {
    specific <- data.frame(
      criterion = c(
        "2A", "2B", "2C", "2D", "2E", "2F", "2G", "2H",
        "2I-PVS1", "2I-PVS1_Strong", "2I-NA", "2J", "2K", "2L"
      ),
      default_score = c(1, 0, -1, -1, 0, -1, 0, 0, .90, .45, 0, 0, .45, 0),
      min_score = c(1, 0, -1, -1, 0, -1, 0, 0, .45, .30, 0, 0, .45, 0),
      max_score = c(1, 0, -1, -1, 0, 0, 0, 0, .90, .90, 0, 0, .45, 0),
      stringsAsFactors = FALSE
    )
  }
  specific$cap_group <- specific$criterion
  specific$group_min <- specific$min_score
  specific$group_max <- specific$max_score
  specific$repeatable <- FALSE
  out <- rbind(common, specific)
  out$cnv_type <- cnv_type
  out[, c(
    "cnv_type", "criterion", "default_score", "min_score", "max_score",
    "cap_group", "group_min", "group_max", "repeatable"
  )]
}

.cnv_2019_criteria <- rbind(.cnv_2019_rows("loss"), .cnv_2019_rows("gain"))

#' ACMG/ClinGen 2019 constitutional CNV criterion contract
#'
#' Return the machine-readable criterion identifiers and numeric ranges used by
#' the independently authored scoring lane. This is a scoring contract, not an
#' automatic clinical evidence adjudicator. Dynamic section 5 criteria have an
#' `NA` default and require an explicit score selected by the curator.
#'
#' @param cnv_type One or both of `"loss"` and `"gain"`.
#' @return A data frame with allowed per-evidence ranges and aggregate caps.
#' @export
cnv_2019_criteria <- function(cnv_type = c("loss", "gain")) {
  cnv_type <- match.arg(cnv_type, c("loss", "gain"), several.ok = TRUE)
  .cnv_2019_criteria[.cnv_2019_criteria$cnv_type %in% cnv_type, , drop = FALSE]
}

#' Convert an ACMG/ClinGen 2019 CNV score to a five-tier classification
#'
#' @param score Numeric score vector.
#' @return An ordered factor with the five constitutional CNV classes.
#' @export
classify_cnv_score <- function(score) {
  if (!is.numeric(score) || any(!is.finite(score))) {
    stop("'score' must contain finite numeric values", call. = FALSE)
  }
  label <- ifelse(
    score >= 0.99, "pathogenic",
    ifelse(
      score >= 0.90, "likely_pathogenic",
      ifelse(score <= -0.99, "benign", ifelse(score <= -0.90, "likely_benign", "uncertain_significance"))
    )
  )
  factor(
    label,
    levels = c("benign", "likely_benign", "uncertain_significance", "likely_pathogenic", "pathogenic"),
    ordered = TRUE
  )
}

.normalize_cnv_evidence <- function(evidence) {
  required <- c(
    "cnv_id", "cnv_type", "criterion", "score", "source_id",
    "source_release", "evidence_id"
  )
  if (!is.data.frame(evidence)) stop("'evidence' must be a data frame", call. = FALSE)
  missing <- setdiff(required, names(evidence))
  if (length(missing)) stop("evidence is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  out <- as.data.frame(evidence, stringsAsFactors = FALSE)
  for (column in setdiff(required, "score")) out[[column]] <- as.character(out[[column]])
  out$score <- suppressWarnings(as.numeric(out$score))
  if (!"selected" %in% names(out)) out$selected <- TRUE
  if (!is.logical(out$selected) || anyNA(out$selected)) stop("evidence selected must be TRUE or FALSE", call. = FALSE)
  if (any(!nzchar(out$cnv_id)) || any(!out$cnv_type %in% c("loss", "gain")) ||
      any(!nzchar(out$criterion)) || any(!is.finite(out$score)) ||
      any(!nzchar(out$source_id)) || any(!nzchar(out$source_release)) ||
      any(!nzchar(out$evidence_id))) {
    stop("evidence identifiers, provenance, CNV type, and score must be explicit", call. = FALSE)
  }
  keys <- paste(out$cnv_type, out$criterion, sep = ":")
  contract_keys <- paste(.cnv_2019_criteria$cnv_type, .cnv_2019_criteria$criterion, sep = ":")
  contract_row <- match(keys, contract_keys)
  if (anyNA(contract_row)) {
    stop("unknown criterion for CNV type: ", paste(unique(keys[is.na(contract_row)]), collapse = ", "), call. = FALSE)
  }
  contract <- .cnv_2019_criteria[contract_row, , drop = FALSE]
  tolerance <- 1e-12
  invalid_score <- out$score < contract$min_score - tolerance | out$score > contract$max_score + tolerance
  if (any(invalid_score)) {
    stop("criterion score outside its published range: ", paste(unique(keys[invalid_score]), collapse = ", "), call. = FALSE)
  }
  out$cap_group <- contract$cap_group
  out$group_min <- contract$group_min
  out$group_max <- contract$group_max
  out$repeatable <- contract$repeatable
  nonrepeatable <- out$selected & !out$repeatable
  duplicate_key <- paste(out$cnv_id, keys, sep = ":")
  if (anyDuplicated(duplicate_key[nonrepeatable])) {
    stop("non-repeatable criteria may be selected once per CNV", call. = FALSE)
  }
  out
}

.score_cnv_2019_groups <- function(evidence) {
  evidence <- .normalize_cnv_evidence(evidence)
  selected <- evidence[evidence$selected, , drop = FALSE]
  if (!nrow(selected)) {
    return(data.frame(
      cnv_id = character(), cnv_type = character(), cap_group = character(),
      raw_score = numeric(), applied_score = numeric(), source_count = integer()
    ))
  }
  key <- interaction(selected$cnv_id, selected$cnv_type, selected$cap_group, drop = TRUE, lex.order = TRUE)
  split_rows <- split(seq_len(nrow(selected)), key)
  rows <- lapply(split_rows, function(i) {
    raw <- sum(selected$score[i])
    lower <- unique(selected$group_min[i])
    upper <- unique(selected$group_max[i])
    if (length(lower) != 1L || length(upper) != 1L) stop("inconsistent criterion cap group", call. = FALSE)
    data.frame(
      cnv_id = selected$cnv_id[i[[1L]]], cnv_type = selected$cnv_type[i[[1L]]],
      cap_group = selected$cap_group[i[[1L]]], raw_score = raw,
      applied_score = min(upper, max(lower, raw)), source_count = length(i),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$cnv_id, out$cnv_type, out$cap_group), , drop = FALSE]
}

#' Score transparent ACMG/ClinGen 2019 constitutional CNV evidence
#'
#' The input is an auditable evidence relation. Every row names its CNV,
#' deletion/gain context, criterion, numeric score, source, source release, and
#' evidence record. Repeated case evidence is capped by the published criterion
#' group maximum/minimum. Evidence selection remains explicit and is never
#' inferred from free text or phenotype similarity by this function.
#'
#' @param evidence Data frame with `cnv_id`, `cnv_type`, `criterion`, `score`,
#'   `source_id`, `source_release`, and `evidence_id`; optional logical
#'   `selected` defaults to `TRUE`.
#' @param detail `"summary"` for one row per CNV or `"groups"` for the capped
#'   criterion-group audit relation.
#' @return A data frame.
#' @export
score_cnv_2019 <- function(evidence, detail = c("summary", "groups")) {
  detail <- match.arg(detail)
  groups <- .score_cnv_2019_groups(evidence)
  if (identical(detail, "groups")) return(groups)
  if (!nrow(groups)) {
    return(data.frame(
      cnv_id = character(), cnv_type = character(), score = numeric(),
      classification = classify_cnv_score(numeric())
    ))
  }
  key <- interaction(groups$cnv_id, groups$cnv_type, drop = TRUE, lex.order = TRUE)
  split_rows <- split(seq_len(nrow(groups)), key)
  rows <- lapply(split_rows, function(i) {
    score <- sum(groups$applied_score[i])
    data.frame(
      cnv_id = groups$cnv_id[i[[1L]]], cnv_type = groups$cnv_type[i[[1L]]],
      score = score, classification = as.character(classify_cnv_score(score)),
      criterion_groups = length(i), stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$cnv_id, out$cnv_type), , drop = FALSE]
}

.one_provenance_value <- function(x, column, label) {
  if (!column %in% names(x)) stop(label, " is missing provenance column ", column, call. = FALSE)
  value <- unique(as.character(x[[column]]))
  value <- value[!is.na(value) & nzchar(value)]
  if (length(value) != 1L) stop(label, " must contain one explicit ", column, call. = FALSE)
  value
}

#' Derive objective ACMG/ClinGen 2019 CNV content evidence
#'
#' Derive only the criteria determined mechanically by gene content and
#' established ClinGen dosage records: section 1 content, section 2A complete
#' dosage-feature containment, and section 3 protein-coding gene count. Case,
#' phenotype, segregation, population, partial-gene, and inheritance evidence
#' are intentionally not inferred. Inputs are BED half-open relations and must
#' carry source identifiers and releases.
#'
#' @param cnvs Data frame with `cnv_id`, `chr`, `start`, `end`, and `cnv_type`
#'   (`"loss"` or `"gain"`).
#' @param genes Gene relation with `chr`, `start`, `end`, `gene_id`, `biotype`,
#'   `source_id`, and `source_release`.
#' @param dosage ClinGen dosage relation with `chr`, `start`, `end`,
#'   `record_id`, `hi_score`, `ts_score`, `source_id`, and `source_release`.
#' @param functional_biotypes Biotypes counted as functionally important for
#'   sections 1 and 3; defaults to `"protein_coding"`.
#' @param established_score ClinGen score treated as established; defaults to 3.
#' @return A provenance-complete evidence data frame accepted by
#'   `score_cnv_2019()`.
#' @export
derive_cnv_2019_evidence <- function(
  cnvs, genes, dosage, functional_biotypes = "protein_coding",
  established_score = 3
) {
  cnv_required <- c("cnv_id", "chr", "start", "end", "cnv_type")
  gene_required <- c("chr", "start", "end", "gene_id", "biotype", "source_id", "source_release")
  dosage_required <- c(
    "chr", "start", "end", "record_id", "hi_score", "ts_score",
    "source_id", "source_release"
  )
  requirements <- list(cnvs = cnv_required, genes = gene_required, dosage = dosage_required)
  objects <- list(cnvs = cnvs, genes = genes, dosage = dosage)
  for (label in names(requirements)) {
    missing <- setdiff(requirements[[label]], names(objects[[label]]))
    if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  cnvs <- as.data.frame(cnvs, stringsAsFactors = FALSE)
  genes <- as.data.frame(genes, stringsAsFactors = FALSE)
  dosage <- as.data.frame(dosage, stringsAsFactors = FALSE)
  if (anyDuplicated(cnvs$cnv_id) || any(!cnvs$cnv_type %in% c("loss", "gain"))) {
    stop("cnv_id must be unique and cnv_type must be loss or gain", call. = FALSE)
  }
  gene_source <- .one_provenance_value(genes, "source_id", "genes")
  gene_release <- .one_provenance_value(genes, "source_release", "genes")
  dosage_source <- .one_provenance_value(dosage, "source_id", "dosage")
  dosage_release <- .one_provenance_value(dosage, "source_release", "dosage")
  functional <- genes[genes$biotype %in% functional_biotypes, , drop = FALSE]
  context <- .new_overlap_context("duckdb")
  on.exit(.close_overlap_context(context), add = TRUE)
  gene_pairs <- .bed_overlap_pairs(cnvs, functional, context = context)
  dosage_pairs <- .bed_overlap_pairs(cnvs, dosage, context = context)
  rows <- list()
  add <- function(cnv, criterion, score, source_id, source_release, evidence_id) {
    rows[[length(rows) + 1L]] <<- data.frame(
      cnv_id = as.character(cnv$cnv_id), cnv_type = as.character(cnv$cnv_type),
      criterion = criterion, score = score, source_id = source_id,
      source_release = source_release, evidence_id = evidence_id,
      selected = TRUE, stringsAsFactors = FALSE
    )
  }
  for (i in seq_len(nrow(cnvs))) {
    cnv <- cnvs[i, , drop = FALSE]
    gene_index <- unique(gene_pairs$b_id[gene_pairs$a_id == i & !is.na(gene_pairs$b_id)])
    gene_ids <- sort(unique(as.character(functional$gene_id[gene_index])))
    gene_count <- length(gene_ids)
    if (gene_count) {
      add(cnv, "1A", 0, gene_source, gene_release, paste0("genes:", paste(gene_ids, collapse = ",")))
    } else {
      add(cnv, "1B", -0.60, gene_source, gene_release, "genes:none")
    }
    gene_criterion <- if (identical(cnv$cnv_type[[1L]], "loss")) {
      if (gene_count >= 35L) "3C" else if (gene_count >= 25L) "3B" else "3A"
    } else {
      if (gene_count >= 50L) "3C" else if (gene_count >= 35L) "3B" else "3A"
    }
    gene_score <- c(`3A` = 0, `3B` = .45, `3C` = .90)[[gene_criterion]]
    add(cnv, gene_criterion, gene_score, gene_source, gene_release, paste0("gene-count:", gene_count))

    dosage_index <- unique(dosage_pairs$b_id[dosage_pairs$a_id == i & !is.na(dosage_pairs$b_id)])
    if (length(dosage_index)) {
      score_column <- if (identical(cnv$cnv_type[[1L]], "loss")) "hi_score" else "ts_score"
      established <- suppressWarnings(as.numeric(dosage[[score_column]][dosage_index])) == established_score
      fully_contained <-
        as.numeric(cnv$start[[1L]]) <= as.numeric(dosage$start[dosage_index]) &
        as.numeric(cnv$end[[1L]]) >= as.numeric(dosage$end[dosage_index])
      records <- dosage_index[established & fully_contained]
      if (length(records)) {
        ids <- sort(unique(as.character(dosage$record_id[records])))
        add(cnv, "2A", 1, dosage_source, dosage_release, paste0("dosage:", paste(ids, collapse = ",")))
      }
    }
  }
  if (!length(rows)) {
    return(data.frame(
      cnv_id = character(), cnv_type = character(), criterion = character(),
      score = numeric(), source_id = character(), source_release = character(),
      evidence_id = character(), selected = logical()
    ))
  }
  do.call(rbind, rows)
}
