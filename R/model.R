.xcnv_model_features <- c(
  "SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred",
  "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred",
  "RadialSVM_pred", "LR_pred", "VEST3_score", "CADD_phred", "GERP++_RS",
  "phyloP46way_placental", "phyloP100way_vertebrate", "SiPhy_29way_logOdds",
  "CDTS_1st", "CDTS_5th", "gain.freq", "loss.freq", "pELS", "CTCF-bound",
  "PLS", "dELS", "CTCF-only", "DNase-H3K4me3", "pLI", "Episcore", "GHIS",
  "Type", "Length"
)

.portable_model_metadata <- function(path) {
  comments <- readLines(path, warn = FALSE)
  comments <- comments[grepl("^#", comments)]
  fields <- strsplit(sub("^#\\s*", "", comments), "=", fixed = TRUE)
  fields <- fields[lengths(fields) == 2L]
  out <- vapply(fields, `[[`, character(1), 2L)
  names(out) <- vapply(fields, `[[`, character(1), 1L)
  out
}

.read_portable_model <- function(path) {
  metadata <- .portable_model_metadata(path)
  required_metadata <- c("format", "objective", "base_score", "feature_names")
  missing_metadata <- setdiff(required_metadata, names(metadata))
  if (length(missing_metadata)) {
    stop("portable model is missing metadata: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
  }
  if (!identical(metadata[["format"]], "xcnv-tree-v1")) {
    stop("unsupported portable model format: ", metadata[["format"]], call. = FALSE)
  }
  if (!identical(metadata[["objective"]], "reg:logistic")) {
    stop("unsupported portable model objective: ", metadata[["objective"]], call. = FALSE)
  }
  base_score <- suppressWarnings(as.numeric(metadata[["base_score"]]))
  if (length(base_score) != 1L || is.na(base_score) || base_score <= 0 || base_score >= 1) {
    stop("portable model base_score must be strictly between zero and one", call. = FALSE)
  }
  feature_names <- strsplit(metadata[["feature_names"]], ",", fixed = TRUE)[[1L]]
  feature_names <- trimws(feature_names)
  if (!length(feature_names) || any(!nzchar(feature_names)) || anyDuplicated(feature_names)) {
    stop("portable model feature_names are invalid", call. = FALSE)
  }
  nodes <- data.table::fread(
    path, sep = "\t", header = TRUE, comment.char = "#", na.strings = c("", "NA"),
    data.table = FALSE, showProgress = FALSE
  )
  required_columns <- c(
    "tree_id", "node_id", "feature", "threshold", "yes_id", "no_id",
    "missing_id", "leaf_value"
  )
  if (!identical(names(nodes), required_columns)) {
    stop("portable model columns must be: ", paste(required_columns, collapse = ", "), call. = FALSE)
  }
  nodes$tree_id <- suppressWarnings(as.integer(nodes$tree_id))
  nodes$node_id <- suppressWarnings(as.integer(nodes$node_id))
  nodes$yes_id <- suppressWarnings(as.integer(nodes$yes_id))
  nodes$no_id <- suppressWarnings(as.integer(nodes$no_id))
  nodes$missing_id <- suppressWarnings(as.integer(nodes$missing_id))
  nodes$threshold <- suppressWarnings(as.numeric(nodes$threshold))
  nodes$leaf_value <- suppressWarnings(as.numeric(nodes$leaf_value))
  if (!nrow(nodes) || anyNA(nodes$tree_id) || anyNA(nodes$node_id) ||
      anyDuplicated(paste(nodes$tree_id, nodes$node_id, sep = ":"))) {
    stop("portable model tree/node identifiers are invalid", call. = FALSE)
  }
  leaf <- !is.na(nodes$leaf_value)
  branch <- !leaf
  if (any(!is.na(nodes$feature[leaf])) ||
      any(!is.na(nodes$threshold[leaf])) ||
      any(!is.na(nodes$yes_id[leaf])) ||
      any(!is.na(nodes$no_id[leaf])) ||
      any(!is.na(nodes$missing_id[leaf]))) {
    stop("portable model leaf rows must contain only leaf_value", call. = FALSE)
  }
  if (any(is.na(nodes$feature[branch])) || any(!nodes$feature[branch] %in% feature_names) ||
      anyNA(nodes$threshold[branch]) || anyNA(nodes$yes_id[branch]) ||
      anyNA(nodes$no_id[branch]) || anyNA(nodes$missing_id[branch])) {
    stop("portable model branch rows are invalid", call. = FALSE)
  }
  for (tree in unique(nodes$tree_id)) {
    tree_nodes <- nodes[nodes$tree_id == tree, , drop = FALSE]
    if (!0L %in% tree_nodes$node_id) stop("every tree must contain root node 0", call. = FALSE)
    child_ids <- c(tree_nodes$yes_id[!is.na(tree_nodes$yes_id)],
                   tree_nodes$no_id[!is.na(tree_nodes$no_id)],
                   tree_nodes$missing_id[!is.na(tree_nodes$missing_id)])
    if (any(!child_ids %in% tree_nodes$node_id)) {
      stop("portable model branch references a missing child node", call. = FALSE)
    }
  }
  structure(
    list(
      model_type = "xcnv-tree-v1", objective = metadata[["objective"]],
      base_score = base_score, feature_names = feature_names, nodes = nodes,
      metadata = metadata
    ),
    class = "xcnv_portable_model"
  )
}

#' Load an X-CNV model
#'
#' Load a portable, tabular X-CNV tree model. The portable artifact is exported
#' from the pinned model during resource staging and is evaluated without an
#' XGBoost runtime. Serialized `.Rdata`, `.rds`, and native XGBoost files are
#' deliberately rejected; use `tools/stage_xcnv_model.R` to convert them in an
#' environment where XGBoost is installed.
#'
#' @param model A portable model object or path to an `xcnv-tree-v1` TSV file.
#'   The bundled published model is used when omitted or `NULL`.
#' @return A model object suitable for `predict_cnv()`.
#' @export
load_xcnv_model <- function(model = NULL) {
  if (is.null(model)) return(xcnv_bundled_model())
  if (!is.character(model) || length(model) != 1L || is.na(model) || !file.exists(model)) {
    if (is.list(model)) return(model)
    stop("'model' must be an existing model path or model object", call. = FALSE)
  }
  if (!grepl("\\.tsv$", model, ignore.case = TRUE)) {
    stop(
      "native XGBoost/R serialized models are staging inputs, not runtime models; convert with tools/stage_xcnv_model.R",
      call. = FALSE
    )
  }
  .read_portable_model(model)
}

#' Return the bundled published X-CNV model
#'
#' @return A validated `xcnv_portable_model` object.
#' @export
xcnv_bundled_model <- function() {
  data_environment <- new.env(parent = emptyenv())
  loaded <- utils::data(
    "xcnv_portable_model",
    package = "XCNV",
    envir = data_environment
  )
  if (!"xcnv_portable_model" %in% loaded) {
    stop("the bundled X-CNV model is not available", call. = FALSE)
  }
  model <- data_environment$xcnv_portable_model
  if (!inherits(model, "xcnv_portable_model")) {
    stop("the bundled X-CNV model is invalid", call. = FALSE)
  }
  model
}

.model_feature_names <- function(model) {
  if (inherits(model, "xcnv_portable_model")) return(model$feature_names)
  if (is.list(model) && !is.null(model$feature_names)) return(model$feature_names)
  character()
}

.predict_portable_tree <- function(model, matrix) {
  margin <- rep.int(stats::qlogis(model$base_score), nrow(matrix))
  for (tree in sort(unique(model$nodes$tree_id))) {
    nodes <- model$nodes[model$nodes$tree_id == tree, , drop = FALSE]
    rownames(nodes) <- as.character(nodes$node_id)
    current <- rep.int(0L, nrow(matrix))
    active <- rep.int(TRUE, nrow(matrix))
    tree_value <- numeric(nrow(matrix))
    steps <- 0L
    while (any(active)) {
      steps <- steps + 1L
      if (steps > nrow(nodes) + 1L) stop("portable model tree contains a cycle", call. = FALSE)
      for (node_id in unique(current[active])) {
        rows <- which(active & current == node_id)
        node <- nodes[as.character(node_id), , drop = FALSE]
        if (!nrow(node)) stop("portable model reached an unknown node", call. = FALSE)
        if (!is.na(node$leaf_value[[1L]])) {
          tree_value[rows] <- node$leaf_value[[1L]]
          active[rows] <- FALSE
        } else {
          values <- matrix[rows, node$feature[[1L]]]
          next_id <- ifelse(
            is.na(values), node$missing_id[[1L]],
            ifelse(values < node$threshold[[1L]], node$yes_id[[1L]], node$no_id[[1L]])
          )
          current[rows] <- as.integer(next_id)
        }
      }
    }
    margin <- margin + tree_value
  }
  stats::plogis(margin)
}

.feature_matrix <- function(features, model) {
  wanted <- .model_feature_names(model)
  if (!length(wanted)) {
    wanted <- .xcnv_model_features
  }
  values <- vector("list", length(wanted))
  names(values) <- wanted
  for (i in seq_along(wanted)) {
    key <- wanted[[i]]
    source <- if (identical(key, "Type")) "Type.1" else key
    if (!source %in% names(features)) {
      stop("annotation output is missing model feature: ", key, call. = FALSE)
    }
    values[[i]] <- suppressWarnings(as.numeric(features[[source]]))
  }
  matrix(unlist(values, use.names = FALSE), ncol = length(wanted),
         dimnames = list(NULL, wanted))
}

#' Predict X-CNV MVP scores
#'
#' Annotate CNVs and compute the X-CNV MVP score. With an upstream model this
#' evaluates a staged `xcnv-tree-v1` artifact in ordinary R. Installing and
#' running the package does not require XGBoost.
#'
#' @param x A path, data frame, or matrix containing CNV records.
#' @param resources An `xcnv_resources` object or an explicit resource directory.
#' @param model A model object/path. If `NULL`, a model in `resources` is used
#'   when present, otherwise the bundled published model is used.
#' @param overlap_backend Passed to `annotate_cnv()`.
#' @return A data frame containing the legacy annotation columns and a final
#'   `MVP_score` column.
#' @export
predict_cnv <- function(
  x, resources, model = NULL, overlap_backend = c("duckdb", "reference")
) {
  overlap_backend <- match.arg(overlap_backend)
  resources <- if (inherits(resources, "xcnv_resources")) {
    resources
  } else {
    xcnv_resources(resources, require = if (is.null(model)) "all" else "annotations")
  }
  if (is.null(model)) {
    model_path <- resources$files$model
    model <- if (!is.null(model_path) && !is.na(model_path)) {
      load_xcnv_model(model_path)
    } else {
      xcnv_bundled_model()
    }
  }
  model <- if (is.character(model)) load_xcnv_model(model) else model
  features <- annotate_cnv(x, resources, overlap_backend = overlap_backend)
  matrix <- .feature_matrix(features, model)
  if (inherits(model, "xcnv_portable_model")) {
    score <- .predict_portable_tree(model, matrix)
  } else {
    stop("model must be a validated xcnv-tree-v1 portable model", call. = FALSE)
  }
  score <- as.numeric(score)
  if (length(score) != nrow(features)) {
    stop("model returned one score per CNV is required", call. = FALSE)
  }
  features$MVP_score <- score
  features
}
