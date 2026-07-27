.xcnv_compatibility_commit <- "2237e765fa9c31555dd8f4591b7f2693cff1ff8a"
.xcnv_compatibility_release <- paste0("upstream-", substr(.xcnv_compatibility_commit, 1L, 12L))
.xcnv_model_id <- "zhang-2021-xcnv-tree-v1"
.xcnv_policy_id <- "acmg-clingen-2019-constitutional-cnv"
.xcnv_provider_id <- "xcnv"
.xcnv_acmg_provider_id <- "xcnv-acmg-clingen-2019"

.xcnv_canonical_assembly <- function(assembly) {
  if (length(assembly) != 1L || is.na(assembly) || !nzchar(as.character(assembly))) {
    return(NA_character_)
  }
  value <- trimws(as.character(assembly))
  aliases <- c(hg19 = "GRCh37", `GRCh38.p14` = "GRCh38")
  if (value %in% names(aliases)) aliases[[value]] else value
}

.xcnv_relation_assembly <- function(x, label) {
  if (!nrow(x) || !"assembly" %in% names(x)) return(NA_character_)
  values <- unique(vapply(as.character(x$assembly), .xcnv_canonical_assembly, character(1)))
  if (anyNA(values) || length(values) != 1L) {
    stop(label, " must declare one non-missing assembly", call. = FALSE)
  }
  values[[1L]]
}

.xcnv_validate_relation_assemblies <- function(objects, labels) {
  assemblies <- vapply(
    seq_along(objects),
    function(i) .xcnv_relation_assembly(objects[[i]], labels[[i]]),
    character(1)
  )
  nonempty <- vapply(objects, nrow, integer(1)) > 0L
  declared <- !is.na(assemblies)
  if (any(declared) && any(nonempty & !declared)) {
    stop("nonempty CNV, gene, and dosage relations must all declare assembly", call. = FALSE)
  }
  if (sum(declared) > 1L && length(unique(assemblies[declared])) != 1L) {
    stop("CNV, gene, and dosage relations have conflicting assemblies", call. = FALSE)
  }
  assembly <- if (any(declared)) assemblies[[which(declared)[[1L]]]] else NA_character_
  status <- if (is.na(assembly)) {
    "unspecified"
  } else if (assembly %in% c("GRCh37", "GRCh38")) {
    "supported"
  } else {
    "unsupported"
  }
  list(assembly = assembly, status = status)
}

.xcnv_resource_contract <- function(assembly = "GRCh37") {
  assembly <- .xcnv_canonical_assembly(assembly)
  status <- if (identical(assembly, "GRCh37")) "supported" else "unsupported"
  list(
    provider_id = .xcnv_provider_id,
    assembly = assembly,
    resource_release = .xcnv_compatibility_release,
    model_id = .xcnv_model_id,
    policy_id = NA_character_,
    status = status
  )
}
