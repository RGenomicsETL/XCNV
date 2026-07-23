#' @export
print.xcnv_resources <- function(x, ...) {
  cat("XCNV resource bundle\n")
  cat("root:", x$root, "\n")
  present <- vapply(x$files, function(path) !is.na(path), logical(1))
  cat("files:", sum(present), "/", length(present), "resolved\n")
  invisible(x)
}
