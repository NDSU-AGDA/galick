# Regenerate website signatures from the package namespace.
pkgload::load_all(".", quiet = TRUE)
exports <- sort(getNamespaceExports("galick"))
exports <- exports[!grepl(".", exports, fixed = TRUE)]
lines <- c(
  "# API signatures", "",
  "Generated from the package source by `tools/build-reference.R`.",
  "Use `?function_name` in R for complete parameter documentation.",
  "Signatures below are declarations; required inputs are not example data.", ""
)
for (name in exports) {
  parameters <- formals(getExportedValue("galick", name))
  fields <- vapply(seq_along(parameters), function(i) {
    value <- paste(deparse(parameters[[i]], width.cutoff = 80L), collapse = " ")
    if (!nzchar(value)) names(parameters)[i] else paste(names(parameters)[i], "=", value)
  }, character(1))
  signature <- if (length(fields) <= 2L) {
    paste0(name, "(", paste(fields, collapse = ", "), ")")
  } else {
    c(paste0(name, "("), paste0("  ", fields, c(rep(",", length(fields) - 1L), "")), ")")
  }
  lines <- c(lines, paste0("## `", name, "()`"), "", "```r", signature, "```", "")
}
writeLines(lines[-length(lines)], "docs/reference/api.md")
