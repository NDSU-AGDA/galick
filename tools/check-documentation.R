# Parse every published R fence and check named arguments against the public API.
pkgload::load_all(".", quiet = TRUE)
files <- c("README.md", list.files("docs", pattern = "[.]md$", recursive = TRUE, full.names = TRUE))
exports <- getNamespaceExports("galick")
checked <- 0L
snippets <- list()
inspect <- function(expression, location) {
  if (missing(expression)) {
    return(invisible(NULL))
  }
  if (!is.call(expression)) {
    return(invisible(NULL))
  }
  head <- expression[[1]]
  name <- if (is.symbol(head)) as.character(head) else ""
  if (name %in% exports) {
    fun <- getExportedValue("galick", name)
    supplied <- names(as.list(expression)[-1])
    supplied <- supplied[nzchar(supplied)]
    accepted <- names(formals(fun))
    if (!"..." %in% accepted && any(!supplied %in% accepted)) {
      stop(location, ": unsupported arguments to ", name, ": ", paste(setdiff(supplied, accepted), collapse = ", "))
    }
  }
  for (part in as.list(expression)[-1]) {
    inspect(part, location)
  }
  invisible(NULL)
}
for (file in files) {
  lines <- readLines(file, warn = FALSE)
  inside <- FALSE
  code <- character()
  snippets[[file]] <- list()
  for (i in seq_along(lines)) {
    line <- trimws(lines[i])
    if (!inside && identical(line, "```r")) {
      inside <- TRUE
      start <- i
      code <- character()
    } else if (inside && identical(line, "```")) {
      location <- paste0(file, ":", start)
      expressions <- tryCatch(parse(text = code), error = function(e) stop(location, ": ", conditionMessage(e)))
      for (expression in expressions) inspect(expression, location)
      snippets[[file]][[length(snippets[[file]]) + 1L]] <- expressions
      checked <- checked + 1L
      inside <- FALSE
    } else if (inside) {
      code <- c(code, lines[i])
    }
  }
  if (inside) stop(file, ": unclosed R fence")
}
cat("Parsed and checked", checked, "R examples across", length(files), "documents.\n")

if ("--run" %in% commandArgs(trailingOnly = TRUE)) {
  has_beagle <- nzchar(Sys.getenv("BEAGLE_PATH")) && file.exists(Sys.getenv("BEAGLE_PATH")) &&
    nzchar(Sys.which("java"))
  plans <- list(
    "README.md" = if (has_beagle) c(2L, 3L) else 2L,
    "docs/index.md" = 1L,
    "docs/getting-started/quick-tour.md" = if (has_beagle) 1:3 else 1:2,
    "docs/getting-started/result-classes.md" = 1:2,
    "docs/examples/csv-to-vcf.md" = seq_along(snippets[["docs/examples/csv-to-vcf.md"]]),
    "docs/examples/hapmap-to-vcf.md" = seq_along(snippets[["docs/examples/hapmap-to-vcf.md"]])
  )
  if (has_beagle) {
    plans[["docs/examples/benchmark-accuracy.md"]] <- 1:2
    plans[["docs/examples/run-beagle.md"]] <- c(1L, 2L, 4L)
  }
  executed <- 0L
  for (file in names(plans)) {
    environment <- new.env(parent = globalenv())
    for (index in plans[[file]]) {
      tryCatch(eval(snippets[[file]][[index]], envir = environment),
        error = function(error) stop(file, " block ", index, ": ", conditionMessage(error))
      )
      executed <- executed + 1L
    }
  }
  cat("Executed", executed, "bundled-data example blocks; Beagle available:", has_beagle, "\n")
}
