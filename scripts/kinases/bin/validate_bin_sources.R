#!/usr/bin/env Rscript
# Validator: check `scripts/kinases/bin/*.R` for `source()` targets and executable bits.
# Outputs a short report to stdout.
suppressPackageStartupMessages({})

files <- list.files("scripts/kinases/bin", pattern = "\\.R$", full.names = TRUE)
if (length(files) == 0) {
  cat("No bin R scripts found.\n")
  quit(status = 0)
}

for (f in files) {
  cat("\nFile:", f, "\n")
  lines <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
  src_idx <- grep("source\\s*\\(", lines)
  if (length(src_idx) == 0) {
    cat("  (no source() lines)\n")
  } else {
    for (i in src_idx) cat("  ", i, ": ", lines[i], "\n")
    quoted <- unlist(lapply(lines[src_idx], function(l) {
      m <- gregexpr("['\"][^'\"]+['\"]", l)
      s <- regmatches(l, m)[[1]]
      s <- gsub("^['\"]|['\"]$", "", s)
      s
    }))
    if (length(quoted) == 0) {
      cat("  (no quoted paths found)\n")
    } else {
      for (q in quoted) {
        if (file.exists(q)) {
          cat("  OK:", q, "(absolute)\n")
        } else if (file.exists(file.path("scripts/kinases", q))) {
          cat("  OK: scripts/kinases/", q, "\n")
        } else if (file.exists(file.path("scripts/kinases/lib", basename(q)))) {
          cat("  OK: scripts/kinases/lib/", basename(q), "\n")
        } else {
          cat("  MISSING:", q, "\n")
        }
      }
    }
  }
  exec <- ifelse(file.access(f, mode = 1) == 0, "yes", "no")
  cat("  Executable:", exec, "\n")
}

cat("\nValidator finished.\n")
