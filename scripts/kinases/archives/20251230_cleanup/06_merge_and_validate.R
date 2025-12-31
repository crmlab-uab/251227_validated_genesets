#!/usr/bin/env Rscript
# 05_merge_and_validate.R
# Runner that merges validation sources and runs validations for the kinases pipeline.

suppressWarnings(suppressMessages({
	library(data.table)
	library(optparse)
}))

cat("Starting merge_and_validate runner...\n", file=stderr())

# Prefer canonical bin entrypoint if present, otherwise fallback to lib script
if (file.exists("scripts/kinases/bin/05_merge_and_validate.R")) {
	# avoid self-sourcing; source lib instead
	if (file.exists("scripts/kinases/lib/fetch_val_sources_and_merge.R")) {
		source("scripts/kinases/lib/fetch_val_sources_and_merge.R")
	} else if (file.exists("kinases/fetch_val_sources_and_merge.R")) {
		source("kinases/fetch_val_sources_and_merge.R")
	} else {
		stop("No merge/validate implementation found (expected scripts/kinases/lib/fetch_val_sources_and_merge.R)")
	}
} else {
	# fallback
	if (file.exists("scripts/kinases/lib/fetch_val_sources_and_merge.R")) source("scripts/kinases/lib/fetch_val_sources_and_merge.R") else stop("No merge/validate implementation found")
}

cat("merge_and_validate runner finished.\n", file=stderr())

# Ensure logs for this run are recorded in sessions/
## compute repo_root so logs/sessions are placed at repo root
cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArgIdx <- grep("--file=", cmdArgs)
fileArg <- if (length(fileArgIdx) > 0) sub("--file=", "", cmdArgs[fileArgIdx[1]]) else NA_character_
script_dir <- if (!is.na(fileArg) && nzchar(fileArg)) dirname(normalizePath(fileArg)) else normalizePath(".")
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))

sess_dir <- file.path(repo_root, "sessions", format(Sys.time(), "%y%m%d_%H%M%S"))
dir.create(sess_dir, recursive=TRUE, showWarnings=FALSE)
logf <- file.path(sess_dir, "05_merge_and_validate.log")
cat(sprintf("Log: %s\n", logf), file=stderr())

# Append a timestamped copy of console output
sink(logf, append=TRUE, split=TRUE)
sink(logf, type="message", append=TRUE)
on.exit({
	sink(type="message")
	sink()
}, add=TRUE)
