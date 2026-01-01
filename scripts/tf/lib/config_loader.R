# config_loader.R - Configuration and path management for TF pipeline
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)
# Created: 2025-12-31
# Purpose: Provide consistent config loading and path resolution for TF scripts
#
# Usage:
#   source("lib/config_loader.R")  # from scripts/tf/
#
# After sourcing, these variables are available:
#   - repo_root: Absolute path to repository root
#   - cfg: Parsed YAML config (list)
#   - paths$input_dir: Absolute path to inputs directory
#   - paths$output_dir: Absolute path to outputs directory
#   - paths$cache_dir: Absolute path to cache directory
#
# Helper functions:
#   - cfg_get(path, default): Get nested config value with fallback
#   - input_path(filename): Get absolute path to file in inputs dir
#   - output_path(filename): Get absolute path to file in outputs dir
#   - ensure_dir(path): Create directory if it doesn't exist

suppressPackageStartupMessages({
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required for config loading. Install with: install.packages('yaml')")
  }
  library(yaml)
})

# -----------------------------------------------------------------------------
# Repository Root Detection
# -----------------------------------------------------------------------------

# Find repo root by searching upward for genesets_config.yaml
.find_repo_root <- function(start_dir, max_levels = 10) {
  cur <- normalizePath(start_dir, mustWork = FALSE)
  for (i in seq_len(max_levels)) {
    config_path <- file.path(cur, "genesets_config.yaml")
    if (file.exists(config_path)) {
      return(cur)
    }
    parent <- dirname(cur)
    if (parent == cur) break
    cur <- parent
  }
  stop("Could not find repository root (genesets_config.yaml) starting from: ", start_dir)
}

# Determine script directory robustly (works for Rscript, source(), and interactive)
.get_script_dir <- function() {
  # Try --file= argument (Rscript command line)
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg_idx <- grep("--file=", cmd_args)
  if (length(file_arg_idx) > 0) {
    script_path <- sub("--file=", "", cmd_args[file_arg_idx[1]])
    return(dirname(normalizePath(script_path, mustWork = FALSE)))
  }

  # Try source() context
  if (sys.nframe() > 0) {
    for (i in seq_len(sys.nframe())) {
      if (!is.null(sys.call(i)) &&
          as.character(sys.call(i)[[1]]) %in% c("source", "sys.source")) {
        src_file <- sys.call(i)[[2]]
        if (is.character(src_file) && file.exists(src_file)) {
          return(dirname(normalizePath(src_file, mustWork = FALSE)))
        }
      }
    }
  }

  # Fallback to current working directory
  return(normalizePath(getwd(), mustWork = FALSE))
}

# Set repo_root
repo_root <- .find_repo_root(.get_script_dir())

# -----------------------------------------------------------------------------
# Configuration Loading
# -----------------------------------------------------------------------------

# Load config file
.config_file <- file.path(repo_root, "genesets_config.yaml")
if (!file.exists(.config_file)) {
  warning("Config file not found at: ", .config_file, ". Using defaults.")
  cfg <- list()
} else {
  cfg <- yaml::read_yaml(.config_file)
}

#' Get nested config value with fallback
#' @param path Dot-separated path (e.g., "input_files.manning")
#' @param default Default value if path not found
#' @return Config value or default
cfg_get <- function(path, default = NULL) {
  parts <- strsplit(path, "\\.")[[1]]
  cur <- cfg
  for (p in parts) {
    if (is.null(cur) || is.null(cur[[p]])) return(default)
    cur <- cur[[p]]
  }
  return(cur)
}

# -----------------------------------------------------------------------------
# Path Configuration
# -----------------------------------------------------------------------------

# Define TF-specific paths (override kinase defaults)
paths <- list(
  input_dir = file.path(repo_root, cfg_get("tf.input_dir", "curated/tf/inputs")),
  output_dir = file.path(repo_root, cfg_get("tf.output_dir", "curated/tf/outputs")),
  cache_dir = file.path(repo_root, "cache")
)

# Normalize all paths
paths <- lapply(paths, function(p) normalizePath(p, mustWork = FALSE))

#' Ensure directory exists (create if needed)
#' @param dir_path Path to directory
#' @return TRUE if directory exists or was created
ensure_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(dir_path)) {
    stop("Failed to create directory: ", dir_path)
  }
  return(TRUE)
}

#' Get absolute path to file in inputs directory
#' @param filename Filename (without directory)
#' @return Absolute path
input_path <- function(filename) {
  file.path(paths$input_dir, filename)
}

#' Get absolute path to file in outputs directory
#' @param filename Filename (without directory)
#' @return Absolute path
output_path <- function(filename) {
  file.path(paths$output_dir, filename)
}

# -----------------------------------------------------------------------------
# Gene Exclusion List
# -----------------------------------------------------------------------------

#' Get list of genes to exclude from validation
#' @return Character vector of gene symbols to exclude
get_exclude_genes <- function() {

  exclude_list <- cfg_get("validation.exclude_genes", character(0))
  if (is.null(exclude_list)) return(character(0))
  return(as.character(exclude_list))
}

# -----------------------------------------------------------------------------
# File Validation Helpers
# -----------------------------------------------------------------------------

#' Check that file exists, stop with message if not
#' @param f File path
#' @param msg Optional custom error message
check_file_exists <- function(f, msg = NULL) {
  if (!file.exists(f)) {
    stop(if (is.null(msg)) paste0("Missing required file: ", f) else msg, call. = FALSE)
  }
}

#' Check that file exists and is non-empty, stop with message if not
#' @param f File path
#' @param msg Optional custom error message
check_file_nonempty <- function(f, msg = NULL) {
  if (!file.exists(f) || file.info(f)$size == 0) {
    stop(if (is.null(msg)) paste0("File missing or empty: ", f) else msg, call. = FALSE)
  }
}

# -----------------------------------------------------------------------------
# Initialization Output
# -----------------------------------------------------------------------------

# Print configuration summary (can be suppressed with TF_QUIET=1)
if (Sys.getenv("TF_QUIET", "0") != "1") {
  cat(sprintf("[config_loader] repo_root: %s\n", repo_root))
  cat(sprintf("[config_loader] input_dir: %s\n", paths$input_dir))
  cat(sprintf("[config_loader] output_dir: %s\n", paths$output_dir))
  cat(sprintf("[config_loader] species: %s\n", cfg_get("species", "human")))
}

# Ensure output directory exists
ensure_dir(paths$output_dir)
