#!/usr/bin/env Rscript
# 07_generate_coral_tree.R
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)
# Created: 2025-12-31
# Purpose: Generate CORAL kinome tree SVG figures from expression data
# Usage: Rscript 07_generate_coral_tree.R --log2fc=<DEG_all.csv> --comparison=<name>
#
# This script uses CORAL (https://github.com/dphansti/CORAL) rendering code
# to generate kinome tree SVG visualizations directly in R.
# Configuration is read from genesets_config.yaml (coral: section)

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# Load centralized config
.script_dir <- (function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg[1]))))
  }
  return(normalizePath("."))
})()
source(file.path(.script_dir, "lib", "config_loader.R"))

# Load CORAL settings from YAML config (with fallback defaults)
.cfg_get_safe <- function(key, default) {
  val <- tryCatch(cfg_get(key, default), error = function(e) default)
  if (is.null(val)) default else val
}

# Resolve CORAL path relative to repo root
.coral_path_raw <- .cfg_get_safe("coral.repo_path", "external/CORAL")
CORAL_PATH <- if (startsWith(.coral_path_raw, "/")) .coral_path_raw else file.path(repo_root, .coral_path_raw)
CORAL_URL <- .cfg_get_safe("coral.repo_url", "https://github.com/dphansti/CORAL.git")
CORAL_FC_LIMIT <- as.numeric(.cfg_get_safe("coral.fc_limit", 3.0))
CORAL_COLOR_LOW <- .cfg_get_safe("coral.colors.low", "blue")
CORAL_COLOR_MID <- .cfg_get_safe("coral.colors.mid", "white")
CORAL_COLOR_HIGH <- .cfg_get_safe("coral.colors.high", "red")
CORAL_COLOR_BG <- .cfg_get_safe("coral.colors.background", "#D3D3D3")
CORAL_NODE_RADIUS <- as.numeric(.cfg_get_safe("coral.node_radius", 5))
CORAL_NODE_RADIUS_MATCHED <- as.numeric(.cfg_get_safe("coral.node_radius_matched", 6))
CORAL_OUTPUT_SVG <- as.logical(.cfg_get_safe("coral.output_svg", TRUE))
CORAL_OUTPUT_PDF <- as.logical(.cfg_get_safe("coral.output_pdf", TRUE))

# Parse command-line arguments (defaults from YAML config)
option_list <- list(
  make_option(c("-e", "--expression"), type = "character", default = NULL,
              help = "Path to expression data CSV (gene_id, gene_name, sample columns)"),
  make_option(c("-c", "--comparison"), type = "character", default = NULL,
              help = "Comparison name for differential expression (e.g., 'Treatment_vs_Control')"),
  make_option(c("-s", "--samples"), type = "character", default = NULL,
              help = "Comma-separated sample names to compare (e.g., 'TK028,TK029,TK030:TK049,TK050,TK051')"),
  make_option(c("-l", "--log2fc"), type = "character", default = NULL,
              help = "Path to pre-computed log2 fold-change file (gene_name, log2FC columns)"),
  make_option(c("-o", "--output-prefix"), type = "character", default = "kinome",
              help = "Output file prefix [default: %default]"),
  make_option(c("--min-expression"), type = "double", default = 10,
              help = "Minimum mean expression to include [default: %default]"),
  make_option(c("--fc-limit"), type = "double", default = CORAL_FC_LIMIT,
              help = "Symmetric fold-change limit for color scale [default: %default]"),
  make_option(c("--coral-path"), type = "character", default = CORAL_PATH,
              help = "Path to CORAL repository [default: %default]"),
  make_option(c("--title"), type = "character", default = NULL,
              help = "Title for the kinome tree figure")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
CORAL_PATH <- opt$`coral-path`
fc_limit <- opt$`fc-limit`

cat("=== CORAL Kinome Tree Generator ===\n\n")

# Check for CORAL repository
if (!dir.exists(CORAL_PATH)) {
  cat("CORAL repository not found at:", CORAL_PATH, "\n")
  cat("Cloning from:", CORAL_URL, "\n")
  system2("git", c("clone", CORAL_URL, CORAL_PATH),
          stdout = TRUE, stderr = TRUE)
}

# Source CORAL functions
coral_files <- c("colorby.R", "map2color.R", "writekinasetree.R", "legendfunctions.R")
for (f in coral_files) {
  fpath <- file.path(CORAL_PATH, "R", f)
  if (file.exists(fpath)) {
    source(fpath)
  } else {
    stop("Missing CORAL file: ", fpath)
  }
}

# Load CORAL dataframe (contains kinase tree structure)
coral_df_file <- file.path(CORAL_PATH, "Data", "coral_dataframe.tsv")
if (!file.exists(coral_df_file)) {
  stop("Missing CORAL dataframe: ", coral_df_file)
}
coral_df <- fread(coral_df_file, header = TRUE)
cat("Loaded CORAL tree with", nrow(coral_df), "kinases\n")

# Load kinase list with mouse-to-human mapping
cat("Loading kinase annotations...\n")
candidates <- list.files(paths$input_dir, pattern = "201006_composite_kinases_curated.*\\.csv$",
                         full.names = TRUE, ignore.case = TRUE)
if (length(candidates) == 0) {
  stop("Missing kinase annotation file in ", paths$input_dir)
}
kinases_file <- sort(candidates, decreasing = TRUE)[1]
kinases <- fread(kinases_file, header = TRUE)
cat("  Loaded", nrow(kinases), "kinases from", basename(kinases_file), "\n")

# Standardize column names
if ("Mouse_symbol" %in% names(kinases)) {
  setnames(kinases, "Mouse_symbol", "Mouse_Symbol")
}

# Filter out excluded genes
exclude_genes <- get_exclude_genes()
if (length(exclude_genes) > 0) {
  n_before <- nrow(kinases)
  kinases <- kinases[!Mouse_Symbol %in% exclude_genes]
  n_excluded <- n_before - nrow(kinases)
  if (n_excluded > 0) {
    cat("  Excluded", n_excluded, "genes per config\n")
  }
}

cat("  Kinases with HGNC mapping:", sum(!is.na(kinases$HGNC) & kinases$HGNC != ""), "\n\n")

# Determine input mode and load expression data
if (!is.null(opt$log2fc)) {
  # Mode 1: Pre-computed log2 fold-change
  cat("Loading pre-computed log2FC from:", opt$log2fc, "\n")
  fc_data <- fread(opt$log2fc, header = TRUE)

  # Find gene name and fold-change columns
  gene_col <- intersect(c("gene_name", "Gene_Name", "symbol", "Symbol", "gene"), names(fc_data))[1]
  fc_col <- intersect(c("log2FC", "log2FoldChange", "logFC", "log2fc", "FC"), names(fc_data))[1]

  if (is.na(gene_col) || is.na(fc_col)) {
    stop("Could not find gene name or fold-change columns. Available: ", paste(names(fc_data), collapse = ", "))
  }

  expression_values <- data.table(
    gene_name = fc_data[[gene_col]],
    value = fc_data[[fc_col]]
  )
  comparison_name <- if (!is.null(opt$comparison)) opt$comparison else "log2FC"

} else if (!is.null(opt$expression) && !is.null(opt$samples)) {
  # Mode 2: Calculate fold-change from expression matrix
  cat("Loading expression data from:", opt$expression, "\n")
  expr_data <- fread(opt$expression, header = TRUE)

  # Parse sample groups
  sample_groups <- strsplit(opt$samples, ":")[[1]]
  if (length(sample_groups) != 2) {
    stop("--samples must specify two groups separated by ':' (e.g., 'S1,S2:S3,S4')")
  }
  group1_samples <- strsplit(sample_groups[1], ",")[[1]]
  group2_samples <- strsplit(sample_groups[2], ",")[[1]]

  cat("  Group 1 samples:", paste(group1_samples, collapse = ", "), "\n")
  cat("  Group 2 samples:", paste(group2_samples, collapse = ", "), "\n")

  # Find gene name column
  gene_col <- intersect(c("gene_name", "Gene_Name", "symbol", "external_gene_name"), names(expr_data))[1]
  if (is.na(gene_col)) {
    stop("Could not find gene name column. Available: ", paste(names(expr_data), collapse = ", "))
  }

  # Calculate mean expression per group
  group1_mean <- rowMeans(as.matrix(expr_data[, ..group1_samples]), na.rm = TRUE)
  group2_mean <- rowMeans(as.matrix(expr_data[, ..group2_samples]), na.rm = TRUE)

  # Calculate log2 fold-change (group2 vs group1)
  pseudocount <- 1
  log2fc <- log2((group2_mean + pseudocount) / (group1_mean + pseudocount))

  # Filter by minimum expression
  mean_expr <- (group1_mean + group2_mean) / 2
  keep <- mean_expr >= opt$`min-expression`

  expression_values <- data.table(
    gene_name = expr_data[[gene_col]][keep],
    value = log2fc[keep]
  )

  comparison_name <- if (!is.null(opt$comparison)) opt$comparison else paste0(sample_groups[2], "_vs_", sample_groups[1])
  cat("  Calculated log2FC for", nrow(expression_values), "genes (min expr >=", opt$`min-expression`, ")\n")

} else if (!is.null(opt$expression)) {
  # Mode 3: Use mean expression across all samples
  cat("Loading expression data from:", opt$expression, "\n")
  cat("No comparison specified, using mean expression across all samples\n")

  expr_data <- fread(opt$expression, header = TRUE)

  # Find gene name column
  gene_col <- intersect(c("gene_name", "Gene_Name", "symbol", "external_gene_name"), names(expr_data))[1]
  if (is.na(gene_col)) {
    stop("Could not find gene name column. Available: ", paste(names(expr_data), collapse = ", "))
  }

  # Find numeric columns (expression values)
  numeric_cols <- names(expr_data)[sapply(expr_data, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, c("gene_id", "entrez_id", "Entrez_ID"))

  if (length(numeric_cols) == 0) {
    stop("No numeric sample columns found")
  }

  # Calculate mean expression
  mean_expr <- rowMeans(as.matrix(expr_data[, ..numeric_cols]), na.rm = TRUE)

  # Log2 transform
  log2_expr <- log2(mean_expr + 1)

  expression_values <- data.table(
    gene_name = expr_data[[gene_col]],
    value = log2_expr
  )

  comparison_name <- "mean_log2_expression"
  cat("  Calculated mean log2 expression for", nrow(expression_values), "genes\n")

} else {
  cat("\nUsage examples:\n")
  cat("  # From pre-computed log2FC:\n")
  cat("  Rscript 07_generate_coral_tree.R --log2fc=de_results.csv\n\n")
  cat("  # From expression matrix with sample comparison:\n")
  cat("  Rscript 07_generate_coral_tree.R --expression=counts.csv --samples='TK028,TK029:TK049,TK050'\n\n")
  cat("  # Mean expression across all samples:\n")
  cat("  Rscript 07_generate_coral_tree.R --expression=counts.csv\n\n")
  stop("Must provide --expression or --log2fc")
}

cat("\n")

# Map mouse gene names to human HGNC symbols
cat("Mapping mouse genes to HGNC symbols...\n")

# Create lookup: mouse symbol -> HGNC
mouse_to_hgnc <- kinases[!is.na(HGNC) & HGNC != "", .(Mouse_Symbol, HGNC)]
setnames(mouse_to_hgnc, "Mouse_Symbol", "gene_name")

# Merge expression with kinase mapping
coral_input <- merge(expression_values, mouse_to_hgnc, by = "gene_name", all = FALSE)
cat("  Matched", nrow(coral_input), "kinases with expression data\n")

if (nrow(coral_input) == 0) {
  # Try direct HGNC matching (maybe expression data is already human)
  cat("  No mouse matches. Trying direct HGNC matching...\n")
  expression_values[, HGNC := gene_name]
  coral_input <- merge(expression_values, kinases[, .(HGNC)], by = "HGNC", all = FALSE)
  cat("  Matched", nrow(coral_input), "kinases\n")
}

# Remove duplicates (keep first)
coral_input <- coral_input[!duplicated(HGNC)]

# Calculate summary statistics
cat("\n=== Expression Statistics ===\n")
cat("  Kinases matched:", nrow(coral_input), "/", nrow(kinases[!is.na(HGNC) & HGNC != ""]), "\n")
cat("  Value range:", round(min(coral_input$value, na.rm = TRUE), 3), "to",
    round(max(coral_input$value, na.rm = TRUE), 3), "\n")
cat("  Mean:", round(mean(coral_input$value, na.rm = TRUE), 3), "\n")
cat("  Median:", round(median(coral_input$value, na.rm = TRUE), 3), "\n")

# Show top up/down regulated kinases
if (grepl("log2|FC|fold", comparison_name, ignore.case = TRUE)) {
  setorder(coral_input, -value)

  cat("\n=== Top 10 Upregulated Kinases ===\n")
  top_up <- head(coral_input[value > 0], 10)
  if (nrow(top_up) > 0) {
    print(top_up[, .(HGNC, value = round(value, 3))])
  }

  cat("\n=== Top 10 Downregulated Kinases ===\n")
  top_down <- head(coral_input[value < 0][order(value)], 10)
  if (nrow(top_down) > 0) {
    print(top_down[, .(HGNC, value = round(value, 3))])
  }
}

# Create output directory
coral_output_dir <- file.path(paths$output_dir, "coral")
if (!dir.exists(coral_output_dir)) {
  dir.create(coral_output_dir, recursive = TRUE)
}

# Generate CORAL input file (for reference / web app use)
output_file <- file.path(coral_output_dir, paste0(opt$`output-prefix`, "_", comparison_name, ".txt"))
coral_output <- coral_input[, .(HGNC, value)]
fwrite(coral_output, output_file, sep = "\t", col.names = FALSE)
cat("\n✓ CORAL input saved to:", output_file, "\n")

# Save metadata
metadata_file <- file.path(coral_output_dir, paste0(opt$`output-prefix`, "_", comparison_name, "_metadata.csv"))
coral_metadata <- merge(coral_input, kinases[, .(HGNC, Mouse_Symbol, Group, Family, SubFamily, gene_name)],
                        by = "HGNC", all.x = TRUE)
fwrite(coral_metadata, metadata_file)
cat("✓ Metadata saved to:", metadata_file, "\n")

# === Generate SVG using CORAL rendering ===
cat("\n=== Generating Kinome Tree SVG ===\n")

# Prepare recolor dataframe for CORAL (id, value format)
recolordf <- data.frame(
  kinase = coral_input$HGNC,
  userinfo = coral_input$value,
  stringsAsFactors = FALSE
)

# Convert HGNC to CORAL IDs
coral_df_copy <- as.data.frame(coral_df)

# Match our data to CORAL kinases
matched_idx <- match(recolordf$kinase, coral_df_copy$id.HGNC)
valid_matches <- !is.na(matched_idx)
cat("  Matched to CORAL tree:", sum(valid_matches), "/", nrow(recolordf), "kinases\n")

# Create recolor dataframe with CORAL IDs
recolordf_coral <- data.frame(
  kinase = coral_df_copy$id.coral[matched_idx[valid_matches]],
  userinfo = recolordf$userinfo[valid_matches],
  stringsAsFactors = FALSE
)

# Set color scale limits (from config or CLI)
heatrange <- c(-fc_limit, fc_limit)

# Define diverging color palette from config
colors <- colorRampPalette(c(CORAL_COLOR_LOW, CORAL_COLOR_MID, CORAL_COLOR_HIGH))(100)

# Apply colors to nodes using CORAL's color.by.value function
color_result <- color.by.value(
  df = coral_df_copy,
  recolordf = recolordf_coral,
  colors = colors,
  heatrange = heatrange,
  bg.col = CORAL_COLOR_BG
)

# Update coral dataframe with new colors
coral_df_copy$node.col <- color_result[[1]]
coral_df_copy$node.val.col <- color_result[[2]]

# Set node radius for matched kinases (from config)
coral_df_copy$node.radius <- CORAL_NODE_RADIUS
coral_df_copy$node.radius[coral_df_copy$id.coral %in% recolordf_coral$kinase] <- CORAL_NODE_RADIUS_MATCHED

# Build SVG info structure
svginfo <- list()

# Read base SVG header
base_svg <- file.path(CORAL_PATH, "Data", "basetree.svg")
if (file.exists(base_svg)) {
  svglines <- readLines(base_svg)
  # Extract header (first few lines until defs close)
  header_end <- grep("</defs>", svglines)[1]
  if (!is.na(header_end)) {
    svginfo$header <- svglines[1:header_end]
  } else {
    svginfo$header <- "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">"
  }
} else {
  svginfo$header <- "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" viewBox=\"0 0 850 600\">"
}

# Set title
fig_title <- if (!is.null(opt$title)) opt$title else comparison_name
svginfo$title <- fig_title

# Build color legend using config colors
legend_text <- sprintf("log2FC: [%.1f, %.1f]", -fc_limit, fc_limit)
svginfo$legend <- paste0(
  "<rect x=\"10\" y=\"10\" width=\"120\" height=\"20\" fill=\"url(#legendGradient)\" stroke=\"black\" stroke-width=\"0.5\"/>",
  "<defs><linearGradient id=\"legendGradient\">",
  "<stop offset=\"0%\" stop-color=\"", CORAL_COLOR_LOW, "\"/>",
  "<stop offset=\"50%\" stop-color=\"", CORAL_COLOR_MID, "\"/>",
  "<stop offset=\"100%\" stop-color=\"", CORAL_COLOR_HIGH, "\"/>",
  "</linearGradient></defs>",
  "<text x=\"65\" y=\"45\" text-anchor=\"middle\" font-size=\"8\" font-family=\"sans-serif\">", legend_text, "</text>"
)

# Add branch order and node order columns if not present
if (!"branchorder" %in% names(coral_df_copy)) {
  coral_df_copy$branchorder <- 1:nrow(coral_df_copy)
}
if (!"nodeorder" %in% names(coral_df_copy)) {
  coral_df_copy$nodeorder <- 1:nrow(coral_df_copy)
}
if (!"node.selected" %in% names(coral_df_copy)) {
  coral_df_copy$node.selected <- -1
}
if (!"node.opacity" %in% names(coral_df_copy)) {
  coral_df_copy$node.opacity <- 1
}

svginfo$dataframe <- coral_df_copy

# Extract group labels from base SVG
if (file.exists(base_svg)) {
  groups_start <- grep("<g id=\"GROUPS\">", svglines)
  groups_end <- grep("</g>", svglines)
  groups_end <- groups_end[groups_end > groups_start[length(groups_start)]][1]
  if (!is.na(groups_start[1]) && !is.na(groups_end)) {
    svginfo$groups <- svglines[(groups_start[length(groups_start)] + 1):(groups_end - 1)]
  } else {
    svginfo$groups <- character(0)
  }
} else {
  svginfo$groups <- character(0)
}

# Write SVG file (if enabled in config)
svg_file <- file.path(coral_output_dir, paste0(opt$`output-prefix`, "_", comparison_name, ".svg"))

if (CORAL_OUTPUT_SVG) {
  writekinasetree(
    svginfo = svginfo,
    destination = svg_file,
    font = "'Roboto', sans-serif",
    labelselect = "HGNC",
    groupcolor = "black"
  )
  cat("✓ SVG saved to:", svg_file, "\n")
} else {
  # Still write SVG temporarily for PDF conversion
  writekinasetree(
    svginfo = svginfo,
    destination = svg_file,
    font = "'Roboto', sans-serif",
    labelselect = "HGNC",
    groupcolor = "black"
  )
}

# Convert to PDF if enabled in config and rsvg is available
if (CORAL_OUTPUT_PDF) {
  pdf_file <- file.path(coral_output_dir, paste0(opt$`output-prefix`, "_", comparison_name, ".pdf"))
  if (requireNamespace("rsvg", quietly = TRUE)) {
    tryCatch({
      rsvg::rsvg_pdf(svg_file, pdf_file)
      cat("✓ PDF saved to:", pdf_file, "\n")
    }, error = function(e) {
      cat("  Note: PDF conversion failed (", e$message, ")\n")
    })
  } else {
    cat("  Note: Install 'rsvg' package for PDF output\n")
  }
}

# Clean up SVG if only PDF was requested
if (!CORAL_OUTPUT_SVG && CORAL_OUTPUT_PDF && file.exists(svg_file)) {
  unlink(svg_file)
}

cat("\nCORAL tree generation complete.\n")
cat("Configuration loaded from: genesets_config.yaml (coral: section)\n")
