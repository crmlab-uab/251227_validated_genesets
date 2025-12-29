# Map human UniProt IDs to mouse UniProt IDs for kinase validation
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)

library(data.table)
library(UniProt.ws)
library(parallel)

# Load kinase list
kinases <- fread('genesets/curated/kinases/201006_composite_kinases_curated.csv', header=TRUE)



# Initialize UniProt.ws for both human and mouse
up_human <- UniProt.ws(taxId=9606)
# up_mouse will be initialized in each parallel worker

# Prepare output columns
kinases[, UniProt_Mouse := NA_character_]
kinases[, Mouse_UniProt_MatchType := NA_character_]

cat("Querying UniProt.ws for mouse UniProt IDs by gene symbol (Gene_Name) in parallel...\n")

# Function to query UniProt for a single kinase row
query_mouse_uniprot <- function(idx) {
  mouse_symbol <- kinases$Mouse_Symbol[idx]
  human_symbol <- kinases$Human_Symbol[idx]
  up_mouse <- UniProt.ws(taxId=10090) # Each process needs its own instance
  # Try mouse gene symbol first
  mouse_res <- tryCatch({
    select(up_mouse,
      keys = mouse_symbol,
      columns = c("Entry", "Entry.Name", "Protein.names"),
      keytype = "Gene_Name"
    )
  }, error = function(e) NULL)
  if (!is.null(mouse_res) && nrow(mouse_res) > 0) {
    mouse_row <- mouse_res[grepl("_MOUSE$", mouse_res$Entry.Name), ]
    if (nrow(mouse_row) > 0) {
      return(list(UniProt_Mouse=mouse_row$Entry[1], MatchType="mouse_symbol"))
    }
  }
  # If not found, try human gene symbol in mouse
  if (!is.na(human_symbol) && human_symbol != "") {
    mouse_res2 <- tryCatch({
      select(up_mouse,
        keys = human_symbol,
        columns = c("Entry", "Entry.Name", "Protein.names"),
        keytype = "Gene_Name"
      )
    }, error = function(e) NULL)
    if (!is.null(mouse_res2) && nrow(mouse_res2) > 0) {
      mouse_row2 <- mouse_res2[grepl("_MOUSE$", mouse_res2$Entry.Name), ]
      if (nrow(mouse_row2) > 0) {
        return(list(UniProt_Mouse=mouse_row2$Entry[1], MatchType="human_symbol_in_mouse"))
      }
    }
  }
  return(list(UniProt_Mouse=NA_character_, MatchType="not_found"))
}

# Use mclapply for parallel processing (Linux/macOS only)
num_cores <- min(4, parallel::detectCores())
results <- parallel::mclapply(1:nrow(kinases), query_mouse_uniprot, mc.cores=num_cores)

# Assign results back to kinases table
kinases$UniProt_Mouse <- sapply(results, function(x) x$UniProt_Mouse)
kinases$Mouse_UniProt_MatchType <- sapply(results, function(x) x$MatchType)


# Save updated kinase list with today's date
output_file <- "251227_curated_kinases.csv"
fwrite(kinases, output_file)
cat(paste0("✓ Mouse UniProt mapping complete. Output: ", output_file, "\n"))

# To run in the background, use:
#   nohup Rscript map_human_to_mouse_uniprot.R &
