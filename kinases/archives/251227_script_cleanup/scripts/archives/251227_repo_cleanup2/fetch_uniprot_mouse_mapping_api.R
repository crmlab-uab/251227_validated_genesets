# Fetch UniProt mouse gene name to UniProtKB mapping using UniProt REST API
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)

library(httr)
library(data.table)

# Set up API parameters
base_url <- "https://rest.uniprot.org/uniprotkb/search"
query <- "(organism_id:10090) AND (gene:* )"
fields <- "accession,gene_names"
size <- 500

# Prepare output file
output_file <- "../uniprot_mouse_idmapping_selected.tab"

cat("Fetching UniProt mouse gene name to UniProtKB mapping via API...\n")

# Paginate through results and write in chunks
offset <- 0
first_chunk <- TRUE
repeat {
  url <- paste0(base_url,
                "?query=", URLencode(query),
                "&fields=", fields,
                "&format=tsv",
                "&size=", size,
                "&offset=", offset)
  resp <- GET(url)
  if (status_code(resp) != 200) break
  dat <- fread(text=content(resp, as="text", encoding="UTF-8"), sep="\t", header=TRUE)
  if (nrow(dat) == 0) break
  setnames(dat, c("Entry", "Gene Names"), c("UniProtKB", "Gene_Name"))
  # Write chunk to file (append after first chunk)
  fwrite(dat, output_file, sep="\t", append=!first_chunk, col.names=first_chunk)
  first_chunk <- FALSE
  offset <- offset + size
  cat("Fetched ", offset, " records...\n")
  if (nrow(dat) < size) break
}
if (file.exists(output_file) && file.info(output_file)$size > 0) {
  cat("✓ Partial mapping file saved to ", output_file, "\n")
} else {
  cat("No mapping data retrieved.\n")
}
