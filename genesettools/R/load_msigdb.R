#' Load MSigDB GMT files (simple loader)
#'
#' Downloads a GMT file for a given collection (if not cached) and returns a named
#' list of gene vectors.
#'
#' @param collection character vector of MSigDB collections (e.g. "H", "C2")
#' @param species character, species code suffix used in filenames (e.g. "Hs", "Mm").
#' @param cache_dir character path to store downloaded GMTs. Defaults to "inst/gmt".
#' @return named list where names are gene set names and values are character vectors
#' @export
load_msigdb <- function(collection = "H", species = "Hs", cache_dir = "inst/gmt") {
  stopifnot(is.character(collection))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  all_sets <- list()
  for (c in collection) {
    # construct filename (example format used by msigdb: m<collection>.all.v2024.1.Hs.symbols.gmt)
    fname <- file.path(cache_dir, paste0("m", tolower(c), ".all.v2024.1.", species, ".symbols.gmt"))
    if (!file.exists(fname)) {
      # try to download from GSEA MSigDB (user may need credential or manual placement)
      url <- sprintf("https://data.gsea-msigdb.org/msigdb/release/2024.1/%s", basename(fname))
      try(
        httr::GET(url, httr::write_disk(fname, overwrite = TRUE), httr::timeout(60)),
        silent = TRUE
      )
    }
    if (file.exists(fname)) {
      lines <- readr::read_lines(fname)
      for (ln in lines) {
        parts <- strsplit(ln, "\t")[[1]]
        if (length(parts) >= 3) {
          name <- parts[1]
          genes <- parts[-c(1,2)]
          all_sets[[name]] <- genes
        }
      }
    }
  }
  return(all_sets)
}
