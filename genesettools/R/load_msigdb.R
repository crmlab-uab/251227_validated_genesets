#' Load MSigDB GMT files (simple loader)
#'
#' Download (if missing) and load MSigDB GMT files into R as a named list of
#' character vectors. This is a lightweight helper intended for small-scale
#' integration tests and pipelines.
#'
#' @param collection Character vector of MSigDB collections (e.g. "H", "C2").
#' @param species Character, species code suffix used in filenames (e.g. "Hs", "Mm").
#' @param cache_dir Character path to store downloaded GMTs. Defaults to "inst/gmt".
#'
#' @details
#' The function constructs a filename using the convention used by MSigDB
#' release GMT files (example: `mh.all.v2024.1.Hs.symbols.gmt`). If the file is
#' missing it will attempt an HTTP GET to the canonical GSEA URL and write the
#' file to `cache_dir`. Download attempts are wrapped in `try(..., silent=TRUE)`
#' so callers can rely on graceful failure (an empty list) when network access is
#' unavailable.
#'
#' @return A named list where each element is a character vector of gene symbols.
#' @examples
#' \dontrun{
#'   # load local cached collection (if present)
#'   load_msigdb("H", "Hs", cache_dir = tempdir())
#' }
#' @export
#' @import httr
#' @importFrom readr read_lines
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
