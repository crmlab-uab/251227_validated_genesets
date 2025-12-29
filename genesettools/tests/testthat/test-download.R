library(testthat)

test_that("load_msigdb reads a local GMT file correctly", {
  tmp <- tempfile(pattern = "gmtdir")
  dir.create(tmp)
  # construct expected filename used by the loader
  fname <- file.path(tmp, paste0("m", tolower("H"), ".all.v2024.1.", "Hs", ".symbols.gmt"))
  sample <- c(
    "GS1\tdesc\tGENE1\tGENE2",
    "GS2\tdesc\tGENE3"
  )
  writeLines(sample, fname)

  res <- genesettools::load_msigdb(collection = "H", species = "Hs", cache_dir = tmp)
  expect_true("GS1" %in% names(res))
  expect_equal(res$GS1, c("GENE1", "GENE2"))
  expect_true("GS2" %in% names(res))
  expect_equal(res$GS2, c("GENE3"))
})

test_that("load_msigdb handles missing file (download failure) gracefully", {
  tmp <- tempfile(pattern = "gmtdir")
  dir.create(tmp)
  # no file present; function should attempt download but ultimately return empty list
  res <- genesettools::load_msigdb(collection = "ZNOTEXIST", species = "Hs", cache_dir = tmp)
  expect_type(res, "list")
  expect_equal(length(res), 0)
})
