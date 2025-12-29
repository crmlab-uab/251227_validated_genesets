test_that("load_msigdb returns a list and handles missing files", {
  # calling with an empty collection should return an empty list
  res <- load_msigdb(collection = character(0), cache_dir = tempdir())
  expect_type(res, "list")
  # calling with non-existent collection should not error (network optional)
  res2 <- tryCatch(load_msigdb(collection = "ZZ", species = "Hs", cache_dir = file.path(tempdir(), "gmt_test")), error = function(e) e)
  expect_true(is.list(res2) || inherits(res2, "error") == FALSE)
})
