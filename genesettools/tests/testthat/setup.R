## Test setup: load package in-place so tests can call package functions
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cran.rstudio.com")
}
devtools::load_all()
