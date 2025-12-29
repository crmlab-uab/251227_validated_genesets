#!/usr/bin/env Rscript
# Simple pipeline runner for kinases workflow (fetch -> generate -> validate)
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

steps <- c('fetch','generate','validate')
option_list <- list(
  make_option(c('-s','--steps'), type='character', default=paste(steps,collapse=','), help='Comma-separated steps: fetch,generate,validate')
)
opt <- parse_args(OptionParser(option_list=option_list))
run_steps <- strsplit(opt$steps,',')[[1]]

if ('fetch' %in% run_steps) {
  message('Running fetcher index...')
  source('scripts/kinases/fetchers_index.R')
}
if ('generate' %in% run_steps) {
  message('Generating kinases from BioMart (GO + InterPro)')
  system('Rscript scripts/kinases/generate_kinases_from_biomart.R', intern=FALSE)
}
if ('validate' %in% run_steps) {
  message('Running validations index...')
  source('scripts/kinases/validations_index.R')
}

message('Pipeline completed.')
