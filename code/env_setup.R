library(renv)

renv::init()

renv::install("tidyverse")
renv::install("bioc::GenomicRanges")
renv::install("Matrix")
renv::install("furrr")

renv::snapshot()
