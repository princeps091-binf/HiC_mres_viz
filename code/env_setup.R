library(renv)

renv::init()

renv::install("tidyverse")
renv::install("bioc::GenomicRanges")
renv::install("Matrix")
renv::install("furrr")
renv::install("data.tree")
renv::install("viridis")
renv::install("glue")

renv::install("caret")
renv::install("igraph")

renv::snapshot()
