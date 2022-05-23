#!/usr/bin/env Rscript

# Chris Sansam
# version 01
# May 23, 2022

# load libraries
library(cowplot)
library(knitr)
library(rmarkdown)

# get arguments from Snakemake S4 object
samples <- snakemake@params[["samples"]]
samples <- strsplit(snakemake@params[["samples"]], split = ",")[[1]]
rmarkdown::render("resources/Report.Rmd",output_dir="results/",params=list(samples=samples))
