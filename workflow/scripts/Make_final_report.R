#!/usr/bin/env Rscript

# Chris Sansam
# version 01
# May 23, 2022

# load libraries
library(cowplot)
library(knitr)
library(rmarkdown)
library(magick)

# get arguments from Snakemake S4 object
samples <- snakemake@params[["samples"]]
rmarkdown::render("resources/Report.Rmd",output_dir="results/",params=list(samples=samples))
