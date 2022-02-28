library(dplyr)
library(fuzzyjoin)
library(readxl)
library(IRanges)
library(stringr)
library(tidyr)
library(tibble)
library(data.table)
library(foreach)
library(doParallel)

samples <- snakemake@params[["samples"]]
all_types <- snakemake@params[["all_types"]]


samples
all_types



