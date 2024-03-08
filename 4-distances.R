#!/usr/bin/Rscript

library(glue)
library(stringr)
library(tidyverse)
library(segmenTools)

# args <- commandArgs(trailingOnly = TRUE)

# Globals -----------------------------------------------------------------

# Input
HITS <- "GCF_000699465.1_Test0 .tsv"

get_genome <- function(path) {
  str_extract(path, "GC[FA]_[0-9]+\\.[0-9]")
}

GENOME <- get_genome(HITS)

GENOME

# Output
BASE <- "distance"
OUT <- glue("{GENOME}_{BASE}.tsv")


# Read the Data -----------------------------------------------------------

hits <- read_tsv(HITS)

# Carterian product
# Copy from in_bet.R