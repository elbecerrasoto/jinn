#!/usr/bin/Rscript
# TODO: change name to pairs

library(glue)
library(stringr)
library(tidyverse)

# args <- commandArgs(trailingOnly = TRUE)

# Globals -----------------------------------------------------------------

# Two genes
TARGETS <- c("WP_003243987.1", "WP_003243213.1")

# A distance is between 2 things
stopifnot(length(TARGETS) == 2)

# Input
HITS <- "GCF_000699465.1_Test0.tsv"

get_genome <- function(path) {
  str_extract(path, "GC[FA]_[0-9]+\\.[0-9]")
}

GENOME <- get_genome(HITS)

# Output
BASE <- "distance"
OUT <- glue("{GENOME}_{BASE}.tsv")

graceful_exit <- function(assertion) {
  if (!assertion) {
    write_tsv(tibble(), OUT)
    quit(status = 0)
  }
}

# genome
# distance
# genes_inbet
# contig
# strand
# q1 p1 o1 s1 e1 # sort
# q2 p2 o2 s2 e2 # sort


# Read the Data -----------------------------------------------------------

hits <- read_tsv(HITS)
stopifnot(nrow(hits) > 0)

hits <- hits |>
  filter(query %in% TARGETS)

# Query to factor
# Used on counting pairs
hits <- hits |>
  mutate(
    query = as_factor(query),
    query = `levels<-`(query, TARGETS)
  )

CONTIGS <- hits |>
  pull(contig) |>
  unique()

# Calculate the number of distance operations
count_pairs <- function(hits) {
  n <- 0

  contig_query <- hits |>
    count(contig, query, .drop = FALSE)

  for (contig in CONTIGS) {
    x <- contig_query |>
      filter(contig == {{ contig }}) |>
      pull(n) |>
      reduce(`*`)
    n <- n + x
  }
  n
}


# for(contig in contigs) {
#  for(query in query)
#
# }
# for (i in seq_len(nrow(genes_q1))) {
#   for (j in seq_len(nrow(genes_q2))) {
#     n <- n + 1
#     gene1 <- genes_q1[i, ]
#     gene2 <- genes_q2[j, ]
#     lresults[[n]] <- calc(gene1, gene2)
#   }
# }
#
