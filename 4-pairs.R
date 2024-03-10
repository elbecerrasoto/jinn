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


# Calc --------------------------------------------------------------------


# genome
# distance
# genes_inbet
# contig
# strand
# q1 p1 o1 s1 e1 # sort
# q2 p2 o2 s2 e2 # sort


calc <- function(gene1, gene2) {
  stopifnot(is_tibble(gene1), is_tibble(gene2))
  stopifnot(gene1$order != gene2$order)
  stopifnot(gene1$contig == gene2$contig)

  first <- ifelse(gene1$order < gene2$order, gene1, gene2)
  second <- ifelse(gene1$order > gene2$order, gene1, gene2)

  print(glue("first \n{first}"))
  print(glue("second \n{second}"))

  print(glue("names {names(first)}"))

  distance <- second$start - first$end
  genes_inbet <- second$order - first$order

  tibble(
    genome = GENOME,
    distance = distance,
    genes_inbet = genes_inbet,
    contig = first$contig,
    query_1 = first$query,
    pid_1 = first$order,
    order_1 = first$order,
    start_1 = first$start,
    end_1 = first$end,
    strand_1 = first$strand,
    locustag_1 = first$locus_tag,
    query_2 = second$query,
    pid_2 = second$order,
    order_2 = second$order,
    start_2 = second$start,
    end_2 = second$end,
    strand_2 = second$strand,
    locustag_2 = second$locus_tag,
  )
}


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

N_PAIRS <- count_pairs(hits)
results <- vector(mode = "list", length = N_PAIRS)

i <- 0
for (contig in CONTIGS) {
  same_contig <- hits |>
    filter(contig == {{ contig }}) %>%
    split(.$query)

  print(glue("same contig {same_contig}"))

  genes_q1 <- same_contig[[1]]
  genes_q2 <- same_contig[[2]]

  for (g1_idx in seq_len(nrow(genes_q1))) {
    for (g2_idx in seq_len(nrow(genes_q2))) {
      i <- i + 1
      gene1 <- genes_q1[g1_idx, ]
      gene2 <- genes_q2[g2_idx, ]
      results[[i]] <- calc(gene1, gene2)
    }
  }
}

stopifnot(N_PAIRS == i)

hits

calc(hits[1, ], hits[2, ])

hits[1, ]
results
