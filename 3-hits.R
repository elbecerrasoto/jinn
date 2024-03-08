#!/usr/bin/Rscript

library(glue)
library(stringr)
library(tidyverse)
library(segmenTools)

# args <- commandArgs(trailingOnly = TRUE)

# Globals -----------------------------------------------------------------

# Input
GFF <- "results/genomes/GCF_000699465.1/GCF_000699465.1.gff"
MAPPINGS <- "mappings_filtered.tsv"


get_genome <- function(path) {
  str_replace(path, ".*(GC[FA]_[0-9]+\\.[0-9])\\.gff+", "\\1")
}

GENOME <- get_genome(GFF)

# Output
BASE <- "hits"
OUT <- glue("{GENOME}_{BASE}.tsv")

OUT_COLS <- c(
  "genome",
  "pid",
  "gene",
  "q_alias",
  "order",
  "start",
  "end",
  "contig",
  "strand",
  "query",
  "domains",
  "product"
)


# Read the Data -----------------------------------------------------------


mappings <- read_tsv(MAPPINGS)

gff <- segmenTools::gff2tab(GFF) |>
  tibble() |>
  filter(feature == "CDS") |> # only CDS
  select_if({
    \(x) !(all(is.na(x)) | all(x == ""))
  }) # exclude empty cols

# Remove pseudogenes
if ("pseudo" %in% names(gff)) {
  gff <- gff |>
    filter(is.na(pseudo))
}


# definition of neighbor
# same contig, order by start position
gff <- gff |>
  group_by(seqname) |>
  arrange(start) |>
  mutate(order = seq_along(start)) |>
  relocate(order)


# Policy ------------------------------------------------------------------


hits <- inner_join(mappings, gff, join_by(pid == protein_id)) |>
  mutate(genome = GENOME, contig = seqname) |>
  select(all_of(OUT_COLS))


hits |>
  write_tsv(OUT)
