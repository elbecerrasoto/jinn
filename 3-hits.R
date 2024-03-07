#!/usr/bin/Rscript
library(tidyverse)
library(stringr)
library(segmenTools)

# Out table ---------------------------------------------------------------


# 1. genome
# 2. pid
# 3. q_alias
# 4. query
# 5. order
# 6. start
# 7. end
# 8. contig
# 9. strand
# 10. pid_txt
# 11. interpro

# Globals -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# Tables
GFF <- "results/genomes/GCF_000699465.1/GCF_000699465.1.gff"
BLAST <- "blasts_filtered.tsv"
ISCAN <- "results/iscan.tsv"

QUERIES <- c("WP_003243987.1", "WP_003243213.1")
QUERIES_ALIAS <- c("YwqJ", "YwqL") |>
  `names<-`(QUERIES)

get_genome <- function(path) {
  str_replace(path, ".*(GC[FA]_[0-9]+\\.[0-9])\\.gff+", "\\1")
}

GENOME <- get_genome(GFF)


# Read Tabs --------------------------------------------------------------------

# The mapping data is generated in previous file

blast <- read_tsv(BLAST)
iscan <- read_tsv(ISCAN, na = c("NA", "-", ""))

gff <- segmenTools::gff2tab(GFF) |>
  tibble() |>
  filter(feature == "CDS") |> # only CDS
  select_if({
    \(x) !(all(is.na(x)) | all(x == ""))
  }) # exclude empty cols

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
