#!/usr/bin/Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

# Filter blasts by domain
BLASTS <- "results/blasts.tsv"
ISCAN <- "results/iscan.tsv"

# TODO: the filtering data into table
# Something like
# q_alias, query, IP_filter, IP_txt 
QUERIES <- c("WP_003243987.1", "WP_003243213.1")
QUERIES_ALIAS <- c("YwqJ", "YwqL")
names(QUERIES_ALIAS) <- QUERIES
FILTER_DOMAINS <- list(c("IPR027797", "IPR025968"), c("IPR007581"))
names(FILTER_DOMAINS) <- QUERIES


blasts <- read_tsv(BLASTS)
iscan <- read_tsv(ISCAN)

# Mappings
# qs -> pids -> doms
# Map each query to pid (is a 1-to-many mapping)
q_pids <- blast |>
  group_by(qseqid) |>
  summarise(pids = list(unique(sseqid)))
# Map each pid to domains (is a 1-to-many mapping)
pid_doms <- iscan |>
  group_by(protein) |>
  summarise(doms = list(unique(interpro[interpro != "-"])))


pid_doms
