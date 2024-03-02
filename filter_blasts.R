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
QUERIES_ALIAS <- list("YwqJ", "YwqL") |>
  `names<-`(QUERIES)

FILTER_DOMAINS <- list(c("IPR027797", "IPR025968"), c("IPR007581")) |>
  `names<-`(QUERIES)


blasts <- read_tsv(BLASTS)
iscan <- read_tsv(ISCAN)



# Mappings
# qs -> pids -> doms
# Use summarise + list to get a different shape
# Map each query to pid (is a 1-to-many mapping)
q_pids <- blasts |>
  group_by(qseqid) |>
  reframe(pid = unique(sseqid)) |>
  rename(query = qseqid)


# Map each pid to domains (is a 1-to-many mapping)
pid_doms <- iscan |>
  group_by(protein) |>
  reframe(domain = unique(interpro[interpro != "-"])) |>
  rename(pid = protein)


switch_vect <- function(v_chr, switch_list) {
  f <- function(v_chr) {
    do.call(
      switch,
      c(v_chr, switch_list, NA)
    )
  }

  map_chr(v_chr, f)
}

q_pids <- q_pids |>
  mutate(q_alias = switch_vect(query, QUERIES_ALIAS)) |>
  relocate(q_alias)


x <- left_join(q_pids, pid_doms)

view(x)
