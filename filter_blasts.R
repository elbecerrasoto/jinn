#!/usr/bin/Rscript
library(tidyverse)

OUT <- "blasts_filtered.tsv"

args <- commandArgs(trailingOnly = TRUE)

# Filter blasts by domain
BLASTS <- "results/blasts.tsv"
ISCAN <- "results/iscan.tsv"

# TODO: the filtering data into table
# Something like
# q_alias, query, IP_filter, IP_txt

# YwqJ
# LXG - IPR006829 TOSEARCH
# PT-TG - IPR027797
# YwqJ-like - IPR025968 TOSEARCH

# YwqL
# Endonuclease-V - IPR007581 TOSEARCH


QUERIES <- c("WP_003243987.1", "WP_003243213.1")
QUERIES_ALIAS <- list("YwqJ", "YwqL") |>
  `names<-`(QUERIES)

FILTER_DOMAINS <- list(c("IPR006829", "IPR025968"), c("IPR007581")) |>
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


# Add the aliases
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

# Join both tables
q_pids_domains <- left_join(q_pids, pid_doms) |>
  group_by(q_alias, query, pid) |>
  summarise(domains = list(domain))

# Filter by domains
filter_by_domain <- function(pids, qs, doms, filtering_doms) {
  names(doms) <- pids
  names(qs) <- pids
  return_vec <- vector(mode = "logical", length = length(pids))

  for (i in seq_along(pids)) {
    pid <- pids[i]
    return_current <- FALSE

    for (query in names(filtering_doms)) {
      if (qs[pid] == query) {
        return_current <- all(filtering_doms[[query]] %in% doms[[pid]])
        break()
      }
    }
    return_vec[i] <- return_current
  }
  return_vec
}

q_pids_domains_filtered <- q_pids_domains |>
  filter(filter_by_domain(pid, query, domains, FILTER_DOMAINS))

# Filter blasts
blasts_filtered <- semi_join(blasts,
  q_pids_domains_filtered,
  by = join_by(sseqid == pid)
)

blasts_filtered |>
  write_tsv(OUT)
