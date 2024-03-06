#!/usr/bin/Rscript
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)


# Globals -----------------------------------------------------------------


OUT_BLAST <- "blasts_filtered.tsv"
OUT_MAPPINGS <- "mappings_filtered.tsv"

# Filter blasts by domain
BLASTS <- "results/blasts.tsv"
MAPPINGS <- "mappings.tsv"

# YwqJ
# LXG - IPR006829 TO_SEARCH
# PT-TG - IPR027797
# YwqJ-like - IPR025968 TO_SEARCH

# YwqL
# Endonuclease-V - IPR007581 TO_SEARCH

FILTER_DOMAINS <- list(
  WP_003243987.1 = c("IPR006829", "IPR025968"),
  WP_003243213.1 = c("IPR007581")
)


# Reading Data ------------------------------------------------------------


blasts <- read_tsv(BLASTS) # only to filter it
mappings <- read_tsv(MAPPINGS) # operate on this table


# Code --------------------------------------------------------------------

# vectorized boolean function
# to be used on dplyr::filter steps
# tbl var | l = n | query: chr
# tbl var | l = n | domains: list chr
# domains_to_check: list chr
# OUT | l = n | lgl


check_domains <- function(query, domains, domains_to_check) {
  n <- length(domains)
  filter_lgl <- rep(TRUE, n) # default value

  if (length(domains_to_check) == 0) {
    return(filter_lgl)
  } else {
    for (i in seq_along(domains)) {
      for (Q in names(domains_to_check)) {
        if (query[i] == Q) {
          filter_lgl[i] <- domains_to_check[[Q]] %in% domains[[i]] |>
            all()

          break()
        }
      }
    }
  }
  filter_lgl
}
