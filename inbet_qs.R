#!/usr/bin/Rscript
library(tidyverse)
library(stringr)
library(segmenTools)

args <- commandArgs(trailingOnly = TRUE)

# Tables
GFF <- "results/genomes/GCF_000699465.1/GCF_000699465.1.gff"
BLAST <- "results/blasts.tsv" 
ISCAN <- "results/iscan.tsv"
  
# Queries
# YwqJ: WP_003243987.1
# YwqL: WP_003243213.1

QUERIES <- c("WP_003243987.1", "WP_003243213.1")
QUERIES_ALIAS <- c("YwqJ", "YwqL")
names(QUERIES_ALIAS) <- QUERIES  

# Filter Domains
# J: IPR027797, IPR025968
# J: Pre-toxin, YwqJ-like
# L: IPR007581
# L: Endonuclease V
  
FILTER_DOMAINS <- list(c("IPR027797", "IPR025968"), c("IPR007581"))
names(FILTER_DOMAINS) <- QUERIES

get_genome <- function(path) {
  str_replace(path, ".*(GC[FA]_[0-9]+\\.[0-9])\\.gff+", "\\1")
}

GENOME <- get_genome(GFF)

OUT <- paste0(GENOME, "_", QUERIES_ALIAS[QUERIES[1]], "_", QUERIES_ALIAS[QUERIES[2]], ".tsv")

graceful_exit <- function() {
  write_tsv(tibble(), OUT)
  quit(status = 0)
}

Rgff::check_gff(BLAST)

blast <- read_tsv(BLAST)
iscan <- read_tsv(ISCAN)

genome <- get_genome(GFF)

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


# Select query 1 and query 2 
q1_blast <- blast |> filter(qseqid == QUERIES[1])
q2_blast <- blast |> filter(qseqid == QUERIES[2])


mtcars %>%
  group_by(cyl) %>%
  summarise(mean = mean(disp), n = n())

filter_PIDs_by_doms <- function(PID, query) {
  
}

pid <- iscan$protein[1]
pid
names(iscan)

iscan |> group_by(protein) 


# query-subject mappings
PID_q1 <- q1_blast |> pull(sseqid) |> unique()
PID_q2 <- q2_blast |> pull(sseqid) |> unique()









q1_gff <- gff |>
  filter(protein_id %in% q1_map_PID)
q2_gff <- gff |>
  filter(protein_id %in% q2_map_PID)

q1_gff
q2_gff

calc <- function(gene1, gene2) {
  contig1 <- gene1$seqname
  contig2 <- gene2$seqname

  if (identical(contig1, contig2)) {
    contig <- contigs |>
      filter(seqname == contig1)

    contig_gff <- gff |>
      filter(seqname == contig1)

    s2 <- max(gene1$start, gene2$start)
    e1 <- min(gene1$end, gene2$end)
    distance <- s2 - e1
    gene_count <- abs(gene1$order - gene2$order)
    same_strand <- identical(gene1$strand, gene2$strand)
    gene1_first <- gene1$order < gene2$order
    return(list(
      genome = GENOME,
      gene1 = gene1$protein_id,
      gene2 = gene2$protein_id,
      distance = distance,
      gene_count = gene_count,
      same_strand = same_strand,
      gene1_order = gene1$order,
      gene2_order = gene2$order,
      gene1_first = gene1_first,
      contig = contig1,
      contig_ostart = min(contig_gff$order),
      contig_oend = max(contig_gff$order),
      contig_start = contig$start,
      contig_end = contig$end,
      gene1_start = gene1$start,
      gene1_end = gene1$end,
      gene2_start = gene2$start,
      gene2_end = gene2$end,
      refseq = REFSEQ
    ))
  } else {
    return(list())
  }
}

lresults <- vector("list", length = nrow(genes_q1) * nrow(genes_q2))
n <- 0

for (i in seq_len(nrow(genes_q1))) {
  for (j in seq_len(nrow(genes_q2))) {
    n <- n + 1
    gene1 <- genes_q1[i, ]
    gene2 <- genes_q2[j, ]
    lresults[[n]] <- calc(gene1, gene2)
  }
}

results <- bind_rows(map(lresults, as_tibble))
write_tsv(results, OUT, col_names = FALSE)
