#!/usr/bin/env Rscript
#use: generate sequence similarity matrix per cluster: R pdb_id cluster_0 --no-save < sequence_muscle.R
#use: R output_dir --no-save < A3-assam_sequence.R

options(warn = -1)
library("Biostrings")
query = readAAStringSet("query.fasta")
hits = readAAStringSet("hits.fasta")
for (i in seq_along(query)){
    for (j in seq_along(hits)) {
    write(paste(names(query[i]), names(hits[j]), pid(pairwiseAlignment( query[i], hits[j]))), file="sequence_identity.txt", append=TRUE)
    }}