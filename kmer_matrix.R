# load libraries 
library(Biostrings)

setwd("C:/Users/danie/OneDrive - Colostate/NIFA_PROJECT/Obj1/analysis/promoter_analysis/genome_gff/")

# load promoters fasta file generated from extract_upstream_features.R
promoters = readDNAStringSet("rpadi_upstream_promoters.fasta")

# lets define a k-value to pass to kmer function. Trying a value of 8 to start
k = 8

# making matrix of kmers with len 8 
kmer_matrix = oligonucleotideFrequency(promoters, width = k)
head(kmer_matrix)
