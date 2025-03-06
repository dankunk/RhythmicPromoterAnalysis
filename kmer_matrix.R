# load libraries 
library(Biostrings)
library(dplyr)

setwd("C:/Users/danie/OneDrive - Colostate/NIFA_PROJECT/Obj1/analysis/promoter_analysis/genome_gff/")

# load promoters fasta file generated from extract_upstream_features.R
promoters = readDNAStringSet("rpadi_upstream_promoters.fasta")

# Lets load a file with all the rhythmic data (example file - filtered for rhythmicity)
rhythmic_df = read.csv("example_rhythmic_dataset_jtk.csv")

# lets grab all the rhythmic and non-rhythmic genes sorting by p-value
rhythmic_promoter = rhythmic_df %>%
                        filter(ADJ.P <= 0.05)
  
nonrhythmic_promoter = rhythmic_df %>%
                        filter(ADJ.P > 0.05)

# lets define a k-value to pass to kmer function. Trying a value of 8 to start
k = 6

# making matrix of kmers with len 8 
rhythmic_kmer = oligonucleotideFrequency(rhythm)


















