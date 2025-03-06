# load packages
#BiocManager::install("memes")

library(memes)
library(Biostrings)

# setting dir to the project data dir
setwd("C:/Users/danie/OneDrive - Colostate/NIFA_PROJECT/Obj1/analysis/promoter_analysis/genome_gff/")

# setting fasta location 
rhythmic_fa = readDNAStringSet("rpadi_upstream_promoters_DREME_rhythmic_5-3.fasta")
# and for the non rhythmic
nonrhythmic_fa = system.file("", package = "memes")



# Running dreme 
runDreme(rhythmic_fa, "shuffle", outdir = "meme_output_rhythmic")
