# load libraries
#library(txdbmaker)
library(GenomicFeatures)
#BiocManager::install("Rsamtools")
library(Rsamtools)
library(Biostrings)
library(dplyr)

# sql database can be loaded, its good to reset your dir here. 
setwd('C:/Users/danie/OneDrive - Colostate/NIFA_PROJECT/Obj1/analysis/promoter_analysis/genome_gff/')

# set path to indexed genome 
rpadi_genome = "C:/Users/danie/OneDrive - Colostate/NIFA_PROJECT/Obj1/analysis/promoter_analysis/genome_gff/R_padi_v2.fasta"

# make FaFile if not already loaded
rpadi_fafile = FaFile(rpadi_genome)

# set path to database (it should have info for chromosome starts and ends, see create_TxDb.R)
rpadi_txdb = loadDb("rpadi_V2_txdb.sqlite")

# extract genes 
rpadi_genes = genes(rpadi_txdb)

# view extracted genes
head(rpadi_genes)

# extract promoter regions around each gene TSS:
rpadi_promoters = promoters(rpadi_genes, upstream=1500, downstream=0)

# trim them so they don't exceed contig boundaries:
rpadi_promoters = trim(rpadi_promoters)

# now extract these (possibly shortened) promoter ranges
rpadi_upseqs = getSeq(rpadi_fafile, rpadi_promoters)

## if you want them oriented in the gene’s 5'->3' direction,
neg_idx = strand(rpadi_promoters) == "-"
rpadi_upseqs[neg_idx] = reverseComplement(rpadi_upseqs[neg_idx])

# rename the sequences using your gene IDs
names(rpadi_upseqs) = rpadi_genes$gene_id

# # grabbing meta info 
# prom_chr    = as.character(seqnames(rpadi_promoters))
# prom_start  = start(rpadi_promoters)
# prom_end    = end(rpadi_promoters)
# prom_strand = as.character(strand(rpadi_promoters))
# gene_ids    = rpadi_genes$gene_id
# 
# # Build a descriptive FASTA name per gene promoter
# # e.g. "g1|scf7180000007785:2235-6402(-)"
# fancy_names = paste0(
#   gene_ids, "|",
#   prom_chr, ":",
#   prom_start, "-",
#   prom_end, "(",
#   prom_strand, ")"
# )
# 
# # Make a named vector that maps from gene_id -> fancy_name
# fancy_map = setNames(fancy_names, gene_ids)



rhythmic_df = read.csv("example_rhythmic_dataset_jtk.csv")

rhythmic_IDs = rhythmic_df %>%
  dplyr::filter(ADJ.P <= 0.05)
rhythmic_IDs$gene_ID = sub("\\.t.*", "", rhythmic_IDs$CycID)

non_rhythmic_IDs = rhythmic_df %>%
  dplyr::filter(ADJ.P > 0.05)
non_rhythmic_IDs$gene_ID = sub("\\.t.*", "", non_rhythmic_IDs$CycID)

rhythmic_seqs = rpadi_upseqs[names(rpadi_upseqs) %in% rhythmic_IDs$gene_ID]
non_rhythmic_seqs = rpadi_upseqs[names(rpadi_upseqs) %in% non_rhythmic_IDs$gene_ID]


# names(rhythmic_seqs)      = fancy_map[names(rhythmic_seqs)]
# names(non_rhythmic_seqs)  = fancy_map[names(non_rhythmic_seqs)]


head(rhythmic_seqs)
head(non_rhythmic_seqs)


writeXStringSet(rhythmic_seqs,
                "rpadi_upstream_promoters_DREME_rhythmic_5-3.fasta")
writeXStringSet(non_rhythmic_seqs,
                "rpadi_upstream_promoters_DREME_non-rhythmic_5-3.fasta")

