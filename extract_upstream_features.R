# load libraries
#library(txdbmaker)
library(GenomicFeatures)
#BiocManager::install("Rsamtools")
library(Rsamtools)
library(Biostrings)

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

## if you want them oriented in the gene’s 5'->3' direction, most of the time this is not what we want
# neg_idx = strand(rpadi_promoters) == "-"
# rpadi_upseqs[neg_idx] = reverseComplement(rpadi_upseqs[neg_idx])

# rename the sequences using your gene IDs
names(rpadi_upseqs) = rpadi_genes$gene_id

## write to a FASTA file with gene names. Blast the seq and view in jbrowse with v2 annotations labelled to see if things line up
writeXStringSet(rpadi_upseqs, "rpadi_upstream_promoters.fasta")
