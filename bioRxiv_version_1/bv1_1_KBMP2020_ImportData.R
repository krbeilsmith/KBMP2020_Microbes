# KBMP2020_ImportData
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Autumn 2019

# The commands below broadly follow the phyloseq tutorial here: https://joey711.github.io/phyloseq/import-data.html

# The sequencing data used here were collected by Matt Perisin as part of his 2016 dissertation "THE DYNAMICS OF BACTERIAL COMMUNITIES 
# ASSOCIATED WITH ARABIDOPSIS THALIANA" for the Committee on Microbiology at the University of Chicago. The raw data consisted of 1427 sets of
# paired 250bp Illumina MiSeq reads. These reads are deposited in the NCBI SRA under BioProject PRJNA607544: Bacterial assemblages in field-grown 
# Arabidopsis thaliana from southeast Michigan.

# Data pre-processing was performed with in Qiime2. Primers were removed with cutadapt and DADA2 was used to model errors and find ASVs. 
# The initial table included 1421 samples with 10,987 sequence variants. Sequences were filtered for PhiX, mitochondrial, and chloroplast 
# sequences (removed 35 sequence variants and 102 samples). Samples with notes in the metadata indicating any irregularities in the collection 
# were excluded (removed 45 samples). Sequence variants with a frequency lower than 2 counts and samples with fewer than 10 reads were
# excluded (removed 19 features and 8 samples). This left a table with 10,803 taxa and 1,272 samples.

# Load tools

# 16S data handling, normalization, phylogenetics
library(phyloseq)
library(ape)
library(picante)
library(DESeq2)
library(edgeR)

# Finding indicator taxa
library(indicspecies)

# File reading and writing
library(readr)
library(MASS)

# Working with dataframes
library(reshape)
library(reshape2)
library(dplyr)
library(tidyr)

# Making plots and tables
library(ggplot2)
library(ggpubr)
library(directlabels)
library(ggrepel)
library(ztable)
library(venn)
library(VennDiagram)
library(eulerr)
library(gridExtra)
library(ggstance)
library(viridis)
library(gplots)

# Parallelizing
library(BiocParallel)
register(MulticoreParam(4))
library(doParallel)
registerDoParallel(cores=3)
library(foreach)

# 1. Import Phyloseq

# Import OTU TABLE from .txt file and make matrix
# The table has counts for DADA2 amplicon sequence variants (ASVs) in each sample. In phyloseq it is still called the "otu" table.
otu_table <- read.table("~/Documents/KBMP2020_Microbes/Dataset/ASV_table.txt", header=TRUE, sep="\t", quote="",row.names=1)
otu_table = as.matrix(otu_table)

# Import METADATA
# In the original metadata file, a collection date was recorded for soil samples but the "stage" column only said "soil."
# In the recoded file, the stage of the plants on the date the soil was collected is recorded so that side-by-side plots of soil and plant 
# samples across developmental time can be made. In the original file, plant positions were recorded. In the recoded file, these positions
# were used to assign plant IDs so that samples from different tissues in the same plant could be grouped.
sample_data <- read.table("~/Documents/KBMP2020_Microbes/Dataset/July18_Mapping_recoded.txt", header=TRUE, sep="\t", quote="",row.names=1)

# Import TREE (rooted)
phy_tree = read_tree("~/Documents/KBMP2020_Microbes/Dataset/ASV_tree.nwk")

# Import TAXONOMY from csv file and make matrix
tax_table <- read.table("~/Documents/KBMP2020_Microbes/Dataset/ASV_taxonomy.tsv",header=TRUE,sep=",",quote="",row.names=1,stringsAsFactors = FALSE)

# Modify TAXONOMY
# In the taxonomy table, there are a lot of ambiguous assignments. 
# unique(tax_table$Rank6)
# To simplify grouping for plots, the ASVs which were unassigned or given a generic assignment of uncultured or ambiguous were recoded as 
# "aggregate unclassified [rank]" in the table. 
# This allowed the ambiguous ASVs to remain in the dataset for calculations of relative abundances and diversity metrics 
# but to be easily excluded from plots.

AmbiguousGenera <- c("uncultured beta proteobacterium","uncultured alpha proteobacterium","g__","uncultured bacterium","uncultured soil bacterium","uncultured Chlorobi bacterium",
                     "uncultured","uncultured Bacteroidetes bacterium","Ambiguous_taxa","possible genus 04","uncultured Acidobacteria bacterium","uncultured forest soil bacterium",
                     "uncultured delta proteobacterium","uncultured Gemmatimonadetes bacterium","uncultured actinobacterium","uncultured Verrucomicrobia subdivision 3 bacterium",
                     "uncultured Chloroflexi bacterium","uncultured Ktedonobacter sp.","uncultured Acidobacteriales bacterium","uncultured Firmicutes bacterium",
                     "uncultured proteobacterium","uncultured organism","uncultured Candidatus Saccharibacteria bacterium","uncultured Armatimonadetes bacterium",
                     "uncultured Candidatus Saccharibacteria bacterium", "uncultured Desulfuromonadales bacterium")

AmbiguousGroups <- c("p__","c__","o__","f__")

tax_table$Rank6<-if_else(tax_table$Rank6 %in% AmbiguousGenera, "aggregate unclassified genus",tax_table$Rank6)
tax_table$Rank5<-if_else(tax_table$Rank5 %in% AmbiguousGroups | grepl("uncultur",tax_table$Rank5) | grepl("Ambig",tax_table$Rank5) | grepl("Unkn",tax_table$Rank5), "aggregate unclassified family",tax_table$Rank5)
tax_table$Rank4<-if_else(tax_table$Rank4 %in% AmbiguousGroups | grepl("uncultur",tax_table$Rank4) | grepl("Ambig",tax_table$Rank4) | grepl("Unkn",tax_table$Rank4), "aggregate unclassified order",tax_table$Rank4)
tax_table$Rank3<-if_else(tax_table$Rank3 %in% AmbiguousGroups | grepl("uncultur",tax_table$Rank3) | grepl("Ambig",tax_table$Rank3) | grepl("Unkn",tax_table$Rank3), "aggregate unclassified class",tax_table$Rank3)
tax_table$Rank2<-if_else(tax_table$Rank2 %in% AmbiguousGroups | grepl("uncultur",tax_table$Rank2) | grepl("Ambig",tax_table$Rank2) | grepl("Unkn",tax_table$Rank2), "aggregate unclassified phylum",tax_table$Rank2)

tax_table = as.matrix(tax_table)

# Use phyloseq functions to turn the imported data into phyloseq components
OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(tax_table)
META = sample_data(sample_data)

# Merge phyloseq components into a class 
physeq = phyloseq(TAX,OTU,META,phy_tree)

# Assess missing taxonomy information
physeq # 10803 ASVs
subset_taxa(physeq, Rank6=="aggregate unclassified genus") # 4747 ASVs, 43.94%
subset_taxa(physeq, Rank5=="aggregate unclassified family") # 1914 ASVs, 17.71%
subset_taxa(physeq, Rank4=="aggregate unclassified order") # 883 ASVs, 8.17%
subset_taxa(physeq, Rank3=="aggregate unclassified class") # 321 ASVs, 2.97%
subset_taxa(physeq, Rank2=="aggregate unclassified phylum") # 62 tASVs, 0.57%
