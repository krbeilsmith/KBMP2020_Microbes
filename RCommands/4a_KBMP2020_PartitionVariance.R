# KBMP2020_PartitionVariance
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Summer 2020

##########################################################################################################################################################

# The commands for PERMANOVA and ordination below broadly follow the tutorials here: http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
# and here: https://bioconductor.riken.jp/packages/3.0/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html

# When using the distance and dissimilarities below, note that 16S data is compositional: https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full
# The data below is rarefied, or randomly subsampled, to even read counts for each sample in order to compare abundances between samples.

##########################################################################################################################################################

# The first goal is to find the study variables significantly associated with bacterial community variation.
# In this case, we're considering variation in community composition, specifically differences in the presence and relative abundance of ASVs between 
# the communities in each sample.

# Presence-absence variation can be represented by the Raup-Crick dissimilarity between samples. This is a probability of samples having different
# ASV composition. It is assessed by simulation given a null model in which the probability of selecting ASVs is proportional to ASV frequency.
# Documentation: https://rdrr.io/rforge/vegan/man/raupcrick.html
# Reference: Chase, J.M., Kraft, N.J.B., Smith, K.G., Vellend, M. and Inouye, B.D. (2011). Ecosphere. [doi:10.1890/ES10-00117.1]

# Alternatively, the Bray-Curtis dissimilarity quantifies the differences between the bacterial communities based on the abundance differences between ASVs
# in each sample. Documentation: https://rdrr.io/rforge/vegan/man/vegdist.html

# Finally, the UniFrac distance incorporates presence-absence variation as well as phylogenetic relatedness between the ASVs present in samples based on the
# 16S gene tree. The weighted version of UniFrac incorporates abundance differences too. 
# Documentation: https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/UniFrac

##########################################################################################################################################################

# The matrix of Raup-Crick dissimilarities, Bray-Curtis dissimilarities, or UniFrac distances between samples can be used to attribute additive 
# proportions of the total variance to variables. 
# Ecology 82, no. 1 (2001): 290-297.

# An F ratio statistic can be calculated to compare the sum of squared differences within groups defined by these variables and between them. A 
# non-parametric p-value can be assigned to this statstic with a permutational approach. (PERMANOVA)
# Austral ecology 26, no. 1 (2001): 32-46.
# Documentation: https://rdrr.io/rforge/vegan/man/adonis.html

##########################################################################################################################################################

# OUTLINE:
# Encode sample data as factors
# Rarefy data
# Calculate pairwise sample dissimilarity or distance (Raup-Crick, Bray-Curtis, UniFrac, Weighted UniFrac)
# Use matrix of sample dissimilarities and sample data as input for a permutational ANOVA on each study variable
# Select the study variables significantly associated with sample dissimilarities
# Nest terms based on study design and rank nested terms by cumulative partial regression coefficients, R^2
# Using the full model, partition variance in community composition across model terms

##########################################################################################################################################################

# Make a copy of the phyloseq object and encode sample data variables as factors.
physeq_factored <- physeq

sample_data(physeq_factored)$Year <- factor(sample_data(physeq_factored)$Year, levels=c("1", "2"))
sample_data(physeq_factored)$Stage <- factor(sample_data(physeq_factored)$Stage, levels=c("Soil", "TwoLeaf","FourLeaf","SixLeaf","EightLeaf","Flowering","Senescent"))
sample_data(physeq_factored)$PlantPart <- factor(sample_data(physeq_factored)$PlantPart, levels=c("Soil", "Roots","RosLeaves","Stems","CauLeaves","Flowers","Siliques"))
sample_data(physeq_factored)$Ecotype <- factor(sample_data(physeq_factored)$Ecotype)
sample_data(physeq_factored)$Site <- factor(sample_data(physeq_factored)$Site, levels=c("ME", "WW"))
sample_data(physeq_factored)$SamplePlate <- factor(sample_data(physeq_factored)$SamplePlate)
sample_data(physeq_factored)$MiSeqRun <- factor(as.character(as.numeric(sample_data(physeq_factored)$MiSeqRun)))
sample_data(physeq_factored)$PlantID <- factor(as.character(sample_data(physeq_factored)$PlantID))
sample_data(physeq_factored)$Soil <- factor(sample_data(physeq_factored)$Soil,levels=c("Y","N"))

# Focus on plant samples
physeq_factored_plants <- subset_samples(physeq_factored,PlantPart!="Soil")
physeq_factored_plants <- prune_taxa(taxa_sums(physeq_factored_plants)>0,physeq_factored_plants)

# Rarefy samples to even depth (will cause samples to drop from set)
R_physeq <- rarefy_even_depth(physeq_factored_plants,1000)

# Raup Crick index for Rarefied Set
raup_Rare <- phyloseq::distance(R_physeq, method = "raup")

# Bray Curtis dissimilarity for Rarefied Set
bc_Rare <- phyloseq::distance(R_physeq, method = "bray")

# Calculate UniFrac distance for Rarefied Set
uni_Rare <- UniFrac(R_physeq, weighted=FALSE, normalized=TRUE, parallel=TRUE, fast=TRUE)

# Calculate Weighted UniFrac distance for Rarefied Set
wuni_Rare <- UniFrac(R_physeq, weighted=TRUE, normalized=TRUE, parallel=TRUE, fast=TRUE)

# These are the variables to test in the permutational ANOVA to find those that are significant at alpha=0.001 and should be included in the model
TestTerms = c("Ecotype","PlantPart","Stage","Site","Year","SamplePlate","MiSeqRun")

# These are the dissimilarities or distances used to quantify community variation
Distances = c("raup_Rare","bc_Rare","uni_Rare","wuni_Rare")

# This creates a dataframe for the ANOVA results.
ResultsDF <- NULL
ResultsDF <- data.frame(matrix(ncol = 8, nrow = 0),stringsAsFactors = FALSE)
ResultsDF <- setNames(data.frame(ResultsDF), c("Distance","Var","Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)"))

for(b in Distances){ # For each measure of community variation
  # For each variable in the study
  for(a in TestTerms){
    # ANOVA on the matrix of sample-to-sample variation
    permanova_result <- adonis(as.formula(paste(b,a,sep='~')),data.frame(sample_data(R_physeq)),permutations=999)
    # Store results in dataframe and in a file
    write.table(data.frame(permanova_result$aov.tab),file=paste(paste("~/Documents/KBMP2020_Microbes/PERMANOVAs/Plant_ASV_",b,sep=""),a,sep="_"), sep="\t", row.names=TRUE,quote=FALSE)
    ResultsDF[nrow(ResultsDF)+1,] <- c(b,a,as.character(permanova_result$aov.tab[a,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")]))
  }
}

##########################################################################################################################################################

# This will filter the dataframe of ANOVA results for individual variables to those significant at alpha=0.001
Vars_for_Model <- ResultsDF[ResultsDF[,'Pr(>F)']<=0.001,]
Vars_for_Model <- Vars_for_Model[with(Vars_for_Model, order(Vars_for_Model$R2,decreasing=T)),]
# Store the table of model terms
write.table(Vars_for_Model,file="~/Documents/KBMP2020_Microbes/PERMANOVAs/ModelTerms_Alpha001", sep="\t", row.names=TRUE,quote=FALSE)
# Display the table
z = ztable(Vars_for_Model)
z

##########################################################################################################################################################

# To see whether a PlantID term helps explain sample variance, the data was pruned to plants which had more than one sample after quality control.

# Find plant IDs with more than one sample in the dataset
Retain <- c()
for(SampID in sample_names(physeq_factored_plants)){
  SampPlantID <- as.character(sample_data(physeq_factored_plants)[SampID,"PlantID"][[1]]) # get the plant ID of the focal sample
  # prune the data to just samples sharing the plant ID
  samps_to_keep = as.character(get_variable(physeq_factored_plants, "PlantID")) == SampPlantID & as.character(rownames(sample_data(physeq_factored_plants))) != SampID
  if(table(samps_to_keep)["FALSE"] == length(sample_names(physeq_factored_plants))){ # if there are no samples with matching plant IDs:
    Retain <- Retain # add nothing to the list
  }else{
    Retain <- c(Retain, SampID) # add sample ID to the list
  }
}

# Prune data to plants with multiple samples in order to test effect of plant ID in combination with other variables:
physeq_factored_plants_ID <-prune_samples(sample_names(physeq_factored_plants) %in% Retain, physeq_factored_plants)
physeq_factored_plants_ID <- prune_taxa(taxa_sums(physeq_factored_plants_ID)>0,physeq_factored_plants_ID)

permanova_results_ID <- adonis(t(otu_table(physeq_factored_plants_ID)) ~PlantID,
                               data.frame(sample_data(physeq_factored_plants_ID)),method="raup",permutations=999)
permanova_results_ID

sink(file = "~/Documents/KBMP2020_Microbes/PERMANOVAs/PlantID_Raup", append = TRUE, type = c("output"), split = FALSE)
print(permanova_results_ID)
closeAllConnections()

# Not significant at alpha=0.001 for Raup-Crick
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# PlantID   375    63.849 0.17026  1.2534 0.44763  0.012 *

permanova_results_ID <- adonis(t(otu_table(physeq_factored_plants_ID)) ~PlantID,
                               data.frame(sample_data(physeq_factored_plants_ID)),method="bray",permutations=999)
permanova_results_ID

sink(file = "~/Documents/KBMP2020_Microbes/PERMANOVAs/PlantID_Bray", append = TRUE, type = c("output"), split = FALSE)
print(permanova_results_ID)
closeAllConnections()

# Significant at alpha=0.001 for Bray-Curtis
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# PlantID   375    189.42 0.50511  1.1724 0.43118  0.001 ***

uni_ID <- UniFrac(physeq_factored_plants_ID, weighted=FALSE, normalized=TRUE, parallel=TRUE, fast=TRUE)

permanova_results_ID <- adonis(uni_ID~PlantID,
                               data.frame(sample_data(physeq_factored_plants_ID)),permutations=999)
permanova_results_ID

sink(file = "~/Documents/KBMP2020_Microbes/PERMANOVAs/PlantID_Uni", append = TRUE, type = c("output"), split = FALSE)
print(permanova_results_ID)
closeAllConnections()

# Not significant at alpha=0.001 for UniFrac
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# PlantID   375    100.44 0.26785  1.0674 0.40832  0.003 **

##########################################################################################################################################################

# Plant Part is nested within Stage because some tissues were only available to sample at some stages. 
# Stages are nested within Years because there were some differences in the sampling times between years.
# Sample Plate is nested within MiSeqRun because multiple sample plates were included on each sequencing run.
# Site is treated as an independent term. 

# Partition variance to terms in the full model with PERMANOVA
for(b in Distances){
  print(b)
  permanova_result <- adonis(as.formula(paste(b,"Year/Stage/PlantPart + Site + MiSeqRun/SamplePlate",sep='~')),data.frame(sample_data(R_physeq)),permutations=999)
  print(permanova_result$aov.tab["Residuals","R2"])
  # Store results in a file
  write.table(data.frame(permanova_result$aov.tab[,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")]),file=paste("~/Documents/KBMP2020_Microbes/PERMANOVAs/Plant_FullModel",b,sep="_"), sep="\t", row.names=TRUE,quote=FALSE)
}

####################################################################################################################################################

# Given that tissue type and stage are both associated with community composition, the next goal is to tease out whether the variance explained by
# tissue type changes at different stages. The R2 of tissue type in PERMANOVAs on root and rosette samples was compared between six leaf and
# and flowering plants since these stages were both sampled in all sites and years. Then, the R2 of tissue type in PERMANOVAs on root and stem 
# samples were compared between flowering and senescent stages.

# Tissue effect over time
# 1. Vegetative vs Flowering Roots and Rosettes

Vegetative_PhylloRhizo <- subset_samples(R_physeq, PlantPart %in% c("Roots","RosLeaves") & Stage %in% c("SixLeaf"))
Vegetative_PhylloRhizo <- prune_taxa(taxa_sums(Vegetative_PhylloRhizo)>0,Vegetative_PhylloRhizo)
Vegetative_PhylloRhizo_permanova_results <- adonis(t(otu_table(Vegetative_PhylloRhizo)) ~PlantPart,data.frame(sample_data(Vegetative_PhylloRhizo)),method="raup",permutations=999)
Vegetative_PhylloRhizo_permanova_results

#            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# PlantPart  1   0.20217 0.202173  15.397 0.18687  0.013 *

Flowering_PhylloRhizo <- subset_samples(R_physeq, PlantPart %in% c("Roots","RosLeaves") & Stage %in% c("Flowering"))
Flowering_PhylloRhizo <- prune_taxa(taxa_sums(Flowering_PhylloRhizo)>0,Flowering_PhylloRhizo)
Flowering_PhylloRhizo_permanova_results <- adonis(t(otu_table(Flowering_PhylloRhizo)) ~PlantPart,data.frame(sample_data(Flowering_PhylloRhizo)),method="raup",permutations=999)
Flowering_PhylloRhizo_permanova_results

#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# PlantPart   1   0.90458 0.90458  263.32 0.62057  0.001 ***

# 2. Flowering vs Senescent Roots and Stems

Flowering_PhylloRhizo <- subset_samples(R_physeq, PlantPart %in% c("Roots","Stems") & Stage %in% c("Flowering"))
Flowering_PhylloRhizo <- prune_taxa(taxa_sums(Flowering_PhylloRhizo)>0,Flowering_PhylloRhizo)
Flowering_PhylloRhizo_permanova_results <- adonis(t(otu_table(Flowering_PhylloRhizo)) ~PlantPart,data.frame(sample_data(Flowering_PhylloRhizo)),method="raup",permutations=999)
Flowering_PhylloRhizo_permanova_results

#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# PlantPart  1   0.46692 0.46692  59.057 0.39096  0.001 ***

Senescent_PhylloRhizo <- subset_samples(R_physeq, PlantPart %in% c("Roots","Stems") & Stage %in% c("Senescent"))
Senescent_PhylloRhizo <- prune_taxa(taxa_sums(Senescent_PhylloRhizo)>0,Senescent_PhylloRhizo)
Senescent_PhylloRhizo_permanova_results <- adonis(t(otu_table(Senescent_PhylloRhizo)) ~PlantPart,data.frame(sample_data(Senescent_PhylloRhizo)),method="raup",permutations=999)
Senescent_PhylloRhizo_permanova_results

#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# PlantPart   1    7.0536  7.0536  1582.1 0.88085  0.001 ***
