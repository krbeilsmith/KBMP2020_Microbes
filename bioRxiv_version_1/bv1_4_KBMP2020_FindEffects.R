# KBMP2020_FindEffects
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Autumn 2019

# These are the commands used to find the variables associated with variation in the plant bacterial community composition at different levels of taxonomic
# grouping from fine (amplicon sequence variants) to coarse (phylum).

# The commands for ordination and PERMANOVA below broadly follow the tutorials here: http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
# and here: https://bioconductor.riken.jp/packages/3.0/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html

####################################################################################################################################################
# The first goal is to find the study variables significantly associated with bacterial community variation.
# In this case, we're considering variation in community composition, specifically differences in the presence of ASVs between the communities in 
# each sample.
# This kind of variation can be represented by the Raup-Crick dissimilarity between samples. This is a probability of samples having different
# ASV composition. It is assessed by simulation given a null model in which the probability of selecting ASVs is proportional to ASV frequency.
# Documentation: https://rdrr.io/rforge/vegan/man/raupcrick.html
# Reference: Chase, J.M., Kraft, N.J.B., Smith, K.G., Vellend, M. and Inouye, B.D. (2011). Ecosphere. [doi:10.1890/ES10-00117.1]

# The matrix of Raup-Crick dissimilarities between samples can be used to attribute additive proportions of the total variance to variables. 
# Ecology 82, no. 1 (2001): 290-297.

# An F ratio statistic can be calculated to compare the sum of squared differences within groups defined by these variables and between them. A 
# non-parametric p-value can be assigned to this statstic with a permutational approach.
# Austral ecology 26, no. 1 (2001): 32-46.
# (PERMANOVA)
# Documentation: https://rdrr.io/rforge/vegan/man/adonis.html

# Prepare data for testing
# Encode sample data as factors
# Select plant samples
# Test for associations
# Calculate pairwise sample dissimilarity (Raup-Crick)
# Use matrix of sample dissimilarities and sample data as input for a permutational ANOVA on each study variable
# Select terms for model
# Select the study variables significantly associated with sample dissimilarities
# Nest terms based on study design
# PERMANOVA with nested terms
# Rank nested terms by cumulative partial regression coefficients, R^2, to create final model
# Partition variance in community composition across model terms
# Run PERMANOVA with the full model and the matrix of sample dissimilarities based on ASV presence

####################################################################################################################################################

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

# Calculate dissimilarity
raup_physeq_factored_plants <- phyloseq::distance(physeq_factored_plants, method = "raup")

# These are the variables to test one by one in the permutational ANOVA to find those that are significant at p=0.001 and should be included in the model.
TestTerms = c("PlantID","Ecotype","PlantPart","Stage","Site","Year","SamplePlate","MiSeqRun")

# This creates a dataframe for the ANOVA results.
ResultsDF <- NULL
ResultsDF <- data.frame(matrix(ncol = 7, nrow = 0),stringsAsFactors = FALSE)
ResultsDF <- setNames(data.frame(ResultsDF), c("Var","Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)"))

# For each variable in the study, repeat the ANOVA on the matrix of sample-to-sample Raup-Crick distances. Print the full results and store.
for(a in TestTerms){
  permanova_result <- adonis(as.formula(paste('raup_physeq_factored_plants',a,sep='~')),data.frame(sample_data(physeq_factored_plants)),permutations=999)
  write.table(data.frame(permanova_result$aov.tab),file=paste("~/Documents/KBMP2020_Microbes/PERMANOVAs/Plant_ASV_",a,sep=""), sep="\t", row.names=TRUE,quote=FALSE)
  ResultsDF[nrow(ResultsDF)+1,] <- c(a,as.character(permanova_result$aov.tab[a,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")]))
}

# This will filter the dataframe of ANOVA results for individual variables to those significant at p=0.001 and then order them by descending R2.
Vars_for_Model<-ResultsDF[ResultsDF[,'Pr(>F)']<=0.001,]
Vars_for_Model<-Vars_for_Model[with(Vars_for_Model, order(Vars_for_Model$R2,decreasing=T)),]
write.table(Vars_for_Model,file="~/Documents/KBMP2020_Microbes/PERMANOVAs/Plant_ASV_VarsforModel", sep="\t", row.names=TRUE,quote=FALSE)

# Display the table
z = ztable(Vars_for_Model)
z

# The significant terms are Plant Part, Stage, Sample Plate, Site, Year, and MiSeq Run. 

# The term for PlantID was significant at a lower threshold (p = 0.02). To see whether adding this term helps explain sample variance, the data was pruned to
# plants which had more than one sample after quality control and a PERMANOVA was run with the terms for Plant Part, Stage, Site, Year, MiSeq Run, Sample Plate,
# and PlantID.

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

# The term for PlantID is nested within Stage because a plant could only be harvested once. PlantPart was left independent because samples from the same plant
# could be from different tissue types. 

# If the nested Stage/PlantID term is put first in the model, it is insignificant (p=0.490). In addition, the order of R^2 values for the other terms does not
# change and their magnitudes do not change greatly.

permanova_results_ID <- adonis(t(otu_table(physeq_factored_plants_ID)) ~Stage/PlantID + PlantPart + SamplePlate + Site + Year + MiSeqRun,
                               data.frame(sample_data(physeq_factored_plants_ID)),method="raup",permutations=999)
permanova_results_ID

# Plant Part is nested within Stage because some tissues were only available to sample at some stages. 
# Stages are nested within Years because there were some differences in the sampling times between years.
# Sample Plate is nested within MiSeqRun because multiple sample plates were included on each sequencing run.
# Site is treated as an independent term. 
# The two nested terms and Site were used in PERMANOVAs below and ranked by total R2 to obtain the final model for the PERMANOVA.

# These are the variables (with nesting) to test
TestTerms = c("Year/Stage/PlantPart","MiSeqRun/SamplePlate","Site")

# This creates a dataframe for the PERMANOVA results.
ResultsDF <- NULL
ResultsDF <- data.frame(matrix(ncol = 7, nrow = 0),stringsAsFactors = FALSE)
ResultsDF <- setNames(data.frame(ResultsDF), c("Var","Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)"))

# This creates a dataframe for the total R2 of each nested term.
NestedDF <- NULL
NestedDF <- data.frame(matrix(ncol = 2, nrow = 0),stringsAsFactors = FALSE)
NestedDF <- setNames(data.frame(NestedDF), c("Term","Total R2"))

# For each set of nested variables, repeat the PERMANOVA on the matrix of sample-to-sample Raup-Crick dissimilarities. Print the full results and store.
for(a in TestTerms){
  permanova_result <- adonis(as.formula(paste('raup_physeq_factored_plants',a,sep='~')),data.frame(sample_data(physeq_factored_plants)),permutations=999)
  write.table(data.frame(permanova_result$aov.tab),file=paste("~/Documents/KBMP2020_Microbes/PERMANOVAs/Nested_Plant_ASV_",gsub("/", "_", a),sep=""), sep="\t", row.names=TRUE,quote=FALSE)
  #print(paste(a,sum(permanova_result$aov.tab[1:length(strsplit(a,"/")[[1]]),"R2"]),sep=": "))
  NestedDF[nrow(NestedDF)+1,] <- c(a,sum(permanova_result$aov.tab[1:length(strsplit(a,"/")[[1]]),"R2"]))
  rowcount=1
  for(row in strsplit(a,"/")[[1]]){
    ResultsDF[nrow(ResultsDF)+1,] <- c(row,as.character(permanova_result$aov.tab[rowcount,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")]))
    rowcount=rowcount+1
  }
}

# This will filter the dataframe of PERMANOVA results for individual variables to those significant at p=0.001 and then order them by descending R^2.
Vars_for_Model<-ResultsDF[ResultsDF[,'Pr(>F)']<=0.001,]
Vars_for_Model<-Vars_for_Model[with(Vars_for_Model, order(Vars_for_Model$R2,decreasing=T)),]
write.table(Vars_for_Model,file="~/Documents/KBMP2020_Microbes/PERMANOVAs/Plant_ASV_VarsforModel_WithNesting", sep="\t", row.names=TRUE,quote=FALSE)

# This will store the total R^2 for each nested term and print the table.
write.table(NestedDF,file="~/Documents/KBMP2020_Microbes/PERMANOVAs/Plant_ASV_VarsforModel_NestedTermOrder", sep="\t", row.names=TRUE,quote=FALSE)
z = ztable(NestedDF)
z

# Based on the sum of the R^2 for the nested terms, they are ordered in the model as ~Year/Stage/PlantPart + MiSeqRun/SamplePlate + Site
# The order of terms in the model will affect the variation in sample dissimilarities explained by each term.
# This model is used to run a PERMANOVA on the Raup-Crick matrix for all plant samples.
# This time, the otu_table is provided so the coefficients for individual ASVs can be extracted from the PERMANOVA results.

full_permanova_results <- adonis(t(otu_table(physeq_factored_plants)) ~Year/Stage/PlantPart + MiSeqRun/SamplePlate + Site,data.frame(sample_data(physeq_factored_plants)),method="raup",permutations=999)
full_permanova_results

write.table(data.frame(full_permanova_results$aov.tab[,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")]),file="~/Documents/KBMP2020_Microbes/PERMANOVAs/Plant_ASV_FullModelResult_Raup", sep="\t", row.names=TRUE,quote=FALSE)

# This prints the results of the PERMANOVA to show how much variation in sample dissimilarities is explained by each study variable.
z = ztable(data.frame(perma_raup_physeq_factored_plants$aov.tab[,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")]))
z

# To see how types of samples vary in composition, eigendecomposition of the dissimilarity matrix is used to produce distances in a coordinate space. The
# first two axes of these space have the most distance between samples and their relative eigenvalues (eigenvalues divided by sum of all eigenvalues) are
# reported in the axes labels.

# Ordinate samples with Raup distance

#PCoA_raup  <- ordinate(physeq_factored_plants, "PCoA", distance=raup_physeq_factored_plants)

# PCoA_raup$values$Eigenvalues[1] / sum(PCoA_raup$values$Eigenvalues)
# PCoA_raup$values$Eigenvalues[2] / sum(PCoA_raup$values$Eigenvalues)

scree <- plot_scree(PCoA_raup)+coord_cartesian(xlim=c(0,25))+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(x="\nPCoA Axis", y="Fraction of Eigenvalues\n")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12))

ggsave("~/Documents/KBMP2020_Microbes/Figures/Raup_PCoA_Scree.tiff", plot = scree, device = NULL, path = NULL,
       scale = 1.3, width = 5, height = 4, units = c("in"),
       dpi = 600, limitsize = TRUE)

# Look for negative eigenvalues
# sum(abs(PCoA_raup$values$Eigenvalues[PCoA_raup$values$Eigenvalues < 0])) / sum(abs(PCoA_raup$values$Eigenvalues))

# Perform PCoA with correction for negative eigenvalues:
PCoA_raup_corrected <- ape::pcoa(raup_physeq_factored_plants,correction="cailliez")

# AXES LABELS (Eigenvalues (corrected) for PCoA vector / sum of eigenvalues)
# PCoA_raup_corrected$values$Eigenvalues[1] / sum(PCoA_raup_corrected$values$Eigenvalues)
# 45.5%
# PCoA_raup_corrected$values$Eigenvalues[2] / sum(PCoA_raup_corrected$values$Eigenvalues)
# 27.7%

# Plot Ordination with Raup distance
# !!! REMOVE LABELS !!! when changing the plot input. They have the relative eigenvalues for each PCoA axis hard-coded in for the publication version!
Palette<-c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")

pcoa <- plot_ordination(physeq_factored_plants, PCoA_raup_corrected, color="PlantPart") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  geom_point(size=1.5)+
  scale_colour_manual(values=Palette,labels=c("Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+ 
  #scale_shape_manual(values=c("\u25AC","\u25D6","\u25A0","\u25B2","\u25AE","\u25C6"),labels=c('Two Leaf','Four Leaf','Six Leaf','Eight Leaf','Flowering','Senescent'),name='Development Stage')+
  #scale_shape_manual(values=c("\u25CF","\u25CF","\u25CF","\u25CF","\u25FC","\u25C6"),labels=c('Two Leaf','Four Leaf','Six Leaf','Eight Leaf','Flowering','Senescent'),name='Development Stage')+
  theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),legend.direction="horizontal")+
  theme(legend.background = element_rect(color="black", size=.5),legend.position="top")+
  guides(color=guide_legend(override.aes = list(size=5)))+
  coord_cartesian(ylim=c(-0.75,0.75),xlim=c(-0.75,0.75))

pcoa <- pcoa  + labs(x="\nPCoA Axis 1 (45.5%)", y="PCoA Axis 2 (27.7%)\n")

ggsave("~/Documents/KBMP2020_Microbes/Figures/Raup_PCoA.tiff", plot = pcoa, device = NULL, path = NULL,
       scale = 1.3, width = 5, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)

# Focus on plant samples for flowering phyllosphere tissues
physeq_factored_plants_flowering <- subset_samples(physeq_factored,!(PlantPart %in% c("Soil","Roots")) & Stage%in%c("Flowering"))
physeq_factored_plants_flowering <- prune_taxa(taxa_sums(physeq_factored_plants_flowering)>0,physeq_factored_plants_flowering)

flowering_permanova_results <- adonis(t(otu_table(physeq_factored_plants_flowering)) ~Year/PlantPart + MiSeqRun/SamplePlate + Site,data.frame(sample_data(physeq_factored_plants_flowering)),method="raup",permutations=999)
flowering_permanova_results

# Tissue still has an effect just within the phyllosphere:
# Year:PlantPart         8     3.054 0.38174  1.9659 0.04085  0.010 **

# Calculate dissimilarity for flowering phyllosphere tissues
raup_physeq_factored_plants_flowering <- phyloseq::distance(physeq_factored_plants_flowering, method = "raup")

# Perform PCoA with correction for negative eigenvalues for flowering phyllosphere tissues:
PCoA_raup_corrected <- ape::pcoa(raup_physeq_factored_plants_flowering)
# correction="cailliez" unnecessary here

# AXES LABELS (Eigenvalues (corrected) for PCoA vector / sum of eigenvalues)
# PCoA_raup_corrected$values$Eigenvalues[1] / sum(PCoA_raup_corrected$values$Eigenvalues)
# 30.5%
# PCoA_raup_corrected$values$Eigenvalues[2] / sum(PCoA_raup_corrected$values$Eigenvalues)
# 18.0%

# Plot Ordination with Raup distance for flowering phyllosphere tissues
# !!! REMOVE LABELS !!! when changing the plot input. They have the relative eigenvalues for each PCoA axis hard-coded in for the publication version!
PaletteF<-c("#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")

pcoa_F <- plot_ordination(physeq_factored_plants_flowering, PCoA_raup_corrected, color="PlantPart") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  geom_point(size=1.5)+
  scale_colour_manual(values=PaletteF,labels=c("Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+ 
  #scale_shape_manual(values=c("\u25AC","\u25D6","\u25A0","\u25B2","\u25AE","\u25C6"),labels=c('Two Leaf','Four Leaf','Six Leaf','Eight Leaf','Flowering','Senescent'),name='Development Stage')+
  #scale_shape_manual(values=c("\u25CF","\u25CF","\u25CF","\u25CF","\u25FC","\u25C6"),labels=c('Two Leaf','Four Leaf','Six Leaf','Eight Leaf','Flowering','Senescent'),name='Development Stage')+
  theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),legend.direction="horizontal")+
  theme(legend.background = element_rect(color="black", size=.5),legend.position="top")+
  guides(color=guide_legend(override.aes = list(size=5),nrow=2))+
  coord_cartesian(ylim=c(-0.75,0.75),xlim=c(-0.75,0.75))

pcoa_F <- pcoa_F  + labs(x="\nPCoA Axis 1 (30.5%)", y="PCoA Axis 2 (18.0%%)\n")

ggsave("~/Documents/KBMP2020_Microbes/Figures/Raup_PCoA_Flowering.tiff", plot = pcoa_F, device = NULL, path = NULL,
       scale = 1.3, width = 5, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)

# Plot together
ggsave("~/Documents/KBMP2020_Microbes/Figures/TissueEffect.tiff", plot = grid.arrange(
  freq+ggtitle("A\n")+theme(plot.title = element_text(face="bold",size=20)),
  pcoa+ggtitle("B\n")+theme(plot.title = element_text(face="bold",size=20)),
  pcoa_F+ggtitle("C\n")+theme(plot.title = element_text(face="bold",size=20)),
  layout_matrix = rbind(c(1,2,3))), device = NULL, path = NULL,
  scale = 2, width = 10, height = 3, units = c("in"),
  dpi = 600, limitsize = TRUE)

###################################################################################################################################################

# The next goal is to find the levels of taxonomic groupings at which the model terms are significantly associated with the community composition.

# For each taxonomic level above ASV:
# Remove the unassigned taxa.
# Group the ASV counts at the desired level
# Calculate the matrix of Raup-Crick dissimilarities
# Run PERMANOVA with the model from above and the new Raup-Crick distances
# Store and print the results

####################################################################################################################################################
RankList = c("Rank6","Rank5", "Rank4", "Rank3", "Rank2")
for(b in RankList){
  mphyseq <- prune_taxa(row.names(tax_table(physeq_factored_plants))[!(grepl("aggregate",tax_table(physeq_factored_plants)[,b]))],physeq_factored_plants) # remove unknowns
  mphyseq <-prune_samples(sample_sums(mphyseq)>0,mphyseq) # remove any samples with 0 counts
  mphyseq <- tax_glom(mphyseq, taxrank=b) # group ASVs
  raup_mphyseq <- phyloseq::distance(mphyseq, method = "raup") # Raup-Crick dissimilarity
  perma_raup_mphyseq <- adonis(raup_mphyseq ~Year/Stage/PlantPart + MiSeqRun/SamplePlate + Site,data.frame(sample_data(mphyseq)),permutations=999) # PERMANOVA
  PERMAResults <- data.frame(perma_raup_mphyseq$aov.tab[,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")])
  PERMAResults
  write.table(PERMAResults,file=paste("~/Documents/KBMP2020_Microbes/PERMANOVAs/Plant_ASV_",b,sep=""), sep="\t", row.names=TRUE,quote=FALSE)
}

# Here, we want to explore whether the PERMANOVA terms for plant part and stage (which explain most of the variation in sample dissimilarities based on ASV presence)
# are significant and large when ASVs are randomly sorted into groups the same size as the families, orders, classes, and phyla in the dataset.
# If randomly assigning the ASVs to asymmetrically sized groups is sufficient to produce effects, then the observed effects in the original PERMANOVAs might be due
# to artifacts rather than real biological differences between the taxonomic groups related to their prevalence in different tissues or developmental stages.

# Remove ASVs with ambiguous classification at any level between family and phylum
physeq_classified <- prune_taxa(taxa_names(physeq) %in% row.names(tax_table(physeq))[tax_table(physeq)[,"Rank2"]!="aggregate unclassified phylum"] &
                                  taxa_names(physeq) %in% row.names(tax_table(physeq))[tax_table(physeq)[,"Rank3"]!="aggregate unclassified class"] &
                                  taxa_names(physeq) %in% row.names(tax_table(physeq))[tax_table(physeq)[,"Rank4"]!="aggregate unclassified order"] &
                                  taxa_names(physeq) %in% row.names(tax_table(physeq))[tax_table(physeq)[,"Rank5"]!= "aggregate unclassified family"], physeq)

# Remove any samples that no longer have ASVs                          
physeq_classified <-prune_samples(sample_sums(physeq_classified)>0,physeq_classified)
# 8889 taxa and 1272 samples

# Encode sample data variables as factors.
sample_data(physeq_classified)$Year <- factor(sample_data(physeq_classified)$Year, levels=c("1", "2"))
sample_data(physeq_classified)$Stage <- factor(sample_data(physeq_classified)$Stage, levels=c("Soil", "TwoLeaf","FourLeaf","SixLeaf","EightLeaf","Flowering","Senescent"))
sample_data(physeq_classified)$PlantPart <- factor(sample_data(physeq_classified)$PlantPart, levels=c("Soil", "Roots","RosLeaves","Stems","CauLeaves","Flowers","Siliques"))
sample_data(physeq_classified)$Ecotype <- factor(sample_data(physeq_classified)$Ecotype)
sample_data(physeq_classified)$Site <- factor(sample_data(physeq_classified)$Site, levels=c("ME", "WW"))
sample_data(physeq_classified)$SamplePlate <- factor(sample_data(physeq_classified)$SamplePlate)
sample_data(physeq_classified)$MiSeqRun <- factor(as.character(as.numeric(sample_data(physeq_classified)$MiSeqRun)))
sample_data(physeq_classified)$PlantID <- factor(as.character(sample_data(physeq_classified)$PlantID))
sample_data(physeq_classified)$Soil<-factor(sample_data(physeq_classified)$Soil,levels=c("Y","N")) # this will get used later for testing soil vs. plant samples

# Focus on plant samples
physeq_classified <- subset_samples(physeq_classified,PlantPart!="Soil")
physeq_classified <- prune_taxa(taxa_sums(physeq_classified)>0,physeq_classified)
# 7114 taxa and 1195 samples

# Create dataframe for storing results
RandomTaxaDF <- NULL
RandomTaxaDF <- data.frame(matrix(ncol = 5, nrow = 0),stringsAsFactors = FALSE)
RandomTaxaDF <- setNames(data.frame(RandomTaxaDF), c("Rank","Tissue_p","Tissue_r2","Stage_p","Stage_r2"))
# Recreate the taxonomic group structure with random taxonomy assignments
RankList = c("Rank5", "Rank4", "Rank3", "Rank2")
for(b in RankList){ # for each level of grouping
  Tab<-as.data.frame(table(tax_table(physeq_classified)[,b])) # get the sizes of the groups
  BinSizeList<-Tab[,"Freq"] # get the sizes of the groups
  for(r in 1:10){ # repeat ten times
    mphyseq<-physeq_classified # make a copy of the phyloseq object for storing randomized groups
    iphyseq<-physeq_classified # make a copy of the phyloseq object for pruning between random assignments
    count <- 0 
    for(k in BinSizeList){ # for each group size at the level of grouping
      count<-count+1
      print(paste(b,r,k,sep="; "))
      SelectedTaxa<-sample(taxa_names(iphyseq),k,replace=FALSE) # randomly sample taxa to create a group of the same size
      if(k != tail(BinSizeList,n=1)){ # if not the end of the list, remove those taxa from the list before the next random group is made
        iphyseq<-prune_taxa(!(taxa_names(iphyseq)%in%SelectedTaxa),iphyseq)
      }
      tax_table(mphyseq)[,b][row.names(tax_table(mphyseq))%in% SelectedTaxa] <- paste(b,"_iTaxon_",count) # replace assignment with random group
    }
    mphyseq <- tax_glom(mphyseq, taxrank=b) # collapse the ASVs into the random groups
    raup_mphyseq <- phyloseq::distance(mphyseq, method = "raup") # calculate distance
    perma_raup_mphyseq <- adonis(raup_mphyseq ~Year/Stage/PlantPart + MiSeqRun/SamplePlate + Site,data.frame(sample_data(mphyseq)),permutations=999) # PERMANOVA
    RandomTaxaResults <- data.frame(perma_raup_mphyseq$aov.tab[,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")])
    write.table(RandomTaxaResults,file=paste("~/Documents/KBMP2020_Microbes/PERMANOVAs/RandomTaxControl/Plant",paste(b,r,sep="_"),sep=""), sep="\t", row.names=TRUE,quote=FALSE)
    #assign(paste(b,r,sep="_"),perma_raup_mphyseq$aov.tab) # rename PERMANOVA results
    RandomTaxaDF[nrow(RandomTaxaDF)+1,] <- c(b,as.numeric(as.character(perma_raup_mphyseq$aov.tab["Year:Stage:PlantPart",])[6]),
                                             as.numeric(as.character(perma_raup_mphyseq$aov.tab["Year:Stage:PlantPart",])[5]),
                                             as.numeric(as.character(perma_raup_mphyseq$aov.tab["Year:Stage",])[6]),
                                             as.numeric(as.character(perma_raup_mphyseq$aov.tab["Year:Stage",])[5]))
  }
}

write.table(RandomTaxaDF,file="~/Documents/KBMP2020_Microbes/PERMANOVAs/Plant_ASV_RandomizedTaxaResults", sep="\t", row.names=TRUE,quote=FALSE)
z = ztable(RandomTaxaDF)
z
