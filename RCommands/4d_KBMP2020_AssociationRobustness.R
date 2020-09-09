# KBMP2020_AssociationRobustness
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Autumn 2019

# The associations between plant bacterial community composition and study variables vary with the normalizations applied to the data and the distance or
# dissimilarity used to quantify composition shfts between samples. These commands reproduce the PERMANOVA and PCoA with several approaches to show that
# the same variables are significantly associated with community composition regardless of how counts are processed or composition compared.

# The normalization approaches used are discussed here: Weisset al. Microbiome (2017) 5:27 DOI 10.1186/s40168-017-0237-y
# The DESeq2 normalization approach follows the tutorial here: https://joey711.github.io/phyloseq-extensions/DESeq2.html
# The TMM normalization approach follows the tutorial here: https://joey711.github.io/phyloseq-extensions/edgeR.html

######################################################################################################################################
# Approach 1: Use Bray-Curtis dissimilarity on relative abundance matrix
######################################################################################################################################

# Convert sample data to relative abundances
physeq_RA <- transform_sample_counts(physeq_factored_plants,  function(x) x / sum(x)) 

# Bray Curtis dissimilarity
bc_RA <- phyloseq::distance(physeq_RA, method = "bray")

######################################################################################################################################
# Approach 2: Rarefy data and use weighted UNIFRAC distance on abundance matrix
######################################################################################################################################

# Rarefy samples to even depth (will cause samples to drop from set)
Rare_physeq <- rarefy_even_depth(physeq_factored_plants,1000)

# Weighted UniFrac
wuni_Rare <- UniFrac(Rare_physeq, weighted=TRUE, normalized=TRUE, parallel=TRUE, fast=TRUE)

######################################################################################################################################
# Approach 3: Variance-stabilized count transformation and Bray-Curtis dissimilarity
######################################################################################################################################

# The commands below make a DESeq object and get the variance stabilization transformed counts for the OTU table.
# These are written to a file.

# DESeq_physeq <- phyloseq_to_deseq2(physeq_factored,~PlantPart+Stage+Site+Year)

# Transform counts with geometric mean to avoid zeros
# gm_mean = function(x, na.rm=TRUE){
# exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
# }
# geoMeans = apply(counts(DESeq_physeq), 1, gm_mean)

# Get variance stabilized OTU table and make negative counts zero
# DESeq_physeq = estimateSizeFactors(DESeq_physeq,geoMeans = geoMeans)
# DESeq_physeq = estimateDispersions(DESeq_physeq)
# VST_table = getVarianceStabilizedData(DESeq_physeq)
# dim(VST_table)
# VST_table = as.matrix(VST_table)
# VST_table[VST_table<0.0] <- 0.0

# write.table(otu_table(N_physeq),file="~/Documents/KBMP2020_Microbes/Outputs/VST_Physeq_ASV", sep="\t", row.names=TRUE,quote=FALSE)

N_physeq <- physeq_factored

# Replace OTU table with variance stabilized counts
VST_table <- read.table("~/Documents/KBMP2020_Microbes/Outputs/VST_Physeq_ASV", header=TRUE, sep="\t", quote="",row.names=1)
VST_table = as.matrix(VST_table)
otu_table(N_physeq) <- otu_table(VST_table, taxa_are_rows = TRUE)
N_physeq <- prune_taxa(taxa_names(N_physeq)%in%taxa_names(physeq_factored_plants), N_physeq)
N_physeq <- prune_samples(sample_names(N_physeq)%in%sample_names(physeq_factored_plants), N_physeq)

# Bray-Curtis dissimilarity
bc_VST <- phyloseq::distance(N_physeq, "bray")

######################################################################################################################################
# Approach 4: Prune rare taxa (present in fewer than three samples) and Raup-Crick dissimilarity
######################################################################################################################################

# Transform to presence-absence table and get the names of ASVs in at least 3 plant samples
physeq_factored_plants_PA <- transform_sample_counts(physeq_factored_plants, function(x) ifelse(x>0,1,0))
physeq_factored_plants_PA <- prune_taxa(taxa_sums(physeq_factored_plants_PA)>2, physeq_factored_plants_PA)
TaxaToKeep <- taxa_names(physeq_factored_plants_PA)

# Now we prune the data to just the ASVs present in >3 samples.
PrevalentTaxa <- prune_taxa(taxa_names(physeq_factored_plants) %in% TaxaToKeep, physeq_factored_plants)
PrevalentTaxa <- prune_samples(sample_sums(PrevalentTaxa)>0,PrevalentTaxa)

# Raup-Crick index
raup_PrevalentTaxa <- phyloseq::distance(PrevalentTaxa, method = "raup")

######################################################################################################################################
# Approach 5: Prune scarce taxa (fewer than 1000 total observations) and Raup-Crick dissimilarity
######################################################################################################################################

# Prune to just ASVs counted >1000 times in the plant samples
AbundantTaxa <- prune_taxa(taxa_sums(physeq_factored_plants)>1000, physeq_factored_plants)
AbundantTaxa <- prune_samples(sample_sums(AbundantTaxa)>0,AbundantTaxa)

# Raup-Crick index
raup_AbundantTaxa <- phyloseq::distance(AbundantTaxa, method = "raup")

######################################################################################################################################
# Approach 6: TMM normalization and Bray-Curtis dissimilarity
######################################################################################################################################

# The commands below make a normalized OTU count table, which is written to a file.

TMMphyseq <- physeq_factored_plants

# Add pseudocount to avoid zeros
# otu_table(TMMphyseq) <- otu_table(TMMphyseq) + 0.01

# TMM normalization
# NormFactors <- calcNormFactors(otu_table(TMMphyseq))
# NormTable<-(otu_table(TMMphyseq)%*%diag(NormFactors))
# colnames(NormTable)<-sample_names(TMMphyseq)

# Replace OTU table with normalized counts
# otu_table(TMMphyseq) <- otu_table(NormTable, taxa_are_rows = TRUE)

# write.table(otu_table(TMMphyseq),file="~/Documents/KBMP2020_Microbes/Outputs/EdgeR_Physeq_ASV", sep="\t", row.names=TRUE,quote=FALSE)

NormTable <- read.table("~/Documents/KBMP2020_Microbes/Outputs/EdgeR_Physeq_ASV", header=TRUE, sep="\t", quote="",row.names=1)

# Replace OTU table with normalized counts
otu_table(TMMphyseq) <- otu_table(NormTable, taxa_are_rows = TRUE)

bc_TMM <- phyloseq::distance(TMMphyseq, "bray")

######################################################################################################################################

# Run PERMANOVA for each version of the data
# Loop over distances based on the same phyloseq object.

distance_matrices <- c('bc_RA', 'bc_VST', 'bc_TMM')

for(b in distance_matrices){ # For each measure of community variation
    # ANOVA on the matrix of sample-to-sample variation
    permanova_result <- adonis(as.formula(paste(b,"Year/Stage/PlantPart + Site + MiSeqRun/SamplePlate",sep='~')),data.frame(sample_data(physeq_factored_plants)),permutations=999)
    # Save result
    sink(file = "~/Documents/KBMP2020_Microbes/PERMANOVAs/Robustness", append = TRUE, type = c("output"), split = FALSE)
    print(b)
    print(permanova_result)
    closeAllConnections()
}

distance_matrices <- c('wuni_Rare')

for(b in distance_matrices){ # For each measure of community variation
  # ANOVA on the matrix of sample-to-sample variation
  permanova_result <- adonis(as.formula(paste(b,"Year/Stage/PlantPart + Site + MiSeqRun/SamplePlate",sep='~')),data.frame(sample_data(Rare_physeq)),permutations=999)
  # Save result
  sink(file = "~/Documents/KBMP2020_Microbes/PERMANOVAs/Robustness", append = TRUE, type = c("output"), split = FALSE)
  print(b)
  print(permanova_result)
  closeAllConnections()
}

distance_matrices <- c('raup_PrevalentTaxa','raup_AbundantTaxa')

for(b in distance_matrices){ # For each measure of community variation
  po <- strsplit(b,"_")[[1]][2]
  phyloseq_object <- get(po)
  # ANOVA on the matrix of sample-to-sample variation
  permanova_result <- adonis(as.formula(paste(b,"Year/Stage/PlantPart + Site + MiSeqRun/SamplePlate",sep='~')),data.frame(sample_data(phyloseq_object)),permutations=999)
  # Save result
  sink(file = "~/Documents/KBMP2020_Microbes/PERMANOVAs/Robustness", append = TRUE, type = c("output"), split = FALSE)
  print(b)
  print(permanova_result)
  closeAllConnections()
}
 
######################################################################################################################################

# Assemble figure with PCoA plots for each treatment of the data.
# Loop over distances based on the same phyloseq object.

# !!! The ape pcoa function throws warnings about NaNs for some distance matrices. Function will run outside of loop in these cases, 
# but need to look into why it is still trying to take the sqaure root of negative eigenvalues after the correction in these cases.

distance_matrices <- c('bc_RA', 'bc_TMM', 'bc_VST')

for(d in distance_matrices){
  distance_matrix <- get(d)
  print(d)
  # Ordinate samples with distance
  PCoA_dist  <- ordinate(physeq_factored_plants, "PCoA", distance=distance_matrix)
  # Look for negative eigenvalues
  neg_check <- sum(abs(PCoA_dist$values$Eigenvalues[PCoA_dist$values$Eigenvalues > 0])) / sum(abs(PCoA_dist$values$Eigenvalues))
  if(neg_check < 1){
    print("correcting")
    # Perform PCoA with correction for negative eigenvalues:
    PCoA_dist <- ape::pcoa(distance_matrix,correction="cailliez")
  }
  assign(paste("PCoA",d,sep="_"),PCoA_dist)
  
  # AXES LABELS (Eigenvalues for PCoA vector / sum of eigenvalues)
  eig1 <- (PCoA_dist$values$Rel_corr_eig[1] / sum(PCoA_dist$values$Rel_corr_eig))*100
  eig2 <- (PCoA_dist$values$Rel_corr_eig[2] / sum(PCoA_dist$values$Rel_corr_eig))*100
  
  # Plot Ordination with distance
  Palette<-c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
  
  pcoa <- plot_ordination(physeq_factored_plants, PCoA_dist, color="PlantPart") + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))+
    font("ylab", size = 16, color = "black",face = "bold")+
    font("xlab", size = 16, color = "black",face = "bold")+
    font("xy.text", size = 12, color = "black", face = "bold")+
    geom_point(size=1.5)+
    scale_colour_manual(values=Palette,labels=c("Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+ 
    theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),
          legend.direction="horizontal")+
    theme(legend.background = element_rect(color="black", size=.5),legend.position="top")+
    guides(color=guide_legend(override.aes = list(size=5)))+
    coord_cartesian(ylim=c(-0.75,0.75),xlim=c(-0.75,0.75))
  
  # Modify axes title formats
  pcoa <- pcoa  + labs(x=paste(paste("\nPCoA Axis 1 (",round(eig1,1),sep=""),"%)"), y=paste(paste("PCoA Axis 2 (",round(eig2,1),sep=""), "%)\n"))
  
  # Assign plot unique name for making full figure later
  assign(paste("PCoA_Plot",d,sep="_"),pcoa)
}

distance_matrices <- c('wuni_Rare')

for(d in distance_matrices){
  distance_matrix <- get(d)
  print(d)
  # Ordinate samples with distance
  PCoA_dist  <- ordinate(Rare_physeq, "PCoA", distance=distance_matrix)
  # Look for negative eigenvalues
  neg_check <- sum(abs(PCoA_dist$values$Eigenvalues[PCoA_dist$values$Eigenvalues > 0])) / sum(abs(PCoA_dist$values$Eigenvalues))
  if(neg_check < 1){
    print("correcting")
    # Perform PCoA with correction for negative eigenvalues:
    PCoA_dist <- ape::pcoa(distance_matrix,correction="cailliez")
  }
  assign(paste("PCoA",d,sep="_"),PCoA_dist)
  
  # AXES LABELS (Eigenvalues for PCoA vector / sum of eigenvalues)
  eig1 <- (PCoA_dist$values$Rel_corr_eig[1] / sum(PCoA_dist$values$Rel_corr_eig))*100
  eig2 <- (PCoA_dist$values$Rel_corr_eig[2] / sum(PCoA_dist$values$Rel_corr_eig))*100
  
  # Plot Ordination with distance
  Palette<-c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
  
  pcoa <- plot_ordination(Rare_physeq, PCoA_dist, color="PlantPart") + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))+
    font("ylab", size = 16, color = "black",face = "bold")+
    font("xlab", size = 16, color = "black",face = "bold")+
    font("xy.text", size = 12, color = "black", face = "bold")+
    geom_point(size=1.5)+
    scale_colour_manual(values=Palette,labels=c("Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+ 
    theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),
          legend.direction="horizontal")+
    theme(legend.background = element_rect(color="black", size=.5),legend.position="top")+
    guides(color=guide_legend(override.aes = list(size=5)))+
    coord_cartesian(ylim=c(-0.75,0.75),xlim=c(-0.75,0.75))
  
  # Modify axes title formats
  pcoa <- pcoa  + labs(x=paste(paste("\nPCoA Axis 1 (",round(eig1,1),sep=""),"%)"), y=paste(paste("PCoA Axis 2 (",round(eig2,1),sep=""), "%)\n"))
  
  # Assign plot unique name for making full figure later
  assign(paste("PCoA_Plot",d,sep="_"),pcoa)
}

distance_matrices <- c('raup_AbundantTaxa','raup_PrevalentTaxa')

d='raup_AbundantTaxa'

for(d in distance_matrices){
  # Split string d to get the name of the phyloseq object
  po <- strsplit(d,"_")[[1]][2]
  phyloseq_object <- get(po)
  distance_matrix <- get(d)
  print(d)
  # Ordinate samples with distance
  PCoA_dist  <- ordinate(phyloseq_object, "PCoA", distance=distance_matrix)
  # Look for negative eigenvalues
  neg_check <- sum(abs(PCoA_dist$values$Eigenvalues[PCoA_dist$values$Eigenvalues > 0])) / sum(abs(PCoA_dist$values$Eigenvalues))
  if(neg_check < 1){
    print("correcting")
    # Perform PCoA with correction for negative eigenvalues:
    PCoA_dist <- ape::pcoa(distance_matrix,correction="cailliez")
  }
  assign(paste("PCoA",d,sep="_"),PCoA_dist)
  
  # AXES LABELS (Eigenvalues for PCoA vector / sum of eigenvalues)
  eig1 <- (PCoA_dist$values$Rel_corr_eig[1] / sum(PCoA_dist$values$Rel_corr_eig))*100
  eig2 <- (PCoA_dist$values$Rel_corr_eig[2] / sum(PCoA_dist$values$Rel_corr_eig))*100
  
  # Plot Ordination with distance
  Palette<-c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
  
  pcoa <- plot_ordination(phyloseq_object, PCoA_dist, color="PlantPart") + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))+
    font("ylab", size = 16, color = "black",face = "bold")+
    font("xlab", size = 16, color = "black",face = "bold")+
    font("xy.text", size = 12, color = "black", face = "bold")+
    geom_point(size=1.5)+
    scale_colour_manual(values=Palette,labels=c("Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+ 
    theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),
          legend.direction="horizontal")+
    theme(legend.background = element_rect(color="black", size=.5),legend.position="top")+
    guides(color=guide_legend(override.aes = list(size=5)))+
    coord_cartesian(ylim=c(-0.75,0.75),xlim=c(-0.75,0.75))
  
  # Modify axes title formats
  pcoa <- pcoa  + labs(x=paste(paste("\nPCoA Axis 1 (",round(eig1,1),sep=""),"%)"), y=paste(paste("PCoA Axis 2 (",round(eig2,1),sep=""), "%)\n"))
  
  # Assign plot unique name for making full figure later
  assign(paste("PCoA_Plot",d,sep="_"),pcoa)
}

##########################################################################################################################################

# Assemble plot

legendCol <- cowplot::get_legend(pcoa)

tiff(file = "~/Documents/KBMP2020_Microbes/Figures/PCoA_Robustness.tiff", width = 7000, height = 8800, units = "px", res = 800)
grid.arrange(
  legendCol,
  PCoA_Plot_wuni_Rare + theme(legend.position = "none",plot.title = element_text(face="bold",size=12))+ggtitle("A: Rarefied, Weighted UniFrac\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")), 
  PCoA_Plot_bc_RA+ theme(legend.position = "none",plot.title = element_text(face="bold",size=12))+ggtitle("B: Relative Abundance, Bray-Curtis\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_bc_VST+ theme(legend.position = "none",plot.title = element_text(face="bold",size=12))+ggtitle("C: VST, Bray-Curtis\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_bc_TMM+ theme(legend.position = "none",plot.title = element_text(face="bold",size=12))+ggtitle("D: TMM, Bray-Curtis\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_raup_PrevalentTaxa + theme(legend.position = "none",plot.title = element_text(face="bold",size=12))+ggtitle("E: Prevalence >3, Raup-Crick\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_raup_AbundantTaxa + theme(legend.position = "none",plot.title = element_text(face="bold",size=12))+ggtitle("F: Count >1000, Raup-Crick\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  widths = c(4, 4),
  heights = c(1,4,4,4),
  layout_matrix = rbind(c(1),
                        c(2, 3),
                        c(4, 5),
                        c(6, 7))
)
dev.off()

png(file = "~/Documents/KBMP2020_Microbes/Figures/PCoA_Robustness.png", width = 8800, height = 6800, units = "px", res = 800)
grid.arrange(
  legendCol,
  PCoA_Plot_wuni_Rare + theme(legend.position = "none",plot.title = element_text(face="bold",size=12))+ggtitle("A: Rarefied, Weighted UniFrac\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")), 
  PCoA_Plot_bc_RA+ theme(legend.position = "none",plot.title = element_text(face="bold",size=12))+ggtitle("B: Relative Abundance, Bray-Curtis\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_bc_VST+ theme(legend.position = "none",plot.title = element_text(face="bold",size=12))+ggtitle("C: VST, Bray-Curtis\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_bc_TMM+ theme(legend.position = "none",plot.title = element_text(face="bold",size=12))+ggtitle("D: TMM, Bray-Curtis\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_raup_PrevalentTaxa + theme(legend.position = "none",plot.title = element_text(face="bold",size=12))+ggtitle("E: Prevalence >3, Raup-Crick\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_raup_AbundantTaxa + theme(legend.position = "none",plot.title = element_text(face="bold",size=12))+ggtitle("F: Count >1000, Raup-Crick\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  widths = c(4, 4, 4),
  heights = c(1,4,4),
  layout_matrix = rbind(c(1),
                        c(2, 6, 7),
                        c(3, 4, 5))
)
dev.off()