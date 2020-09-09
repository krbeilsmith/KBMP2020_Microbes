# KBMP2020_Ordination
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Summer 2020

##########################################################################################################################################################

# The commands for PERMANOVA and ordination below broadly follow the tutorials here: http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
# and here: https://bioconductor.riken.jp/packages/3.0/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html

# When using the distance and dissimilarities below, note that 16S data is compositional: https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full

# To see how types of samples vary in composition, eigendecomposition of the dissimilarity matrix is used to produce distances in a coordinate space. The
# first two axes of these space have the most distance between samples and their relative eigenvalues (eigenvalues divided by sum of all eigenvalues) are
# reported in the axes labels.

distance_matrices <- c('raup_Rare', 'bc_Rare', 'uni_Rare', 'wuni_Rare')

for(d in distance_matrices){
  distance_matrix <- get(d)
  print(d)
  # Ordinate samples with distance
  PCoA_dist  <- ordinate(R_physeq, "PCoA", distance=distance_matrix)
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
  
  pcoa <- plot_ordination(R_physeq, PCoA_dist, color="PlantPart") + 
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
  
  # Save plot to a file
  ggsave(paste(paste("~/Documents/KBMP2020_Microbes/Figures/PCoA_Plot_",d,sep=""),".tiff",sep=""), plot = pcoa, device = NULL, path = NULL,
         scale = 1.3, width = 5, height = 5, units = c("in"),
         dpi = 600, limitsize = TRUE)
  
  # Assign plot unique name for making full figure later
  assign(paste("PCoA_Plot",d,sep="_"),pcoa)
}

##########################################################################################################################################################

# Repeat the PERMANOVA for just samples at flowering in the phyllosphere

# Focus on plant samples for flowering phyllosphere tissues
R_physeq_flowering <- subset_samples(R_physeq,!(PlantPart %in% c("Soil","Roots")) & Stage%in%c("Flowering"))
R_physeq_flowering <- prune_taxa(taxa_sums(R_physeq_flowering)>0,R_physeq_flowering)

# PERMANOVA
flowering_permanova_results <- adonis(t(otu_table(R_physeq_flowering)) ~Year/PlantPart + Site + MiSeqRun/SamplePlate,data.frame(sample_data(R_physeq_flowering)),method="raup",permutations=999)
flowering_permanova_results

# Store results
write.table(flowering_permanova_results$aov.tab,file="~/Documents/KBMP2020_Microbes/PERMANOVAs/Flowering_raup_ASV", sep="\t", row.names=TRUE,quote=FALSE)

# Tissue still has an association with composition just within the phyllosphere at flowering:
#                       Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#Year:PlantPart         7   0.56329  0.080470  11.385  0.38022  0.010 ** 
  
##########################################################################################################################################################

# Calculate dissimilarity / distance

# Raup Crick
raup_Rare_f <- phyloseq::distance(R_physeq_flowering, method = "raup")

# Bray Curtis
bc_Rare_f <- phyloseq::distance(R_physeq_flowering, method = "bray")

# UniFrac
uni_Rare_f <- UniFrac(R_physeq_flowering, weighted=FALSE, normalized=TRUE, parallel=TRUE, fast=TRUE)

# Weighted UniFrac
wuni_Rare_f <- UniFrac(R_physeq_flowering, weighted=TRUE, normalized=TRUE, parallel=TRUE, fast=TRUE)

distance_matrices <- c('raup_Rare_f', 'bc_Rare_f', 'uni_Rare_f', 'wuni_Rare_f')

for(d in distance_matrices){
  distance_matrix <- get(d)
  print(d)
  # Ordinate samples with distance
  PCoA_dist  <- ordinate(R_physeq_flowering, "PCoA", distance=distance_matrix)
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
  Palette<-c("#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
  
  pcoa <- plot_ordination(R_physeq_flowering, PCoA_dist, color="PlantPart") + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))+
    font("ylab", size = 16, color = "black",face = "bold")+
    font("xlab", size = 16, color = "black",face = "bold")+
    font("xy.text", size = 12, color = "black", face = "bold")+
    geom_point(size=1.5)+
    scale_colour_manual(values=Palette,labels=c("Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+ 
    theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),
          legend.direction="horizontal")+
    theme(legend.background = element_rect(color="black", size=.5),legend.position="top")+
    guides(color=guide_legend(override.aes = list(size=5)))+
    coord_cartesian(ylim=c(-0.75,0.75),xlim=c(-0.75,0.75))
  
  pcoa <- pcoa  + labs(x=paste(paste("\nPCoA Axis 1 (",round(eig1,1),sep=""),"%)"), y=paste(paste("PCoA Axis 2 (",round(eig2,1),sep=""), "%)\n"))
  
  ggsave(paste(paste("~/Documents/KBMP2020_Microbes/Figures/f_PCoA_Plot_",d,sep=""),".tiff",sep=""), plot = pcoa, device = NULL, path = NULL,
         scale = 1.3, width = 5, height = 5, units = c("in"),
         dpi = 600, limitsize = TRUE)
  
  assign(paste("PCoA_Plot",d,sep="_"),pcoa)
}

##########################################################################################################################################################

# Make figure

legendCol <- cowplot::get_legend(PCoA_Plot_raup_Rare+guides(color=guide_legend(override.aes = list(size=5),nrow=1)))

ggsave("~/Documents/KBMP2020_Microbes/Figures/TissuePCoA.tiff", plot = grid.arrange(
  legendCol,
  PCoA_Plot_raup_Rare+ggtitle("A: Raup-Crick, all stages\n")+theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_bc_Rare+ggtitle("B: Bray-Curtis, all stages\n")+theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_uni_Rare+ggtitle("C: UniFrac, all stages\n")+theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_raup_Rare_f+ggtitle("D: Raup-Crick, flowering\n")+theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_bc_Rare_f+ggtitle("E: Bray-Curtis, flowering\n")+theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_uni_Rare_f+ggtitle("F: UniFrac, flowering\n")+theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  layout_matrix = rbind(c(NA,1,NA),
                        c(2,3,4),
                        c(5,6,7)),
  heights=c(1,4,4),widths=c(4,4,4)), device = NULL, path = NULL,
  scale = 1.5, width = 7.5, height = 5, units = c("in"),
  dpi = 600, limitsize = TRUE)

ggsave("~/Documents/KBMP2020_Microbes/Figures/TissuePCoA.png", plot = grid.arrange(
  legendCol,
  PCoA_Plot_raup_Rare+ggtitle("A: Raup-Crick, all stages\n")+theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_bc_Rare+ggtitle("B: Bray-Curtis, all stages\n")+theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_uni_Rare+ggtitle("C: UniFrac, all stages\n")+theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_raup_Rare_f+ggtitle("D: Raup-Crick, flowering\n")+theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_bc_Rare_f+ggtitle("E: Bray-Curtis, flowering\n")+theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  PCoA_Plot_uni_Rare_f+ggtitle("F: UniFrac, flowering\n")+theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  layout_matrix = rbind(c(NA,1,NA),
                        c(2,3,4),
                        c(5,6,7)),
  heights=c(1,4,4),widths=c(4,4,4)), device = NULL, path = NULL,
  scale = 1.5, width = 7.5, height = 5, units = c("in"),
  dpi = 600, limitsize = TRUE)
