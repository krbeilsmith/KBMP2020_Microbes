# KBMP2020_PartitionVariance
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Summer 2020

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

# Calculate dissimilarity / distance

# Raup Crick
raup_Rare <- phyloseq::distance(R_physeq, method = "raup")

# Bray Curtis
bc_Rare <- phyloseq::distance(R_physeq, method = "bray")

# UniFrac
uni_Rare <- UniFrac(R_physeq, weighted=FALSE, normalized=TRUE, parallel=TRUE, fast=TRUE)

# Weighted UniFrac
wuni_Rare <- UniFrac(R_physeq, weighted=TRUE, normalized=TRUE, parallel=TRUE, fast=TRUE)

distance_matrices <- c('raup_Rare', 'bc_Rare', 'uni_Rare', 'wuni_Rare')

for(d in distance_matrices){
  distance_matrix <- get(d)
  full_permanova_results <- adonis(distance_matrix ~ Year/Stage/PlantPart + MiSeqRun/SamplePlate + Site,data.frame(sample_data(R_physeq)),permutations=999)
  # print(permanova_results)
  write.table(data.frame(full_permanova_results$aov.tab[,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")]),file=paste("~/Documents/KBMP2020_Microbes/PERMANOVAs/Plant_ASV_FullModelResult",d,sep="_"), sep="\t", row.names=TRUE,quote=FALSE)
}

# To see how types of samples vary in composition, eigendecomposition of the dissimilarity matrix is used to produce distances in a coordinate space. The
# first two axes of these space have the most distance between samples and their relative eigenvalues (eigenvalues divided by sum of all eigenvalues) are
# reported in the axes labels.

distance_matrices <- c('raup_Rare', 'bc_Rare', 'uni_Rare', 'wuni_Rare')

d='raup_Rare'

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
    #scale_shape_manual(values=c("\u25AC","\u25D6","\u25A0","\u25B2","\u25AE","\u25C6"),labels=c('Two Leaf','Four Leaf','Six Leaf','Eight Leaf','Flowering','Senescent'),name='Development Stage')+
    #scale_shape_manual(values=c("\u25CF","\u25CF","\u25CF","\u25CF","\u25FC","\u25C6"),labels=c('Two Leaf','Four Leaf','Six Leaf','Eight Leaf','Flowering','Senescent'),name='Development Stage')+
    theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),legend.direction="horizontal")+
    theme(legend.background = element_rect(color="black", size=.5),legend.position="top")+
    guides(color=guide_legend(override.aes = list(size=5)))+
    coord_cartesian(ylim=c(-0.75,0.75),xlim=c(-0.75,0.75))
  
  pcoa <- pcoa  + labs(x=paste(paste("\nPCoA Axis 1 (",round(eig1,1),sep=""),"%)"), y=paste(paste("PCoA Axis 2 (",round(eig2,1),sep=""), "%)\n"))
  
  ggsave(paste(paste("~/Documents/KBMP2020_Microbes/Figures/PCoA_Plot_",d,sep=""),".tiff",sep=""), plot = pcoa, device = NULL, path = NULL,
         scale = 1.3, width = 5, height = 5, units = c("in"),
         dpi = 600, limitsize = TRUE)
  
  assign(paste("PCoA_Plot",d,sep="_"),pcoa)
}

PCoAs <- c('PCoA_raup_Rare','PCoA_bc_Rare','PCoA_uni_Rare','PCoA_wuni_Rare')

ggsave("~/Documents/KBMP2020_Microbes/Figures/TissuePCoA.tiff", plot = grid.arrange(
  PCoA_Plot_raup_Rare+ggtitle("A\n")+theme(plot.title = element_text(face="bold",size=20)),
  PCoA_Plot_bc_Rare+ggtitle("B\n")+theme(plot.title = element_text(face="bold",size=20)),
  PCoA_Plot_uni_Rare+ggtitle("C\n")+theme(plot.title = element_text(face="bold",size=20)),
  layout_matrix = rbind(c(1,2,3))), device = NULL, path = NULL,
  scale = 2, width = 10, height = 3, units = c("in"),
  dpi = 600, limitsize = TRUE)

# Focus on plant samples for flowering phyllosphere tissues
R_physeq_flowering <- subset_samples(R_physeq,!(PlantPart %in% c("Soil","Roots")) & Stage%in%c("Flowering"))
R_physeq_flowering <- prune_taxa(taxa_sums(R_physeq_flowering)>0,R_physeq_flowering)

flowering_permanova_results <- adonis(t(otu_table(R_physeq_flowering)) ~Year/PlantPart + MiSeqRun/SamplePlate + Site,data.frame(sample_data(R_physeq_flowering)),method="raup",permutations=999)
flowering_permanova_results

# Tissue still has an effect just within the phyllosphere:
# Year:PlantPart         8     3.054 0.38174  1.9659 0.04085  0.010 **

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
    #scale_shape_manual(values=c("\u25AC","\u25D6","\u25A0","\u25B2","\u25AE","\u25C6"),labels=c('Two Leaf','Four Leaf','Six Leaf','Eight Leaf','Flowering','Senescent'),name='Development Stage')+
    #scale_shape_manual(values=c("\u25CF","\u25CF","\u25CF","\u25CF","\u25FC","\u25C6"),labels=c('Two Leaf','Four Leaf','Six Leaf','Eight Leaf','Flowering','Senescent'),name='Development Stage')+
    theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),legend.direction="horizontal")+
    theme(legend.background = element_rect(color="black", size=.5),legend.position="top")+
    guides(color=guide_legend(override.aes = list(size=5)))+
    coord_cartesian(ylim=c(-0.75,0.75),xlim=c(-0.75,0.75))
  
  pcoa <- pcoa  + labs(x=paste(paste("\nPCoA Axis 1 (",round(eig1,1),sep=""),"%)"), y=paste(paste("PCoA Axis 2 (",round(eig2,1),sep=""), "%)\n"))
  
  ggsave(paste(paste("~/Documents/KBMP2020_Microbes/Figures/f_PCoA_Plot_",d,sep=""),".tiff",sep=""), plot = pcoa, device = NULL, path = NULL,
         scale = 1.3, width = 5, height = 5, units = c("in"),
         dpi = 600, limitsize = TRUE)
  
  assign(paste("PCoA_Plot",d,sep="_"),pcoa)
}

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
  scale = 1.5, width = 9, height = 6, units = c("in"),
  dpi = 600, limitsize = TRUE)
