# KBMP2020_EffectRobustness
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
# Approach 1: Rarefy data and use Bray-Curtis dissimilarity on abundance matrix
######################################################################################################################################

# Rarefy samples to even depth (will cause samples to drop from set)
R_physeq <- rarefy_even_depth(physeq_factored_plants,1000)

# Bray distance for Rarefied Set
bc_Rare <- phyloseq::distance(R_physeq, method = "bray")

# PERMANOVA with Bray distance for Rarefied Set
perma_bc_Rare <- adonis(bc_Rare ~ Year/Stage/PlantPart + MiSeqRun/SamplePlate + Site,data.frame(sample_data(R_physeq)),permutations=999)
perma_bc_Rare

write.table(data.frame(perma_bc_Rare$aov.tab[,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")]),file="~/Documents/KBMP2020_Microbes/PERMANOVAs/Plant_ASV_FullModelResult_Bray", sep="\t", row.names=TRUE,quote=FALSE)

# Ordinate with Bray distance for for Rarefied Set
PCoA_bc_Rare  <- ape::pcoa(bc_Rare,correction="cailliez")

# Variance explained
PCoA_bc_Rare$values$Eigenvalues[1]/sum(PCoA_bc_Rare$values$Eigenvalues)
PCoA_bc_Rare$values$Eigenvalues[2]/sum(PCoA_bc_Rare$values$Eigenvalues)

p <- plot_ordination(R_physeq, PCoA_bc_Rare, color="PlantPart") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  geom_point(size=1.5)+
  scale_colour_manual(values=Palette,labels=c("Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+ 
  #scale_shape_manual(values=c("\u25CF","\u25AC","\u25D6","\u25A0","\u25B2","\u25AE","\u25C6"),labels=c('No Plant','Two Leaf','Four Leaf','Six Leaf','Eight Leaf','Flowering','Senescent'),name='Development Stage')
  theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),legend.direction="horizontal")+
  theme(legend.background = element_rect(color="black", size=.5),legend.position="top")+
  guides(color=guide_legend(override.aes = list(size=5)))+
  coord_cartesian(ylim=c(-0.75,0.75),xlim=c(-0.75,0.75))

RareBrayPCoA <- p + labs(x="\nPCoA Axis 1 (12.0%)", y="PCoA Axis 2 (6.5%)\n")

ggsave("~/Documents/KBMP2020_Microbes/Figures/Bray_PCoA.tiff", plot = RareBrayPCoA, device = NULL, path = NULL,
       scale = 1.3, width = 5, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)

######################################################################################################################################
# Approach 2: Rarefy data and use UNIFRAC distance on abundance matrix
######################################################################################################################################

uni_Rare <- UniFrac(R_physeq, weighted=FALSE, normalized=TRUE, parallel=TRUE, fast=TRUE)

perma_uni_Rare <- adonis(uni_Rare ~ Year/Stage/PlantPart + MiSeqRun/SamplePlate + Site,data.frame(sample_data(R_physeq)),permutations=999)
perma_uni_Rare

write.table(data.frame(perma_uni_Rare$aov.tab[,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")]),file="~/Documents/KBMP2020_Microbes/PERMANOVAs/Plant_ASV_FullModelResult_Uni", sep="\t", row.names=TRUE,quote=FALSE)

PCoA_uni_Rare  <- ape::pcoa(uni_Rare,correction="cailliez")

# Variance explained
PCoA_uni_Rare$values$Eigenvalues[1]/sum(PCoA_uni_Rare$values$Eigenvalues)
PCoA_uni_Rare$values$Eigenvalues[2]/sum(PCoA_uni_Rare$values$Eigenvalues)

p <- plot_ordination(R_physeq, PCoA_uni_Rare, color="PlantPart") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  geom_point(size=1.5)+
  scale_colour_manual(values=Palette,labels=c("Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+ 
  #scale_shape_manual(values=c("\u25CF","\u25AC","\u25D6","\u25A0","\u25B2","\u25AE","\u25C6"),labels=c('No Plant','Two Leaf','Four Leaf','Six Leaf','Eight Leaf','Flowering','Senescent'),name='Development Stage')
  theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),legend.direction="horizontal")+
  theme(legend.background = element_rect(color="black", size=.5),legend.position="top")+
  guides(color=guide_legend(override.aes = list(size=5)))+
  coord_cartesian(ylim=c(-0.75,0.75),xlim=c(-0.75,0.75))

RareUniPCoA <- p + labs(x="\nPCoA Axis 1 (17.8%)", y="PCoA Axis 2 (8.4%)\n")

ggsave("~/Documents/KBMP2020_Microbes/Figures/Uni_PCoA.tiff", plot = RareUniPCoA, device = NULL, path = NULL,
       scale = 1.3, width = 5, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)

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

N_physeq <- prune_samples(sample_names(N_physeq)%in%sample_names(physeq_factored_plants), N_physeq)

bc_VST <- phyloseq::distance(N_physeq, "bray")

PCoA_bc_VST  <- ordinate(N_physeq, "PCoA", distance=bc_VST)

p <- plot_ordination(N_physeq, PCoA_bc_VST, color="PlantPart") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  geom_point(size=1.5)+
  scale_colour_manual(values=Palette,labels=c("Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+ 
  #scale_shape_manual(values=c("\u25CF","\u25AC","\u25D6","\u25A0","\u25B2","\u25AE","\u25C6"),labels=c('No Plant','Two Leaf','Four Leaf','Six Leaf','Eight Leaf','Flowering','Senescent'),name='Development Stage')
  theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),legend.direction="horizontal")+
  theme(legend.background = element_rect(color="black", size=.5),legend.position="top")+
  guides(color=guide_legend(override.aes = list(size=5)))+
  coord_cartesian(ylim=c(-0.75,0.75),xlim=c(-0.75,0.75))

VSTBrayPCoA <-p + labs(x="\nPCoA Axis 1 (7.9%)", y="PCoA Axis 2 (5.6%)\n")

######################################################################################################################################
# Approach 4: Prune rare taxa (present in fewer than three samples) and Raup-Crick dissimilarity
######################################################################################################################################

physeq_factored_plants_PA <- transform_sample_counts(physeq_factored_plants, function(x) ifelse(x>0,1,0))
physeq_factored_plants_PA <- prune_taxa(taxa_sums(physeq_factored_plants_PA)>2, physeq_factored_plants_PA)
TaxaToKeep <- taxa_names(physeq_factored_plants_PA)

# Now we prune the data to just the taxa present in >3 samples.
PrevalentTaxa <- prune_taxa(taxa_names(physeq_factored_plants) %in% TaxaToKeep, physeq_factored_plants)
PrevalentTaxa <- prune_samples(sample_sums(PrevalentTaxa)>0,PrevalentTaxa)

# Raup distance for All Samples
raup_PrevalentTaxa <- phyloseq::distance(PrevalentTaxa, method = "raup")

# Ordinate 
PCoA_raup_PrevalentTaxa  <- ape::pcoa(raup_PrevalentTaxa,correction="lingoes")

# Variance explained
PCoA_raup_PrevalentTaxa$values$Eigenvalues[1]/sum(PCoA_raup_PrevalentTaxa$values$Eigenvalues)
PCoA_raup_PrevalentTaxa$values$Eigenvalues[2]/sum(PCoA_raup_PrevalentTaxa$values$Eigenvalues)

p <- plot_ordination(PrevalentTaxa, PCoA_raup_PrevalentTaxa, color="PlantPart") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  geom_point(size=1.5)+
  scale_colour_manual(values=Palette,labels=c("Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+ 
  #scale_shape_manual(values=c("\u25CF","\u25AC","\u25D6","\u25A0","\u25B2","\u25AE","\u25C6"),labels=c('No Plant','Two Leaf','Four Leaf','Six Leaf','Eight Leaf','Flowering','Senescent'),name='Development Stage')
  theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),legend.direction="horizontal")+
  theme(legend.background = element_rect(color="black", size=.5),legend.position="top")+
  guides(color=guide_legend(override.aes = list(size=5)))+
  coord_cartesian(ylim=c(-0.75,0.75),xlim=c(-0.75,0.75))

Raup_Prev_PCoA <- p + labs(x="\nPCoA Axis 1 (47.1%)", y="PCoA Axis 2 (27.0%)\n")

######################################################################################################################################
# Approach 5: Prune scarce taxa (fewer than 1000 total observations) and Raup-Crick dissimilarity
######################################################################################################################################

AbundantTaxa <- prune_taxa(taxa_sums(physeq_factored_plants)>1000, physeq_factored_plants)
AbundantTaxa <- prune_samples(sample_sums(AbundantTaxa)>0,AbundantTaxa)

# Raup distance for All Samples
raup_AbundantTaxa <- phyloseq::distance(AbundantTaxa, method = "raup")

# Ordinate 
PCoA_raup_AbundantTaxa  <- ape::pcoa(raup_AbundantTaxa,correction="cailliez")

# Variance explained
PCoA_raup_AbundantTaxa$values$Eigenvalues[1]/sum(PCoA_raup_AbundantTaxa$values$Eigenvalues)
PCoA_raup_AbundantTaxa$values$Eigenvalues[2]/sum(PCoA_raup_AbundantTaxa$values$Eigenvalues)

p <- plot_ordination(AbundantTaxa, PCoA_raup_AbundantTaxa, color="PlantPart") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  geom_point(size=1.5)+
  scale_colour_manual(values=Palette,labels=c("Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+ 
  #scale_shape_manual(values=c("\u25CF","\u25AC","\u25D6","\u25A0","\u25B2","\u25AE","\u25C6"),labels=c('No Plant','Two Leaf','Four Leaf','Six Leaf','Eight Leaf','Flowering','Senescent'),name='Development Stage')
  theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),legend.direction="horizontal")+
  theme(legend.background = element_rect(color="black", size=.5),legend.position="top")+
  guides(color=guide_legend(override.aes = list(size=5)))+
  coord_cartesian(ylim=c(-0.75,0.75),xlim=c(-0.75,0.75))

Raup_Abun_PCoA <- p + labs(x="\nPCoA Axis 1 (48.3%)", y="PCoA Axis 2 (24.8%)\n")

######################################################################################################################################
# Approach 6: TMM normalization and Bray-Curtis dissimilarity
######################################################################################################################################

Lphyseq <- physeq_factored_plants

# Add pseudocount to avoid zeros
otu_table(Lphyseq) <- otu_table(Lphyseq) + 0.01

# Replace OTU table with normalized counts
NormFactors <- calcNormFactors(otu_table(Lphyseq))
NormTable<-(otu_table(Lphyseq)%*%diag(NormFactors))
colnames(NormTable)<-sample_names(Lphyseq)
otu_table(Lphyseq) <- otu_table(NormTable, taxa_are_rows = TRUE)

# write.table(otu_table(Lphyseq),file="~/Documents/KBMP2020_Microbes/Outputs/EdgeR_Physeq_ASV", sep="\t", row.names=TRUE,quote=FALSE)

bc_Lphyseq <- phyloseq::distance(Lphyseq, "bray")

PCoA_bc_Lphyseq  <- ape::pcoa(bc_Lphyseq,correction="cailliez")

# Variance explained
PCoA_bc_Lphyseq$values$Eigenvalues[1]/sum(PCoA_bc_Lphyseq$values$Eigenvalues)
PCoA_bc_Lphyseq$values$Eigenvalues[2]/sum(PCoA_bc_Lphyseq$values$Eigenvalues)

p <- plot_ordination(Lphyseq, PCoA_bc_Lphyseq, color="PlantPart") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  geom_point(size=1.5)+
  scale_colour_manual(values=Palette,labels=c("Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+ 
  #scale_shape_manual(values=c("\u25CF","\u25AC","\u25D6","\u25A0","\u25B2","\u25AE","\u25C6"),labels=c('No Plant','Two Leaf','Four Leaf','Six Leaf','Eight Leaf','Flowering','Senescent'),name='Development Stage')
  theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),legend.direction="horizontal")+
  theme(legend.background = element_rect(color="black", size=.5),legend.position="top")+
  guides(color=guide_legend(override.aes = list(size=5)))+
  coord_cartesian(ylim=c(-0.75,0.75),xlim=c(-0.75,0.75))

TMMBrayPCoA <-p + labs(x="\nPCoA Axis 1 (13.0%)", y="PCoA Axis 2 (8.5%)\n")

##########################################################################################################################################
legendCol <- cowplot::get_legend(p)
tiff(file = "~/Documents/KBMP2020_Microbes/Figures/PCoA_Robustness.tiff", width = 6800, height = 8800, units = "px", res = 800)
grid.arrange(
  legendCol,
  RareUniPCoA+ theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("A: Rarefied, UniFrac\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")), 
  RareBrayPCoA+ theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("B: Rarefied, Bray-Curtis\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  VSTBrayPCoA+ theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("C: VST, Bray-Curtis\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  TMMBrayPCoA+ theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("D: TMM, Bray-Curtis\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  Raup_Prev_PCoA+ theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("E: Prevalence >3, Raup-Crick\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  Raup_Abun_PCoA+ theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("F: Count >1000, Raup-Crick\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
  widths = c(4, 4),
  heights = c(1,4,4,4),
  layout_matrix = rbind(c(1),
                        c(2, 3),
                        c(4, 5),
                        c(6, 7))
)
dev.off()
