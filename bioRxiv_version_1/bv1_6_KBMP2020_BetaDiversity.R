# KBMP2020_BetaDiversity
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Winter 2020

# Select only the plant samples with counts over 100
physeq_forBetaDiv <-subset_samples(physeq, PlantPart!="Soil")
physeq_forBetaDiv <- prune_samples(sample_sums(physeq_forBetaDiv)>100,physeq_forBetaDiv)

# Retain only the samples with at least 20 ASVs
Retain <- c()
for(SampID in sample_names(physeq_forBetaDiv)){
  print(SampID)
  SampTaxa <- subset_samples(physeq_forBetaDiv,sample_names(physeq_forBetaDiv) %in% SampID) # prune the data to just the focal sample
  SampTaxa <- prune_taxa(taxa_sums(SampTaxa)>0,SampTaxa) # prune the data to just the ASVs in the focal sample
  if(length(taxa_names(SampTaxa))>20){ # if more than 20 ASVs present
    Retain <- c(Retain, SampID) # add to list
  }
}
physeq_forBetaDiv <-prune_samples(sample_names(physeq_forBetaDiv) %in% Retain, physeq_forBetaDiv)

# Prune to roots and rosette leaves only
physeq_forBetaDiv <-subset_samples(physeq_forBetaDiv, PlantPart%in%c("Roots","RosLeaves"))
# Prune to ASVs remaining with nonzero counts
physeq_forBetaDiv <- prune_taxa(taxa_sums(physeq_forBetaDiv)>0,physeq_forBetaDiv)

# Set variables for looping
SampleStage <- c(levels(sample_data(physeq_forBetaDiv)$Stage)) # each stage by which the samples can be grouped
SampleStage <- SampleStage[!(SampleStage%in%c("Senescent"))] # remove senescent because no rosette leaves at this stage
SampleTissue <- c("Roots","RosLeaves") # each tissue type
SampleSite <- c("ME","WW") # each site

# Rarefy and perform UniFrac
R_physeq <- rarefy_even_depth(physeq_forBetaDiv,500)
uni_Rare <- UniFrac(R_physeq, weighted=FALSE, normalized=FALSE, parallel=TRUE, fast=TRUE)
UniFracDF <- as.data.frame(as.matrix(uni_Rare))

# Make a dataframe for the results
BetaDiv <- NULL
BetaDiv <- data.frame(matrix(ncol = 5, nrow = 0),stringsAsFactors = FALSE)
BetaDiv <- setNames(data.frame(BetaDiv), c("Stage","Tissue","Site","Distance","Group"))
NextRow <- NULL

for(stage in SampleStage){ # for each stage
  print(stage)
  for(tissue in SampleTissue){ # for each tissue
    print(tissue)
    for(site in SampleSite){ # for each site
      print(site)
      
      # Find the names of samples collected in the given stage, tissue and SAME site
      print("Finding samples...")
      samps_to_keep = as.character(get_variable(physeq_forBetaDiv, "PlantPart")) == tissue & as.character(get_variable(physeq_forBetaDiv, "Stage")) == stage & as.character(get_variable(physeq_forBetaDiv, "Site")) == site
      Sharing_Property = prune_samples(samps_to_keep, physeq_forBetaDiv)
      SelectSamples <- sample_names(Sharing_Property)
      
      # Prune the samples in the UniFrac results to those present in the SAME site
      print("Pruning to samples from same tissue...")
      SelectedDF <- UniFracDF[(row.names(UniFracDF) %in% SelectSamples),(row.names(UniFracDF) %in% SelectSamples)]
      UniFracMat <- as.matrix(SelectedDF)
      UniFracMat <- UniFracMat[upper.tri(UniFracMat,diag=FALSE)] # take the values above the diagonal, exclude the diagonal (0s for same samples)
      print("Adding to dataframe...")
      for(row in UniFracMat){
        for(i in row){
          NextRow <- NULL
          NextRow <- c(stage,tissue,site,i,"in")
          BetaDiv[nrow(BetaDiv)+1,] <- NextRow
        }
      }
      
      # Find the names of samples collected in the given stage, tissue and DIFFERENT site
      print("Finding samples...")
      samps_to_keep = as.character(get_variable(physeq_forBetaDiv, "Stage")) == stage & as.character(get_variable(physeq_forBetaDiv, "PlantPart")) == tissue
      Sharing_Stage = prune_samples(samps_to_keep, physeq_forBetaDiv)
      SelectSamples_Stage <- sample_names(Sharing_Stage)
      samps_to_keep = as.character(get_variable(physeq_forBetaDiv, "PlantPart")) == tissue & as.character(get_variable(physeq_forBetaDiv, "Stage")) == stage & as.character(get_variable(physeq_forBetaDiv, "Site")) == site
      Sharing_Property = prune_samples(samps_to_keep, physeq_forBetaDiv)
      SelectSamples <- sample_names(Sharing_Property)
      
      # Prune the samples in the UniFrac results to those present in a DIFFERENT site
      print("Pruning to samples from different tissue...")
      SelectedDF <- UniFracDF[(row.names(UniFracDF) %in% SelectSamples_Stage),(row.names(UniFracDF) %in% SelectSamples_Stage)]
      SelectedDF <- SelectedDF[(row.names(SelectedDF) %in% SelectSamples),!(row.names(SelectedDF) %in% SelectSamples)]
      UniFracMat <- as.matrix(SelectedDF)
      print("Adding to dataframe...")
      for(row in UniFracMat){
        for(i in row){
          NextRow <- NULL
          NextRow <- c(stage,tissue,site,i,"out")
          BetaDiv[nrow(BetaDiv)+1,] <- NextRow
        }
      }
    }
    
  }
}

# save results
write.table(BetaDiv,file="~/Documents/KBMP2020_Microbes/Outputs/BetaDiv_Site", sep="\t", row.names=TRUE,quote=FALSE)

# BetaDiv <- read.table("~/Documents/KBMP2020_Microbes/Outputs/BetaDiv_Site", sep="\t", header=TRUE)

# Make columns numeric or organize factor levels for plotting
BetaDiv$Distance <- as.numeric(BetaDiv$Distance)
BetaDiv$Stage <- factor(BetaDiv$Stage, levels = c("TwoLeaf","FourLeaf","SixLeaf","EightLeaf","Flowering"))
BetaDiv$Group <- factor(BetaDiv$Group, levels=c("in","out"))
BetaDiv$Tissue <- factor(BetaDiv$Tissue, levels=c("RosLeaves","Roots"))

Stage.labs <- c("Two Leaf","Four Leaf","Six Leaf","Eight Leaf","Flowering")
names(Stage.labs) <- c("TwoLeaf","FourLeaf","SixLeaf","EightLeaf","Flowering")

Tissue.labs <- c("Roots","Rosettes")
names(Tissue.labs) <- c("Roots","RosLeaves")

# Frequency plots
freqB <- ggplot(BetaDiv[BetaDiv$Stage!="TwoLeaf",], aes(Distance,stat(ndensity),fill=interaction(Group,Tissue),color=Tissue))+geom_density(alpha=0.2)+
  facet_grid(cols=vars(Stage),rows=vars(Tissue),switch="y",labeller=labeller(Stage=Stage.labs, Tissue=Tissue.labs))+
  coord_cartesian(xlim=c(0.1,1))+
  theme_bw()+
  scale_fill_manual(values=c("#56B4E9","white","#E69F00","white"),
                    labels=c("Rosettes within Site","Rosettes between Sites", "Roots within Site", "Roots between Sites"),name='Samples Compared')+
  scale_color_manual(values=c("#56B4E9","#E69F00"),labels=c("Rosettes","Roots"),name='Tissue',guide=FALSE)+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+ 
  font("xlab", size = 16, color = "black",face = "bold")+ 
  theme(strip.background=element_blank(),strip.text.x = element_text(size = 16, color = "black",face="bold"), 
        strip.text.y = element_text(size = 16, color = "black",face="bold",margin = margin(0,0,0,0.3, "cm"),vjust=3),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(strip.placement = "outside")+
  labs(x="\nUniFrac Distance (rarefied 500 counts)",y="Frequency (in 469 plant samples)\n")+ 
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5))+ guides(linetype=guide_legend(keywidth = 5, keyheight = 1))

############################################################################################################################################################
############################################################################################################################################################

# Make a dataframe for the results.
BetaDiv <- NULL
BetaDiv <- data.frame(matrix(ncol = 4, nrow = 0),stringsAsFactors = FALSE)
BetaDiv <- setNames(data.frame(BetaDiv), c("Stage","Tissue","Distance","Group"))
NextRow <- NULL

for(stage in SampleStage){ # for each stage
  print(stage)
  for(tissue in SampleTissue){ # for each tissue
    print(tissue)
    
    # Find the names of samples collected in the given stage and SAME tissue type
    print("Finding samples...")
    samps_to_keep = as.character(get_variable(physeq_forBetaDiv, "PlantPart")) == tissue & as.character(get_variable(physeq_forBetaDiv, "Stage")) == stage
    Sharing_Property = prune_samples(samps_to_keep, physeq_forBetaDiv)
    SelectSamples <- sample_names(Sharing_Property)
    
    # Prune the samples in the UniFrac results to those present in the SAME tissue at the same stage
    print("Pruning to samples from same tissue...")
    SelectedDF <- UniFracDF[(row.names(UniFracDF) %in% SelectSamples),(row.names(UniFracDF) %in% SelectSamples)]
    UniFracMat <- as.matrix(SelectedDF)
    UniFracMat <- UniFracMat[upper.tri(UniFracMat,diag=FALSE)] # take the values above the diagonal, exclude the diagonal (0s for same samples)
    print("Adding to dataframe...")
    for(row in UniFracMat){
      for(i in row){
        NextRow <- NULL
        NextRow <- c(stage,tissue,i,"in")
        BetaDiv[nrow(BetaDiv)+1,] <- NextRow
      }
    }
    
    # Find the names of samples collected in the given stage and DIFFERENT tissue type
    print("Finding samples...")
    samps_to_keep = as.character(get_variable(physeq_forBetaDiv, "Stage")) == stage
    Sharing_Stage = prune_samples(samps_to_keep, physeq_forBetaDiv)
    SelectSamples_Stage <- sample_names(Sharing_Stage)
    samps_to_keep = as.character(get_variable(physeq_forBetaDiv, "PlantPart")) == tissue & as.character(get_variable(physeq_forBetaDiv, "Stage")) == stage
    Sharing_Property = prune_samples(samps_to_keep, physeq_forBetaDiv)
    SelectSamples <- sample_names(Sharing_Property)
    
    # Prune the samples in the UniFrac results to those present in a DIFFERENT tissue at the same stage
    print("Pruning to samples from different tissue...")
    SelectedDF <- UniFracDF[(row.names(UniFracDF) %in% SelectSamples_Stage),(row.names(UniFracDF) %in% SelectSamples_Stage)]
    SelectedDF <- SelectedDF[(row.names(SelectedDF) %in% SelectSamples),!(row.names(SelectedDF) %in% SelectSamples)]
    UniFracMat <- as.matrix(SelectedDF)
    print("Adding to dataframe...")
    for(row in UniFracMat){
      for(i in row){
        NextRow <- NULL
        NextRow <- c(stage,tissue,i,"out")
        BetaDiv[nrow(BetaDiv)+1,] <- NextRow
      }
    }
  }
}

# save results
write.table(BetaDiv,file="~/Documents/KBMP2020_Microbes/Outputs/BetaDiv_Tissue", sep="\t", row.names=TRUE,quote=FALSE)

# BetaDiv <- read.table("~/Documents/KBMP2020_Microbes/Outputs/BetaDiv_Tissue", sep="\t", header=TRUE)

# Make columns numeric or organize factor levels for plotting
BetaDiv$Distance <- as.numeric(BetaDiv$Distance)
BetaDiv$Stage <- factor(BetaDiv$Stage, levels = c("TwoLeaf","FourLeaf","SixLeaf","EightLeaf","Flowering"))
BetaDiv$Group <- factor(BetaDiv$Group, levels=c("in","out"))
BetaDiv$Tissue <- factor(BetaDiv$Tissue, levels=c("RosLeaves","Roots"))


# Frequency plots
freqA <- ggplot(BetaDiv[BetaDiv$Stage!="TwoLeaf",], aes(Distance,stat(ndensity),fill=Group))+geom_density(alpha=0.2)+
  facet_grid(cols=vars(Stage),rows=vars(Tissue),labeller=labeller(Stage=Stage.labs, Tissue=Tissue.labs),switch="y")+
  coord_cartesian(xlim=c(0.25,1))+
  theme_bw()+
  theme(strip.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=c("black","white"),labels=c("Same Tissue","Different Tissue"),name='Samples Compared')+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+ 
  font("xlab", size = 16, color = "black",face = "bold")+ 
  theme(strip.background=element_blank(),strip.text.x = element_text(size = 16, color = "black",face="bold"), 
        strip.text.y = element_text(size = 16, color = "black",face="bold",margin = margin(0,0,0,0.3, "cm"),vjust=3),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(strip.placement = "outside")+
  labs(x="\nUniFrac Distance (rarefied 500 counts)",y="Frequency (in 469 plant samples)\n")+ 
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5))+ guides(linetype=guide_legend(keywidth = 5, keyheight = 1))

# Plot together
ggsave("~/Documents/KBMP2020_Microbes/Figures/BetaDiv.tiff", plot = grid.arrange(
  freqA+ggtitle("A\n")+theme(plot.title = element_text(face="bold",size=20)),
  freqB+ggtitle("B\n")+theme(plot.title = element_text(face="bold",size=20)),
  layout_matrix = rbind(c(1),
                        c(2))), device = NULL, path = NULL,
  scale = 2, width = 5, height = 5, units = c("in"),
  dpi = 600, limitsize = TRUE)

########################################################################################################################################################

# Tissue effect over time
# 1. Vegetative vs Flowering Roots and Rosettes
# 2. Flowering vs Senescent Roots and Stems

Roots_Rosettes_Veg <- subset_samples(physeq_factored_plants, PlantPart %in% c("Roots","RosLeaves") & Stage %in% c("SixLeaf","EightLeaf"))

PruneList_Veg <- prune_taxa(taxa_sums(Roots_Rosettes_Veg)>0,Roots_Rosettes_Veg)
PruneList_Veg <- taxa_names(PruneList_Veg)

Roots_Rosettes_Flo <- subset_samples(physeq_factored_plants, PlantPart %in% c("Roots","RosLeaves") & Stage %in% c("Flowering"))

PruneList_Flo <- prune_taxa(taxa_sums(Roots_Rosettes_Flo)>0,Roots_Rosettes_Flo)
PruneList_Flo <- taxa_names(PruneList_Flo)

Roots_Rosettes_Veg <- prune_taxa(taxa_names(Roots_Rosettes_Veg) %in% PruneList_Flo, Roots_Rosettes_Veg)
Roots_Rosettes_Veg <- prune_samples(sample_sums(Roots_Rosettes_Veg)>0, Roots_Rosettes_Veg)
Roots_Rosettes_Veg <- prune_taxa(taxa_sums(Roots_Rosettes_Veg)>0,Roots_Rosettes_Veg)

Roots_Rosettes_Flo <- prune_taxa(taxa_names(Roots_Rosettes_Flo) %in% PruneList_Veg, Roots_Rosettes_Flo)
Roots_Rosettes_Flo <- prune_samples(sample_sums(Roots_Rosettes_Flo)>0, Roots_Rosettes_Flo)
Roots_Rosettes_Flo <- prune_taxa(taxa_sums(Roots_Rosettes_Flo)>0,Roots_Rosettes_Flo)

Roots_Rosettes_Veg
Roots_Rosettes_Veg_permanova_results <- adonis(t(otu_table(Roots_Rosettes_Veg)) ~Year/PlantPart + MiSeqRun/SamplePlate + Site,data.frame(sample_data(Roots_Rosettes_Veg)),method="raup",permutations=999)
Roots_Rosettes_Veg_permanova_results

write.table(data.frame(Roots_Rosettes_Veg_permanova_results$aov.tab[,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")]),file="~/Documents/KBMP2020_Microbes/PERMANOVAs/Vegetative_RootsRosettes_Raup", sep="\t", row.names=TRUE,quote=FALSE)

Roots_Rosettes_Flo
Roots_Rosettes_Flo_permanova_results <- adonis(t(otu_table(Roots_Rosettes_Flo)) ~Year/PlantPart + MiSeqRun/SamplePlate + Site,data.frame(sample_data(Roots_Rosettes_Flo)),method="raup",permutations=999)
Roots_Rosettes_Flo_permanova_results

write.table(data.frame(Roots_Rosettes_Flo_permanova_results$aov.tab[,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")]),file="~/Documents/KBMP2020_Microbes/PERMANOVAs/Flowering_RootsRosettes_Raup", sep="\t", row.names=TRUE,quote=FALSE)

####################################################################################################################################################

Roots_Stems_Flo <- subset_samples(physeq_factored_plants, PlantPart %in% c("Roots","Stems") & Stage %in% c("Flowering"))

PruneList_Flo <- prune_taxa(taxa_sums(Roots_Stems_Flo)>0,Roots_Stems_Flo)
PruneList_Flo <- taxa_names(PruneList_Flo)

Roots_Stems_Sen <- subset_samples(physeq_factored_plants, PlantPart %in% c("Roots","Stems") & Stage %in% c("Senescent"))

PruneList_Sen <- prune_taxa(taxa_sums(Roots_Stems_Sen)>0,Roots_Stems_Sen)
PruneList_Sen <- taxa_names(PruneList_Sen)

Roots_Stems_Flo <- prune_taxa(taxa_names(Roots_Stems_Flo) %in% PruneList_Sen, Roots_Stems_Flo)
Roots_Stems_Flo <- prune_samples(sample_sums(Roots_Stems_Flo)>0, Roots_Stems_Flo)
Roots_Stems_Flo <- prune_taxa(taxa_sums(Roots_Stems_Flo)>0, Roots_Stems_Flo)

Roots_Stems_Sen <- prune_taxa(taxa_names(Roots_Stems_Sen) %in% PruneList_Flo, Roots_Stems_Sen)
Roots_Stems_Sen <- prune_samples(sample_sums(Roots_Stems_Sen)>0, Roots_Stems_Sen)
Roots_Stems_Sen <- prune_taxa(taxa_sums(Roots_Stems_Sen)>0, Roots_Stems_Sen)

Roots_Stems_Flo
Roots_Stems_Flo_permanova_results <- adonis(t(otu_table(Roots_Stems_Flo)) ~Year/PlantPart + MiSeqRun/SamplePlate + Site,data.frame(sample_data(Roots_Stems_Flo)),method="raup",permutations=999)
Roots_Stems_Flo_permanova_results

write.table(data.frame(Roots_Stems_Flo_permanova_results$aov.tab[,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")]),file="~/Documents/KBMP2020_Microbes/PERMANOVAs/Flowering_RootsStems_Raup", sep="\t", row.names=TRUE,quote=FALSE)

Roots_Stems_Sen
Roots_Stems_Sen_permanova_results <- adonis(t(otu_table(Roots_Stems_Sen)) ~Year/PlantPart + MiSeqRun/SamplePlate + Site,data.frame(sample_data(Roots_Stems_Sen)),method="raup",permutations=999)
Roots_Stems_Sen_permanova_results

write.table(data.frame(Roots_Stems_Sen_permanova_results$aov.tab[,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")]),file="~/Documents/KBMP2020_Microbes/PERMANOVAs/Senescent_RootsStems_Raup", sep="\t", row.names=TRUE,quote=FALSE)
