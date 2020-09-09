# KBMP2020_BetaDiversity_revised
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Summer 2020

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
SampleStage <- SampleStage[!(SampleStage%in%c("TwoLeaf"))] # remove 2L because not enough samples after rarefaction

SampleTissue <- c("Roots","RosLeaves") # each tissue type
SampleSite <- c("ME","WW") # each site

# Rarefy and get distance matrices
R_physeq <- rarefy_even_depth(physeq_forBetaDiv,500)
uni_Rare <- UniFrac(R_physeq, weighted=FALSE, normalized=FALSE, parallel=TRUE, fast=TRUE)
raup_Rare <- phyloseq::distance(R_physeq, method = "raup")
bray_Rare <- phyloseq::distance(R_physeq, method = "bray")

distance_list <- c("uni_Rare","raup_Rare","bray_Rare")

# Make a dataframe for the results.
BetaDiv <- NULL
BetaDiv <- data.frame(matrix(ncol = 5, nrow = 0),stringsAsFactors = FALSE)
BetaDiv <- setNames(data.frame(BetaDiv), c("Stage","Tissue","Distance","Group","Metric"))
NextRow <- NULL

# For each distance matrix...
for(d in distance_list){
  print(d)
  dist_Rare <- get(d)
  DistDF <- as.data.frame(as.matrix(dist_Rare))
  
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
      SelectedDF <- DistDF[(row.names(DistDF) %in% SelectSamples),(row.names(DistDF) %in% SelectSamples)]
      UniFracMat <- as.matrix(SelectedDF)
      UniFracMat <- UniFracMat[upper.tri(UniFracMat,diag=FALSE)] # take the values above the diagonal, exclude the diagonal (0s for same samples)
      print("Adding to dataframe...")
      for(row in UniFracMat){
        for(i in row){
          NextRow <- NULL
          NextRow <- c(stage,tissue,i,"in",d)
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
      SelectedDF <- DistDF[(row.names(DistDF) %in% SelectSamples_Stage),(row.names(DistDF) %in% SelectSamples_Stage)]
      SelectedDF <- SelectedDF[(row.names(SelectedDF) %in% SelectSamples),!(row.names(SelectedDF) %in% SelectSamples)]
      UniFracMat <- as.matrix(SelectedDF)
      print("Adding to dataframe...")
      for(row in UniFracMat){
        for(i in row){
          NextRow <- NULL
          NextRow <- c(stage,tissue,i,"out",d)
          BetaDiv[nrow(BetaDiv)+1,] <- NextRow
        }
      }
    }
  }
  
  # save results
  write.table(BetaDiv,file="~/Documents/KBMP2020_Microbes/Outputs/BetaDiv_Tissue_MultiDist", sep="\t", row.names=TRUE,quote=FALSE)
}

# Make columns numeric or organize factor levels for plotting
BetaDiv$Distance <- as.numeric(BetaDiv$Distance)
BetaDiv$Stage <- factor(BetaDiv$Stage, levels = c("FourLeaf","SixLeaf","EightLeaf","Flowering"))
BetaDiv$Group <- factor(BetaDiv$Group, levels=c("in","out"))
BetaDiv$Tissue <- factor(BetaDiv$Tissue, levels=c("RosLeaves","Roots"))

Stage.labs <- c("Two Leaf","Four Leaf","Six Leaf","Eight Leaf","Flowering")
names(Stage.labs) <- c("TwoLeaf","FourLeaf","SixLeaf","EightLeaf","Flowering")

Tissue.labs <- c("Roots","Rosettes")
names(Tissue.labs) <- c("Roots","RosLeaves")

# Frequency plots
raup_plot <- ggplot(BetaDiv[BetaDiv$Tissue=="RosLeaves" & BetaDiv$Metric=="raup_Rare",], 
                    aes(log10(Distance),stat(ndensity),fill=Group))+geom_density(alpha=0.2)+
  facet_grid(cols=vars(Stage),rows=vars(Tissue),labeller=labeller(Stage=Stage.labs, Tissue=Tissue.labs),switch="y")+
  coord_cartesian(xlim=c(-15,0))+
  theme_bw()+
  theme(strip.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=c("black","white"),labels=c("Rosettes vs. Rosettes","Rosettes vs. Roots"),name='Samples Compared')+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+ 
  font("xlab", size = 16, color = "black",face = "bold")+ 
  theme(strip.background=element_blank(),strip.text.x = element_text(size = 16, color = "black",face="bold"), 
        strip.text.y = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(strip.placement = "outside")+
  labs(x=bquote(bold(atop("\n",~log[10]~"Raup-Crick Index (rarefied 500 counts)"))),y=" \n")+ 
  theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5))+ guides(linetype=guide_legend(keywidth = 5, keyheight = 1))

bray_plot <- ggplot(BetaDiv[BetaDiv$Tissue=="RosLeaves" & BetaDiv$Metric=="bray_Rare",], 
                    aes(Distance,stat(ndensity),fill=Group))+geom_density(alpha=0.2)+
  facet_grid(cols=vars(Stage),rows=vars(Tissue),labeller=labeller(Stage=Stage.labs, Tissue=Tissue.labs),switch="y")+
  coord_cartesian(xlim=c(0,1))+
  theme_bw()+
  theme(strip.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=c("black","white"),labels=c("Rosettes vs. Rosettes","Rosettes vs. Roots"),name='Samples Compared',guide=FALSE)+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+ 
  font("xlab", size = 16, color = "black",face = "bold")+ 
  theme(strip.background=element_blank(),strip.text.x = element_blank(), 
        strip.text.y = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(strip.placement = "outside")+
  labs(x="\n Bray-Curtis Dissimilarity (rarefied 500 counts)",y="Frequency (in 469 samples)\n")+ 
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5))+ guides(linetype=guide_legend(keywidth = 5, keyheight = 1))

uni_plot <- ggplot(BetaDiv[BetaDiv$Tissue=="RosLeaves" & BetaDiv$Metric=="uni_Rare",], 
                    aes(Distance,stat(ndensity),fill=Group))+geom_density(alpha=0.2)+
  facet_grid(cols=vars(Stage),rows=vars(Tissue),labeller=labeller(Stage=Stage.labs, Tissue=Tissue.labs),switch="y")+
  coord_cartesian(xlim=c(0,1))+
  theme_bw()+
  theme(strip.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=c("black","white"),labels=c("Rosettes vs. Rosettes","Rosettes vs. Roots"),name='Samples Compared',guide=FALSE)+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+ 
  font("xlab", size = 16, color = "black",face = "bold")+ 
  theme(strip.background=element_blank(),strip.text.x = element_blank(), 
        strip.text.y = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(strip.placement = "outside")+
  labs(x="\n UniFrac Distance (rarefied 500 counts)",y=" \n")+ 
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5))+ guides(linetype=guide_legend(keywidth = 5, keyheight = 1))

# Plot together
legendCol <- cowplot::get_legend(raup_plot)
ggsave("~/Documents/KBMP2020_Microbes/Figures/BetaDiv_MultiDist.tiff", plot = grid.arrange(
  legendCol,
  raup_plot+ggtitle("A\n")+theme(legend.position = "none",plot.title = element_text(face="bold",size=20)),
  bray_plot+ggtitle("B\n")+theme(plot.title = element_text(face="bold",size=20)),
  uni_plot+ggtitle("C\n")+theme(plot.title = element_text(face="bold",size=20)),
  layout_matrix = rbind(c(1),
                        c(2),
                        c(3),
                        c(4)),
  heights = c(1,4,4,4)), device = NULL, path = NULL,
  scale = 2, width = 7, height = 5, units = c("in"),
  dpi = 600, limitsize = TRUE)

ggsave("~/Documents/KBMP2020_Microbes/Figures/BetaDiv_MultiDist.png", plot = grid.arrange(
  legendCol,
  raup_plot+ggtitle("A\n")+theme(legend.position = "none",plot.title = element_text(face="bold",size=20)),
  bray_plot+ggtitle("B\n")+theme(plot.title = element_text(face="bold",size=20)),
  uni_plot+ggtitle("C\n")+theme(plot.title = element_text(face="bold",size=20)),
  layout_matrix = rbind(c(1),
                        c(2),
                        c(3),
                        c(4)),
  heights = c(1,4,4,4)), device = NULL, path = NULL,
  scale = 2, width = 7, height = 5, units = c("in"),
  dpi = 600, limitsize = TRUE)

