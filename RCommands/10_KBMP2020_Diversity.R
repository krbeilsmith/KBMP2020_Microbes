# KBMP2020_Diversity
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Autumn 2019

# This loops through all the samples. For each one it grabs the relevant metadata, sequence counts and mass for correcting 
# diversity measures, and the diversity measures themselves and puts them in a dataframe.

# Weighted Faith's Phylogenetic Diversity uses the function here: https://github.com/NGSwenson/lefse_0.5/blob/master/R/weighted.faith.R
# with 
library(geiger)

# prepare dataframe
RowVector <- NULL
DiversityMat <- NULL
DiversityMat <- data.frame(matrix(ncol = 20, nrow = 0),stringsAsFactors = FALSE)
DiversityMat <- setNames(data.frame(DiversityMat), c("Mass", "Part", "Site", "Year", "Stage", "SampleSum", "NumTaxa", "NumTaxaCorrect", "LogNumTaxa", "PhyDiv", "PhyDivCorrect", "PhyDivWeighted", "Shannon", "ShannonRare", "MaxRelAbun", "LogMaxRelAbun", "LargestTax", "LargestClass", "LargestFamily", "LargestOrder"))

for(i in sample_names(physeq)){ # for each sample
  print(i)
  
  # 1. METADATA
  Mass <- as.numeric(sample_data[i,"Full_mass"])
  Part <- as.character(sample_data[i,"PlantPart"])
  Site <- as.character(sample_data[i,"Site"])
  Year <- as.character(sample_data[i,"Year"])
  Stage <- as.character(sample_data(physeq)[i,"Stage"]$Stage)
  
  # 2. ALPHA DIVERSITY / SPECIES RICHNESS
  # get the OTU table for the sample and prune it to only the taxa present
  s <- otu_table(physeq)[,i]
  p <- prune_taxa(taxa_sums(s)>0,s)
  # get the total sequence counts, needed for correcting alpha diversity measures for sampling effort
  SampleSum <- as.numeric(as.character(sample_sums(p))) # total counts in sample
  # Alpha Diversity (number of ASVs present)
  NumTaxa <- as.numeric(as.character(length(p))) # alpha diversity
  NumTaxaCorrect <- as.numeric(as.character((length(p)/sum(p))/Mass)) # corrected alpha diversity
  LogNumTaxa <- log(NumTaxaCorrect) # log alpha diversity
  
  # 3. PHYLOGENETIC DIVERSITY AND EVENNESS
  # If there is only one ASV in the sample, set the phylogenetic diversity and evenness to 0
  if(length(row.names(p))==1){
    PhyDiv=0
    PhyDivCorrect=0
    PhyDivWeighted=0
    Shannon=0
    ShannonRare=0
    # Otherwise, make a phyloseq object merging the sample OTU table with metadata and the tree. 
  }else{
    m <- merge_phyloseq(p,tax_table(physeq),sample_data(physeq),phy_tree)
    m <- prune_taxa(taxa_sums(m)>0,m)
    phy <- phy_tree(m)
    comm <- as.data.frame(t(otu_table(m))[,phy$tip.label]) # format objects for phylogenetic diversity calculation
    PhyDiv <- as.numeric(as.character(pd(comm,phy)$PD)) # phylogenetic diversity with picante
    PhyDivCorrect <- as.numeric(as.character(PhyDiv / (pd(comm,phy)$SR))) # phylogenetic diversity corrected for alpha diversity
    PhyDivWeighted <- as.numeric(as.character(weighted.faith(phy,comm))) # weighted phylogenetic diversity
    Shannon<-as.numeric(as.character(vegan::diversity(comm,index="shannon"))) # evenness with vegan
    
    if(SampleSum<100){ # if fewer than 100 counts are in the data, set the rarefied evenness to zero
      ShannonRare=NA
    }else{ # otherwise, rarefy the sample to 100 counts and recalculate Shannon evenness
      mR <- rarefy_even_depth(m, sample.size = 100, rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
      phyR <- phy_tree(mR)
      commR <- as.data.frame(t(otu_table(mR))[,phyR$tip.label])
      ShannonRare <- as.numeric(as.character(vegan::diversity(commR,index="shannon"))) # evenness
    }
  }
  
  # 4. MAXIMUM RELATIVE ABUNDANCE
  MaxRelAbun <- as.numeric(as.character(max(p)/sum(p))) # maximum counts for an ASV in the sample divided by total counts in the sample
  LogMaxRelAbun <- as.numeric(as.character(log(max(p)/sum(p)))) # log of the maximum relative abundance
  LargestTax <- as.character(row.names(p[p==max(p),])) # this is the ASV ID of the highest relative abundance feature
  # Taxonomy information for the ASV with maximum relative abundance
  LargestClass <- as.character(tax_table(physeq)[LargestTax,"Rank3"])
  LargestFamily <- as.character(tax_table(physeq)[LargestTax,"Rank4"])
  LargestOrder <- as.character(tax_table(physeq)[LargestTax,"Rank5"])
  
  # Store data
  RowVector <- NULL
  RowVector <- c(Mass, Part, Site, Year, Stage, SampleSum, NumTaxa, NumTaxaCorrect, LogNumTaxa, PhyDiv, PhyDivCorrect, PhyDivWeighted, 
                 Shannon, ShannonRare, MaxRelAbun, LogMaxRelAbun, LargestTax, LargestClass, LargestFamily, LargestOrder)
  DiversityMat[nrow(DiversityMat)+1,] <- RowVector
}

# If needed to set NAs to zero:
# DiversityMat[is.na(DiversityMat)] <- 0

# Select samples from the second year of study when soil and plants were sampled at multiple timepoints

# Plant samples only, timepoints 2L, 6L, Fl, Sen
Plant_DiversityMat <- DiversityMat[DiversityMat$Year=="2" & DiversityMat$Stage!="Soil" & DiversityMat$Part!="Soil",]
Plant_DiversityMat$Stage <- factor(Plant_DiversityMat$Stage, levels=c("TwoLeaf","SixLeaf","Flowering","Senescent"))

# Soil samples only, timepoints 2L, 6L, Fl, Sen
Soil_DiversityMat <- DiversityMat[DiversityMat$Year=="2" & DiversityMat$Part=="Soil" & DiversityMat$Stage!="Soil",]
Soil_DiversityMat$Stage <- factor(Soil_DiversityMat$Stage, levels=c("TwoLeaf","SixLeaf","Flowering","Senescent"))

Soil_DiversityMat <- Soil_DiversityMat %>% mutate(Type="Soil")
Plant_DiversityMat<- Plant_DiversityMat %>% mutate(Type="Plant")

Div_TS <- Soil_DiversityMat  %>% 
  full_join(Plant_DiversityMat)

Div_TS$Type <- factor(Div_TS$Type, levels=c("Plant","Soil"))

Stages <- c("Soil","TwoLeaf","FourLeaf","SixLeaf","EightLeaf","Flowering","Senescent")

write.table(Div_TS,file="~/Documents/KBMP2020_Microbes/Outputs/Diversity_TimeSeries", sep=",", row.names=TRUE,quote=FALSE)

Div_TS <- read.table("~/Documents/KBMP2020_Microbes/Outputs/Diversity_TimeSeries", sep=",", header=TRUE)

############################################## SPECIES RICHNESS
Div_TS$NumTaxaCorrect <- as.numeric(as.character(Div_TS$NumTaxaCorrect))

comparisons <- list( c("TwoLeaf", "SixLeaf"), c("SixLeaf", "Flowering"), c("Flowering", "Senescent"))
SRTime<-NULL
SRTime <- ggline(Div_TS, x = "Stage", y="NumTaxaCorrect",add = "mean_se", add.params=list(group="Type"),
                 size=0.75,linetype="Type",color="Type",palette = c("black", "gray"),
                 order = c("TwoLeaf", "SixLeaf", "Flowering","Senescent")) +  
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(x="\n ",y="(ASVs / Sequence Count) / Mass (g)\n")+theme(legend.position="top")+
  scale_x_discrete(breaks=Stages[c(2,4,6,7)],labels=c("Two Leaf", "Six Leaf", "Flowering","Senescent"))

SRTime <-SRTime +theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(data=Div_TS%>%filter(Type=="Plant"),comparisons=comparisons,method = "wilcox.test",label="p.signif",fontface=2,size=5,label.y=c(0.08,0.075,0.07))+
  guides(color=guide_legend(title="Sample Type",direction="horizontal"),linetype=guide_legend(title="Sample Type",position="none"))+
  theme(legend.text=element_text(size=12,face=2),legend.title=element_text(size=16,face=2))+coord_cartesian(ylim=c(0.02,0.09))

############################################## WEIGHTED PHY DIV
Div_TS$PhyDivWeighted <- as.numeric(as.character(Div_TS$PhyDivWeighted))

my_comparisons <- list( c("TwoLeaf", "SixLeaf"), c("SixLeaf", "Flowering"), c("Flowering", "Senescent"))
PDTime<-NULL
PDTime <- ggline(Div_TS, x = "Stage", y="PhyDivWeighted",add = "mean_se", add.params=list(group="Type"),
                 size=0.75,linetype="Type",color="Type",palette = c("black", "gray"),
                 order = c("TwoLeaf", "SixLeaf", "Flowering","Senescent")) +  
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(x="\n ",y="\n Weighted Faith's PD\n")+
  scale_x_discrete(breaks=Stages[c(2,4,6,7)],labels=c("Two Leaf", "Six Leaf", "Flowering","Senescent"))

PDTime <-PDTime + rremove("legend")+theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(data=Div_TS%>%filter(Type=="Plant"),comparisons=my_comparisons,method = "wilcox.test",label="p.signif",fontface=2,size=5,label.y=c(8,9,10))

############################################## RAREFIED SHANNON DIV
Div_TS$ShannonRare <- as.numeric(as.character(Div_TS$ShannonRare))

my_comparisons <- list( c("TwoLeaf", "SixLeaf"), c("SixLeaf", "Flowering"), c("Flowering", "Senescent"))
RShTime<-NULL
RShTime <- ggline(subset(Div_TS, !is.na(ShannonRare)), x = "Stage", y="ShannonRare",add = "mean_se", add.params=list(group="Type"),
                  size=0.75,linetype="Type",color="Type",palette = c("black", "gray"),
                  order = c("TwoLeaf", "SixLeaf", "Flowering","Senescent")) +  
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(x="\n ",y="Rarefied Shannon Diversity \n")+
  scale_x_discrete(breaks=Stages[c(2,4,6,7)],labels=c("Two Leaf", "Six Leaf", "Flowering","Senescent"))

RShTime <-RShTime + rremove("legend")+theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(data=subset(Div_TS, !is.na(ShannonRare))%>%filter(Type=="Plant"),comparisons=my_comparisons,method = "wilcox.test",label="p.signif",fontface=2,size=5,label.y=c(2.7,2.9,3.1))

############################################## MAX RELATIVE ABUNDANCE
Div_TS$MaxRelAbun <- as.numeric(as.character(Div_TS$MaxRelAbun))

my_comparisons <- list( c("TwoLeaf", "SixLeaf"), c("SixLeaf", "Flowering"), c("Flowering", "Senescent"),c("SixLeaf","Senescent"))
MRATime<-NULL
MRATime <- ggline(Div_TS, x = "Stage", y="MaxRelAbun",add = "mean_se", add.params=list(group="Type"),
                  size=0.75,linetype="Type",color="Type",palette = c("black", "gray"),
                  order = c("TwoLeaf", "SixLeaf", "Flowering","Senescent")) +  
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(x="\n ",y="\n Max Relative Abundance \n")+coord_cartesian(ylim=c(0.1,0.6))+
  scale_x_discrete(breaks=Stages[c(2,4,6,7)],labels=c("Two Leaf", "Six Leaf", "Flowering","Senescent"))

MRATime <-MRATime + rremove("legend")+theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(data=Div_TS%>%filter(Type=="Plant"),comparisons=my_comparisons,method = "wilcox.test",label="p.signif",fontface=2,size=5,label.y=c(0.5,0.5,0.5,0.4))

#########################################################################################################################################
# DIVERSITY BY PART
#########################################################################################################################################
Div_TS$Part <- factor(Div_TS$Part, levels=c("Soil","Roots","RosLeaves","Stems","CauLeaves","Flowers","Siliques"))

ByPart_SRTime<-NULL
ByPart_SRTime <- ggline(Div_TS%>%filter(Part!="Soil"),
                        x = "Stage", y="NumTaxaCorrect",add = "mean_se", facet.by="Part",size=0.75,panel.labs = list(Part=c("Soil","Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques")),
                        order = c("TwoLeaf", "SixLeaf", "Flowering","Senescent")) +  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(x="\n ",y="(ASVs / Sequence Count) / Mass (g)\n")+
  scale_x_discrete(breaks=Stages[c(2,4,6,7)],labels=c("Two Leaf", "Six Leaf", "Flowering","Senescent"))+
  theme(strip.background.x=element_blank(),strip.text.x = element_text(size = 12, color = "black",face = "bold"))
ByPart_SRTime <-ByPart_SRTime + rremove("legend")+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ByPart_PDTime<-NULL
ByPart_PDTime <- ggline(Div_TS%>%filter(Part!="Soil"),
                        x = "Stage", y="PhyDivWeighted",add = "mean_se", facet.by="Part",size=0.75,panel.labs = list(Part=c("Soil","Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques")),
                        order = c("TwoLeaf", "SixLeaf", "Flowering","Senescent")) +  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(x="\n ",y="\n Weighted Faith's PD\n")+
  scale_x_discrete(breaks=Stages[c(2,4,6,7)],labels=c("Two Leaf", "Six Leaf", "Flowering","Senescent"))+
  theme(strip.background.x=element_blank(),strip.text.x = element_text(size = 12, color = "black",face = "bold"))
ByPart_PDTime <-ByPart_PDTime + rremove("legend")+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ByPart_RShTime<-NULL
ByPart_RShTime <- ggline(subset(Div_TS, !is.na(ShannonRare))%>%filter(Part!="Soil"),
                         x = "Stage", y="ShannonRare",add = "mean_se", facet.by="Part",size=0.75,panel.labs = list(Part=c("Soil","Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques")),
                         order = c("TwoLeaf", "SixLeaf", "Flowering","Senescent")) +  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(x="\n ",y="Rarefied Shannon Diversity \n")+
  scale_x_discrete(breaks=Stages[c(2,4,6,7)],labels=c("Two Leaf", "Six Leaf", "Flowering","Senescent"))+
  theme(strip.background.x=element_blank(),strip.text.x = element_text(size = 12, color = "black",face = "bold"))
ByPart_RShTime <-ByPart_RShTime + rremove("legend")+theme(axis.text.x = element_text(angle = 45, hjust = 1))

ByPart_MRATime<-NULL
ByPart_MRATime <- ggline(Div_TS%>%filter(Part!="Soil"),
                         x = "Stage", y="MaxRelAbun",add = "mean_se", facet.by="Part",size=0.75,panel.labs = list(Part=c("Soil","Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques")),
                         order = c("TwoLeaf", "SixLeaf", "Flowering","Senescent")) +  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(x="\n ",y="\n Max Relative Abundance \n")+
  scale_x_discrete(breaks=Stages[c(2,4,6,7)],labels=c("Two Leaf", "Six Leaf", "Flowering","Senescent"))+
  theme(strip.background.x=element_blank(),strip.text.x = element_text(size = 12, color = "black",face = "bold"))
ByPart_MRATime <-ByPart_MRATime + rremove("legend")+theme(axis.text.x = element_text(angle = 45, hjust = 1))

MyLegend<-get_legend(SRTime)

GridDiv<-grid.arrange(
  MyLegend,
  SRTime+ theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("A\n"), 
  PDTime+ theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("B\n"),
  RShTime+ theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("C\n"),
  MRATime+ theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("D\n"),
  heights=c(1,5,5),
  widths = c(5,5),
  layout_matrix = rbind(c(1),
                        c(2,3),
                        c(4,5))
)

ggsave("~/Documents/KBMP2020_Microbes/Figures/DiversityPlot.tiff", plot = GridDiv, device = NULL, path = NULL,
       scale = 1.8, width = 5, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)

ggsave("~/Documents/KBMP2020_Microbes/Figures/DiversityPlot.png", plot = GridDiv, device = NULL, path = NULL,
       scale = 1.8, width = 5, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)

GridDiv_ByPart<-grid.arrange(
  ByPart_SRTime+ theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("A\n"), 
  ByPart_PDTime+ theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("B\n"),
  ByPart_RShTime+ theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("C\n"),
  ByPart_MRATime+ theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("D\n"),
  heights=c(5,5),
  widths = c(5,5),
  layout_matrix = rbind(c(1,2),
                        c(3,4))
)

ggsave("~/Documents/KBMP2020_Microbes/Figures/DiversityPlot_ByPart.tiff", plot = GridDiv_ByPart, device = NULL, path = NULL,
       scale = 1.8, width = 6, height = 6, units = c("in"),
       dpi = 600, limitsize = TRUE)

ggsave("~/Documents/KBMP2020_Microbes/Figures/DiversityPlot_ByPart.png", plot = GridDiv_ByPart, device = NULL, path = NULL,
       scale = 1.8, width = 6, height = 6, units = c("in"),
       dpi = 600, limitsize = TRUE)

#########################################################################################################################################
# RELATIVE ABUNDANCE DISTRIBUTIONS
#########################################################################################################################################

physeq_forDiversity <-subset_samples(physeq, PlantPart!="Soil")
physeq_forDiversity <- prune_samples(sample_sums(physeq_forDiversity)>100,physeq_forDiversity)

SampleStages <- as.character(unique(sample_data(physeq_forDiversity)[,"Stage"])$Stage)

Diversity <- NULL
Diversity <- data.frame(matrix(ncol = 2, nrow = 0),stringsAsFactors = FALSE)
Diversity <- setNames(data.frame(Diversity), c("SampRelAbund","Stage"))

for(stage in SampleStages){
  print(stage)
  physeq_partition <- subset_samples(physeq_forDiversity,Stage == stage)
  print(length(sample_names(physeq_partition)))
  for(SampID in sample_names(physeq_partition)){
    print(SampID)
    SampTaxa <- subset_samples(physeq_partition,sample_names(physeq_partition) %in% SampID)
    SampTaxa <- prune_taxa(taxa_sums(SampTaxa)>0,SampTaxa)
    RA_table  <- transform_sample_counts(SampTaxa,  function(x) x / sum(x))
    Results <- c()
    Results <- cbind(as.data.frame(taxa_sums(RA_table))[,1], stage)
    colnames(Results) <- c("SampRelAbund","Stage")
    Diversity <- rbind(Diversity, Results)
  }
}

Diversity$Stage <- factor(Diversity$Stage, levels = c("TwoLeaf","FourLeaf","SixLeaf","EightLeaf","Flowering","Senescent"))
Diversity$SampRelAbund <- as.numeric(as.character(Diversity$SampRelAbund))

Stage.labs <- c("Two Leaf","Four Leaf","Six Leaf","Eight Leaf","Flowering","Senescent")
names(Stage.labs) <- c("TwoLeaf","FourLeaf","SixLeaf","EightLeaf","Flowering","Senescent")

library(plyr)
medians <- ddply(Diversity, "Stage", summarise, grp.median=median(SampRelAbund))

# Frequency plots
freq <- ggplot(Diversity, aes(log10(SampRelAbund),stat(ndensity),color=Stage)) + geom_density(alpha=1,size=1.2)+
  # geom_histogram(size=1,binwidth=0.1,fill="gray70",color="gray70")+
  #scale_color_manual(values=c("grey","gray50","gray50","gray50","gray50","black"),
                     #labels=c("Two Leaf","Four Leaf", "Six Leaf","Eight Leaf","Flowering","Senescent"))+
  scale_color_grey(start=0.8, end=0,labels=c("Two Leaf","Four Leaf", "Six Leaf","Eight Leaf","Flowering","Senescent"))+
  theme_bw()+
  theme(strip.background=element_blank(),strip.text.x = element_text(size = 16, color = "black",face = "black"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+ 
  font("xlab", size = 16, color = "black",face = "bold")+ 
  labs(x=bquote(bold(atop("\n",~log[10]~"ASV Relative Abundance"))),y="Frequency in Plant Samples\n")+
  theme(strip.text.y = element_text(size=12, face="bold",margin = margin(0,0,0,1, "cm"),vjust=5))+
  theme(strip.placement = "outside")+
  theme(panel.spacing = unit(1, "lines"))+
  #geom_vline(data=medians, aes(xintercept=log10(grp.median)),linetype="dashed", color="black",size=1)+rremove("legend")
  theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),legend.position="right")+
  theme(legend.background = element_rect(color="black", size=.5))+ guides(color=guide_legend(keywidth = 0.75, keyheight = 0.75, reverse=FALSE))

ggsave("~/Documents/KBMP2020_Microbes/Figures/RelativeAbundancesByStage.tiff", plot = freq, device = NULL, path = NULL,
       scale = 1.5, width = 6, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)

ggsave("~/Documents/KBMP2020_Microbes/Figures/RelativeAbundancesByStage.png", plot = freq, device = NULL, path = NULL,
       scale = 1.5, width = 6, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)
