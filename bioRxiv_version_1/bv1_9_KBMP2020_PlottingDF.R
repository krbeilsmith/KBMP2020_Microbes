# KBMP2020_PlottingDF
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Autumn 2019

# This will collect the data for making time series plots and the dotplots to represent spatial patterns for the ASVs that drive the associations of community
# dissimilarities with tissue and developmental stage.

# Select plant samples that are among the top drivers of the associations with host variables (threshold based on the summed absolute values of PERMANOVA 
# coefficients for regression with Plant Part). See 7_KBMP2020_ASVDrivers.R for seleciton procedure.

physeq_plants <- subset_samples(physeq,PlantPart!="Soil")
physeq_plants <- prune_taxa(taxa_sums(physeq_plants)>0,physeq_plants)
physeq_select <- subset_taxa(physeq_plants,taxa_names(physeq_plants)%in%unlist(DrivingASVs))
physeq_select <- prune_samples(sample_sums(physeq_select)>0,physeq_select)

# Get all the unique combinations of harvest conditions (tissue, stage, site, year)
ConditionsFrame <- data.frame(unique(sample_data(physeq_select)[,c("Year","Site","Stage","PlantPart")]),row.names=c())
ConditionsFrame <- data.frame(lapply(ConditionsFrame, as.character), stringsAsFactors=FALSE)

# Make dataframe for storing results
PlottingFrame <- NULL
PlottingFrame <- data.frame(matrix(ncol = 13, nrow = 0),stringsAsFactors = FALSE)
PlottingFrame  <- setNames(data.frame(PlottingFrame), c("ASV","Tissue","Stage","Site","Year","Prevalence","MeanRA","sdRA",
                                                        "Genus","Family","Order","Class","Phylum"))

for(m in unique(unlist(DrivingASVs))){ # for each ASV (referred to as focal ASV below)
  print(m)
  for(row in rownames(ConditionsFrame)){ # for each set of harvest conditions
    pp <- ConditionsFrame[row,"PlantPart"] # get metadata for the harvest conditions
    st <- ConditionsFrame[row,"Stage"]
    si <- ConditionsFrame[row,"Site"]
    yr <- ConditionsFrame[row,"Year"]
    physeq_partition <-subset_samples(physeq_select, PlantPart==pp & Stage==st & Site==si & Year==yr) # prune to data for the harvest conditions
    physeq_partition <- prune_taxa(taxa_sums(physeq_partition)>0, physeq_partition) # prune to ASVs present in harvest conditions
    if(m %in% taxa_names(physeq_partition)){ # if focal ASV is present in harvest conditions
      PA_table  <- transform_sample_counts(physeq_partition, function(x) ifelse(x>0,1,0)) # convert sample data to presence-absence binary
      Prevalence <- taxa_sums(PA_table)[m]/length(sample_sums(PA_table)) # calculate prevalence of the ASV in samples from the harvest conditions
      RA_table  <- transform_sample_counts(physeq_partition,  function(x) x / sum(x)) # convert sample data to relative abundances
      RelAbund <- data.frame(otu_table(RA_table)[m]) 
      MeanRelAbund <- mean(RelAbund[RelAbund>0]) # calculate mean of the relative abundances of the ASV for samples in the harvest conditions
      SDRelAbund <- sd(RelAbund[RelAbund>0]) # calculate standard deviation of the relative abundances of the ASV for samples in the harvest conditions
    }else{ # if the focal ASV is not present, set prevalence and relative abundance measures to zero
      Prevalence <- 0
      MeanRelAbund <- 0
      SDRelAbund <- 0
    }
    Genus <- as.character(tax_table(physeq)[m,"Rank6"]) # get taxonomy data for focal ASV
    Family <- as.character(tax_table(physeq)[m,"Rank5"])
    Order <- as.character(tax_table(physeq)[m,"Rank4"])
    Class <- as.character(tax_table(physeq)[m,"Rank3"])
    Phylum <- as.character(tax_table(physeq)[m,"Rank2"])
    PlottingFrame[nrow(PlottingFrame)+1,] <- c(m,pp,st,si,yr,Prevalence,MeanRelAbund,SDRelAbund,Genus,Family,Order,Class,Phylum) # add values to dataframe
  }
}

# Store dataframe
write.table(PlottingFrame,file="~/Documents/KBMP2020_Microbes/Outputs/PlottingFrame", sep="\t", row.names=TRUE,quote=FALSE)

# PlottingFrame <- read.table("~/Documents/KBMP2020_Microbes/Outputs/PlottingFrame", header=TRUE, sep="\t",row.names=1)

# Set factor levels for plotting
PlottingFrame$Tissue <-factor(PlottingFrame$Tissue, levels=c("Roots","RosLeaves","Stems","CauLeaves","Flowers","Siliques"))
PlottingFrame$Stage <-factor(PlottingFrame$Stage, levels=c("TwoLeaf","FourLeaf","SixLeaf","EightLeaf","Flowering","Senescent"))
PlottingFrame$Site <-factor(PlottingFrame$Site, levels=c("ME","WW"))
PlottingFrame$Year <-factor(PlottingFrame$Year, levels=c("1","2"))
PlottingFrame$Genus <-factor(PlottingFrame$Genus)
PlottingFrame$Family <-factor(PlottingFrame$Family)
PlottingFrame$Order <-factor(PlottingFrame$Order)
PlottingFrame$Class <-factor(PlottingFrame$Class)
PlottingFrame$Phylum <-factor(PlottingFrame$Phylum)

# Make prevalence and relative abundance data numeric for plotting
PlottingFrame$Prevalence <- as.numeric(PlottingFrame$Prevalence)
PlottingFrame$MeanRA <- as.numeric(PlottingFrame$MeanRA)
PlottingFrame$sdRA <- as.numeric(PlottingFrame$sdRA)

#################################################################################################################################################
# Make dot plots for the ASVs driving the tissue effects to show where in the plant they peak in prevalence and mean relative abundance
# in both sites and years.
#################################################################################################################################################

detach(package:plyr)    
library(dplyr)

# For dotplot, take ASVs that consistently reach peak prevalence in the same tissue and have over 5% prevalence in all sites/years.

# In what tissue does ASV prevalence peak?
PeakLocations <- PlottingFrame %>% filter(ASV %in% rownames(PlantTissueEffects)) %>%
  group_by(ASV,Site,Year) %>% dplyr::slice(which.max(Prevalence))

# 2. For each tissue, look at the ASVs peaking there and determine how many times they peaked there (site/year).
# Take ASVs that peak in the same tissue in all 4 replicate plantings across sites and years.

Tissues <- c("Roots","RosLeaves","Stems","Siliques")
for(t in Tissues){ # for each tissue
  Taxa <- data.frame(table(PeakLocations[PeakLocations$Tissue==t,]$ASV)) # find ASVs that peak in that tissue and get frequency (1-4)
  colnames(Taxa) <- c("ASV","Freq")
  Taxa$ASV <- as.character(Taxa$ASV)
  Taxa <- unique(Taxa[Taxa$Freq==4,]$ASV) # if observed in all 4 site/year combinations, add to list
  assign(paste(t,"Taxa",sep=""),Taxa) # name list
  print(t)
}

ConsistentASVs <- c(RootsTaxa, RosLeavesTaxa, StemsTaxa, SiliquesTaxa)

# 3. Option A:
# This will select all the consistently peaking ASVs that have >5% prevalence in each site and year combination and that have PERMANOVA coefficients in the top 
# 5% for the regression with plant tissue type. See 7_KB_A19_ASV_Drivers.R for seleciton procedure. 
# The plot made below will have a dot colored for the tissue type with a size to reflect the peak prevalence (across stages) in that tissue and a transparency 
# to reflect the peak mean relative abundance (across stages) in that tissue. 
DotPlotData <- PlottingFrame %>% filter(ASV %in% rownames(PlantTissueEffects) & ASV %in% ConsistentASVs) %>% group_by(ASV,Site,Year) %>% 
  mutate(SYSum = sum(Prevalence)) %>% filter(all(SYSum > 0.05)) %>% ungroup %>% group_by(ASV,Site,Year,Tissue) %>%
  mutate(MP = max(Prevalence), MRA=max(MeanRA))
DotPlotData$ASV <- factor(DotPlotData$ASV, levels=c(RootsTaxa,RosLeavesTaxa,StemsTaxa,SiliquesTaxa))

# 3. Option B:
# Most of the ASVs do not show a consistent spatial distribution (do not peak in the same tissue across sites and years). 
# This will select all the inconsistent ASVs in the top 50 for driving tissue effects and that have >5% prevalence in each site and year combination.
# The plot made below will have a dot colored for the tissue type with a size to reflect the peak prevalence (across stages) in that tissue and a transparency 
# to reflect the peak mean relative abundance (across stages) in that tissue. 

TopTissueDrivers <- PlantTissueEffects %>% mutate(ASV = rownames(PlantTissueEffects)) %>% arrange(EffectSums) %>% top_n(50,EffectSums)

DotPlotData <- PlottingFrame %>% filter(ASV %in% rownames(PlantTissueEffects) & !(ASV %in% ConsistentASVs) &
                                          ASV %in% TopTissueDrivers$ASV) %>%
  group_by(ASV,Site,Year) %>% mutate(SYSum = sum(Prevalence)) %>% filter(all(SYSum > 0.05)) %>% ungroup %>% 
  group_by(ASV,Site,Year,Tissue) %>% mutate(MP = max(Prevalence), MRA=max(MeanRA))

DotPlotData$Genus <- ifelse(DotPlotData$Genus=="aggregate unclassified genus","unclassified",as.character(DotPlotData$Genus))
DotPlotData$Genus <- ifelse(DotPlotData$Genus=="Burkholderia-Paraburkholderia","Paraburkholderia",as.character(DotPlotData$Genus))

Site.labs<-c("Site: ME", "Site: WW")
names(Site.labs)<-levels(DotPlotData$Site)

Year.labs<-c("Year: 1", "Year: 2")
names(Year.labs)<-levels(DotPlotData$Year)

Plot_Filter_SP<-NULL
Plot_Filter_SP<-ggplot(data= DotPlotData, aes(x=ASV, y=Tissue,group = interaction(ASV,Site))) + 
  geom_point(aes(size = MP, alpha=MRA, color=Tissue))+
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"))+
  theme_bw()+
  theme(strip.background.y=element_blank(),strip.text.x = element_text(size = 20, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 16, color = "black", face = "bold")+
  labs(y="",x="\n16S Amplicon Sequence Variant")+
  scale_alpha_continuous(range=c(0.1,1),limits=c(0,1),breaks = c(0,0.2,0.4,0.6,0.8,1),labels=c("0%","20%","40%","60%","80%","100%"))+
  scale_size_continuous(range=c(0,5),limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1))+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),
        legend.position="right")+
  theme(legend.background = element_rect(color="black", size=.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12))+
  facet_grid(rows=vars(Site,Year),cols=vars(Genus),scales="free_x",space="free_x",switch="y",labeller=labeller(Site=Site.labs, Year=Year.labs))+
  scale_x_discrete(expand=c(0.025,1),label=function(x) abbreviate(x, minlength=8))+
  theme(strip.text.y = element_text(size=16, face="bold", vjust = 2),strip.placement="outside")+
  guides(size = guide_legend(title="Peak Prevalence",order=1,nrow=1),
         color=FALSE,alpha=guide_legend(title="Peak Mean Relative Abundance",override.aes = list(size=5),nrow=1))+
  theme(strip.text.x = element_text(size = 16,angle=90),strip.background.x= element_rect(colour="black",fill=NA,size=1))+
  theme(plot.margin = unit(c(1,3,1,1), "cm"))+
  theme(strip.placement = "outside")+ theme(panel.spacing = unit(0, "lines"))+ theme(panel.spacing.y = unit(1, "lines"))+
  theme(legend.text=element_text(size=14,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),
        legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5))

ggsave("~/Documents/KBMP2020_Microbes/Figures/DotPlot_SpatialPatterns_B.tiff", plot = Plot_Filter_SP, device = NULL, path = NULL,
       scale = 3.0, width = 5.2, height = 4, units = c("in"),
       dpi = 600, limitsize = TRUE)

#################################################################################################################################################
# Make a time series plot for the ASVs with the largest summed coefficients for stage in the PERMANOVA
#################################################################################################################################################

# This selects ASVs among the top 50 for tissue effect that also drove stage effects (based on the threshold set at top 5% of regression coefficients).
# Root and rosette samples are selected and the mean prevalence/abundance of ASVs at each stage, site, and year in each of those tissues is calculated.
MyData <- PlottingFrame %>% filter(ASV %in% rownames(PlantStageEffects)) %>% filter(Tissue%in%c("Roots","RosLeaves")) %>% filter(ASV %in% TopTissueDrivers$ASV) %>%
  group_by(ASV,Site,Year,Stage) %>% mutate(MP = mean(Prevalence),MRA=mean(MeanRA))

MyData$Tissue <- factor(MyData$Tissue, levels=c("RosLeaves","Roots"))

MyPalette <- c("#20A486FF","#EEDD82FF","#404688FF","#000004FF","#F56B5CFF","#75D054FF","#9C2E7FFF")

Site.labs<-c("Site: ME", "Site: WW")
names(Site.labs)<-levels(MyData$Site)

Year.labs<-c("Year: 1", "Year: 2")
names(Year.labs)<-levels(MyData$Year)

Tissue.labs<-c("Rosettes", "Roots")
names(Tissue.labs)<-levels(Select_MyData$Tissue)

TSA <- NULL
TSA <- ggplot(data=MyData, aes(x=as.factor(Stage), y=MP, group = interaction(ASV,Site,Year),color=Class)) +
  geom_line(data=MyData,aes(x=as.factor(Stage), y=MP, group = interaction(ASV,Site,Year),color=Class),size=0.5,alpha=0.5) +
  geom_point(data=MyData,aes_string(x="Stage",y="MP"),alpha=0.5,size=0.5)+
  theme_bw()+
  theme(strip.background=element_blank(),strip.text.x = element_text(size = 16, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+#+rremove("legend")+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(y="Prevalence",x=" ")+
  scale_x_discrete(labels=c("Two Leaf", "Four Leaf", "Six Leaf", "Eight Leaf", "Flowering", "Senescent"))+
  scale_y_continuous(breaks=c(0.5,1),limits=c(0,1))+
  scale_color_manual(name="Genus",values=MyPalette,guide=FALSE)+
  facet_grid(cols = vars(Site,Year),rows=vars(Tissue),labeller=labeller(Site=Site.labs, Year=Year.labs, Tissue=Tissue.labs),switch="y")+
  theme(panel.spacing = unit(2, "lines"))+
  guides(color=guide_legend(title="Class",override.aes = list(size=5)))+
  theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),
        legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.text.y = element_text(size=16, face="bold",margin = margin(0,0,0,1, "cm"),vjust=5))+
  theme(strip.text.x = element_text(size=16, face="bold",margin = margin(0,0,0.1,0, "cm")))+
  theme(strip.placement = "outside")

ggsave("~/Documents/KBMP2020_Microbes/Figures/TimeSeries_MeanAcrossTissues.tiff", plot = TSA, device = NULL, path = NULL,
       scale = 1.7, width = 6, height = 4, units = c("in"),
       dpi = 600, limitsize = TRUE)

#################################################################################################################################################
# Make a time series plot to break down the temporal trends in Alphaproteobacteria and Betaproteobacteria across root and rosette tissues
#################################################################################################################################################

#################################################################################################################################################
# Alphaproteobacteria
#################################################################################################################################################

# Select data from known genera and factor it. The time series for each of these ASVs will be plotted as a transparent line.
MyData <- PlottingFrame %>% filter(ASV %in% rownames(PlantStageEffects)) %>% filter(Genus != "aggregate unclassified genus")

Select_MyData<-MyData[MyData$Genus%in%c("Rhizobium","Methylobacterium","Rhodopseudomonas","Aureimonas"),]
Select_MyData<-Select_MyData[Select_MyData$Tissue%in%c("Roots","RosLeaves"),]

Select_MyData$Genus<-factor(Select_MyData$Genus,levels=c("Rhizobium","Methylobacterium","Rhodopseudomonas","Aureimonas"))
Select_MyData$Tissue<-factor(Select_MyData$Tissue,levels=c("RosLeaves","Roots"))

# Find points with prevalence over 70% and get the names of the ASVs with those points
LabelsList<-Select_MyData
LabelsList<- LabelsList%>%
  group_by(ASV,Site,Year,Tissue) %>%
  filter(Prevalence == max(Prevalence)) %>%
  filter(Prevalence>0.70)
Highlights<-unique(LabelsList$ASV)

# Make a copy of the dataframe for highlighted ASVs (those that reach above 70% prevalence somewhere in the dataset). Each ASV will have a solid line
# with points sized by relative abundance for the time series.
Highlight_MyData<-Select_MyData[Select_MyData$ASV%in%Highlights,]
Highlight_MyData$Genus<-factor(Highlight_MyData$Genus,levels=c("Rhizobium","Methylobacterium","Rhodopseudomonas","Aureimonas"))

# Make a copy of the dataframe with shortened ASV names for labels at the points of maximum nonzero prevalence for highlighted ASVs
MyLabels<-Select_MyData[Select_MyData$ASV%in%Highlights,]
MyLabels<- MyLabels%>%
  group_by(ASV,Site,Year,Tissue) %>%
  filter(Prevalence == max(Prevalence))
MyLabels<-MyLabels[MyLabels$Prevalence>0,]
MyLabels$ASV<-substring(MyLabels$ASV,1,3)

MyPalette <- c("#FD9B6BFF","#9C2E7FFF","#20A486FF","#D5E21AFF")

Site.labs<-c("Site: ME", "Site: WW")
names(Site.labs)<-levels(Select_MyData$Site)

Year.labs<-c("Year: 1", "Year: 2")
names(Year.labs)<-levels(Select_MyData$Year)

Tissue.labs<-c("Rosettes", "Roots")
names(Tissue.labs)<-levels(Select_MyData$Tissue)

Alpha_TS <-NULL
Alpha_TS <- ggplot(data=Select_MyData, aes(x=as.factor(Stage), y=Prevalence, group = interaction(ASV,Site,Tissue),color=Genus)) +
  geom_line(data=Highlight_MyData,aes(x=as.factor(Stage), y=Prevalence, group = interaction(ASV,Site,Tissue),color=Genus),size=.75,alpha=1) +
  geom_line(data=Select_MyData,aes(x=as.factor(Stage), y=Prevalence, group = interaction(ASV,Site,Tissue),color=Genus),size=0.75,alpha=0.25) +
  geom_point(data=Highlight_MyData,aes_string(x="Stage",y="Prevalence",size="MeanRA"),alpha=1)+
  theme_bw()+
  theme(strip.background=element_blank(),strip.text.x = element_text(size = 16, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+#+rremove("legend")+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(y="Prevalence",x=" ")+
  scale_size_continuous(name="Mean Relative Abundance",range=c(1,6),breaks=c(0.0,0.2,0.4,0.6,0.8,1),limits=c(0,1))+
  scale_x_discrete(labels=c("Two Leaf", "Four Leaf", "Six Leaf", "Eight Leaf", "Flowering", "Senescent"))+
  scale_y_continuous(breaks=c(0.5,1),limits=c(0,1))+
  scale_color_manual(name="Genus",values=MyPalette,guide=FALSE)+
  geom_label_repel(data=MyLabels,size=4,aes_string(x="Stage",y="Prevalence",label = "ASV",color="Genus"),fontface=2,alpha=1,hjust=1,vjust=2,direction="both",show.legend = FALSE)+
  facet_grid(cols = vars(Site,Year),rows=vars(Tissue),labeller=labeller(Site=Site.labs, Year=Year.labs, Tissue=Tissue.labs),switch="y")+
  theme(panel.spacing = unit(2, "lines"))+
  guides(size=guide_legend(title="Relative Abundance",nrow=1,order=1),color=guide_legend(title="Genus",override.aes = list(size=5)))+
  theme(legend.text=element_text(size=14,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5),legend.box="vertical")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.text.y = element_text(size=16, face="bold",margin = margin(0,0,0,1, "cm"),vjust=5))+
  theme(strip.text.x = element_text(size=16, face="bold",margin = margin(0,0,0.1,0, "cm")))+
  theme(strip.placement = "outside")

ggsave("~/Documents/KBMP2020_Microbes/Figures/TimeSeries_Alpha_TS.tiff", plot = Alpha_TS, device = NULL, path = NULL,
       scale = 1.7, width = 6, height = 4, units = c("in"),
       dpi = 600, limitsize = TRUE)

#################################################################################################################################################
# Betaproteobacteria
#################################################################################################################################################

MyDataB <- PlottingFrame %>% filter(ASV %in% rownames(PlantStageEffects)) %>% filter(Genus != "aggregate unclassified genus")

Select_MyDataB<-MyDataB[MyDataB$Family%in%c("Oxalobacteraceae","Comamonadaceae"),]
Select_MyDataB<-Select_MyDataB[Select_MyDataB$Tissue%in%c("Roots","RosLeaves"),]
Select_MyDataB$Genus<-factor(Select_MyDataB$Genus,levels=c("Duganella","Massilia","Pelomonas","Janthinobacterium"))
Select_MyDataB$Tissue<-factor(Select_MyDataB$Tissue,levels=c("RosLeaves","Roots"))

LabelsListB<-Select_MyDataB
LabelsListB<- LabelsListB %>%
  group_by(ASV,Site,Year,Tissue) %>%
  filter(Prevalence == max(Prevalence)) %>%
  filter(Prevalence>0.70)
HighlightsB<-unique(LabelsListB$ASV)

Highlight_MyDataB<-Select_MyDataB[Select_MyDataB$ASV%in%HighlightsB,]
Highlight_MyDataB$Genus<-factor(Highlight_MyDataB$Genus,levels=c("Duganella","Massilia","Pelomonas","Janthinobacterium"))

MyLabelsB<-Select_MyDataB[Select_MyDataB$ASV%in%HighlightsB,]
MyLabelsB<- MyLabelsB%>%
  group_by(ASV,Site,Year,Tissue) %>%
  filter(Prevalence == max(Prevalence))
MyLabelsB<-MyLabelsB[MyLabelsB$Prevalence>0,]
MyLabelsB$ASV<-substring(MyLabelsB$ASV,1,3)

MyPaletteB <- c("#000004FF","#F56B5CFF","#404688FF","#75D054FF","#FD9B6BFF","#D8456CFF","#9C2E7FFF","#D5E21AFF")

Site.labs<-c("Site: ME", "Site: WW")
names(Site.labs)<-levels(Select_MyDataB$Site)

Year.labs<-c("Year: 1", "Year: 2")
names(Year.labs)<-levels(Select_MyDataB$Year)

Tissue.labs<-c("Rosettes", "Roots")
names(Tissue.labs)<-levels(Select_MyDataB$Tissue)

Beta_TS <-NULL
Beta_TS <- ggplot(data=Select_MyDataB, aes(x=as.factor(Stage), y=Prevalence, group = interaction(ASV,Site,Tissue),color=Genus)) +
  geom_line(data=Highlight_MyDataB,aes(x=as.factor(Stage), y=Prevalence, group = interaction(ASV,Site,Tissue),color=Genus),size=.75,alpha=1) +
  geom_line(data=Select_MyDataB,aes(x=as.factor(Stage), y=Prevalence, group = interaction(ASV,Site,Tissue),color=Genus),size=0.75,alpha=0.25) +
  geom_point(data=Highlight_MyDataB,aes_string(x="Stage",y="Prevalence",size="MeanRA"),alpha=1)+
  theme_bw()+
  theme(strip.background=element_blank(),strip.text.x = element_text(size = 16, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+#+rremove("legend")+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(y="Prevalence",x=" ")+
  scale_size_continuous(name="Mean Relative Abundance",range=c(1,6),breaks=c(0.0,0.2,0.4,0.6,0.8,1),limits=c(0,1))+
  scale_x_discrete(labels=c("Two Leaf", "Four Leaf", "Six Leaf", "Eight Leaf", "Flowering", "Senescent"))+
  scale_y_continuous(breaks=c(0.5,1),limits=c(0,1))+
  scale_color_manual(name="Genus",values=MyPaletteB,guide=FALSE)+
  geom_label_repel(data=MyLabelsB,size=4,aes_string(x="Stage",y="Prevalence",label = "ASV",color="Genus"),fontface=2,alpha=1,hjust=1,vjust=2,direction="both",show.legend = FALSE)+
  facet_grid(cols = vars(Site,Year),rows=vars(Tissue),labeller=labeller(Site=Site.labs, Year=Year.labs, Tissue=Tissue.labs),switch="y")+
  theme(panel.spacing = unit(2, "lines"))+
  #theme(legend.position = "none")+
  guides(size=guide_legend(title="Relative Abundance",nrow=1,order=1),color=guide_legend(title="Genus",override.aes = list(size=5)))+
  theme(legend.text=element_text(size=14,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5),legend.box="vertical")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.text.y = element_text(size=16, face="bold",margin = margin(0,0,0,1, "cm"),vjust=5))+
  theme(strip.text.x = element_text(size=16, face="bold",margin = margin(0,0,0.1,0, "cm")))+
  theme(strip.placement = "outside")

ggsave("~/Documents/KBMP2020_Microbes/Figures/TimeSeries_Beta_TS.tiff", plot = Beta_TS, device = NULL, path = NULL,
       scale = 1.7, width = 6, height = 4, units = c("in"),
       dpi = 600, limitsize = TRUE)

###########################################################################################################################
# Together
###########################################################################################################################

ggsave("~/Documents/KBMP2020_Microbes/Figures/TimeSeries_PrevAndRA.tiff", plot = grid.arrange(
  Alpha_TS+ggtitle("A: Alphaproteobacteria\n")+theme(plot.title = element_text(face="bold",size=20)),
  Beta_TS+ggtitle("B: Betaproteobacteria\n")+theme(plot.title = element_text(face="bold",size=20)),
  layout_matrix = rbind(c(1),
                        c(2))), device = NULL, path = NULL,
  scale = 1.7, width = 8, height = 11, units = c("in"),
  dpi = 600, limitsize = TRUE)
