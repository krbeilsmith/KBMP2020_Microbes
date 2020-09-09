# KBMP2020_DotPlot
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Summer 2020

####################################################################################################################################################

PlottingFrame <- read.table("~/Documents/KBMP2020_Microbes/Outputs/PlottingFrame_Indicators", header=TRUE, sep="\t",row.names=1)

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

####################################################################################################################################################

# In what tissue does ASV prevalence peak?
PeakLocations <- PlottingFrame %>% filter(ASV %in% rownames(PlantPart_Indicators)) %>%
  group_by(ASV,Site,Year) %>% dplyr::slice(which.max(Prevalence))

# 2. For each tissue, look at the ASVs peaking there and determine how many times they peaked there (site/year).
# Take ASVs that peak in the same tissue in all 4 replicate plantings across sites and years.

Tissues <- c("Roots","RosLeaves","Stems","CauLeaves","Flowers","Siliques")
for(t in Tissues){ # for each tissue
  Taxa <- data.frame(table(PeakLocations[PeakLocations$Tissue==t,]$ASV)) # find ASVs that peak in that tissue and get frequency (1-4)
  colnames(Taxa) <- c("ASV","Freq")
  Taxa$ASV <- as.character(Taxa$ASV)
  Taxa <- unique(Taxa[Taxa$Freq==4,]$ASV) # if observed in all 4 site/year combinations, add to list
  assign(paste(t,"Taxa",sep=""),Taxa) # name list
  print(t)
}

ConsistentASVs <- c(RootsTaxa, RosLeavesTaxa, StemsTaxa, CauLeavesTaxa, FlowersTaxa, SiliquesTaxa)

# % of Plant Part indicator ASVs that consistently peak in one tissue: 21.3%
(length(ConsistentASVs) / length(rownames(PlantPart_Indicators)))*100

# % of consistently peaking ASVs for each tissue
(length(RootsTaxa) / length(ConsistentASVs))*100 # roots 78.6%
(length(RosLeavesTaxa) / length(ConsistentASVs))*100 # leaf 11.2%
(length(StemsTaxa) / length(ConsistentASVs))*100 # stem 5.1%
(length(SiliquesTaxa) / length(ConsistentASVs))*100 # silique 5.1%

DotPlotData <- PlottingFrame %>% filter(ASV %in% rownames(PlantPart_Indicators)) %>%
  group_by(ASV,Site,Year) %>% mutate(SYSum = sum(Prevalence)) %>% filter(all(SYSum > 0.00)) %>% ungroup %>% 
  group_by(ASV,Site,Year,Tissue) %>% mutate(MP = max(Prevalence), MRA=max(MeanRA))

DotPlotData$PeakTissue <- ifelse(DotPlotData$ASV%in%RootsTaxa, "Roots", 
                                 ifelse(DotPlotData$ASV%in%RosLeavesTaxa, "RosLeaves",
                                        ifelse(DotPlotData$ASV%in%StemsTaxa, "Stems",
                                               ifelse(DotPlotData$ASV%in%SiliquesTaxa, "Siliques","Inconsistent"))))

####################################################################################################################################################

DotPlotData$Genus <- ifelse(DotPlotData$Genus=="aggregate unclassified genus","unclassified",as.character(DotPlotData$Genus))
DotPlotData$Genus <- ifelse(DotPlotData$Genus=="Burkholderia-Paraburkholderia","Burkholderia",as.character(DotPlotData$Genus))

Select_DotPlotData <- DotPlotData

Select_DotPlotData$PeakTissue<-factor(Select_DotPlotData$PeakTissue, levels=c("Inconsistent","Roots","RosLeaves","Siliques","Stems"))

LevelList <- Select_DotPlotData[order(Select_DotPlotData$Class, Select_DotPlotData$Genus, Select_DotPlotData$ASV),]$Genus
MyData$Genus <- factor(MyData$Genus, levels=unique(LevelList))

# 362 ASVs peak in different tissues across sites and years
length(as.character(unique(Select_DotPlotData[Select_DotPlotData$PeakTissue=="Inconsistent",]$ASV)))
# 98 ASVs peak in the same tissue across sites and years
length(as.character(unique(Select_DotPlotData[Select_DotPlotData$PeakTissue!="Inconsistent",]$ASV)))

Select_DotPlotData <- Select_DotPlotData[Select_DotPlotData$Genus != "unclassified",]

Site.labs<-c("Site: ME", "Site: WW")
names(Site.labs)<-levels(Select_DotPlotData$Site)

Year.labs<-c("Year: 1", "Year: 2")
names(Year.labs)<-levels(Select_DotPlotData$Year)

DotPlotA <- NULL
DotPlotA <- ggplot(data= Select_DotPlotData[Select_DotPlotData$Phylum=="Proteobacteria" & Select_DotPlotData$PeakTissue!="Inconsistent",], aes(x=ASV, y=Tissue)) + 
  geom_point(aes(size = MP, colour=PeakTissue), fill="black")+
  scale_color_manual(values=c("#E69F00","#56B4E9","#D55E00","#009E73"),labels=c("Roots","Rosettes","Siliques","Stems"),name='Sample Type')+
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=9),axis.text.y = element_text(size=12))+
  facet_grid(rows=vars(Site,Year),cols=vars(Genus),scales="free_x",space="free_x",switch="y",labeller=labeller(Site=Site.labs, Year=Year.labs))+
  scale_x_discrete(expand=c(0.025,1),label=function(x) abbreviate(x, minlength=8))+
  theme(strip.text.y = element_text(size=16, face="bold", vjust = 2),strip.placement="outside")+
  guides(size = guide_legend(title="Peak Prevalence",order=2,nrow=1),
         color=guide_legend(title="Peak Tissue",order=1,nrow=1,override.aes = list(size=3)),
         alpha=guide_legend(title="Peak Mean Relative Abundance",override.aes = list(size=5),nrow=1))+
  theme(strip.text.x = element_text(size = 16,angle=90),strip.background.x= element_rect(colour="black",fill=NA,size=1))+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(strip.placement = "outside")+ theme(panel.spacing = unit(0, "lines"))+ theme(panel.spacing.y = unit(1, "lines"))+
  theme(legend.text=element_text(size=14,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),
        legend.position="top",legend.box="vertical")+
  theme(legend.background = element_rect(color="black", size=.5))

ggsave("~/Documents/KBMP2020_Microbes/Figures/DotPlotA.tiff", plot = grid.arrange(
  DotPlotA,
  layout_matrix = rbind(c(1))), device = NULL, path = NULL,
  scale = 1.7, width = 8, height = 5.5, units = c("in"),
  dpi = 600, limitsize = TRUE)

ggsave("~/Documents/KBMP2020_Microbes/Figures/DotPlotA.png", plot = grid.arrange(
  DotPlotA,
  layout_matrix = rbind(c(1))), device = NULL, path = NULL,
  scale = 1.7, width = 8, height = 5.5, units = c("in"),
  dpi = 600, limitsize = TRUE)

####################################################################################################################################################

Select_DotPlotData <- Select_DotPlotData %>% group_by(ASV,Site,Year) %>% mutate(MaxTissue=Tissue[which.max(Prevalence)])

Select_DotPlotData$MaxTissue
  
select_genera <- c("Pseudomonas","Pantoea","Xanthomonas")

DotPlotB <- NULL
DotPlotB <- ggplot(data= Select_DotPlotData[Select_DotPlotData$Genus %in% select_genera,], aes(x=ASV, y=Tissue, color=MaxTissue)) + 
  geom_point(aes(size = MP, colour=MaxTissue), fill="black")+
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"),labels=c("Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=9),axis.text.y = element_text(size=12))+
  facet_grid(rows=vars(Site,Year),cols=vars(Genus),scales="free_x",space="free_x",switch="y",labeller=labeller(Site=Site.labs, Year=Year.labs))+
  scale_x_discrete(expand=c(0.025,1),label=function(x) abbreviate(x, minlength=8))+
  theme(strip.text.y = element_text(size=16, face="bold", vjust = 2),strip.placement="outside")+
  guides(size = guide_legend(title="Peak Prevalence",order=2,nrow=1),
         color=guide_legend(title="Peak Tissue",order=1,nrow=2,override.aes = list(size=3)),
         alpha=guide_legend(title="Peak Mean Relative Abundance",override.aes = list(size=5),nrow=1))+
  theme(strip.text.x = element_text(size = 16,angle=90),strip.background.x= element_rect(colour="black",fill=NA,size=1))+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(strip.placement = "outside")+ theme(panel.spacing = unit(0, "lines"))+ theme(panel.spacing.y = unit(1, "lines"))+
  theme(legend.text=element_text(size=14,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),
        legend.position="top",legend.box="vertical")+
  theme(legend.background = element_rect(color="black", size=.5))

ggsave("~/Documents/KBMP2020_Microbes/Figures/DotPlotB.tiff", plot = grid.arrange(
  DotPlotB,
  layout_matrix = rbind(c(1))), device = NULL, path = NULL,
  scale = 1.7, width = 4.75, height = 5.5, units = c("in"),
  dpi = 600, limitsize = TRUE)

ggsave("~/Documents/KBMP2020_Microbes/Figures/DotPlotB.png", plot = grid.arrange(
  DotPlotB,
  layout_matrix = rbind(c(1))), device = NULL, path = NULL,
  scale = 1.7, width = 4.75, height = 5.5, units = c("in"),
  dpi = 600, limitsize = TRUE)