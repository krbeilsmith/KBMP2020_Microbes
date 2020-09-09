# In what tissue does ASV prevalence peak?
PeakLocations <- PlottingFrame %>% filter(ASV %in% PhylloRhizoIndicators) %>%
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

DotPlotData <- PlottingFrame %>% filter(ASV %in% PhylloRhizoIndicators) %>%
  group_by(ASV,Site,Year) %>% mutate(SYSum = sum(Prevalence)) %>% filter(all(SYSum > 0.00)) %>% ungroup %>% 
  group_by(ASV,Site,Year,Tissue) %>% mutate(MP = max(Prevalence), MRA=max(MeanRA))

# DotPlotData <- DotPlotData[DotPlotData$Genus%in%c("Massilia"),]

DotPlotData$PeakTissue <- ifelse(DotPlotData$ASV%in%RootsTaxa, "Roots", 
                          ifelse(DotPlotData$ASV%in%RosLeavesTaxa, "RosLeaves",
                          ifelse(DotPlotData$ASV%in%StemsTaxa, "Stems",
                          ifelse(DotPlotData$ASV%in%SiliquesTaxa, "Siliques","Inconsistent"))))

DotPlotData$Genus <- ifelse(DotPlotData$Genus=="aggregate unclassified genus","unclassified",as.character(DotPlotData$Genus))
DotPlotData$Genus <- ifelse(DotPlotData$Genus=="Burkholderia-Paraburkholderia","Paraburkholderia",as.character(DotPlotData$Genus))

MyData <- DotPlotData

#MyData <- MyData[MyData$Class%in%c("Actinobacteria","Alphaproteobacteria","Betaproteobacteria","Deinococcus-Thermus"),]
#DotPlotData <- DotPlotData[DotPlotData$Genus%in%c("Massilia","Rhizobium","Pseudomonas"),]

MyData$PeakTissue<-factor(MyData$PeakTissue, levels=c("Inconsistent","Roots","RosLeaves","Siliques","Stems"))

LevelList <- MyData[order(MyData$Class, MyData$Genus, MyData$ASV),]$Genus

MyData$Genus <- factor(MyData$Genus, levels=unique(LevelList))

MyData <- MyData[MyData$Genus != "unclassified",]

Site.labs<-c("Site: ME", "Site: WW")
names(Site.labs)<-levels(DotPlotData$Site)

Year.labs<-c("Year: 1", "Year: 2")
names(Year.labs)<-levels(DotPlotData$Year)

unique(MyData$Phylum)

Plot_Filter_A <- NULL
Plot_Filter_A <- ggplot(data= MyData[MyData$Phylum!="Proteobacteria",], aes(x=ASV, y=Tissue)) + 
  geom_point(aes(size = MP, colour=PeakTissue), fill="black")+
  scale_color_manual(values=c("grey","#E69F00","#56B4E9","#D55E00","#009E73"),labels=c("Inconsistent","Roots","Rosettes","Siliques","Stems"),name='Sample Type')+
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12),axis.text.y = element_text(size=12))+
  facet_grid(rows=vars(Site,Year),cols=vars(Genus),scales="free_x",space="free_x",switch="y",labeller=labeller(Site=Site.labs, Year=Year.labs))+
  scale_x_discrete(expand=c(0.025,1),label=function(x) abbreviate(x, minlength=8))+
  theme(strip.text.y = element_text(size=16, face="bold", vjust = 2),strip.placement="outside")+
  guides(size = guide_legend(title="Peak Prevalence",order=1,nrow=1),
         color=FALSE,alpha=guide_legend(title="Peak Mean Relative Abundance",override.aes = list(size=5),nrow=1))+
  theme(strip.text.x = element_text(size = 16,angle=90),strip.background.x= element_rect(colour="black",fill=NA,size=1))+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(strip.placement = "outside")+ theme(panel.spacing = unit(0, "lines"))+ theme(panel.spacing.y = unit(1, "lines"))+
  theme(legend.text=element_text(size=14,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),
        legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5))

Plot_Filter_B <- NULL
Plot_Filter_B <- ggplot(data= MyData[MyData$Phylum=="Proteobacteria",], aes(x=ASV, y=Tissue)) + 
  geom_point(aes(size = MP, colour=PeakTissue), fill="black")+
  scale_color_manual(values=c("grey","#E69F00","#56B4E9","#D55E00","#009E73"),labels=c("Inconsistent","Roots","Rosettes","Siliques","Stems"),name='Sample Type')+
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12),axis.text.y = element_text(size=12))+
  facet_grid(rows=vars(Site,Year),cols=vars(Genus),scales="free_x",space="free_x",switch="y",labeller=labeller(Site=Site.labs, Year=Year.labs))+
  scale_x_discrete(expand=c(0.025,1),label=function(x) abbreviate(x, minlength=8))+
  theme(strip.text.y = element_text(size=16, face="bold", vjust = 2),strip.placement="outside")+
  guides(size = guide_legend(title="Peak Prevalence",order=1,nrow=1),
         color=FALSE,alpha=guide_legend(title="Peak Mean Relative Abundance",override.aes = list(size=5),nrow=1))+
  theme(strip.text.x = element_text(size = 16,angle=90),strip.background.x= element_rect(colour="black",fill=NA,size=1))+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(strip.placement = "outside")+ theme(panel.spacing = unit(0, "lines"))+ theme(panel.spacing.y = unit(1, "lines"))+
  theme(legend.text=element_text(size=14,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),
        legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5))


legendCol <- cowplot::get_legend(Plot_Filter_A)
ggsave("~/Documents/KBMP2020_Microbes/Figures/DotPlot_SpatialPatterns_Actino.png", plot = grid.arrange(
  legendCol,
  Plot_Filter_A + ggtitle("A: Actinobacteria, Firmicutes, & Deinococcus-Thermus\n") + theme(legend.position = "none",plot.title = element_text(face="bold",size=16)),
  Plot_Filter_B + ggtitle("B: Proteobacteria \n")+theme(legend.position = "none",plot.title = element_text(face="bold",size=16)),
  layout_matrix = rbind(c(1),
                        c(2),
                        c(3)),
  heights = c(1,7,7),widths=c(5)), device = NULL, path = NULL,
  scale = 2, width = 5.5, height=11, units = c("in"),
  dpi = 600, limitsize = TRUE)

