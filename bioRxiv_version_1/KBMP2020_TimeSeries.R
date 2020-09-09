

MyData <- PlottingFrame %>% filter(ASV %in% StageIndicators) %>%
  group_by(ASV,Site,Year) %>% mutate(SYSum = sum(Prevalence)) %>% filter(all(SYSum > 0.00)) %>% ungroup %>% 
  group_by(ASV,Site,Year,Tissue) %>% mutate(MP = max(Prevalence), MRA=max(MeanRA))

MyData$Genus <- ifelse(MyData$Genus=="aggregate unclassified genus","unclassified",as.character(MyData$Genus))
MyData$Genus <- ifelse(MyData$Genus=="Burkholderia-Paraburkholderia","Paraburkholderia",as.character(MyData$Genus))

Site.labs<-c("Site: ME", "Site: WW")
names(Site.labs)<-levels(MyData$Site)

Year.labs<-c("Year: 1", "Year: 2")
names(Year.labs)<-levels(MyData$Year)

#MyData <- MyData[MyData$Phylum%in%c("Proteobacteria"),]
#DotPlotData <- DotPlotData[DotPlotData$Genus%in%c("Massilia","Rhizobium","Pseudomonas"),]
#MyData<-MyData[MyData$Genus != "unclassified",]

LevelList <- MyData[order(MyData$Class, MyData$Genus, MyData$ASV),]$Genus

MyData$Genus <- factor(MyData$Genus, levels=unique(LevelList))

LevelListASV <- MyData[order(MyData$Consistent, MyData$Genus, MyData$ASV),]$ASV

MyData$ASV <- factor(MyData$ASV, levels=unique(LevelListASV))


# This selects ASVs among the top 50 for tissue effect that also drove stage effects (based on the threshold set at top 5% of regression coefficients).
# Root and rosette samples are selected and the mean prevalence/abundance of ASVs at each stage, site, and year in each of those tissues is calculated.
#MyData <- PlottingFrame %>% filter(ASV %in% rownames(PlantStageEffects)) %>% filter(Tissue%in%c("Roots","RosLeaves")) %>% filter(ASV %in% TopTissueDrivers$ASV) %>%
#group_by(ASV,Site,Year,Stage) %>% mutate(MP = mean(Prevalence),MRA=mean(MeanRA))

MyData <- MyData[MyData$Tissue%in%c("Roots","RosLeaves"),]
MyData$Tissue <- factor(MyData$Tissue, levels=c("RosLeaves","Roots"))

#Select_MyData<-MyData[MyData$Genus%in%c("Rhizobium","Methylobacterium","Rhodopseudomonas","Aureimonas"),]
#Select_MyData<-Select_MyData[Select_MyData$Tissue%in%c("Roots","RosLeaves"),]

#Select_MyData$Genus<-factor(Select_MyData$Genus,levels=c("Rhizobium","Methylobacterium","Rhodopseudomonas","Aureimonas"))
#Select_MyData$Tissue<-factor(Select_MyData$Tissue,levels=c("RosLeaves","Roots"))

MyData <- MyData[MyData$Genus != "unclassified",]

Select_MyData <- MyData

Select_MyData <- Select_MyData[Select_MyData$Phylum=="Proteobacteria",]

# Find points with prevalence over 70% and get the names of the ASVs with those points
LabelsList<-Select_MyData
LabelsList<- LabelsList%>%
  group_by(ASV,Site,Year,Tissue) %>%
  filter(Prevalence == max(Prevalence)) %>%
  filter(Prevalence>0.7)

Highlights<-unique(LabelsList$ASV)

Highlights<-Highlights[c(1,3)]

# Make a copy of the dataframe for highlighted ASVs (those that reach above 70% prevalence somewhere in the dataset). Each ASV will have a solid line
# with points sized by relative abundance for the time series.
Highlight_MyData<-Select_MyData[Select_MyData$ASV%in%Highlights,]
#Highlight_MyData$Genus<-factor(Highlight_MyData$Genus,levels=c("Rhizobium","Methylobacterium","Rhodopseudomonas","Aureimonas"))

# Make a copy of the dataframe with shortened ASV names for labels at the points of maximum nonzero prevalence for highlighted ASVs
MyLabels<-Select_MyData[Select_MyData$ASV%in%Highlights,]
MyLabels<- MyLabels%>%
  group_by(ASV,Site,Year,Tissue) %>%
  filter(Prevalence == max(Prevalence))
MyLabels<-MyLabels[MyLabels$Prevalence>0,]
MyLabels$ASV<-abbreviate(MyLabels$ASV, minlength=8)


MyPalette <- c("#D5E21AFF","#75D054FF","#20A446FF","#000004FF","#7A0A0AFF","#004C99FF","#FF8000FF","#404688FF","#F56B5CFF","#9C2E7FFF","#20A496FF","#778899FF","grey")

Site.labs<-c("Site: ME", "Site: WW")
names(Site.labs)<-levels(Select_MyData$Site)

Year.labs<-c("Year: 1", "Year: 2")
names(Year.labs)<-levels(Select_MyData$Year)

Tissue.labs<-c("Rosettes", "Roots")
names(Tissue.labs)<-levels(Select_MyData$Tissue)

P_TS <-NULL
P_TS <- ggplot(data=Select_MyData, aes(x=as.factor(Stage), y=Prevalence, group = interaction(ASV,Site,Tissue),color=Genus)) +
  geom_line(data=Highlight_MyData,aes(x=as.factor(Stage), y=Prevalence, group = interaction(ASV,Site,Tissue),color=Genus),size=.75,alpha=1) +
  geom_line(data=Select_MyData,aes(x=as.factor(Stage), y=Prevalence, group = interaction(ASV,Site,Tissue),color=Genus),size=0.75,alpha=0.25) +
  geom_point(data=Highlight_MyData,aes_string(x="Stage",y="Prevalence"),alpha=1)+
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
  scale_color_manual(name="Genus",values=MyPalette,guide=FALSE, breaks = unique(levels(Select_MyData$Genus)))+
  scale_shape_manual(values=c(21,1),guide=FALSE)+
  geom_label_repel(data=MyLabels,size=3,aes_string(x="Stage",y="Prevalence",label = "ASV",color="Genus"),fontface=2,alpha=1,hjust=1,vjust=2,direction="both",show.legend = FALSE)+
  facet_grid(cols = vars(Site,Year),rows=vars(Tissue),labeller=labeller(Site=Site.labs, Year=Year.labs, Tissue=Tissue.labs),switch="y")+
  theme(panel.spacing = unit(2, "lines"))+
  guides(size=guide_legend(title="Relative Abundance",nrow=1,order=1),color=guide_legend(title="Genus",override.aes = list(size=3)))+
  theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5),legend.box="vertical")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.text.y = element_text(size=16, face="bold",margin = margin(0,0,0,1, "cm"),vjust=5))+
  theme(strip.text.x = element_text(size=16, face="bold",margin = margin(0,0,0.1,0, "cm")))+
  theme(strip.placement = "outside")

ggsave("~/Documents/KBMP2020_Microbes/Figures/TimeSeries_Final.tiff", plot = Alpha_TS, device = NULL, path = NULL,
       scale = 1.7, width = 7, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)

ggsave("~/Documents/KBMP2020_Microbes/Figures/TimeSeries_New_TS.png", plot = Alpha_TS, device = NULL, path = NULL,
       scale = 1.7, width = 6, height = 4, units = c("in"),
       dpi = 600, limitsize = TRUE)

Alpha_TS <-NULL
Alpha_TS <- ggplot(data=Select_MyData, aes(x=as.factor(Stage), y=Prevalence, group = interaction(ASV,Site,Tissue),color=Consistent)) +
  #geom_line(data=Highlight_MyData,aes(x=as.factor(Stage), y=Prevalence, group = interaction(ASV,Site,Tissue),color=Genus),size=.75,alpha=1) +
  geom_line(data=Select_MyData,aes(x=as.factor(Stage), y=Prevalence, group = interaction(ASV,Site,Tissue),color=Consistent),size=0.75,alpha=1) +
  geom_point(data=Select_MyData,aes_string(x="Stage",y="Prevalence"),alpha=1,size=2)+
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
  scale_color_manual(values=c("red3","black"),guide=FALSE)+
  #scale_color_manual(name="Genus",values=MyPalette,guide=FALSE)+
  #scale_shape_manual(values=c(21,1),guide=FALSE)+
  #geom_label_repel(data=MyLabels,size=4,aes_string(x="Stage",y="Prevalence",label = "ASV",color="Consistent"),fontface=2,alpha=1,hjust=1,vjust=2,direction="both",show.legend = FALSE)+
  facet_grid(cols = vars(Site,Year),rows=vars(Tissue),labeller=labeller(Site=Site.labs, Year=Year.labs, Tissue=Tissue.labs),switch="y")+
  theme(panel.spacing = unit(2, "lines"))+
  guides(size=guide_legend(title="Relative Abundance",nrow=1,order=1))+
  theme(legend.text=element_text(size=14,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5),legend.box="vertical")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.text.y = element_text(size=16, face="bold",margin = margin(0,0,0,1, "cm"),vjust=5))+
  theme(strip.text.x = element_text(size=16, face="bold",margin = margin(0,0,0.1,0, "cm")))+
  theme(strip.placement = "outside")

ggsave("~/Documents/KBMP2020_Microbes/Figures/TimeSeries_New2_TS.tiff", plot = Alpha_TS, device = NULL, path = NULL,
       scale = 1.7, width = 6, height = 4, units = c("in"),
       dpi = 600, limitsize = TRUE)

ggsave("~/Documents/KBMP2020_Microbes/Figures/DotPlot_TemporalPatterns_Actino.png", plot = grid.arrange(
  A_TS + ggtitle("A: Actinobacteria, Firmicutes, Deinococcus-Thermus, Bacteroidetes, & Chloroflexi\n") + theme(legend.position="none",plot.title = element_text(face="bold",size=16)),
  P_TS + ggtitle("B: Proteobacteria \n")+theme(legend.position="none",plot.title = element_text(face="bold",size=16)),
  layout_matrix = rbind(c(1),
                        c(2)),
  heights = c(10,10),widths=c(10)), device = NULL, path = NULL,
  scale = 2, width = 5.5, height=7, units = c("in"),
  dpi = 600, limitsize = TRUE)


