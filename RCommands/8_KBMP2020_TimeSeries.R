# KBMP2020_TimeSeries
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

GenusFrame <- read.table("~/Documents/KBMP2020_Microbes/Outputs/PlottingFrame_Genera", header=TRUE, sep="\t",row.names=1)

# Set factor levels for plotting
GenusFrame$Tissue <-factor(GenusFrame$Tissue, levels=c("Roots","RosLeaves","Stems","CauLeaves","Flowers","Siliques"))
GenusFrame$Stage <-factor(GenusFrame$Stage, levels=c("TwoLeaf","FourLeaf","SixLeaf","EightLeaf","Flowering","Senescent"))
GenusFrame$Site <-factor(GenusFrame$Site, levels=c("ME","WW"))
GenusFrame$Year <-factor(GenusFrame$Year, levels=c("1","2"))
GenusFrame$Family <-factor(GenusFrame$Family)
GenusFrame$Order <-factor(GenusFrame$Order)
GenusFrame$Class <-factor(GenusFrame$Class)
GenusFrame$Phylum <-factor(GenusFrame$Phylum)

# Make prevalence and relative abundance data numeric for plotting
GenusFrame$Prevalence <- as.numeric(GenusFrame$Prevalence)
GenusFrame$MeanRA <- as.numeric(GenusFrame$MeanRA)
GenusFrame$sdRA <- as.numeric(GenusFrame$sdRA)

####################################################################################################################################################

# In what stage does ASV prevalence peak?
PeakStages <- PlottingFrame %>% filter(ASV %in% rownames(Stage_Indicators)) %>%
  group_by(ASV,Site,Year) %>% dplyr::slice(which.max(Prevalence))

Stages <- c("TwoLeaf","FourLeaf","SixLeaf","EightLeaf","Flowering","Senescent")
for(s in Stages){ # for each stage
  Taxa <- data.frame(table(PeakStages[PeakStages$Stage==s,]$ASV)) # find ASVs that peak in that stage and get frequency (1-4)
  colnames(Taxa) <- c("ASV","Freq")
  Taxa$ASV <- as.character(Taxa$ASV)
  Taxa <- unique(Taxa[Taxa$Freq==4,]$ASV) # if observed in all 4 site/year combinations, add to list
  assign(paste(s,"Taxa",sep=""),Taxa) # name list
  print(s)
}

ConsistentASVs <- c(TwoLeafTaxa, FourLeafTaxa, SixLeafTaxa, EightLeafTaxa, FloweringTaxa, SenescentTaxa)

# % of Stage indicator ASVs that consistently peak at one point: 23.1%
(length(ConsistentASVs) / length(rownames(Stage_Indicators)))*100

# % of consistently peaking ASVs for each stage
(length(FloweringTaxa) / length(ConsistentASVs))*100 # flowering 29.03%
(length(SenescentTaxa) / length(ConsistentASVs))*100 # senescent 70.9%

####################################################################################################################################################

ASV_series <- PlottingFrame %>% filter(ASV %in% rownames(Stage_Indicators))

ASV_series$Genus <- ifelse(ASV_series$Genus=="aggregate unclassified genus","unclassified",as.character(ASV_series$Genus))
ASV_series <- ASV_series[ASV_series$Genus != "unclassified",]

LevelList <- ASV_series[order(ASV_series$Class, ASV_series$Genus, ASV_series$ASV),]$Genus
ASV_series$Genus <- factor(ASV_series$Genus, levels=unique(LevelList))

ASV_series <- ASV_series[ASV_series$Tissue%in%c("Roots","RosLeaves"),]
ASV_series$Tissue <- factor(ASV_series$Tissue, levels=c("RosLeaves","Roots"))

Select_ASV_series <- ASV_series

Select_ASV_series <- Select_ASV_series[Select_ASV_series$Genus=="Massilia" | Select_ASV_series$Genus=="Methylobacterium",]

####################################################################################################################################################

Genus_series <- GenusFrame

Genus_series$Genus <- ifelse(Genus_series$Genus=="aggregate unclassified genus","unclassified",as.character(Genus_series$Genus))

LevelList <- Genus_series[order(Genus_series$Class, Genus_series$Genus),]$Genus
Genus_series$Genus <- factor(Genus_series$Genus, levels=unique(LevelList))

Genus_series <- Genus_series[Genus_series$Tissue%in%c("Roots","RosLeaves"),]
Genus_series$Tissue <- factor(Genus_series$Tissue, levels=c("RosLeaves","Roots"))

Select_Genus_series <- Genus_series[Genus_series$Genus=="Massilia" | Genus_series$Genus=="Methylobacterium",]

####################################################################################################################################################

Site.labs<-c("Site: ME", "Site: WW")
names(Site.labs)<-levels(Select_ASV_series$Site)

Year.labs<-c("Year: 1", "Year: 2")
names(Year.labs)<-levels(Select_ASV_series$Year)

Tissue.labs<-c("Rosettes", "Roots")
names(Tissue.labs)<-levels(Select_ASV_series$Tissue)

MyPalette <- c("#D5E21AFF","#75D054FF","#20A446FF","#000004FF","#7A0A0AFF","#004C99FF","#FF8000FF","#404688FF","#F56B5CFF","#9C2E7FFF","#20A496FF","#778899FF")

TimeSeriesA <-NULL
TimeSeriesA <- ggplot(data=Select_Genus_series, aes(x=as.factor(Stage), y=Prevalence, group = interaction(Genus,Site,Year,Tissue),color=Genus)) +
  geom_line(data=Select_ASV_series,aes(x=as.factor(Stage), y=Prevalence, group = interaction(ASV,Site,Year,Tissue),color=Genus),alpha=0.15,size=1) +
  geom_line(data=Select_Genus_series,aes(x=as.factor(Stage), y=Prevalence, group = interaction(Genus,Site,Tissue),color=Genus),alpha=1,linetype=3,size=1) +
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
  scale_color_manual(name="Genus",values=MyPalette[c(4,9)],guide=FALSE, breaks = unique(levels(Select_Genus_series$Genus)))+
  scale_shape_manual(values=c(21,1),guide=FALSE)+
  facet_grid(cols = vars(Site,Year),rows=vars(Tissue),labeller=labeller(Site=Site.labs, Year=Year.labs, Tissue=Tissue.labs),switch="y")+
  theme(panel.spacing = unit(2, "lines"))+
  guides(color=guide_legend(title="Genus",override.aes = list(size=3)))+
  theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5),legend.box="vertical")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.text.y = element_text(size=16, face="bold",margin = margin(0,0,0,1, "cm"),vjust=5))+
  theme(strip.text.x = element_text(size=16, face="bold",margin = margin(0,0,0.1,0, "cm")))+
  theme(strip.placement = "outside")

ggsave("~/Documents/KBMP2020_Microbes/Figures/TimeSeriesA.tiff", plot = TimeSeriesA, device = NULL, path = NULL,
       scale = 1.7, width = 7, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)

ggsave("~/Documents/KBMP2020_Microbes/Figures/TimeSeriesA.png", plot = TimeSeriesA, device = NULL, path = NULL,
       scale = 1.7, width = 7, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)

####################################################################################################################################################

# Find points with prevalence over 70% and get the names of the ASVs with those points
LabelsList <- Select_ASV_series
LabelsList <- LabelsList %>%
  group_by(ASV,Site,Year,Tissue) %>%
  filter(Prevalence == max(Prevalence)) %>%
  filter(Prevalence > 0.7)

# Select ASVs to highlight in plot
Highlights <- unique(LabelsList$ASV)
Highlights <- Highlights[c(1,4,6,7)]
Highlight_ASV_series <- Select_ASV_series[Select_ASV_series$ASV%in%Highlights,]

# Make a copy of the dataframe for highlighted ASVs (those that reach above 70% prevalence somewhere in the dataset). Each ASV will have a solid line
# Points can be sized by relative abundance for the time series. ASV time series lines can be labelled.

# Shorten ASV identifiers for figure legend
ASV_Names <- unique(paste(as.character(Highlight_ASV_series$Genus),as.character(abbreviate(Highlight_ASV_series$ASV, minlength=8)),sep=" "))
# Order names for legend
ASV_Names <- ASV_Names[c(2,3,4,1)]
# Line type ordering for legend
override.linetype <- c(1, 3, 2, 4)

TimeSeriesB <-NULL
TimeSeriesB <- ggplot(data=Highlight_ASV_series, aes(x=as.factor(Stage), y=Prevalence, group = interaction(ASV,Site,Year,Tissue),color=interaction(Genus,ASV),linetype=interaction(Genus,ASV))) +
  geom_line(data=Highlight_ASV_series,aes(x=as.factor(Stage), y=Prevalence, group = interaction(ASV,Site,Tissue),color=interaction(Genus,ASV),linetype=interaction(Genus,ASV)),size=.75,alpha=1) +
  geom_point(data=Highlight_ASV_series,aes_string(x="Stage",y="Prevalence"),alpha=1)+
  theme_bw()+
  theme(strip.background=element_blank(),strip.text.x = element_text(size = 16, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+#+rremove("legend")+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(y="Prevalence",x=" ")+
  scale_linetype_manual(values=c(1,3,2,4),guide=FALSE,labels=Names)+
  scale_x_discrete(labels=c("Two Leaf", "Four Leaf", "Six Leaf", "Eight Leaf", "Flowering", "Senescent"))+
  scale_y_continuous(breaks=c(0.5,1),limits=c(0,1))+
  scale_color_manual(name="Genus",values=MyPalette[c(4,9,9,4)],guide=FALSE,labels=ASV_Names)+
  facet_grid(cols = vars(Site,Year),rows=vars(Tissue),labeller=labeller(Site=Site.labs, Year=Year.labs, Tissue=Tissue.labs),switch="y")+
  theme(panel.spacing = unit(2, "lines"))+
  guides(color=guide_legend(nrow=2,title="ASV",override.aes = list(size=1.5,linetype = override.linetype)))+
  theme(legend.key.width = unit(5, "line"))+
  theme(legend.text=element_text(size=12,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5),legend.box="vertical")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.text.y = element_text(size=16, face="bold",margin = margin(0,0,0,1, "cm"),vjust=5))+
  theme(strip.text.x = element_text(size=16, face="bold",margin = margin(0,0,0.1,0, "cm")))+
  theme(strip.placement = "outside")

ggsave("~/Documents/KBMP2020_Microbes/Figures/TimeSeriesB.tiff", plot = TimeSeriesB, device = NULL, path = NULL,
       scale = 1.7, width = 7, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)

ggsave("~/Documents/KBMP2020_Microbes/Figures/TimeSeriesB.png", plot = TimeSeriesB, device = NULL, path = NULL,
       scale = 1.7, width = 7, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)

####################################################################################################################################################

ggsave("~/Documents/KBMP2020_Microbes/Figures/TimeSeries_Full.tiff", plot = grid.arrange(
  TimeSeriesA+ggtitle("A")+theme(plot.title = element_text(face="bold",size=20)),
  TimeSeriesB+ggtitle("B")+theme(plot.title = element_text(face="bold",size=20)),
  layout_matrix = rbind(c(1),
                        c(2))), device = NULL, path = NULL,
  scale = 1.7, width = 5.5, height = 6.5, units = c("in"),
  dpi = 600, limitsize = TRUE)

ggsave("~/Documents/KBMP2020_Microbes/Figures/TimeSeries_Full.png", plot = grid.arrange(
  TimeSeriesA+ggtitle("A")+theme(plot.title = element_text(face="bold",size=20)),
  TimeSeriesB+ggtitle("B")+theme(plot.title = element_text(face="bold",size=20)),
  layout_matrix = rbind(c(1),
                        c(2))), device = NULL, path = NULL,
  scale = 1.7, width = 5.5, height = 6.5, units = c("in"),
  dpi = 600, limitsize = TRUE)