# KBMP2020_CommunitySurvey
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Autumn 2019

# The goal here is to compare the soil and plant samples taken at each timepoint in the second year of the study to see what clades distinguish the
# bacterial assemblages in plant tissues from those in the surrounding soil.

# These commands look at the most counted ASVs in the dataset (>1000 counts total; 680 ASVs in 1258 samples).

# For each harvest of a tissue at a timepoint and site, the abundance and relative abundance of ASVs in those samples is calculated and stored with the 
# condition metadata (tissue, stage, site, and year) and taxonomy assignments in a dataframe for plotting.

# The mean relative abundances of ASVs in each class are plotted for soil and plant samples with a 95% confidence interval. 
# A Kruskal-Wallis test indicates significant differences at alpha < 0.0001 (***), 0.001 (**), and 0.01 (*).
# Since roots were the only tissue sampled throughout the whole course of development, they are used to represent plant assemblages in this comparison.

# The total number of sequencing counts for each ASV was also calculated within each harvest condition. The mean of these raw abundances in samples from
# any plant tissue is plotted for each class in order to survey the most common lineages colonizing plants.

####################################################################################################################################################

# Select ASVs with >1000 counts total in the dataset.
AbundantTaxa <- prune_taxa(taxa_sums(physeq)>1000, physeq)
AbundantTaxa <- prune_samples(sample_sums(AbundantTaxa)>0,AbundantTaxa)

# subset_samples(AbundantTaxa,PlantPart=="Soil")
# subset_samples(AbundantTaxa,PlantPart!="Soil")

# Find all the combinations of sampling material (soil or plant tissues), developmental stage, site, and year in the dataset.
ConditionsFrame <- data.frame(unique(sample_data(AbundantTaxa)[,c("Year","Site","Stage","PlantPart")]),row.names=c())
ConditionsFrame <- data.frame(lapply(ConditionsFrame, as.character), stringsAsFactors=FALSE)

# Create empty dataframe for storing the data.
PlottingFrame_wSoil <- NULL
PlottingFrame_wSoil <- data.frame(matrix(ncol = 12, nrow = 0),stringsAsFactors = FALSE)
PlottingFrame_wSoil  <- setNames(data.frame(PlottingFrame_wSoil), c("ASV","Tissue","Stage","Site","Year","RelAbund","Abund",
                                                                    "Genus","Family","Order","Class","Phylum"))
for(m in taxa_names(AbundantTaxa)){ # for each ASV 
  print(m)
  for(row in rownames(ConditionsFrame)){ # for each condition
    pp <- ConditionsFrame[row,"PlantPart"] # get the condition metadata
    st <- ConditionsFrame[row,"Stage"]
    si <- ConditionsFrame[row,"Site"]
    yr <- ConditionsFrame[row,"Year"]
    physeq_partition <-subset_samples(AbundantTaxa, PlantPart==pp & Stage==st & Site==si & Year==yr) # select only the samples in this condition
    physeq_partition <- prune_taxa(taxa_sums(physeq_partition)>0, physeq_partition) # select only the ASVs in samples from this condition
    if(m %in% taxa_names(physeq_partition)){ # if the focal ASV is among them:
      RelAbund <- taxa_sums(physeq_partition)[m]/sum(taxa_sums(physeq_partition)) # find the ASV fractional abundance over all samples in the condition
      Abund <- taxa_sums(physeq_partition)[m] # get the ASV raw abundance (counts of ASV in condition)
    }else{ # if the ASV is not in this condition, set the relative and raw abundances to zero:
      RelAbund <- 0
      Abund <- 0
    }
    Genus <- as.character(tax_table(physeq)[m,"Rank6"]) # get the ASV taxonomy information
    Family <- as.character(tax_table(physeq)[m,"Rank5"])
    Order <- as.character(tax_table(physeq)[m,"Rank4"])
    Class <- as.character(tax_table(physeq)[m,"Rank3"])
    Phylum <- as.character(tax_table(physeq)[m,"Rank2"])
    PlottingFrame_wSoil[nrow(PlottingFrame_wSoil)+1,] <- c(m,pp,st,si,yr,RelAbund,Abund,Genus,Family,Order,Class,Phylum) # store data on this ASV in the frame
  }
}

# write.table(PlottingFrame_wSoil,file="~/Documents/KBMP2020_Microbes/Outputs/PlottingFrame_wSoil", sep="\t", row.names=TRUE,quote=FALSE)

# PlottingFrame_wSoil <- read.table("~/Documents/KBMP2020_Microbes/Outputs/PlottingFrame_wSoil", header=TRUE, sep="\t",row.names=1)

# Set factor levels for plotting.
PlottingFrame_wSoil$Tissue <-factor(PlottingFrame_wSoil$Tissue, levels=c("Soil","Roots","RosLeaves","Stems","CauLeaves","Flowers","Siliques"))
PlottingFrame_wSoil$Stage <-factor(PlottingFrame_wSoil$Stage, levels=c("Soil","TwoLeaf","FourLeaf","SixLeaf","EightLeaf","Flowering","Senescent"))
PlottingFrame_wSoil$Site <-factor(PlottingFrame_wSoil$Site, levels=c("ME","WW"))
PlottingFrame_wSoil$Year <-factor(PlottingFrame_wSoil$Year, levels=c("1","2"))
PlottingFrame_wSoil$Genus <-factor(PlottingFrame_wSoil$Genus)
PlottingFrame_wSoil$Family <-factor(PlottingFrame_wSoil$Family)
PlottingFrame_wSoil$Order <-factor(PlottingFrame_wSoil$Order)
PlottingFrame_wSoil$Class <-factor(PlottingFrame_wSoil$Class)
PlottingFrame_wSoil$Phylum <-factor(PlottingFrame_wSoil$Phylum)

# Set values to numeric for plotting.
PlottingFrame_wSoil$RelAbund <- as.numeric(PlottingFrame_wSoil$RelAbund)
PlottingFrame_wSoil$Abund <- as.numeric(PlottingFrame_wSoil$Abund)

# Assign samples as either Soil or Plant tissue for plotting.
PlottingFrame_wSoil$IsSoil <- ifelse(PlottingFrame_wSoil$Tissue=="Soil","Soil","Plant")
PlottingFrame_wSoil$IsSoil <- factor(PlottingFrame_wSoil$IsSoil,levels=c("Soil","Plant"))

# Rename the "unclassified" class category. 
PlottingFrame_wSoil$Class <- ifelse(PlottingFrame_wSoil$Class=="aggregate unclassified class","unclassified",as.character(PlottingFrame_wSoil$Class))

# check data
# PlottingFrame_wSoil[PlottingFrame_wSoil$Class=="Bacilli" & PlottingFrame_wSoil$Present>0,]

# Set significance indicators for plotting.
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 1), symbols = c("***", "**", "*", "ns"))

# Filter data for known classes and Year 2, during which both plant and soil were sampled at each timepoint. 
# Plot the mean and bootstrapped 95% confidence intervals for the relative abundance of ASVs in each class across conditions in which they appeared.

CS1 <- PlottingFrame_wSoil %>% filter(Year=="2") %>% filter(Class!="unclassified") %>% filter(Tissue%in%c("Soil","Roots")) %>%
  group_by(ASV,IsSoil) %>% mutate(MeanRA = mean(RelAbund))

# examine the distributions of relative abundances
hist(CS1[CS1$IsSoil=="Soil",]$MeanRA)
hist(CS1[CS1$IsSoil=="Plant",]$MeanRA)
gghistogram(CS1,"MeanRA","..count..")+facet_grid(rows=vars(IsSoil),cols=vars(Class))

CommunitySurvey1 <- ggplot(data=CS1,aes(x=Class,y=MeanRA,color=IsSoil,fill=IsSoil)) +
  stat_summary(fun.data="mean_cl_boot", geom="pointrange",size=0.5,alpha=1)+
  #stat_summary(fun.y="mean", geom="bar",size=0.5,alpha=1,position=position_dodge())+
  stat_compare_means(method = "kruskal.test",label="p.signif",show.legend=FALSE,label.y = 0.02, na.rm = TRUE,
                     fontface=2,size=3,hide.ns = TRUE,symnum.args=symnum.args)+
  facet_grid(cols=vars(Phylum),scales="free_x",space="free_x")+
  theme_bw()+
  theme(strip.background=element_blank(),strip.text.x = element_text(size = 20, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+rremove("legend")+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(y="ASV Relative Abundance\n",x=" ")+
  theme(strip.text.y = element_text(size=16, face="bold",vjust=2))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8))+
  scale_x_discrete(expand=c(0.025, 0.5))+
  scale_color_manual(name="Sample Type",values=c("gray","black"))+
  scale_fill_manual(name="Sample Type",values=c("gray","black"))+
  theme(strip.placement = "outside")+theme(plot.margin = unit(c(1,1,1,1), "cm"))+ theme(panel.spacing = unit(3, "lines"))+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5))+
  theme(strip.text.x = element_blank())+
  theme(plot.margin = unit(c(1,3,1,1), "cm"))+ theme(panel.spacing = unit(0.2, "lines"))+ coord_cartesian(ylim=c(0,0.02))

# stat test results in dataframe
# compare_means(data=PlottingFrame_wSoil %>% filter(Class!="unclassified") %>% filter(Year=="2"),formula= RelAbund~IsSoil,method = "kruskal.test",group.by=c("Class"))

# Filter data for known classes and Year 2 plants
# Plot the mean and bootstrapped 95% confidence intervals for the raw abundance of ASVs in each class across conditions in which they appeared.

CS2 <- PlottingFrame_wSoil %>% filter(Year=="2") %>% filter(Class!="unclassified") %>% filter(Tissue!="Soil") %>%
  group_by(ASV) %>% mutate(SumAbund = sum(Abund))

CommunitySurvey2 <- ggplot(data=CS2, aes(x=Class,y=SumAbund,color=IsSoil,fill=IsSoil)) +
  stat_summary(fun.data="mean_cl_boot", geom="pointrange",size=0.5,alpha=1)+
  facet_grid(cols=vars(Phylum),scales="free_x",space="free_x")+
  theme_bw()+
  theme(strip.background=element_blank(),strip.text.x = element_text(size = 20, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+rremove("legend")+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(y="ASV Raw Abundance\n",x=" ")+
  theme(strip.text.y = element_text(size=16, face="bold",vjust=2))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8))+
  scale_x_discrete(expand=c(0.025, 0.5))+
  scale_color_manual(name="Sample Type",values=c("black","black"))+
  scale_fill_manual(name="Sample Type",values=c("black","black"))+
  theme(strip.placement = "outside")+theme(plot.margin = unit(c(1,1,1,1), "cm"))+ theme(panel.spacing = unit(3, "lines"))+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5))+
  theme(strip.text.x = element_blank())+
  theme(plot.margin = unit(c(1,3,1,1), "cm"))+ theme(panel.spacing = unit(0.2, "lines"))+ coord_cartesian(ylim=c(0,4000))

# Plot together
ggsave("~/Documents/KBMP2020_Microbes/Figures/CommunitySurvey.tiff", plot = grid.arrange(
  CommunitySurvey1+ggtitle("A\n")+theme(plot.title = element_text(face="bold",size=20)),
  CommunitySurvey2+ggtitle("B\n")+theme(plot.title = element_text(face="bold",size=20)),
  layout_matrix = rbind(c(1,2))), device = NULL, path = NULL,
  scale = 2, width = 6, height = 3, units = c("in"),
  dpi = 600, limitsize = TRUE)
