# KBMP2020_ASVProperties
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Autumn 2019

# Here we explore some properties of the ASVs found in plants:

# PERVASIVE VS. TISSUE-SPECIFIC
# Pervasive ASVs are observed in multiple tissues from at least one site/year while tissue-specific ASVs are observed only in one tissue.

# RECURRENT VS. EPHEMERAL
# Recurrent ASVs are observed at multiple developmental stages from at least one site/year while ephemeral ASVs are observed only at one stage.

# PREVALENCE
# Prevalent ASVs are observed in at least x% of the samples for a condition (site/year/stage/tissue) while rare ASVs are below this threshold.
# Multiple thresholds between 25-100% were tested.

# ABUNDANCE
# High-abundance ASVs have a total abundance for all samples in a condition that is above the median abundance plus one standard deviation for all ASVs 
# in that condition (site/year/stage/tissue), while low-abundance ASVs have an abundance that falls below this threshold.

# CONSISTENCY OF SPATIAL PATTERNS
# ASVs that reach maximum prevalence in the same tissue in each replication of the experiment (site/year) in which they were observed have consistent 
# spatial patterns, whereas those that peak in different tissues across replicates are classified as inconsistent.

# Then we look at which of these properties describe the ASVs driving the association between tissue type and assemblage composition. 

# Finally, some groups of ASVs are removed to assess whether the set of assemblage members driving the tissue pattern has been captured.

# Focus on plant samples
physeq_plants <- subset_samples(physeq,PlantPart!="Soil")
physeq_plants <- prune_taxa(taxa_sums(physeq_plants)>0,physeq_plants)

# Get the columns of the sample data corresponding to the variables of interest and find all the combinations of values that occur in the dataset.
ConditionsFrame <- data.frame(unique(sample_data(physeq_plants)[,c("Year","Site","Stage","PlantPart")]),row.names=c())
ConditionsFrame <- data.frame(lapply(ConditionsFrame, as.character), stringsAsFactors=FALSE)

#################################################################################################################################################
# Pervasive vs. Tissue-specific ASVs
#################################################################################################################################################

# For each year and site, iterate through all the plant parts sampled and get the names of ASVs found in at least one sample of that part.
# Add these names to a list and then make a frequency table based on that list to show how many plant tissues each ASV was found in.
# If the ASV appeared in more than one plant part in at least one site/year, then it will be added to the list of pervasive ASVs.

ASVs_Pervasive <- c()
for(year in unique(ConditionsFrame$Year)){
  for(site in unique(ConditionsFrame$Site)){
    RecurrenceList <- c()
    for(part in unique(ConditionsFrame$PlantPart)){
      physeq_partition <-subset_samples(physeq_plants, PlantPart==part & Site==site & Year==year)
      physeq_partition <- prune_taxa(taxa_sums(physeq_partition)>0, physeq_partition)
      RecurrenceList <- c(taxa_names(physeq_partition),RecurrenceList)
    }
    assign(paste(site,year,"spatially","pervasive",sep="_"),names(table(RecurrenceList)[table(RecurrenceList)>1]))
    print(paste(site,year,"spatially","pervasive",sep="_"))
    ASVs_Pervasive  <- c(names(table(RecurrenceList)[table(RecurrenceList)>1]),ASVs_Pervasive)
    ASVs_Pervasive <- unique(ASVs_Pervasive)
  }
}

# This shows 1766 ASVs that, in at least one site/year combination, occur in more than one of the plant tissues sampled.
# These are the "pervasive" ASVs and the complement are the "tissue-specific" ASVs.

#################################################################################################################################################
# Recurrent vs. Ephemeral ASVs
#################################################################################################################################################

# For each year and site, iterate through all the stages sampled and get the names of ASVs found in at least one sample of that stage.
# Add these names to a list and then make a frequency table based on that list to show how many stages each ASV was found in.
# If the ASV appeared in more than one plant stage in at least one site/year, then it will be added to the list ASVs_Temporally_Recurrent.

ASVs_Recurrent <- c()
for(year in unique(ConditionsFrame$Year)){
  for(site in unique(ConditionsFrame$Site)){
    RecurrenceList <- c()
    for(stage in unique(ConditionsFrame[ConditionsFrame$Year==year & ConditionsFrame$Site==site,]$Stage)){
      physeq_partition <-subset_samples(physeq_plants, Stage==stage & Site==site & Year==year)
      physeq_partition <- prune_taxa(taxa_sums(physeq_partition)>0, physeq_partition)
      RecurrenceList <- c(taxa_names(physeq_partition),RecurrenceList)
    }
    assign(paste(site,year,"temporally","recurrent",sep="_"),names(table(RecurrenceList)[table(RecurrenceList)>1]))
    print(paste(site,year,"temporally","recurrent",sep="_"))
    ASVs_Recurrent  <- c(names(table(RecurrenceList)[table(RecurrenceList)>1]),ASVs_Recurrent)
    ASVs_Recurrent <- unique(ASVs_Recurrent)
  }
}

# This shows 1968 ASVs that, in at least one site/year combination, occur in more than one of the plant stages sampled.
# These are the "recurrent" ASVs and the complement are the "ephemeral" ASVs.

# There is a lot of overlap between pervasive and recurrent ASVs:
length(ASVs_Recurrent)
length(ASVs_Recurrent[ASVs_Recurrent %in% ASVs_Pervasive])/length(ASVs_Recurrent) # 76%
length(ASVs_Pervasive[ASVs_Pervasive %in% ASVs_Recurrent])/length(ASVs_Pervasive) # 85%

#################################################################################################################################################
# Prevalent vs. Rare ASVs
#################################################################################################################################################

# For each condition (year, site, plant part, and stage), select the ASVs present in at least X% of the samples for that condition.
# Repeat this for several thresholds: 25%, 50%, 75%, 100%; store the number of ASVs selected for each condition in a dataframe.
# The lists "core_X" will include all the taxa meeting the prevalence threshold in at least one condition.

ThresholdList <- c(0.25,0.5,0.75,1)
CoreDF <- NULL
CoreDF <- data.frame(matrix(ncol = 6, nrow = 0),stringsAsFactors = FALSE)
CoreDF <- setNames(data.frame(CoreDF), c("Threshold","Year","Site","Stage","Tissue","Core ASV Count"))
for(thr in ThresholdList){
  Core_Taxa <- c()
  for(row in rownames(ConditionsFrame)){
    as.character(ConditionsFrame[row,])
    samps_to_keep = as.character(get_variable(physeq_plants, "PlantPart")) == ConditionsFrame[row,"PlantPart"] &
      as.character(get_variable(physeq_plants, "Stage")) == ConditionsFrame[row,"Stage"] &
      as.character(get_variable(physeq_plants, "Site")) == ConditionsFrame[row,"Site"] &
      as.character(get_variable(physeq_plants, "Year")) == ConditionsFrame[row,"Year"]
    physeq_partition = prune_samples(samps_to_keep, physeq_plants)
    physeq_partition <- prune_taxa(taxa_sums(physeq_partition)>0, physeq_partition)
    PA_table  = transform_sample_counts(physeq_partition, function(x) ifelse(x>0,1,0))
    threshold <- (length(sample_names(physeq_partition))*thr)
    if(length(taxa_sums(PA_table)[taxa_sums(PA_table)>threshold])>0){
      print(paste("core",length(taxa_sums(PA_table)[taxa_sums(PA_table)>threshold]),sep=" "))
      physeq_core <- prune_taxa(taxa_sums(PA_table)>threshold,PA_table)
      Core_Taxa <- c(taxa_names(physeq_core),Core_Taxa)
      Core_Taxa <- unique(Core_Taxa)
      CoreDF[nrow(CoreDF)+1,] <- c(thr,as.character(ConditionsFrame[row,]),length(taxa_sums(PA_table)[taxa_sums(PA_table)>threshold]))
    }else{
      print("no core")
      CoreDF[nrow(CoreDF)+1,] <- c(thr,as.character(ConditionsFrame[row,]),0)
    }
  }
  assign(paste("Core",thr,sep="_"),Core_Taxa)
  print(paste(thr,"Core",length(Core_Taxa),sep="_"))
}

#################################################################################################################################################
# High-abundance vs. Low-abundance ASVs
#################################################################################################################################################

# For each condition (year, site, plant part, and stage), convert ASV abundances for each sample to fractional abundances.
# Select the ASVs that are above median ASV fractional abundance plus one standard deviation for samples in the condition.
# The list Abundant_Taxa will include ASVs that are abundant in at least one condition.

ASVs_Abundant <- c()
for(row in rownames(ConditionsFrame)){
  as.character(ConditionsFrame[row,])
  samps_to_keep = as.character(get_variable(physeq_plants, "PlantPart")) == ConditionsFrame[row,"PlantPart"] &
    as.character(get_variable(physeq_plants, "Stage")) == ConditionsFrame[row,"Stage"] &
    as.character(get_variable(physeq_plants, "Site")) == ConditionsFrame[row,"Site"] &
    as.character(get_variable(physeq_plants, "Year")) == ConditionsFrame[row,"Year"]
  physeq_partition = prune_samples(samps_to_keep, physeq_plants)
  physeq_partition <- prune_taxa(taxa_sums(physeq_partition)>0, physeq_partition)
  threshold <- median(taxa_sums(physeq_partition))+sd(taxa_sums(physeq_partition))
  taxa_sums(physeq_partition)[taxa_sums(physeq_partition)>threshold]
  if(length(taxa_sums(physeq_partition)[taxa_sums(physeq_partition)>threshold])>0){
    print(paste("abundant",length(taxa_sums(physeq_partition)[taxa_sums(physeq_partition)>threshold]),sep=" "))
    physeq_big <- prune_taxa(taxa_sums(physeq_partition)>threshold,physeq_partition)
    ASVs_Abundant <- c(taxa_names(physeq_big),ASVs_Abundant)
    ASVs_Abundant <- unique(ASVs_Abundant)
  }
}

#################################################################################################################################################
# Consistency of ASV spatial patterns across sites and years
#################################################################################################################################################

# Use the vector of ASVs driving the tissue effect from KBMP2020_ASVDrivers.R, "PlantTissueEffects"
# Use data from KBMP2020_PlottingDF.R "PlottingFrame"

# Filter for the peak prevalence of each tissue-discriminating ASV within a site and year.
# Add a tally for the number of times this ASV was observed in the four site/year combinations.
PeakLocations <- PlottingFrame %>% filter(ASV %in% rownames(PlantTissueEffects)) %>%
  group_by(ASV,Site,Year) %>% dplyr::slice(which.max(Prevalence)) %>% filter(sum(Prevalence) > 0.00) %>%
  ungroup() %>% group_by(ASV) %>% add_tally()

PeakLocations$ASV <- as.character(PeakLocations$ASV)

# The ASVs driving the tissue effects peak in each type of tissue:
table(PeakLocations$Tissue)

# For each tissue, look at the ASVs peaking there and determine how many times they peaked there (site/year).
# If this number is equal to the number of times observed (tally from above), add it to the list of ASVs with consistent spatial prevalence patterns.

Tissues <- c("Roots","RosLeaves","Stems","Siliques")
for(t in Tissues){ # for each tissue
  Taxa <- data.frame(table(PeakLocations[PeakLocations$Tissue==t,]$ASV)) # find ASVs that peak in that tissue and get frequency
  colnames(Taxa) <- c("ASV","Freq")
  #Taxa$ASV <- as.character(Taxa$ASV)
  Taxa <- Taxa %>% left_join(by="ASV",PeakLocations) # join back to dataframe with the frequency at which ASVs were observed (n)
  Taxa <- unique(Taxa[Taxa$Tissue==t & Taxa$Freq==Taxa$n,]$ASV) # if ASV peaked in same tissue each time observed, add to list
  assign(paste(t,"Taxa",sep=""),Taxa) # name list
  print(t)
}

ConsistentASVs <- c(RootsTaxa, RosLeavesTaxa, StemsTaxa, SiliquesTaxa)

# Fraction of ASVs driving tissue effect that peak consistently in the same tissue for prevalence: 37%
length(ConsistentASVs)/length(TissueVarASVs)

# Repeat with Relative Abundance instead of Prevalence

PeakLocations <- PlottingFrame %>% filter(ASV %in% rownames(PlantTissueEffects)) %>%
  group_by(ASV,Site,Year) %>% dplyr::slice(which.max(MeanRA)) %>% filter(sum(MeanRA) > 0.00) %>%
  ungroup() %>% group_by(ASV) %>% add_tally()

PeakLocations$ASV <- as.character(PeakLocations$ASV)

table(PeakLocations$Tissue)

Tissues <- c("Roots","RosLeaves","Stems","Siliques")
for(t in Tissues){ # for each tissue
  Taxa <- data.frame(table(PeakLocations[PeakLocations$Tissue==t,]$ASV)) # find ASVs that peak in that tissue and get frequency
  colnames(Taxa) <- c("ASV","Freq")
  #Taxa$ASV <- as.character(Taxa$ASV)
  Taxa <- Taxa %>% left_join(by="ASV",PeakLocations) # join back to dataframe with the frequency at which ASVs were observed (n)
  Taxa <- unique(Taxa[Taxa$Tissue==t & Taxa$Freq==Taxa$n,]$ASV) # if ASV peaked in same tissue each time observed, add to list
  assign(paste(t,"Taxa",sep=""),Taxa) # name list
  print(t)
}

ConsistentASVs <- c(RootsTaxa, RosLeavesTaxa, StemsTaxa, SiliquesTaxa)

# Fraction of ASVs driving tissue effect that peak consistently in the same tissue for prevalence: 19%
length(ConsistentASVs)/length(TissueVarASVs)

# Repeat with ASVs driving stage effect that peak consistently at a timepoint: 24% prevalence, 12.5% relative abundance


#################################################################################################################################################
# How much overlap is there among pervasive, recurrent, prevalent, and abundant ASVs?
#################################################################################################################################################

# Code ASVs with tissue effects based on whether they are prevalent, abundant, or recurrent
PlantTissueEffects$Recurrent<-ifelse(rownames(PlantTissueEffects)%in% ASVs_Recurrent,"Y","N")
PlantTissueEffects$Pervasive<-ifelse(rownames(PlantTissueEffects)%in% ASVs_Pervasive,"Y","N")
PlantTissueEffects$Prevalent<-ifelse(rownames(PlantTissueEffects)%in% Core_0.75,"Y","N")
PlantTissueEffects$Abundant<-ifelse(rownames(PlantTissueEffects)%in% ASVs_Abundant,"Y","N")

# Proportion of ASVs with tissue effects that are recurrent: 92%
length(rownames(PlantTissueEffects[PlantTissueEffects$Recurrent=="Y",]))/length(rownames(PlantTissueEffects))

# Proportion of ASVs with tissue effects that are pervasive: 94%
length(rownames(PlantTissueEffects[PlantTissueEffects$Pervasive=="Y",]))/length(rownames(PlantTissueEffects))

# Proportion of ASVs with tissue effects that are prevalent (threshold 75%): 13%
length(rownames(PlantTissueEffects[PlantTissueEffects$Prevalent=="Y",]))/length(rownames(PlantTissueEffects))

# Proportion of ASVs with tissue effects that are abundant (higher than median+1sd relative abundance): 79%
length(rownames(PlantTissueEffects[PlantTissueEffects$Abundant=="Y",]))/length(rownames(PlantTissueEffects))

Properties_PlantTissueEffects <- PlantTissueEffects[,c("Recurrent","Pervasive","Prevalent","Abundant")] 
Properties_PlantTissueEffects <- ifelse(Properties_PlantTissueEffects == "Y",1,0)

venn(Properties_PlantTissueEffects)

palette<-c("#FF0000","#6633CC","#B31A66","lightgray")
Euler <- c("Recurrent" = 0, "Pervasive" = 8, "Prevalent" = 0, "Abundant" = 11, 
           "Recurrent&Pervasive"= 70, "Recurrent&Abundant" = 8, "Recurrent&Prevalent" = 0, "Pervasive&Prevalent" = 0, "Pervasive&Abundant" = 11, "Prevalent&Abundant" = 0,
           "Recurrent&Pervasive&Abundant" = 248, "Recurrent&Pervasive&Prevalent" = 2, "Pervasive&Prevalent&Abundant" = 0, "Abundant&Prevalent&Recurrent" = 0,
           "Recurrent&Pervasive&Prevalent&Abundant" = 50)
EulerPlot <- plot(euler(Euler, shape = "circle"),fontsize=20,quantities = list(fontsize = 14,font=2),fills=palette,lty = c(1,1,2,3))

ggsave("~/Documents/KBMP2020_Microbes/Figures/Venn_TissueDrivingASVs_Properties.tiff", plot = EulerPlot, device = NULL, path = NULL,
       scale = 1.5, width = 5, height = 4, units = c("in"),
       dpi = 600, limitsize = TRUE)

#################################################################################################################################################
# Which ASVs driving the tissue effect are pervasive/recurrent/prevalent/abundant?
#################################################################################################################################################

PlantTissueEffects$Pervasive <- factor(PlantTissueEffects$Pervasive)
PlantTissueEffects$Recurrent <- factor(PlantTissueEffects$Recurrent)
PlantTissueEffects$Abundant <- factor(PlantTissueEffects$Abundant)
PlantTissueEffects$Prevalent <- factor(PlantTissueEffects$Prevalent)

barplot <- ggbarplot(PlantTissueEffects %>% mutate(ASV = rownames(PlantTissueEffects)) %>% arrange(EffectSums) %>% top_n(50,EffectSums), 
                     "ASV", "EffectSums", position = position_dodge(1.5),orientation = "horiz",fill ="Pervasive",color="Prevalent",linetype="Abundant",size=0.75)+
  theme_bw()+
  theme(strip.background=element_blank(),strip.text.x = element_text(size = 16, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(x="Amplicon sequence variants\n",y="\nInfluence on tissue effect\n(sum of PERMANOVA coefficient absolute values)")+
  scale_color_manual(values=c("blue","black"),name="Prevalent\n(75%)")+
  scale_fill_manual(values=c("gray80","gray40"),name="Pervasive\n(> 1 tissue)",drop=FALSE)+
  scale_linetype_manual(values=c(6,1),name="Abundant\n(Median + 1std)",drop=FALSE)+
  #scale_alpha_manual(drop=FALSE,values=c(0.2,1),name="Recurrent\n(> 1 timepoint)")+
  theme(legend.position = "right")+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"))+
  guides(fill = guide_legend(order=1,nrow=2),alpha=guide_legend(override.aes = list(fill="gray40"),order=2),
         linetype = guide_legend(override.aes = list(fill = NA, col = "black"),order=4,nrow=2),
         color=guide_legend(override.aes = list(fill = NA),order=3,nrow=2))+
  scale_x_discrete(expand=c(0.025,1),label=function(x) abbreviate(x, minlength=8))+
  facet_grid(rows=vars(Genus),scales="free_y",space="free_y",switch="y")+
  theme(strip.text.y = element_text(size = 12,angle=180,face="bold"))+
  theme(strip.placement = "outside") + theme(panel.spacing = unit(0.1, "lines"))

ggsave("~/Documents/KBMP2020_Microbes/Figures/Properties_Bar_PERMANOVA_Coefficients.tiff", plot = barplot, device = NULL, path = NULL,
       scale = 1.2, width = 6, height = 9, units = c("in"),
       dpi = 600, limitsize = TRUE)

#################################################################################################################################################
# Which ASVs driving the tissue effect peak consistently?
#################################################################################################################################################

# For barplot, assign the top 50 tissue driving ASVs to the tissue where they consistently peak.
TopTissueDrivers <- PlantTissueEffects %>% mutate(ASV = rownames(PlantTissueEffects)) %>% arrange(EffectSums) %>% top_n(50,EffectSums)
TopTissueDrivers$PeakTissue <- NA
TopTissueDrivers$PeakTissue <- ifelse(TopTissueDrivers$ASV %in% RootsTaxa,"Roots",TopTissueDrivers$PeakTissue)
TopTissueDrivers$PeakTissue <- ifelse(TopTissueDrivers$ASV %in% RosLeavesTaxa,"Rosettes",TopTissueDrivers$PeakTissue)
TopTissueDrivers$PeakTissue <- ifelse(TopTissueDrivers$ASV %in% StemsTaxa,"Stems",TopTissueDrivers$PeakTissue)
TopTissueDrivers$PeakTissue <- ifelse(TopTissueDrivers$ASV %in% SiliquesTaxa,"Siliques",TopTissueDrivers$PeakTissue)
TopTissueDrivers$PeakTissue <- ifelse(is.na(TopTissueDrivers$PeakTissue),"Inconsistent",TopTissueDrivers$PeakTissue)

TopTissueDrivers$PeakTissue <- factor(TopTissueDrivers$PeakTissue,levels=c("Roots","Rosettes","Stems","Siliques","Inconsistent"))
TopTissueDrivers$Genus <- ifelse(TopTissueDrivers$Genus=="aggregate unclassified genus","unclassified",TopTissueDrivers$Genus)

barplot <- ggbarplot(TopTissueDrivers, 
                     "ASV", "EffectSums", position = position_dodge(1.5),orientation = "horiz",fill ="PeakTissue",color=NA,size=1)+
  theme_bw()+
  theme(strip.background=element_blank(),strip.text.x = element_text(size = 16, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(x="Amplicon sequence variants\n",y="\nInfluence on tissue effect\n(sum of PERMANOVA coefficient absolute values)")+
  scale_fill_manual(values=c("#E69F00","#56B4E9","#009E73","#D55E00","black"),name="Tissue with Max Prevalence",labels=c("Roots","Rosettes","Stems","Siliques","Inconsistent"))+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  scale_x_discrete(expand=c(0.025,1),label=function(x) abbreviate(x, minlength=8))+
  facet_grid(rows=vars(Genus),scales="free_y",space="free_y",switch="y")+
  theme(strip.text.y = element_text(size = 12,angle=180,face="bold"))+  #,strip.background= element_rect(colour="black",fill=NA,size=1))+
  theme(strip.placement = "outside") + theme(legend.position="top") + theme(panel.spacing = unit(0.2, "lines"))

ggsave("~/Documents/KBMP2020_Microbes/Figures/PeakTissue_Bar_PERMANOVA_Coefficients.tiff", plot = barplot, device = NULL, path = NULL,
       scale = 1.4, width = 7, height = 9, units = c("in"),
       dpi = 600, limitsize = TRUE)

#################################################################################################################################################
# Repeat PERMANOVA to see how the R^2 for tissue and stage changes with different ASVs removed
#################################################################################################################################################

# Get the lists of samples that remain when datasets are pruned for Pervasive ASVs, Recurrent ASVs, and Driving ASVs ().
ASVs_Driving <- unique(unlist(DrivingASVs))
MeetSomeCriteria <- unique(c(ASVs_Pervasive, ASVs_Recurrent, ASVs_Driving))

SampleSet <- subset_taxa(physeq_factored_plants,!(taxa_names(physeq_factored_plants)%in%MeetSomeCriteria))
SampleSet <- prune_samples(sample_sums(SampleSet)>0,SampleSet)
SampleSet <- sample_names(SampleSet)

# For each list of ASVs, remove it and repeat the PERMANOVA with the constant sample set

SelectionLists <- list(TissueVarASVs, ASVs_Recurrent, ASVs_Pervasive)
SelectionSets <- c("TissueVarASVs","ASVs_Recurrent", "ASVs_Pervasive")

ResultsDF <- NULL
ResultsDF <- data.frame(matrix(ncol = 8, nrow = 0),stringsAsFactors = FALSE)
ResultsDF <- setNames(data.frame(ResultsDF), c("Set","Var","Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)"))

Palette<-c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
for(li in 1:length(SelectionLists)){
  physeq_unselect <- prune_samples(sample_names(physeq_factored_plants)%in%SampleSet,physeq_factored_plants)
  physeq_unselect <- prune_taxa(!(taxa_names(physeq_unselect)%in%SelectionLists[[li]]),physeq_unselect)
  physeq_unselect <- prune_taxa(taxa_sums(physeq_unselect)>0,physeq_unselect)
  raup <- phyloseq::distance(physeq_unselect, method = "raup")
  permanova_result <- adonis(raup ~Year/Stage/PlantPart + MiSeqRun/SamplePlate + Site,data.frame(sample_data(physeq_unselect)),permutations=999)
  write.table(data.frame(permanova_result$aov.tab),file=paste("~/Documents/KBMP2020_Microbes/PERMANOVAs/",SelectionSets[li],sep=""), sep="\t", row.names=TRUE,quote=FALSE)
  for(row in rownames(permanova_result$aov.tab)){
    ResultsDF[nrow(ResultsDF)+1,] <- c(SelectionSets[li],row,as.character(permanova_result$aov.tab[row,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")]))
  }
  #PCoA_raup_All  <- ordinate(physeq_unselect, "PCoA", distance=raup)
  PCoA_raup_All  <- ape::pcoa(raup,correction="cailliez")
  pcoa <- plot_ordination(physeq_unselect, PCoA_raup_All, color="PlantPart",shape="Stage") + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    font("ylab", size = 16, color = "black",face = "bold")+
    font("xlab", size = 16, color = "black",face = "bold")+
    font("xy.text", size = 12, color = "black", face = "bold")+
    geom_point(size=2)+
    scale_colour_manual(values=Palette,labels=c("Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+ 
    scale_shape_manual(values=c("\u25AC","\u25D6","\u25A0","\u25B2","\u25AE","\u25C6"),labels=c('Two Leaf','Four Leaf','Six Leaf','Eight Leaf','Flowering','Senescent'),name='Development Stage')+
    theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),legend.direction="horizontal")+
    theme(legend.background = element_rect(color="black", size=.5))+coord_cartesian(ylim=c(-1,1),xlim=c(-1,1))
  assign(paste("PCoA",SelectionSets[li],sep="_"),pcoa)
}

write.table(ResultsDF,file="~/Documents/KBMP2020_Microbes/PERMANOVAs/Plant_ASV_Removals", sep="\t", row.names=TRUE,quote=FALSE)

PCoA_TissueVarASVs <- PCoA_TissueVarASVs+ guides(color=guide_legend(override.aes = list(size=5)),shape=guide_legend(override.aes = list(size=5)))
PCoA_ASVs_Recurrent
PCoA_ASVs_Pervasive

legendCol <- cowplot::get_legend(PCoA_TissueVarASVs+theme(legend.box="horizontal"))

ggsave("~/Documents/KBMP2020_Microbes/Figures/PCoA_Removals.tiff", plot = 
         grid.arrange(
           legendCol,
           PCoA_TissueVarASVs+ labs(x=" ", y=" ") + theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("A: Without Top 5%\n Tissue-Discriminating ASVs")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")), 
           PCoA_ASVs_Recurrent+ labs(x=" ", y=" ") + theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("B: Without Recurrent ASVs\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
           PCoA_ASVs_Pervasive+ labs(x=" ", y=" ") + theme(legend.position = "none",plot.title = element_text(face="bold",size=16))+ggtitle("C: Without Pervasive ASVs\n")+theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
           widths = c(4, 4, 4),
           heights = c(1,4),
           layout_matrix = rbind(c(1),
                                 c(2,3,4))
         ), 
       device = NULL, path = NULL,
       scale = 1.5, width = 8, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)


