# KBMP2020_Membership
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Autumn 2019

##########################################################################################################################################################

# The ASVs present vary a lot between plant samples, so the work below examines what proportion of the community members are shared between each sample and
# randomly selected plant samples from the dataset and whether that proportion is higher when each sample is compared to a group of samples matched for
# the tissue type, stage, site, or year sampled.

# Select only the plant samples.
physeq_forMembership <-subset_samples(physeq, PlantPart!="Soil")

# Make a dataframe for the results.
Membership <- NULL
Membership <- data.frame(matrix(ncol = 13, nrow = 0),stringsAsFactors = FALSE)
Membership <- setNames(data.frame(Membership), c("Samp","PropShared_PlantID","NumShared_PlantID","PropShared_Stage","NumShared_Stage","PropShared_Site","NumShared_Site","PropShared_Part","NumShared_Part","PropShared_Year","NumShared_Year","PropShared_Random","NumShared_Random"))
NextRow<-NULL

SampleProperties <- c("Stage","Site","Year","PlantPart") # each property by which the samples can be grouped

for(SampID in sample_names(physeq_forMembership)){ # for each plant sample in the dataset
  
  print(SampID)
  SampTaxa <- subset_samples(physeq_forMembership,sample_names(physeq_forMembership) %in% SampID) # prune the data to just the focal sample
  SampTaxa <- prune_taxa(taxa_sums(SampTaxa)>0,SampTaxa) # prune the data to just the ASVs in the focal sample
  
  for(Property in SampleProperties){ # for each property by which the samples can be grouped
    ### Find samples with the same property (stage, for example)
    SampProperty <- as.character(sample_data(physeq_forMembership)[SampID,Property][[1]]) # get the stage of the focal sample
    # prune the data to just samples sharing the same stage
    samps_to_keep = as.character(get_variable(physeq_forMembership, Property)) == SampProperty & as.character(rownames(sample_data(physeq_forMembership))) != SampID
    Sharing_Property = prune_samples(samps_to_keep, physeq_forMembership)
    
    ### Randomly draw from these samples 5x and find the mean number and proportion of ASVs they share with the focal sample
    p<-NULL
    n<-NULL
    for(t in 1:5){
      randomSamp <- sample(sample_names(Sharing_Property), 1, replace=FALSE)
      randomSampTaxa <- prune_samples(randomSamp, Sharing_Property)
      randomSampTaxa <- prune_taxa(taxa_sums(randomSampTaxa)>0,randomSampTaxa)
      n <- c(n,length(intersect(taxa_names(SampTaxa),taxa_names(randomSampTaxa)))) # the number of ASVs shared with the focal sample
      p <- c(p,length(intersect(taxa_names(SampTaxa),taxa_names(randomSampTaxa)))/length(taxa_names(SampTaxa))) # the proportion of ASVs shared with the focal sample
    }
    NumShared_Property <- mean(n) # take the mean of the number of ASVs shared with the focal sample by samples of the same stage
    PropShared_Property <- mean(p) # take the mean of the proportion of ASVs shared with the focal sample by samples of the same stage
    assign(paste("NumShared",Property,sep="_"),NumShared_Property)
    assign(paste("PropShared",Property,sep="_"),PropShared_Property)
  }
  
  ### Randomly draw from the entire sample set and find the mean number and proportion of ASVs the draws share with the focal sample
  p<-NULL
  n<-NULL
  for(t in 1:5){
    random<- sample(sample_names(physeq_forMembership), 1, replace=FALSE)
    randomSet <- prune_samples(random, physeq_forMembership)
    randomSet <- prune_taxa(taxa_sums(randomSet)>0,randomSet)
    p <- c(p,length(intersect(taxa_names(SampTaxa),taxa_names(randomSet)))/length(taxa_names(SampTaxa)))
    n <- c(n,length(intersect(taxa_names(SampTaxa),taxa_names(randomSet))))
  }
  PropShared_Random <- mean(p)
  NumShared_Random <- mean(n)
  
  ### If samples from the same individual plant are present, find the mean number proportion of ASVs that draws share with those samples
  p<-NULL
  n<-NULL
  SampPlantID <- as.character(sample_data(physeq_forMembership)[SampID,"PlantID"][[1]]) # get the plant ID of the focal sample
  # prune the data to just samples sharing the plant ID
  samps_to_keep = as.character(get_variable(physeq_forMembership, "PlantID")) == SampPlantID & as.character(rownames(sample_data(physeq_forMembership))) != SampID
  if(table(samps_to_keep)["FALSE"] == length(sample_names(physeq_forMembership))){ # if there are no samples with matching plant IDs:
    p <- NA # set values to NA
    n <- NA
  }else{ # if there is some number of matching samples:
    Sharing_PlantID = prune_samples(samps_to_keep, physeq_forMembership) # get only the samples with the matching plant IDs
    for(t in sample_names(Sharing_PlantID)){ # for each of these samples:
      nextSampTaxa <- prune_samples(t, Sharing_PlantID) # prune data to the sample drawn
      nextSampTaxa <- prune_taxa(taxa_sums(randomSampTaxa)>0,randomSampTaxa) # prune ASVs to those present in the selected sample with matching plant ID
      n <- c(n,length(intersect(taxa_names(SampTaxa),taxa_names(nextSampTaxa)))) # the number of ASVs shared with the focal sample
      p <- c(p,length(intersect(taxa_names(SampTaxa),taxa_names(nextSampTaxa)))/length(taxa_names(SampTaxa))) # the proportion of ASVs shared with the focal sample
    }
  }
  PropShared_PlantID <- mean(p)
  NumShared_PlantID <- mean(n)
  
  ### Store results
  NextRow <- NULL
  NextRow <- c(SampID,PropShared_PlantID,NumShared_PlantID,PropShared_Stage,NumShared_Stage,PropShared_Site,NumShared_Site,PropShared_PlantPart,NumShared_PlantPart,PropShared_Year,NumShared_Year,PropShared_Random,NumShared_Random)
  Membership[nrow(Membership)+1,] <- NextRow
}

# save results
write.table(Membership,file="~/Documents/KBMP2020_Microbes/Outputs/Membership", sep="\t", row.names=TRUE,quote=FALSE)

##########################################################################################################################################################

# Load data for plotting
Membership <- read.table("~/Documents/KBMP2020_Microbes/Outputs/Membership", sep="\t", header=TRUE)

# Make columns numeric for plotting, select proportions
ForHistogram <- as.data.frame(sapply(Membership[,c(2:13)], as.numeric))
ForHistogram <- ForHistogram %>% select(PropShared_Random,PropShared_Part, PropShared_Stage, PropShared_Site, PropShared_Year,PropShared_PlantID)

# Frequency plots for the proportion of ASV member overlap in samples sharing properties vs. random
freq <- ggplot(melt(ForHistogram), aes(value, linetype = variable,stat(count),color = variable,stat(count))) + geom_freqpoly(binwidth = 0.1, size=1)+
  coord_cartesian(xlim=c(-0.1,1))+
  theme_bw()+
  theme(strip.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_linetype_manual(values=c("solid","solid","twodash","dotdash","dotted","dotted"),labels=c("Random Sample","Same Tissue","Same Stage","Same Site","Same Year","Same Plant"),name='Samples Compared')+
  scale_color_manual(values=c("grey","black","black","black","black","grey"),labels=c("Random Sample","Same Tissue","Same Stage","Same Site","Same Year", "Same Plant"),name='Samples Compared')+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+ 
  font("xlab", size = 16, color = "black",face = "bold")+ 
  labs(x="\nProportion of ASVs Shared",y="Frequency (in 1195 plant samples)\n")+ 
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),legend.position=c(0.65,0.8))+
  theme(legend.background = element_rect(color="black", size=.5))+ guides(linetype=guide_legend(keywidth = 5, keyheight = 1))

# Statistical tests of frequency distributions, null distribution is the proportion of members shared by randomly selected plants
wilcox.test(ForHistogram$PropShared_Random,ForHistogram$PropShared_Stage, alternative = "less")
wilcox.test(ForHistogram$PropShared_Random,ForHistogram$PropShared_Site, alternative = "less")
wilcox.test(ForHistogram$PropShared_Random,ForHistogram$PropShared_Part, alternative = "less")
wilcox.test(ForHistogram$PropShared_Random,ForHistogram$PropShared_Year, alternative = "less")
wilcox.test(ForHistogram$PropShared_Random,ForHistogram$PropShared_PlantID, alternative = "less")

median(ForHistogram$PropShared_Random)
median(ForHistogram$PropShared_Stage)
median(ForHistogram$PropShared_Site)
median(ForHistogram$PropShared_Year)
median(ForHistogram$PropShared_Part)
median(ForHistogram$PropShared_PlantID[!is.na(ForHistogram$PropShared_PlantID)])

ggsave("~/Documents/KBMP2020_Microbes/Figures/Membership.tiff", plot = freq, device = NULL, path = NULL,
       scale = 1, width = 5, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)

ggsave("~/Documents/KBMP2020_Microbes/Figures/Membership.png", plot = freq, device = NULL, path = NULL,
       scale = 1, width = 5, height = 5, units = c("in"),
       dpi = 600, limitsize = TRUE)
