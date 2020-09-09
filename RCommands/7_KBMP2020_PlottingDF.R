# KBMP2020_PlottingDF
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Summer 2020

####################################################################################################################################################

# This will collect the data for making time series plots and the dotplots to represent temporal and spatial patterns for the ASVs significantly 
# associated with tissue and stage.

# Although the PERMANOVAs, ordination, and association tests were performed with rarefied data to dodge differences in read depth, the dataframe
# below will select all the raw data on the ASVs of interest.

# Select plant samples and ASVs with significant associations with tissue type and stage.
physeq_plants <- subset_samples(physeq,PlantPart!="Soil")
physeq_plants <- prune_taxa(taxa_sums(physeq_plants)>0,physeq_plants)
physeq_select <- subset_taxa(physeq_plants,taxa_names(physeq_plants)%in%ASVIndicators)
physeq_select <- prune_samples(sample_sums(physeq_select)>0,physeq_select)

# Get all the unique combinations of harvest conditions (tissue, stage, site, year)
ConditionsFrame <- data.frame(unique(sample_data(physeq_select)[,c("Year","Site","Stage","PlantPart")]),row.names=c())
ConditionsFrame <- data.frame(lapply(ConditionsFrame, as.character), stringsAsFactors=FALSE)

# Make dataframe for storing results
PlottingFrame <- NULL
PlottingFrame <- data.frame(matrix(ncol = 13, nrow = 0),stringsAsFactors = FALSE)
PlottingFrame  <- setNames(data.frame(PlottingFrame), c("ASV","Tissue","Stage","Site","Year","Prevalence","MeanRA","sdRA",
                                                        "Genus","Family","Order","Class","Phylum"))

for(m in ASVIndicators){ # for each ASV (referred to as focal ASV below)
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
write.table(PlottingFrame,file="~/Documents/KBMP2020_Microbes/Outputs/PlottingFrame_Indicators", sep="\t", row.names=TRUE,quote=FALSE)

####################################################################################################################################################

# Now make the plotting dataframe for data aggregated at the genus level.

physeq_glom <- tax_glom(physeq_select,"Rank6")
physeq_glom <- subset_taxa(physeq_glom,Rank6!="aggregate unclassified genus")
taxa_names(physeq_glom) <- as.character(tax_table(physeq_glom)[,"Rank6"])

# Make dataframe for storing results
PlottingFrame_Genera <- NULL
PlottingFrame_Genera <- data.frame(matrix(ncol = 12, nrow = 0),stringsAsFactors = FALSE)
PlottingFrame_Genera  <- setNames(data.frame(PlottingFrame_Genera), c("Genus","Tissue","Stage","Site","Year","Prevalence","MeanRA","sdRA",
                                                                      "Family","Order","Class","Phylum"))

GenusList <- unique(as.character(tax_table(physeq_glom)[rownames(tax_table(physeq_glom)) %in% ASVIndicators,"Rank6"]))
GenusList <- GenusList[GenusList!="aggregate unclassified genus"]

for(m in GenusList){ # for each ASV (referred to as focal ASV below)
  print(m)
  for(row in rownames(ConditionsFrame)){ # for each set of harvest conditions
    pp <- ConditionsFrame[row,"PlantPart"] # get metadata for the harvest conditions
    st <- ConditionsFrame[row,"Stage"]
    si <- ConditionsFrame[row,"Site"]
    yr <- ConditionsFrame[row,"Year"]
    physeq_partition <-subset_samples(physeq_glom, PlantPart==pp & Stage==st & Site==si & Year==yr) # prune to data for the harvest conditions
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
    Family <- as.character(tax_table(physeq_glom)[m,"Rank5"])
    Order <- as.character(tax_table(physeq_glom)[m,"Rank4"])
    Class <- as.character(tax_table(physeq_glom)[m,"Rank3"])
    Phylum <- as.character(tax_table(physeq_glom)[m,"Rank2"])
    PlottingFrame_Genera[nrow(PlottingFrame_Genera)+1,] <- c(m,pp,st,si,yr,Prevalence,MeanRelAbund,SDRelAbund,Family,Order,Class,Phylum) # add values to dataframe
  }
}

# Store dataframe
write.table(PlottingFrame_Genera,file="~/Documents/KBMP2020_Microbes/Outputs/PlottingFrame_Genera", sep="\t", row.names=TRUE,quote=FALSE)

# PlottingFrame_Genera <- read.table("~/Documents/KBMP2020_Microbes/Outputs/PlottingFrame_Genera", header=TRUE, sep="\t",row.names=1)

# Set factor levels for plotting
PlottingFrame_Genera$Tissue <-factor(PlottingFrame_Genera$Tissue, levels=c("Roots","RosLeaves","Stems","CauLeaves","Flowers","Siliques"))
PlottingFrame_Genera$Stage <-factor(PlottingFrame_Genera$Stage, levels=c("TwoLeaf","FourLeaf","SixLeaf","EightLeaf","Flowering","Senescent"))
PlottingFrame_Genera$Site <-factor(PlottingFrame_Genera$Site, levels=c("ME","WW"))
PlottingFrame_Genera$Year <-factor(PlottingFrame_Genera$Year, levels=c("1","2"))
PlottingFrame_Genera$Genus <-factor(PlottingFrame_Genera$Genus)
PlottingFrame_Genera$Family <-factor(PlottingFrame_Genera$Family)
PlottingFrame_Genera$Order <-factor(PlottingFrame_Genera$Order)
PlottingFrame_Genera$Class <-factor(PlottingFrame_Genera$Class)
PlottingFrame_Genera$Phylum <-factor(PlottingFrame_Genera$Phylum)

# Make prevalence and relative abundance data numeric for plotting
PlottingFrame_Genera$Prevalence <- as.numeric(PlottingFrame_Genera$Prevalence)
PlottingFrame_Genera$MeanRA <- as.numeric(PlottingFrame_Genera$MeanRA)
PlottingFrame_Genera$sdRA <- as.numeric(PlottingFrame_Genera$sdRA)