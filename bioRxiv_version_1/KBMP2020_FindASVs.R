# KBMP2020_FindASVs
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Summer 2020

# Make a copy of the phyloseq object and encode sample data variables as factors.
physeq_factored <- physeq

sample_data(physeq_factored)$Year <- factor(sample_data(physeq_factored)$Year, levels=c("1", "2"))
sample_data(physeq_factored)$Stage <- factor(sample_data(physeq_factored)$Stage, levels=c("Soil", "TwoLeaf","FourLeaf","SixLeaf","EightLeaf","Flowering","Senescent"))
sample_data(physeq_factored)$PlantPart <- factor(sample_data(physeq_factored)$PlantPart, levels=c("Soil", "Roots","RosLeaves","Stems","CauLeaves","Flowers","Siliques"))
sample_data(physeq_factored)$Ecotype <- factor(sample_data(physeq_factored)$Ecotype)
sample_data(physeq_factored)$Site <- factor(sample_data(physeq_factored)$Site, levels=c("ME", "WW"))
sample_data(physeq_factored)$SamplePlate <- factor(sample_data(physeq_factored)$SamplePlate)
sample_data(physeq_factored)$MiSeqRun <- factor(as.character(as.numeric(sample_data(physeq_factored)$MiSeqRun)))
sample_data(physeq_factored)$PlantID <- factor(as.character(sample_data(physeq_factored)$PlantID))
sample_data(physeq_factored)$Soil <- factor(sample_data(physeq_factored)$Soil,levels=c("Y","N"))

# Focus on plant samples
physeq_factored_plants <- subset_samples(physeq_factored,PlantPart!="Soil")
physeq_factored_plants <- prune_taxa(taxa_sums(physeq_factored_plants)>0,physeq_factored_plants)

# Rarefy samples to even depth (will cause samples to drop from set)
R_physeq <- rarefy_even_depth(physeq_factored_plants,1000)

# Roots and Rosettes comparison across development
R_physeq_earlyveg <- subset_samples(R_physeq,!(PlantPart %in% c("Soil")) & Stage%in%c("TwoLeaf","FourLeaf"))
R_physeq_earlyveg <- prune_taxa(taxa_sums(R_physeq_earlyveg)>0,R_physeq_earlyveg)

R_physeq_lateveg <- subset_samples(R_physeq,!(PlantPart %in% c("Soil")) & Stage%in%c("SixLeaf","EightLeaf"))
R_physeq_lateveg <- prune_taxa(taxa_sums(R_physeq_lateveg)>0,R_physeq_lateveg)

R_physeq_flowering <- subset_samples(R_physeq,!(PlantPart %in% c("Soil","CauLeaves","Flowers","Siliques","Stems")) & Stage%in%c("Flowering"))
R_physeq_flowering <- prune_taxa(taxa_sums(R_physeq_flowering)>0,R_physeq_flowering)

# Roots and phyllosphere in senescent plants
R_physeq_senescent <- subset_samples(R_physeq,!(PlantPart %in% c("Soil","RosLeaves","CauLeaves","Flowers")) & Stage%in%c("Senescent"))
R_physeq_senescent <- prune_taxa(taxa_sums(R_physeq_senescent)>0,R_physeq_senescent)

options(max.print = 3000)

set_list <- c("R_physeq_earlyveg","R_physeq_lateveg","R_physeq_flowering","R_physeq_senescent")
for(s in set_list){
  sampset <- get(s)
  #print(sampset)
  sample_type <- data.frame(sample_data(sampset))
  sample_type$Sphere <- ifelse(sample_type$PlantPart=="Roots","Rhizo","Phyllo")
  group_list <- as.character(sample_type[order(sample_type$Sphere),]$Sphere)
  group_order <- as.character(rownames(sample_type[order(sample_type$Sphere),]))
  species_matrix <- t(otu_table(sampset))
  sample_type_ordered <- species_matrix[match(group_order,rownames(species_matrix)),]
  indval = multipatt(sample_type_ordered, group_list, func="IndVal.g", control = how(nperm=999))
  # summary(indval,indvalcomp=TRUE)
  sink(paste("~/Documents/KBMP2020_Microbes/ASVList",s,sep="_"))
  print(summary(indval,indvalcomp=TRUE))
  sink() 
}

# Phyllosphere comparison in flowering
R_physeq_flowerphyllo <- subset_samples(R_physeq,!(PlantPart %in% c("Soil","Roots")) & Stage%in%c("Flowering"))
R_physeq_flowerphyllo <- prune_taxa(taxa_sums(R_physeq_flowerphyllo)>0,R_physeq_flowerphyllo)

set_list <- c("R_physeq_flowerphyllo")
for(s in set_list){
  sampset <- get(s)
  #print(sampset)
  sample_type <- data.frame(sample_data(sampset))
  sample_type$Sphere <- ifelse(sample_type$PlantPart=="RosLeaves","LowerPhyllo","UpperPhyllo")
  group_list <- as.character(sample_type[order(sample_type$Sphere),]$Sphere)
  group_order <- as.character(rownames(sample_type[order(sample_type$Sphere),]))
  species_matrix <- t(otu_table(sampset))
  sample_type_ordered <- species_matrix[match(group_order,rownames(species_matrix)),]
  indval = multipatt(sample_type_ordered, group_list, func="IndVal.g", control = how(nperm=999))
  # summary(indval,indvalcomp=TRUE)
  sink(paste("~/Documents/KBMP2020_Microbes/ASVList",s,sep="_"))
  print(summary(indval,indvalcomp=TRUE))
  sink() 
}

summary(indval)

indval$str["ccdbd2c5ce2050852983933be7e2726b","EarlyVeg"]

str(indval)

attr(indval$str,2)
attributes(indval$sign)

# Fix heders, save as CSVs in the Outputs folder

earlyveg_phyllo_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_earlyveg_phyllo.csv", header=TRUE, sep="\t",row.names=1)
earlyveg_rhizo_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_earlyveg_rhizo.csv", header=TRUE, sep="\t",row.names=1)
lateveg_phyllo_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_lateveg_phyllo.csv", header=TRUE, sep="\t",row.names=1)
lateveg_rhizo_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_lateveg_rhizo.csv", header=TRUE, sep="\t",row.names=1)
flowering_phyllo_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_flowering_phyllo.csv", header=TRUE, sep="\t",row.names=1)
flowering_rhizo_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_flowering_rhizo.csv", header=TRUE, sep="\t",row.names=1)
senescent_phyllo_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_senescent_phyllo.csv", header=TRUE, sep="\t",row.names=1)
senescent_rhizo_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_senescent_rhizo.csv", header=TRUE, sep="\t",row.names=1)

pthresh = 0.01
Bthresh = 0.4

PhylloRhizoIndicators <- unique(c(rownames(earlyveg_phyllo_ASVs)[earlyveg_phyllo_ASVs[,"p.value"] <= pthresh & earlyveg_phyllo_ASVs[,"B"] >= Bthresh],
  rownames(earlyveg_rhizo_ASVs)[earlyveg_rhizo_ASVs[,"p.value"] <= pthresh & earlyveg_rhizo_ASVs[,"B"] >= Bthresh],
  rownames(lateveg_phyllo_ASVs)[lateveg_phyllo_ASVs[,"p.value"] <= pthresh & lateveg_phyllo_ASVs[,"B"] >= Bthresh],
  rownames(lateveg_rhizo_ASVs)[lateveg_rhizo_ASVs[,"p.value"] <= pthresh & lateveg_rhizo_ASVs[,"B"] >= Bthresh],
  rownames(flowering_phyllo_ASVs)[flowering_phyllo_ASVs[,"p.value"] <= pthresh & flowering_phyllo_ASVs[,"B"] >= Bthresh],
  rownames(flowering_rhizo_ASVs)[flowering_rhizo_ASVs[,"p.value"] <= pthresh & flowering_rhizo_ASVs[,"B"] >= Bthresh],
  rownames(senescent_phyllo_ASVs)[senescent_phyllo_ASVs[,"p.value"] <= pthresh & senescent_phyllo_ASVs[,"B"] >= Bthresh],
  rownames(senescent_rhizo_ASVs)[senescent_rhizo_ASVs[,"p.value"] <= pthresh & senescent_rhizo_ASVs[,"B"] >= Bthresh]))

PhylloIndicators <- unique(c(rownames(earlyveg_phyllo_ASVs)[earlyveg_phyllo_ASVs[,"p.value"] <= pthresh & earlyveg_phyllo_ASVs[,"B"] >= Bthresh],
                                  rownames(lateveg_phyllo_ASVs)[lateveg_phyllo_ASVs[,"p.value"] <= pthresh & lateveg_phyllo_ASVs[,"B"] >= Bthresh],
                                  rownames(flowering_phyllo_ASVs)[flowering_phyllo_ASVs[,"p.value"] <= pthresh & flowering_phyllo_ASVs[,"B"] >= Bthresh],
                                  rownames(senescent_phyllo_ASVs)[senescent_phyllo_ASVs[,"p.value"] <= pthresh & senescent_phyllo_ASVs[,"B"] >= Bthresh]))

RhizoIndicators <- unique(c(rownames(earlyveg_rhizo_ASVs)[earlyveg_rhizo_ASVs[,"p.value"] <= pthresh & earlyveg_rhizo_ASVs[,"B"] >= Bthresh],
                                  rownames(lateveg_rhizo_ASVs)[lateveg_rhizo_ASVs[,"p.value"] <= pthresh & lateveg_rhizo_ASVs[,"B"] >= Bthresh],
                                  rownames(flowering_rhizo_ASVs)[flowering_rhizo_ASVs[,"p.value"] <= pthresh & flowering_rhizo_ASVs[,"B"] >= Bthresh],
                                  rownames(senescent_rhizo_ASVs)[senescent_rhizo_ASVs[,"p.value"] <= pthresh & senescent_rhizo_ASVs[,"B"] >= Bthresh]))

roots_earlyveg_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_roots_earlyveg.csv", header=TRUE, sep="\t",row.names=1)
roots_lateveg_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_roots_lateveg.csv", header=TRUE, sep="\t",row.names=1)
roots_flowering_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_roots_flowering.csv", header=TRUE, sep="\t",row.names=1)
roots_senescent_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_roots_senescent.csv", header=TRUE, sep="\t",row.names=1)

roots_growing_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_roots_growing.csv", header=TRUE, sep="\t",row.names=1)

roots_prebolt_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_roots_prebolt.csv", header=TRUE, sep="\t",row.names=1)
roots_postbolt_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_roots_postbolt.csv", header=TRUE, sep="\t",row.names=1)

rosettes_flowering_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_rosettes_flowering.csv", header=TRUE, sep="\t",row.names=1)
rosettes_preflowering_ASVs <- read.table("~/Documents/KBMP2020_Microbes/Outputs/ASVList_R_physeq_rosettes_preflowering.csv", header=TRUE, sep="\t",row.names=1)


RootStageIndicators <- unique(c(rownames(roots_earlyveg_ASVs)[roots_earlyveg_ASVs[,"p.value"] <= pthresh & earlyveg_rhizo_ASVs[,"B"] >= Bthresh],
                                  rownames(roots_lateveg_ASVs)[roots_lateveg_ASVs[,"p.value"] <= pthresh & earlyveg_rhizo_ASVs[,"B"] >= Bthresh],
                                  rownames(roots_flowering_ASVs)[roots_flowering_ASVs[,"p.value"] <= pthresh & earlyveg_rhizo_ASVs[,"B"] >= Bthresh],
                                  rownames(roots_senescent_ASVs)[roots_senescent_ASVs[,"p.value"] <= pthresh & earlyveg_rhizo_ASVs[,"B"] >= Bthresh],
                                  rownames(roots_growing_ASVs)[roots_growing_ASVs[,"p.value"] <= pthresh & roots_growing_ASVs[,"B"] >= Bthresh],
                                  rownames(roots_prebolt_ASVs)[roots_prebolt_ASVs[,"p.value"] <= pthresh & roots_prebolt_ASVs[,"B"] >= Bthresh],
                                  rownames(roots_postbolt_ASVs)[roots_postbolt_ASVs[,"p.value"] <= pthresh & roots_postbolt_ASVs[,"B"] >= Bthresh]))

RosetteStageIndicators <- unique(c(rownames(rosettes_preflowering_ASVs)[rosettes_preflowering_ASVs[,"p.value"] <= pthresh & rosettes_preflowering_ASVs[,"B"] >= Bthresh], 
                                   rownames(rosettes_flowering_ASVs[rosettes_flowering_ASVs[,"p.value"] <= pthresh & rosettes_flowering_ASVs[,"B"] >= Bthresh])))

rownames(roots_prebolt_ASVs)

RootStageIndicators %in% RhizoIndicators

RosetteStageIndicators %in% PhylloIndicators

StageIndicators <- c(RootStageIndicators, RosetteStageIndicators)

StageIndicators %in% raup_AllTopASVs

# Roots and Rosettes comparison across development
R_physeq_roots <- subset_samples(R_physeq,PlantPart %in% c("Roots"))
R_physeq_roots <- prune_taxa(taxa_sums(R_physeq_roots)>0,R_physeq_roots)

# Roots and Rosettes comparison across development
R_physeq_rosettes <- subset_samples(R_physeq,PlantPart %in% c("RosLeaves"))
R_physeq_rosettes <- prune_taxa(taxa_sums(R_physeq_rosettes)>0,R_physeq_rosettes)

set_list <- c("R_physeq_roots","R_physeq_rosettes")
for(s in set_list){
  sampset <- get(s)
  #print(sampset)
  sample_type <- data.frame(sample_data(sampset))
  sample_type$DevStage <- ifelse(sample_type$Stage=="Senescent","Senescent",
                                 ifelse(sample_type$Stage=="Flowering","Flowering",
                                        ifelse(sample_type$Stage=="SixLeaf" | sample_type$Stage=="EightLeaf", "LateVeg", "EarlyVeg")))
  group_list <- as.character(sample_type[order(sample_type$DevStage),]$DevStage)
  group_order <- as.character(rownames(sample_type[order(sample_type$DevStage),]))
  species_matrix <- t(otu_table(sampset))
  sample_type_ordered <- species_matrix[match(group_order,rownames(species_matrix)),]
  indval = multipatt(sample_type_ordered, group_list, func="IndVal.g", control = how(nperm=999))
  # summary(indval,indvalcomp=TRUE)
  sink(paste("~/Documents/KBMP2020_Microbes/ASVList",s,sep="_"))
  print(summary(indval,indvalcomp=TRUE))
  sink() 
}
