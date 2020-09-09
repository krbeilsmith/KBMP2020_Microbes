# KBMP2020_VarianceByTaxonomy
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Summer 2020

# These are the commands used to find the variables associated with variation in the plant bacterial community composition at different levels of taxonomic
# grouping from fine (amplicon sequence variants) to coarse (phylum).

# The next goal is to find the levels of taxonomic groupings at which the model terms are significantly associated with the community composition.

####################################################################################################################################################

# Make dataframe for results
TaxResultsDF <- NULL
TaxResultsDF <- data.frame(matrix(ncol = 4, nrow = 0),stringsAsFactors = FALSE)
TaxResultsDF <- setNames(data.frame(TaxResultsDF), c("Rank","R2","Pr(>F)"))

# At ASV level:
# Calculate the matrix of Raup-Crick dissimilarities
# Run ANOVA with the variables for tissue type or stage and the new Raup-Crick distances
# Store results

Terms = c("PlantPart","Stage")
raup_physeq <- phyloseq::distance(R_physeq, method = "raup") # Raup-Crick dissimilarity
for(a in Terms){
  perma_raup_physeq <- adonis(as.formula(paste('raup_physeq',a,sep='~')),data.frame(sample_data(mphyseq)),permutations=999) # PERMANOVA
  PERMAResults <- data.frame(perma_raup_physeq$aov.tab[,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")])
  print(PERMAResults)
  TaxResultsDF[nrow(TaxResultsDF)+1,] <- c("ASV", a, PERMAResults[a,"R2"], PERMAResults[a,"Pr..F."])
}

# For each taxonomic level above ASV:
# Remove the unassigned taxa.
# Group the ASV counts at the desired level
# Calculate the matrix of Raup-Crick dissimilarities
# Run ANOVA with the variables for tissue type or stage and the new Raup-Crick distances
# Store results

RankList = c("Rank6","Rank5", "Rank4", "Rank3", "Rank2")
for(b in RankList){
  mphyseq <- prune_taxa(row.names(tax_table(R_physeq))[!(grepl("aggregate",tax_table(R_physeq)[,b]))],R_physeq) # remove unknowns
  mphyseq <-prune_samples(sample_sums(mphyseq)>0,mphyseq) # remove any samples with 0 counts
  mphyseq <- tax_glom(mphyseq, taxrank=b) # group ASVs
  raup_mphyseq <- phyloseq::distance(mphyseq, method = "raup") # Raup-Crick dissimilarity
  for(a in Terms){
    perma_raup_mphyseq <- adonis(as.formula(paste('raup_mphyseq',a,sep='~')),data.frame(sample_data(mphyseq)),permutations=999) # PERMANOVA
    PERMAResults <- data.frame(perma_raup_mphyseq$aov.tab[,c("Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")])
    print(PERMAResults)
    TaxResultsDF[nrow(TaxResultsDF)+1,] <- c(b, a, PERMAResults[a,"R2"],PERMAResults[a,"Pr..F."])
  }
}

# Display the table
z = ztable(TaxResultsDF,zebra=1)
z

# Save results
write.table(TaxResultsDF,file="~/Documents/KBMP2020_Microbes/PERMANOVAs/VarianceByTaxonomy", sep="\t", row.names=TRUE,quote=FALSE)
