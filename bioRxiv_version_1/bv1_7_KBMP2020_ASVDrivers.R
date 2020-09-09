# KBMP2020_ASVDrivers
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Autumn 2019

# The commands below are used to find the ASVs that are driving the associations between study variables and community composition in the PERMANOVA.

# The approach of taking the coefficients from the PERMANOVA broadly follows the tutorial here: 
# https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/multivariate-comparisons-of-microbial-community-composition.html

####################################################################################################################################################
# Collect coefficients
####################################################################################################################################################

# From PERMANOVA output, collect all the regression coefficients for the variables plant tissue, stage, year, and site.

# rownames(coefficients(full_permanova_results))

PlantTissueEffects <- data.frame(t(coefficients(full_permanova_results)[c(rownames(coefficients(full_permanova_results)[grep(":PlantPart",rownames(coefficients(full_permanova_results))),])),]))

PlantStageEffects <- data.frame(t(coefficients(full_permanova_results)[c(rownames(coefficients(full_permanova_results)[grep(":Stage\\d",rownames(coefficients(full_permanova_results))),])),]))

PlantSiteEffects <- data.frame(coefficients(full_permanova_results)["Site1",])

PlantYearEffects <- data.frame(coefficients(full_permanova_results)["Year1",])

# Histograms of coefficients for each variable
CoeffList <- list(PlantTissueEffects, PlantStageEffects, PlantSiteEffects, PlantYearEffects)
for(coeffs in 1:length(CoeffList)){
  plot <- ggplot(melt(CoeffList[coeffs]), aes(value, color = variable)) + geom_freqpoly(binwidth = 10, size=1)+coord_cartesian(xlim=c(-500,500),ylim=c(0,100))+
    theme_bw()+
    theme(strip.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    font("ylab", size = 16, color = "black",face = "bold")+ # text format
    font("xy.text", size = 12, color = "black", face = "bold")+ # text format
    font("xlab", size = 16, color = "black",face = "bold")+ # text format
    labs(x="\nCoefficients",y="Count of ASVs\n")+ # labels
    rremove("legend")
  assign(paste("CoeffsPlot",coeffs,sep="_"),plot)
}

CoeffsPlot <- grid.arrange(
  CoeffsPlot_1 + ggtitle("Tissue") + theme(plot.title = element_text(face="bold"),axis.title.y = element_blank(),axis.title.x = element_blank(),plot.margin=margin(10,10,10,10)),
  CoeffsPlot_2 + ggtitle("Stage") + theme(plot.title = element_text(face="bold"),axis.title.y = element_blank(),axis.title.x = element_blank(),plot.margin=margin(10,10,10,10)),
  CoeffsPlot_3 + ggtitle("Site") + theme(plot.title = element_text(face="bold"),axis.title.y = element_blank(),axis.title.x = element_blank(),plot.margin=margin(10,10,10,10)),
  CoeffsPlot_4 + ggtitle("Year") + theme(plot.title = element_text(face="bold"),axis.title.y = element_blank(),axis.title.x = element_blank(),plot.margin=margin(10,10,10,10)),
  nrow=2,
  left = textGrob("ASV count\n",rot = 90, gp = gpar(fontface = "bold", cex = 1.3)),
  bottom = textGrob("\nPERMANOVA coefficients", gp = gpar(fontface = "bold", cex = 1.3)),
  vp=viewport(width=0.95, height=0.95)
)

ggsave("~/Documents/KBMP2020_Microbes/Figures/Frequency_PERMANOVA_Coefficients.tiff", plot = CoeffsPlot, device = NULL, path = NULL,
       scale = 1, width = 8, height = 8, units = c("in"),
       dpi = 600, limitsize = TRUE)

####################################################################################################################################################
# Filter coefficients
####################################################################################################################################################

# For each set of coefficients, select the values that are not NAs. 
# For tissue effects (the largest), find the top 5% of coefficient absolute value sums. These sums of absolute coefficient values will represent the
# relative "influence" of ASVs on the observed association with tissue type.
# Take the minimum value of this distribution and use it as the threshold for filtering the coefficients of other effects (stage, site, year).

PlantTissueEffects <- PlantTissueEffects %>% select_if(~ !any(is.na(.)))
PlantTissueEffects$EffectSums<-rowSums(abs(PlantTissueEffects))
TopFivePercent <- PlantTissueEffects %>% mutate(ASV=rownames(PlantTissueEffects)) %>% top_frac(0.05,EffectSums)
PlantTissueEffects <- PlantTissueEffects[rownames(PlantTissueEffects) %in% TopFivePercent$ASV,]

PlantStageEffects <- PlantStageEffects %>% select_if(~ !any(is.na(.)))
PlantStageEffects$EffectSums<-rowSums(abs(PlantStageEffects))
PlantStageEffects<-PlantStageEffects[PlantStageEffects$EffectSums > min(PlantTissueEffects$EffectSums),]

PlantSiteEffects <- PlantSiteEffects %>% select_if(~ !any(is.na(.)))
PlantSiteEffects$EffectSums<-rowSums(abs(PlantSiteEffects))
PlantSiteEffects<-PlantSiteEffects[PlantSiteEffects$EffectSums > min(PlantTissueEffects$EffectSums),]

PlantYearEffects <- PlantYearEffects %>% select_if(~ !any(is.na(.)))
PlantYearEffects$EffectSums<-rowSums(abs(PlantYearEffects))
PlantYearEffects<-PlantYearEffects[PlantYearEffects$EffectSums > min(PlantTissueEffects$EffectSums),]

####################################################################################################################################################
# Compare the ASVs with coefficients passing filter for each variable.
####################################################################################################################################################

# Get the names of the ASVs with coefficients that are above the threshold (mean + 1sd) for each variable.
TissueVarASVs <- rownames(PlantTissueEffects) # 413
StageVarASVs <- rownames(PlantStageEffects) # 168
SiteVarASVs <- rownames(PlantSiteEffects) # 5
YearVarASVs <- rownames(PlantYearEffects) # 6

# Make a venn diagram and use the numbers to produce an Euler diagram
DrivingASVs <- list(TissueVarASVs,StageVarASVs,SiteVarASVs,YearVarASVs)
venn.diagram(DrivingASVs,filename="~/Documents/KBMP2020_Microbes/Figures/Drivers",category.names = c("Tissue" , "Stage " , "Site", "Year") ,output="TRUE")
unique(unlist(DrivingASVs))

# Plot
palette <-c("#f9d50b","#86aa39","#C0C022","#C0C022")
Euler<- c("Stage" = 2, "Tissue" = 247, "Site" = 0, "Year" = 0, 
          "Tissue&Stage"= 158, "Tissue&Year" = 0, "Tissue&Site" = 0, "Stage&Site" = 0, "Stage&Year" = 0, "Site&Year" = 0,
          "Tissue&Stage&Year" = 3, "Tissue&Stage&Site" = 2, "Stage&Site&Year" = 0, "Year&Site&Tissue" = 0,
          "Tissue&Stage&Site&Year" = 3)
EulerPlot <- plot(euler(Euler, shape = "ellipse"),fontsize=20,quantities = list(fontsize = 14,font=2),fills=palette,edges=TRUE,lty = c(1,1,2,3))

ggsave("~/Documents/KBMP2020_Microbes/Figures/Venn_PERMANOVA_Coefficients.tiff", plot = EulerPlot, device = NULL, path = NULL,
       scale = 1.2, width = 5, height = 4, units = c("in"),
       dpi = 600, limitsize = TRUE)

####################################################################################################################################################
# Examine whether the top ASVs driving the tissue effect (largest summed PERMANOVA coefficients) also drive stage, site, and year effects
# Examine the taxonomy of the top ASVs driving the tissue effect
####################################################################################################################################################

# Code ASVs with tissue effects based on whether they are also among those driving stage, site, and year effects.
PlantTissueEffects$StageVar<-ifelse(rownames(PlantTissueEffects)%in% StageVarASVs,"Y","N")
PlantTissueEffects$SiteVar<-ifelse(rownames(PlantTissueEffects)%in% SiteVarASVs,"Y","N")
PlantTissueEffects$YearVar<-ifelse(rownames(PlantTissueEffects)%in% YearVarASVs,"Y","N")

# Proportion of ASVs with top tissue effects that also drive stage effects: 40%
length(rownames(PlantTissueEffects[PlantTissueEffects$StageVar=="Y",]))/length(rownames(PlantTissueEffects))

# Add the summed coefficients for stage if the ASV is one of the 40% that drive stage effects as well as tissue effects
for(r in rownames(PlantTissueEffects)){
  PlantTissueEffects[r,"StageEffectSize"] <- ifelse(r%in%rownames(PlantStageEffects),PlantStageEffects[r,"EffectSums"],0)
}
PlantTissueEffects$StageEffectSize <- as.numeric(as.character(PlantTissueEffects$StageEffectSize))

# Add taxonomy information
for(r in rownames(PlantTissueEffects)){
  PlantTissueEffects[r,"Class"] <- as.character(tax_table(physeq)[r,"Rank3"])
  PlantTissueEffects[r,"Family"] <- as.character(tax_table(physeq)[r,"Rank5"])
  PlantTissueEffects[r,"Genus"] <- as.character(tax_table(physeq)[r,"Rank6"])
}

# Plot the top 50 tissue-driving ASVs and color the bars based on the stage effect. Color and pattern the lines based on whether the ASV
# drives site and year effects.

PlantTissueEffects$Genus <- ifelse(PlantTissueEffects$Genus=="aggregate unclassified genus","unclassified",PlantTissueEffects$Genus)
barplot <- ggbarplot(PlantTissueEffects %>% mutate(ASV = rownames(PlantTissueEffects)) %>% arrange(EffectSums) %>% top_n(50,EffectSums), 
                     "ASV", "EffectSums", position = position_dodge(1.5),orientation = "horiz",fill ="StageEffectSize",color="YearVar",linetype="SiteVar",size=1)+
  theme_bw()+
  theme(strip.background=element_blank(),strip.text.x = element_text(size = 16, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(x="Amplicon sequence variants\n",y="\nInfluence on tissue effect\n(sum of PERMANOVA coefficient absolute values)")+
  scale_color_manual(values=c("black","red"),name="Year-variable",breaks=c("N","Y"))+
  #scale_fill_manual(values=c("grey40","grey80"),name="Stage-variable",breaks=c("N","Y"))+
  scale_fill_gradient(low="black", high="grey90", na.value = "black",name="Influence on\nstage effect")+
  scale_linetype_manual(values=c(1,6),name="Site-variable",breaks=c("N","Y"))+
  theme(legend.position = "right")+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"))+
  guides(linetype = guide_legend(nrow=2,override.aes = list(fill = NA, col = "black")),
         color=guide_legend(nrow=2,override.aes = list(fill = NA)))+
  scale_x_discrete(expand=c(0.025,1),label=function(x) abbreviate(x, minlength=8))+
  facet_grid(rows=vars(Genus),scales="free_y",space="free_y",switch="y")+
  theme(strip.text.y = element_text(size = 12,angle=180,face="bold"))+
  theme(strip.placement = "outside") + theme(panel.spacing = unit(0.1, "lines"))

ggsave("~/Documents/KBMP2020_Microbes/Figures/Bar_PERMANOVA_Coefficients.tiff", plot = barplot, device = NULL, path = NULL,
       scale = 1.2, width = 6, height = 9, units = c("in"),
       dpi = 600, limitsize = TRUE)

####################################################################################################################################################

# Plot the top 50 tissue-driving ASVs and color the bars based on the class of the ASV.

barplot <- ggbarplot(PlantTissueEffects %>% mutate(ASV = rownames(PlantTissueEffects)) %>% arrange(EffectSums) %>% top_n(50,EffectSums), 
                     "ASV", "EffectSums", position = position_dodge(1.5),orientation = "horiz",fill ="Class",color=NA,size=1)+
  theme_bw()+
  theme(strip.background=element_blank(),strip.text.x = element_text(size = 16, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(x="Amplicon sequence variants\n",y="\nInfluence on tissue effect\n(sum of PERMANOVA coefficient absolute values)")+
  scale_fill_manual(values=c("#20A486FF","#EEDD82FF","#000004FF","#F56B5CFF","#75D054FF","#FD9B6BFF","#D8456CFF","#9C2E7FFF","#D5E21AFF"),name="Class")+
  theme(legend.position = "right")+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"))+
  scale_x_discrete(expand=c(0.025,1),label=function(x) abbreviate(x, minlength=8))

ggsave("~/Documents/KBMP2020_Microbes/Figures/Taxonomy_Bar_PERMANOVA_Coefficients.tiff", plot = barplot, device = NULL, path = NULL,
       scale = 1.2, width = 5, height = 8, units = c("in"),
       dpi = 600, limitsize = TRUE)
