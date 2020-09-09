
MyData <- PlottingFrame %>% filter(ASV %in% StageIndicators) %>%
  group_by(ASV,Site,Year) %>% mutate(SYSum = sum(Prevalence)) %>% filter(all(SYSum > 0.00)) %>% ungroup %>% 
  group_by(ASV,Site,Year,Tissue) %>% mutate(MP = max(Prevalence), MRA=max(MeanRA))

DotPlotData <- PlottingFrame %>% filter(!(ASV %in% unique(as.character(MyData$ASV)))) %>%
  group_by(ASV,Site,Year) %>% mutate(SYSum = sum(Prevalence)) %>% filter(all(SYSum > 0.00)) %>% ungroup %>% 
  group_by(ASV,Site,Year,Tissue) %>% mutate(MP = max(Prevalence), MRA=max(MeanRA))

DotPlotData$Genus <- ifelse(DotPlotData$Genus=="aggregate unclassified genus","unclassified",as.character(DotPlotData$Genus))
DotPlotData$Genus <- ifelse(DotPlotData$Genus=="Burkholderia-Paraburkholderia","Paraburkholderia",as.character(DotPlotData$Genus))

Site.labs<-c("Site: ME", "Site: WW")
names(Site.labs)<-levels(DotPlotData$Site)

Year.labs<-c("Year: 1", "Year: 2")
names(Year.labs)<-levels(DotPlotData$Year)

Plot_Filter_C <- NULL
Plot_Filter_C <- ggplot(data= MyData, aes(x=MP, y=MRA)) + 
  geom_point(data = DotPlotData,aes(x=MP, y=MRA),color="grey",alpha=0.2)+
  geom_point(aes(colour=Tissue))+
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"),labels=c("Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+
  theme_bw()+
  theme(strip.background.y=element_blank(),strip.text.x = element_text(size = 20, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 16, color = "black", face = "bold")+
  facet_grid(rows=vars(Site),cols=vars(Year),scales="free_x",space="free_x",switch="y",labeller=labeller(Site=Site.labs, Year=Year.labs))+
  labs(y="Maximum Relative Abundance \n",x="\n Maximum Prevalence")+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),
        legend.position="right")+
  theme(legend.background = element_rect(color="black", size=.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12),axis.text.y = element_text(size=12))+
  theme(strip.text.y = element_text(size=16, face="bold", vjust = 2),strip.placement="outside")+
  guides(size = guide_legend(title="Peak Prevalence",order=1,nrow=1),
         color=FALSE,alpha=guide_legend(title="Peak Mean Relative Abundance",override.aes = list(size=5),nrow=1))+
  theme(strip.text.x = element_text(size = 16),strip.background.x= element_blank())+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(strip.placement = "outside")+ theme(panel.spacing = unit(0, "lines"))+ theme(panel.spacing.y = unit(1, "lines"))+
  theme(legend.text=element_text(size=14,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),
        legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5))

Plot_Filter_D <- NULL
Plot_Filter_D <- ggplot(data= MyData, aes(x=MP, y=stat(ndensity))) + geom_density(alpha=0.5)+
  geom_density(data = DotPlotData,aes(x=MP, y=stat(ndensity)),alpha=0.2,fill="black")+
  scale_fill_manual(values=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"),labels=c("Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+
  theme_bw()+
  theme(strip.background.y=element_blank(),strip.text.x = element_text(size = 20, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 16, color = "black", face = "bold")+
  #facet_grid(rows=vars(forcats::fct_rev(Tissue)),scales="free_x",space="free_x",switch="y",labeller=labeller(Site=Site.labs, Year=Year.labs))+
  labs(y="Frequency \n",x="\n Maximum Prevalence")+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),
        legend.position="right")+
  theme(legend.background = element_rect(color="black", size=.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12),axis.text.y = element_text(size=12))+
  theme(strip.text.y = element_text(size=16, face="bold", vjust = 2),strip.placement="outside")+
  guides(size = guide_legend(title="Peak Prevalence",order=1,nrow=1),
         color=FALSE,alpha=guide_legend(title="Peak Mean Relative Abundance",override.aes = list(size=5),nrow=1))+
  theme(strip.text.x = element_text(size = 16),strip.background.x= element_blank())+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(strip.placement = "outside")+ theme(panel.spacing = unit(0, "lines"))+ theme(panel.spacing.y = unit(1, "lines"))+
  theme(legend.text=element_text(size=14,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),
        legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5))

Plot_Filter_E <- NULL
Plot_Filter_E <- ggplot(data= MyData, aes(x=MRA, y=stat(ndensity))) + geom_density(alpha=0.5)+
  geom_density(data = DotPlotData,aes(x=MRA, y=stat(ndensity)),alpha=0.2,fill="black")+
  scale_fill_manual(values=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"),labels=c("Roots","Rosettes","Stems","Cauline Leaves","Flowers","Siliques"),name='Sample Type')+
  theme_bw()+
  theme(strip.background.y=element_blank(),strip.text.x = element_text(size = 20, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 16, color = "black", face = "bold")+
  #facet_grid(rows=vars(forcats::fct_rev(Tissue)),scales="free_x",space="free_x",switch="y",labeller=labeller(Site=Site.labs, Year=Year.labs))+
  labs(y="Frequency \n",x="\n Maximum Realtive Abundace")+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),legend.title=element_text(size=12,color = "black",face = "bold"),
        legend.position="right")+
  theme(legend.background = element_rect(color="black", size=.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12),axis.text.y = element_text(size=12))+
  theme(strip.text.y = element_text(size=16, face="bold", vjust = 2),strip.placement="outside")+
  guides(size = guide_legend(title="Peak Prevalence",order=1,nrow=1),
         color=FALSE,alpha=guide_legend(title="Peak Mean Relative Abundance",override.aes = list(size=5),nrow=1))+
  theme(strip.text.x = element_text(size = 16),strip.background.x= element_blank())+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(strip.placement = "outside")+ theme(panel.spacing = unit(0, "lines"))+ theme(panel.spacing.y = unit(1, "lines"))+
  theme(legend.text=element_text(size=14,color = "black",face = "bold"),legend.title=element_text(size=14,color = "black",face = "bold"),
        legend.position="top")+
  theme(legend.background = element_rect(color="black", size=.5))






  
