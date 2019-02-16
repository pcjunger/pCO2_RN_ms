##Manucript: CO2 Supersaturation in Tropical Eutrophic Waters
#DOI number:10.1016/j.scitotenv.2019.01.273
###Script used to produce the figures and statistical analysis presented in the main text
### Authors: Pedro C. Junger et al.

#### Loading required libraries ####
library(tidyr)
library(dplyr)
library(ggplot2)
library(pastecs)
library(dunn.test)
library(gtools)
library(ggpubr)
library(visreg)
library(AICcmodavg)
library(car)
library(Rmisc)
library(lme4)
library(lmerTest)
library(ggrepel)
library(grid)
library(scales)
library(gridExtra)
library(lattice)
library(Rcmdr)

library(lavaan) #Package to run the SEM
library(vegan)
library(sem)
library(polycor)
library(MVN)
library(semPLS)
library(BiodiversityR)
BiodiversityRGUI()
require(fBasics)
library(FactoMineR)
library(royston)

#### Data input ####
spatial <- read.csv("spatial_dataset.csv",header = TRUE) #full spatial dataset from 102 ecosystems used in this study
seasonal <- read.csv("seasonal_dataset.csv", sep = ",",header = TRUE) #full temporal dataset used in this study

#for literature comparison
literature<-read.csv("literature_data.csv",header = TRUE)
biomes<-read.csv("biome_final.csv",header = TRUE, sep = " ") #literature biomes dataset

#### Data wrangling ####
##Ordering factors
#for spatial dataset
spatial$trophic_state <- factor(spatial$trophic_state, levels=c("Hypereutrophic","Eutrophic","Mesotrophic","Oligotrophic"))

#for seasonal dataset
seasonal$system <- factor(seasonal$system, levels = c("Gargalheiras", "Cruzeta","Extremoz","Bonfim"))

#for literature comparison
literature$lat_cat_this_study <- factor(literature$lat_cat_this_study, levels = c("Polar", "Subpolar","Temperate","Subtropical","Tropical", "This study"))

##Filtering the spatial dataset by excluding observations with Total Alkalinity < 1meq/L
#Obs.: See detailed reasons in the methods (section 2.4) of the manuscript
spatial_2<- spatial %>%
  filter(lakeID!=67) %>% #excluding outlier "Riacho da Cruz"
  filter(TA>1) #excluding superestimated pCO2 values from the full spatial dataset

seasonal_2 <- seasonal %>%
  filter(TA>1) #excluding superestimated pCO2 values from the full seasonal dataset

literature_2 <- literature %>%
  filter(pH>6.49999) #excluding probably superestimated pCO2 values from the global dataset

biomes_2 <- biomes %>%
  filter(pH>6.49999) #excluding probably superestimated pCO2 values from the Brazilian biome dataset

## Separating seasonal dataset by system 
cruz <-subset(seasonal, system=="Cruzeta")
garg <-subset(seasonal, system=="Gargalheiras")
extr <-subset(seasonal, system=="Extremoz")
bon <-subset(seasonal, system=="Bonfim")

####### Data analysis and visualization ########
### Figure 2: Violin plots showing spatial and seasonal variability of the main variables
#obs.: fine adjustments were made using Adobe Illustrator
##Fig 2A 
A2<-spatial[,-c(1,2,4:7,9:12,14,15,17,19:23,25:27,30:33,35:38,40,43,44,46,48,51:59,60,61)] %>% 
  gather(variable, value,-type) %>% 
  ggplot(aes(y = value, x=variable)) + 
  geom_violin(fill="grey",alpha=0.3)+
  geom_jitter(aes(color=type),alpha=0.5)+
  facet_wrap(~variable, scales = "free",strip.position = "left")+
  scale_colour_manual(values = c("#009E73", "#e79f00"))+
  ggtitle("A")+
  theme_bw()+
  theme(plot.title = element_text(size=18,face="bold"),
        strip.text = element_text(face="bold", size=12),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=14),
        axis.text.y = element_text(size=12, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= " ", y = " ",color= element_blank())

##Fig 2B
B2<-seasonal[,-c(1,2,6:9,11,12,14,16,17,19,22,23,26,30,32,33)] %>% 
  gather(variable, value,-system,-station,-season) %>% 
  ggplot(aes(y = value,x=system)) + 
  geom_violin(fill="grey",alpha=0.3) + 
  geom_jitter(aes(color=season),alpha=0.5)+
  facet_wrap(~variable, scales = "free_y",strip.position = "left")+
  ggtitle("B")+
  theme_bw()+
  theme(plot.title = element_text(size=18,face="bold"),
        strip.text = element_text(face="bold", size=12),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_blank(),
        legend.text = element_text(size=14),
        legend.position = "bottom",
        axis.text.y = element_text(size=12, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= " ", y = " ",color=element_blank())

g1 <- ggplotGrob(A2)
g2 <- ggplotGrob(B2)
g <-rbind(g1, g2, size = "first")
g$widths <- unit.pmax(g1$widths,g2$widths)
grid.newpage()
grid.draw(g) #fine adjustments were made using Adobe Illustrator

#### Prevalance of CO2 saturation ####
###Figure 3: Histogram with relative CO2 saturation
mycolors<-c("green4", "darkolivegreen3", "lightsteelblue3","darkblue") #choosing colors of the trophic state categories
micro <- expression(bold("390 "*mu*"atm"))

ggplot(spatial,aes(co2_rel, fill=trophic_state))+
  geom_histogram(binwidth = 4,alpha = 0.5,color="black")+
  geom_vline(xintercept = 0, linetype="dashed")+
  scale_y_continuous(expand=c(0,0), limits = c(0,18), breaks = seq(0,15,by=3))+
  scale_x_continuous(breaks = seq(-20,120,by=20))+
  scale_fill_manual(values=mycolors)+
  geom_label(aes(x = 0, y = 17, label = as.character(micro)), fill = "white",parse=TRUE)+
  annotate("text",x=-12,y=9,label="Undersaturated",size=6, angle=90)+
  annotate("text",x=70,y=9,label="Supersaturated",size=6, angle=90)+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= expression("Relative CO"[2]*" Saturation"), y = "Frequency",fill="Trophic State")

#### Literature comparisons ####
##Figure 4 + associated statistical analysis
#Fig 4A
fig4a<-ggplot(literature_2, aes(lat_cat_this_study, pCO2))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(outlier.shape = NA)+
  ggtitle("A")+
  theme_bw()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10,10^5))+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        #legend.position=c(.13,.85),
        legend.text = element_text(size=14),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = expression('pCO'[2]*' ('*mu*'atm)'),color = element_blank())

#Kruskal-Wallys test + Dunn's post-hoc test
fit<-kruskal.test(log10(pCO2) ~ lat_cat_this_study,data=literature_2)
dunn.test(literature_2$pCO2, literature_2$lat_cat_this_study)

#Fig 5B
fig5b<-ggplot(biomes, aes(biome,pCO2))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=source), alpha=0.5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10,10^5))+
  ggtitle("B")+
  theme_bw()+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.13,.13),
        legend.text = element_text(size=14),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = expression('pCO'[2]*' ('*mu*'atm)'),color = element_blank())

#Kruskal-Wallys test + Dunn's post-hoc test
fit<-kruskal.test(log10(pCO2) ~ biome,data=biomes_2)
dunn.test(biomes_2$pCO2, biomes_2$biome)

grid.arrange(fig5a, fig5b, ncol=2) #plotting Fig. 5 panel with Fig 5A and 5B

#Descriptive stats of the literature data
##latitudinal regions
stat.desc(literature_2$pCO2)
This <- literature_2 %>%
  filter(lat_cat_sobek_my_simpler=="This study")
stat.desc(This$pCO2)

Tropical <- literature_2 %>%
  filter(lat_cat_this_study=="Tropical")
stat.desc(Tropical$pCO2)

Subtropical <- literature_2 %>%
  filter(lat_cat_this_study=="Subtropical")
stat.desc(Subtropical$pCO2)

Temperado <- literature_2 %>%
  filter(lat_cat_sobek_my_simpler=="Temperate")
stat.desc(Temperado$pCO2)

Subpolar <- literature_2 %>%
  filter(lat_cat_this_study=="Subpolar")
stat.desc(Subpolar$pCO2)
Polar <- literature_2 %>%
  filter(lat_cat_this_study=="Polar")
stat.desc(Polar$pCO2)

##Brazilian biomes
amazon <- biomes_2 %>%
  filter(biome=="Amazonia Forest")
stat.desc(amazon$pCO2)
write.table(stat.desc(amazon$pCO2), "amazon.xls")

pantanal <- biomes_2 %>%
  filter(biome=="Pantanal floodplain")
stat.desc(pantanal$pCO2)
write.table(stat.desc(pantanal$pCO2), "pantanal.xls")

atlantic <- biomes_2 %>%
  filter(biome=="Atlantic Forest")
stat.desc(atlantic$pCO2)
write.table(stat.desc(atlantic$pCO2), "atlantic.xls")

caatinga <- biomes_2 %>%
  filter(biome=="Caatinga")
stat.desc(caatinga$pCO2)
write.table(stat.desc(caatinga$pCO2), "caatinga.xls")

#### Direct and inderect effects of land use, trophic state, carbon allochthony and water volume ####
###PCA scores
##Replacing NAs by the mean if the main variables prior to running PCA:
spatial_2[is.na(spatial_2[,"DOC"]), "DOC"] <- mean(spatial_2[,"DOC"], na.rm = TRUE)
spatial_2[is.na(spatial_2[,"PR"]), "PR"] <- mean(spatial_2[,"PR"], na.rm = TRUE)
spatial_2[is.na(spatial_2[,"a250.a365"]), "a250.a365"] <- mean(spatial_2[,"a250.a365"], na.rm = TRUE)
spatial_2[is.na(spatial_2[,"SUVA250"]), "SUVA250"] <- mean(spatial_2[,"SUVA250"], na.rm = TRUE)

##standardizing values of the main variables representatives of trophic state and C allochthony:
d_std<-decostand(spatial_2[, c("chla","DOC","PR","SUVA250", "TN", "TP", "SOC",
                               "a250.a365","secchi")], method= "standardize", na.rm=TRUE)
###Running PCA
pcaResult<-prcomp(d_std)
biplot(pcaResult) #plotting PCA results
biplot(pcaResult,choices = c(2,3)) #plotting PCA results
summary(pcaResult)

#Generating Scores
pc <- princomp(d_std, cor=TRUE, score=T)
scores<- as.data.frame(pc$scores)
#Inputing the PCA scores to the data.frame
spatial_2$PC1<-scores$Comp.1
spatial_2$PC2<-scores$Comp.2
spatial_2$PC3<-scores$Comp.3

## Structured equation model (SEM)
#standardizing variables in the dataset
dstd_data<-decostand(spatial_2[, c("DO","Dosat","temp","chla","DOC","TP","TN", "DN", "DP","PR", "BR", "BP","BGE","BA","SOC", 
                                   "a365", "a250", "a430", "SUVA250","a250.a365","LA_m2","CA_m2", "CA.LA","precip","vol",
                                   "LP.LA","depth","secchi", "anthropogenic_area", "forested_area")], method= "standardize", na.rm=TRUE)

#Just readding the PCA scores and pH values to the data.frame as they were already statandardized
dstd_data$PC1<-spatial_2$PC1
dstd_data$PC2<-spatial_2$PC2
dstd_data$PC3<-spatial_2$PC3
dstd_data$pH<-spatial_2$pH

#testing multinormality
result <- roystonTest(dstd_data, qqplot = TRUE)
uniPlot(dstd_data, type = "qqplot") 
uniPlot(dstd_data, type = "histogram") 

#running SEM
modelo='
pCO2 ~a*PC1+b*PC2+c*PC3
PC1~d*forested_area+e*anthropogenic_area+f*PC3
PC2~g*forested_area+h*anthropogenic_area+i*PC3

#efeitos indiretos

#Forest on pCO2 through trophic state (PC1)
da:=(d*a)
IndForTSCO2:=d*a
#Anthropic on pCO2 through trophic state (PC1)
ea:=(e*a)
IndAnthTSCO2:=e*a
#Hidrology on pCO2 through trophic state (PC1)
fa:=(f*a)
IndHydroTSCO2:=c+(f*a)
#Forest on pCO2 through allocthont (PC2)
gb:=(g*b)
IndForAlocCO2:=g*b
#Anthropic on pCO2 through allocthony (PC2)
hb:=(h*b)
IndAnthAlocCO2:=h*b
#Hidrology on pCO2 through allocthony (PC2)
ib:=(i*b)
IndHydrohAlocCO2:=c+(i*b)

'

modelo.fit = lavaan::sem(modelo, data = log_data)
modelo.fit = lavaan::sem(modelo, data = log_data, test = "bootstrap", se="bootstrap", verbose = TRUE, bootstrap = 1000) #bootstrap the test statistic
summary(modelo.fit, fit.measures = T, standardized = T, rsquare = T)
parTable(modelo.fit)
parameterEstimates(modelo.fit)
show(modelo.fit)
coef(modelo.fit)
fitted(modelo.fit)
summary(modelo.fit, fit.measures = T, standardized = T, rsquare = T)

###Figure 5 showing the SEM results was built based on the results above using the LibreOffice Impress (Linux)

######### Seasonal analysis #########
##### Running linear models ####
### Table 3: Best-fitting regression models
## Gargalheiras
m1 <-  lm(log10(pCO2) ~ 1, data=garg) #null model
m2 <- lm(log10(pCO2) ~ log10(chla), data=garg)
m3 <- lm(log10(pCO2) ~ log10(TN), data=garg)
m4 <- lm(log10(pCO2) ~ log10(N.P), data=garg)
m5 <-lm(log10(pCO2) ~ log10(vol), data=garg)
m6 <-lm(log10(pCO2) ~ log10(precip+1-min(precip)), data=garg)
m7 <-lm(log10(pCO2) ~ log10(depth), data=garg)
m8 <-lm(log10(pCO2) ~ log10(zeu), data=garg)
m9 <-lm(log10(pCO2) ~ log10(zeu.depth), data=garg)
m10<-  lm (log10(pCO2) ~ log10(chla)+log10(TN), data=garg)
m11<-  lm (log10(pCO2) ~ log10(chla)+log10(vol), data=garg)
m12 <-  lm (log10(pCO2) ~ log10(chla)+log10(depth), data=garg)
m13 <-  lm (log10(pCO2) ~ log10(depth)+log10(TN)+log10(chla), data=garg)
m14 <-  lm (log10(pCO2) ~ log10(depth)+log10(vol)+log10(chla), data=garg)
m15 <-lm(log10(pCO2) ~ log10(precip+1-min(precip))+log10(chla), data=garg)
m16 <-lm(log10(pCO2) ~ log10(precip+1-min(precip))+log10(depth)+log10(chla), data=garg)
m17 <-lm(log10(pCO2) ~ log10(precip+1-min(precip))+log10(vol)+log10(chla), data=garg)

mod_season<- list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17)
aictab(mod_season) #selecting best-fitting models
#selected best models:
summary(m17)
summary(m14)
summary(m11)

## Cruzetas:
m1 <-  lm(log10(pCO2) ~ 1, data=cruz) #null model
m2 <- lm(log10(pCO2) ~ log10(chla), data=cruz)
m3 <- lm(log10(pCO2) ~ log10(BR), data=cruz)
m4 <- lm(log10(pCO2) ~ log10(BR/chla), data=cruz)
m5<-  lm (log10(pCO2) ~ log10(chla)+log10(BR), data=cruz)
m6<-  lm (log10(pCO2) ~ log10(chla)+log10(vol), data=cruz)
m7<-  lm (log10(pCO2) ~ log10(BR)+log10(vol), data=cruz)
m8<-  lm (log10(pCO2) ~ log10(BR/chla)+log10(vol), data=cruz)
m9 <-  lm (log10(pCO2) ~ log10(chla)+log10(depth), data=cruz)
m10 <-  lm (log10(pCO2) ~ log10(BR)+log10(depth), data=cruz)
m11 <-  lm (log10(pCO2) ~ log10(BR/chla)+log10(depth), data=cruz)
m12 <-  lm (log10(pCO2) ~ log10(depth)+log10(BR)+log10(chla), data=cruz)
m13 <-  lm (log10(pCO2) ~ log10(depth)+log10(vol)+log10(chla), data=cruz)
m14 <-  lm (log10(pCO2) ~ log10(depth)+log10(vol)+log10(BR), data=cruz)
m15 <-lm(log10(pCO2) ~ log10(precip+1-min(precip))+log10(chla), data=cruz)
m16 <-lm(log10(pCO2) ~ log10(precip+1-min(precip))+log10(BR), data=cruz)
m17 <-lm(log10(pCO2) ~ log10(precip+1-min(precip))+log10(chla)+log10(BR), data=cruz)
m18 <-lm(log10(pCO2) ~ log10(precip+1-min(precip))+log10(depth)+log10(chla), data=cruz)
m19 <-lm(log10(pCO2) ~ log10(precip+1-min(precip))+log10(vol)+log10(chla), data=cruz)
m20 <-lm(log10(pCO2) ~ log10(precip+1-min(precip))+log10(depth)+log10(chla)+log10(BR), data=cruz)
m21 <-lm(log10(pCO2) ~ log10(precip+1-min(precip))+log10(vol)+log10(chla)+log10(BR), data=cruz)

mod_season<- list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21)
aictab(mod_season) #selecting best-fitting models
#selected best models:
summary(m16)
summary(m17)
summary(m21)
summary(m20)

## Extremoz:
m1 <-  lm(log10(pCO2) ~ 1, data=extr) #null model
m2 <- lm(log10(pCO2) ~ log10(chla), data=extr)
m3 <- lm(log10(pCO2) ~ log10(DOC), data=extr)
m4 <-lm(log10(pCO2) ~ log10(vol), data=extr)
m5 <-lm(log10(pCO2) ~ log10(precip+1-min(precip)), data=extr)
m6 <-lm(log10(pCO2) ~ log10(depth), data=extr)
m7<-  lm (log10(pCO2) ~ log10(chla)+log10(DOC), data=extr)
m8<-  lm (log10(pCO2) ~ log10(chla)+log10(vol), data=extr)
m9 <-  lm (log10(pCO2) ~ log10(chla)+log10(depth), data=extr)
m10 <-  lm (log10(pCO2) ~ log10(chla)+log10(precip+1-min(precip)), data=extr)
m11 <-  lm (log10(pCO2) ~ log10(DOC)+log10(precip+1-min(precip)), data=extr)
m12 <-  lm (log10(pCO2) ~ log10(depth)+log10(precip+1-min(precip))+log10(chla), data=extr)
m13 <-  lm (log10(pCO2) ~ log10(depth)+log10(vol)+log10(chla), data=extr)
m14 <-lm(log10(pCO2) ~ log10(precip+1-min(precip))+log10(vol)+log10(chla), data=extr)

mod_season<- list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14)
aictab(mod_season) #selecting best-fitting models
#selected best models
summary(m5)
summary(m10)
summary(m11)
summary(m1)

# Bonfim:
m1 <-lm(log10(pCO2) ~ 1, data=bon) #null model
m2 <-lm(log10(pCO2) ~ log10(chla), data=bon) 
m3 <-lm(log10(pCO2) ~ log10(TN), data=bon)
m4 <-lm(log10(pCO2) ~ log10(DOC), data=bon) 
m5 <-lm(log10(pCO2) ~ log10(vol), data=bon)
m6 <-lm(log10(pCO2) ~ log10(precip+1-min(precip)), data=bon)
m7 <-lm(log10(pCO2) ~ log10(depth), data=bon)
m8 <-lm(log10(pCO2) ~ log10(zeu.depth), data=bon)
m9<-  lm (log10(pCO2) ~ log10(chla)+log10(TN), data=bon)
m10<-  lm (log10(pCO2) ~ log10(chla)+log10(DOC), data=bon)
m11<-  lm (log10(pCO2) ~ log10(chla)+log10(vol), data=bon)
m12<-  lm (log10(pCO2) ~ log10(DOC)+log10(vol), data=bon)
m13<-  lm (log10(pCO2) ~ log10(DOC)+log10(zeu.depth), data=bon)

mod_season<- list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13)
aictab(mod_season) #selecting best-fitting models
#selected best models
summary(m13)
summary(m8)
summary(m4)

##### Figure 6: Sampling station and Seasonality effect on pCO2 using seasonal dataset ####
#6A) Boxplots with season and stations
A6.1<-ggplot(garg, aes(station,pCO2))+
  geom_boxplot(aes(fill=season),outlier.shape = NA)+
  coord_cartesian(ylim=c(0,8500))+
  #scale_fill_manual(values=c("lightgrey","grey40"))+
  geom_hline(yintercept = 390,color="black",linetype="dashed")+
  ggtitle("Gargalheiras")+
  theme(plot.title = element_text(size=16,face="bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "right",
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = expression('pCO'[2]*' ('*mu*'atm)'),fill = "Season")

B6.1<-ggplot(cruz, aes(station,pCO2))+
  geom_boxplot(aes(fill=season),outlier.shape = NA)+
  #scale_fill_manual(values=c("lightgrey","grey40"))+
  geom_hline(yintercept = 390,color="black",linetype="dashed")+
  coord_cartesian(ylim = c(0,8500))+
  ggtitle("Cruzeta")+
  theme(plot.title = element_text(size=16,face="bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "right",
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = expression('pCO'[2]*' ('*mu*'atm)'),fill = "Season")

C6.1<-ggplot(extr, aes(station,pCO2))+
  geom_boxplot(aes(fill=season),outlier.shape = NA)+
  #scale_fill_manual(values=c("lightgrey","grey40"))+
  geom_hline(yintercept = 390,color="black",linetype="dashed")+
  coord_cartesian(ylim = c(0,10000))+
  ggtitle("Extremoz")+
  theme(plot.title = element_text(size=16,face="bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "right",
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = "",fill = "Season")

D6.1<-ggplot(bon, aes(station,pCO2))+
  geom_boxplot(aes(fill=season),outlier.shape = NA)+
  #scale_fill_manual(values=c("lightgrey","grey40"))+
  geom_hline(yintercept = 390,color="black",linetype="dashed")+
  coord_cartesian(ylim = c(0,30000))+
  ggtitle("Bonfim")+
  theme(plot.title = element_text(size=16,face="bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "right",
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = "",fill = "Season")

A6<-ggarrange(A8.1,C8.1,B8.1,D8.1,ncol=2,nrow=2,widths = 1,heights = 1,common.legend=TRUE,legend="right")

#6B) volume responds to rainfall in the reservoirs, but not in the coastal lagoon
#Delta volume x Precipitation
B6<-ggplot(seasonal, aes(log10(precip+1-min(precip)),log10(d_vol+1-min(d_vol))))+
  stat_smooth(method = "lm",alpha = 0.4, color="black",size = 1) +
  geom_jitter(alpha=0.7,size=2)+
  facet_wrap(~system,scales="free",ncol=2)+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill="white",color = "black"),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(y= expression(Delta*Volume), x = "Precipitation (mm)")

#Obs.:The figures 6A and 6B were gathered using Adobe Illustrator

###Statistic in Fig.6:
#Independent and interactive effects (mean ±1 SE) of seasonality and sampling stations on pCO2 for each ecosystem
##Gargalheiras
m_g <- lm (log10(pCO2) ~ season + station, data = garg)
summary(m_g)#r2=35.16;p<0.0005
m_g <- lm (log10(pCO2) ~ season, data = garg)
summary(m_g)#r2=16.65;p<0.004
m_g <- lm (log10(pCO2) ~ station, data = garg)
summary(m_g)#r2=17.19;p<0.01

##Cruzeta
m_c <- lm (log10(pCO2) ~ season + station, data = cruz)
summary(m_c)
m_c <- lm (log10(pCO2) ~ season, data = cruz)
summary(m_c)
m_c <- lm (log10(pCO2) ~ station, data = cruz)
summary(m_c)

#Extremoz
m_e <- lm (log10(pCO2) ~ season + station, data = extr)
summary(m_e)#r2=16.1;p=0.05
m_e <- lm (log10(pCO2) ~ season, data = extr)
summary(m_e)#r2=14.73;p<0.05
m_e <- lm (log10(pCO2) ~ station, data = extr)
summary(m_e)

#Bonfim
m_b <- lm (log10(pCO2) ~ season + station, data = bon)
summary(m_b)
m_b <- aov(log10(pCO2) ~ season, data = bon)
summary(m_b)
m_b <- lm (log10(pCO2) ~ station, data = bon)
summary(m_b)

#Water volume variation in response to rainfall (Figure 6B)
mdvolprecip_garg<-  lm (log10(d_vol+1-min(d_vol)) ~ log10(precip+1-min(precip)), data=garg) 
summary(mdvolprecip_garg)
mdvolprecip_cruz <-lm(log10(d_vol+1-min(d_vol)) ~ log10(precip+1-min(precip)), data=cruz)
summary(mdvolprecip_cruz)
mdvolprecip_extr <- lm(log10(d_vol+1-min(d_vol)) ~ log10(precip+1-min(precip)), data=extr)
summary(mdvolprecip_extr)
mdvolprecip_bon<-  lm (log10(d_vol+1-min(d_vol)) ~ log10(precip+1-min(precip)), data=bon)
summary(mdvolprecip_bon)

###### Table 4 #####
##General linear mixed-models (GLMM) of the analysis of variance for the effects of seasonality,reservoirs and their interaction on pCO2
model = lmer(log10(pCO2) ~ season * system + (1|station), REML = FALSE,data=seasonal_2)
summary(model)
anova(model, test="F")
lsmeans(model,adjust="tukey")
plot(model)
resid(model)
#post-hoc
TukeyHSD(model) #Tukey test
#Cruzeta-Gargalheiras: p<0.05
#Extremoz-Gargalheiras: p<0.0005
#Bonfim-Gargalheiras: p<0.000005
#Extremoz-Cruzeta: p=0.19
#Bonfim-Cruzeta:p<0.05
#Bonfim-Extremoz: p=0.98
test<-glm(fit_2)
summary(test)
