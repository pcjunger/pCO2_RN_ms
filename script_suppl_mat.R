##Manucript: CO2 Supersaturation in Tropical Eutrophic Waters
#DOI number:10.1016/j.scitotenv.2019.01.273
###Script used to produce the figures and statistical analysis presented in the supplementary material
### Authors: Pedro C. Junger et al.

#### Loading required libraries ####
library(tidyr)
library(dplyr)
library(ggplot2)
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

###### External Functions #####
#function to extract and depict regression parameters in figure
lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}

#funtion to extract min positive value
minpositive = function(x) min(x[x > 0])

#### Data input ####
spatial <- read.csv("spatial_dataset.csv",header = TRUE) #full spatial dataset from 102 ecosystems used in this study
seasonal <- read.csv("seasonal_dataset.csv", sep = ",",header = TRUE) #full temporal dataset used in this study
mean <- read.csv("mean_data.csv", sep=",",header = TRUE) #dataset with mean values used in figures in supplementary material
sur_bot<- read.csv("surf_bottom.csv", header=TRUE) #dataset with surface and bottom pCO2

#### Data wrangling ####

##Ordering factors
#for spatial dataset
spatial$trophic_state <- factor(spatial$trophic_state, levels=c("Hypereutrophic","Eutrophic","Mesotrophic","Oligotrophic"))

#for seasonal dataset
seasonal$system <- factor(seasonal$system, levels = c("Gargalheiras", "Cruzeta","Extremoz","Bonfim"))
seasonal$date <- factor(seasonal$date, levels =c ("10-Dec","11-Jan","11-Feb","11-Mar","11-Apr","11-May",
                                                  "11-Jun","11-Jul","11-Aug","11-Sep","11-Oct","11-Nov",
                                                  "11-Dec","12-Jan","12-Jul","12-Aug","12-Sep","12-Oct","12-Nov","12-Dec",
                                                  "13-Jan","13-Feb","13-Mar","13-Apr","13-May","13-Jun",
                                                  "13-Jul","14-Sep","14-Oct","14-Nov","14-Dec","15-Jan","15-Feb",
                                                  "15-Mar","15-Apr","15-May","15-Jun","15-Jul","15-Aug",
                                                  "15-Sep","15-Oct"))
#for mean values
mean$system <- factor(mean$system, levels = c("Gargalheiras", "Cruzeta","Extremoz","Bonfim"))

## Organizing datasets
seasonal$system <- factor(seasonal$system, levels = c("Gargalheiras", "Cruzeta","Extremoz","Bonfim"))

##Filtering the spatial dataset by excluding observations with Total Alkalinity < 1meq/L
#Obs.: See detailed reasons in the methods (section 2.4) of the manuscript
spatial_2<- spatial %>%
  filter(lakeID!=67) %>% #excluding outlier "Riacho da Cruz"
  filter(TA>1) #excluding superestimated pCO2 values from the full spatial dataset

seasonal_2 <- seasonal %>%
  filter(TA>1) #excluding superestimated pCO2 values from the full seasonal dataset

## Separating seasonal dataset by system 
cruz <-subset(seasonal, system=="Cruzeta")
cruz$date <- factor(cruz$date,levels =c("10-Dec","11-Jan","11-Feb","11-Mar","11-Apr","11-May",
                                        "11-Jun","11-Jul","11-Aug","11-Sep","11-Oct","11-Nov",
                                        "11-Dec","12-Jan"))

garg <-subset(seasonal, system=="Gargalheiras")
garg$date <- factor(garg$date,levels =c("10-Dec","11-Jan","11-Feb","11-Mar","11-Apr","11-May",
                                        "11-Jun","11-Jul","11-Aug","11-Sep","11-Oct","11-Nov",
                                        "11-Dec","12-Jan"))
extr <-subset(seasonal, system=="Extremoz")
extr$date <- factor(extr$date,levels =c("12-Jul","12-Aug","12-Sep","12-Oct","12-Nov","12-Dec",
                                        "13-Jan","13-Feb","13-Mar","13-Apr","13-May","13-Jun",
                                        "13-Jul"))
bon <-subset(seasonal, system=="Bonfim")
bon$date <-factor(bon$date,levels =c("14-Sep","14-Oct","14-Nov","14-Dec","15-Jan","15-Feb",
                                     "15-Mar","15-Apr","15-May","15-Jun","15-Jul","15-Aug",
                                     "15-Sep","15-Oct"))
#creating new required columns
spatial$DOC_chla<- spatial$DOC / spatial$chla*0.04 #creating column with doc:chla ratio
spatial$BR_chla <- spatial$BR / (spatial$chla*40) #creating column with BR:chla ratio
spatial$chla_m3_C <- (spatial$chla*40)/1000
spatial$csoil_chla <- spatial$SOC / spatial$chla_m3_C

#excluding outliers
dochla_data<-subset(spatial_2, DOC_chla<10) #excluding DOC:chla outlier
chla_data<-subset(spatial_2, chla>3.08) #excluding Chl-a outlier from the regression
PR_data<-subset(spatial_2, PR>7) #excluding PR outlier from the regression

#replacing zeros by the lowest value by the minimum value in each category
spatial$fo_ca<- ifelse(spatial$fo_ca==0, minpositive(spatial$fo_ca), spatial$fo_ca)
spatial$anth_ca <- ifelse(spatial$anth_ca==0, minpositive(spatial$anth_ca), spatial$anth_ca)

#replacing 1 by 0.99999 before logit transformation
spatial$fo_ca<- ifelse(spatial$fo_ca==1, 0.99999, spatial$fo_ca)
spatial$anth_ca <- ifelse(spatial$anth_ca==1, 0.99999, spatial$anth_ca)

#logit transformation of land-use data
spatial["fo_logit"]<- logit(spatial$fo_ca)
spatial["anth_logit"]<- logit(spatial$anth_ca)

###### Supplementary figures ######
#Fig S1: Seasonal trends in pCO2, volume, precipitation
## pCO2, chla, rainfall and volume over time in all systems
#precipitação
S1A <- ggplot(mean, aes(date,precip))+
  geom_bar(data=mean[!is.na(mean$precip),],stat="summary", fun.y = "mean", position = "dodge", fill="grey", colour="black",width = 0.7)+
  facet_wrap(~system,scales ="free",ncol = 2)+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill="white",color = "black"),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black",angle=45,hjust=1),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = "Rainfall (mm)")

#Volume
S1B <- ggplot(mean, aes(date,vol/1000))+
  geom_path(group=1)+
  facet_wrap(~system,scales ="free",ncol = 2)+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill="white",color = "black"),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black",angle=45,hjust=1),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = expression("Water volume ("*mm^{3}*")"))

S1C <- ggplot(mean, aes(date,pCO2))+
  geom_point()+
  geom_path(data=mean[!is.na(mean$pCO2),],group=1)+
  stat_boxplot(geom="errorbar",data = seasonal)+
  geom_hline(yintercept = 390,linetype="dashed")+
  facet_wrap(~system,scales ="free",ncol = 2)+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill="white",color = "black"),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black",angle=45,hjust=1),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = expression('pCO'[2]*' ('*mu*'atm)'),shape = element_blank())

S1D <- ggplot(mean, aes(date,chla))+
  geom_point()+
  geom_path(data=mean[!is.na(mean$chla),],group=1)+
  stat_boxplot(geom="errorbar",data = seasonal)+
  facet_wrap(~system,scales ="free",ncol = 2)+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "none",
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill="white",color = "black"),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black",angle=45,hjust=1),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = expression('Chlorophyll-a ('*mu*'g/L)'),shape = element_blank())

grid.arrange(S1A,S1B,S1C,S1D, nrow=2)

#Figure S2: density plots with two datasets for comparison
S2A <- ggplot(spatial, aes(TA,pCO2))+
  geom_jitter(color="black",shape=1)+
  geom_vline(xintercept = 1,linetype="dashed")+
  ggtitle("A")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.85,.90),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16))+
  labs(x= "TA (meq/L)", y = expression('pCO'[2]*' ('*mu*'atm)'),color=element_blank())

SB2 <- ggplot(spatial, aes(x=pCO2))+
  geom_density(aes(fill="All lakes"),alpha=0.5) +
  geom_text(aes(label="N=102",x =11000, y=0.5))+
  geom_density(data=spatial_2, aes(fill="TA > 1 meq/L"),alpha=0.5)+
  geom_text(aes(label="N=44",x = 2000,y=0.5))+
  scale_fill_grey()+
  scale_x_log10()+
  ggtitle("B")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.15,.90),
        legend.text = element_text(size =12),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16))+
  labs(x= expression('pCO'[2]*' ('*mu*'atm)'), y = "Density",fill = element_blank())

S2C <- ggplot(seasonal, aes(TA,pCO2))+
  geom_jitter(color="black",shape=1)+
  geom_vline(xintercept = 1,linetype="dashed")+
  ggtitle("C")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.85,.90),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16))+
  labs(x= "TA (meq/L)", y = expression('pCO'[2]*' ('*mu*'atm)'),color=element_blank())

S2D <- ggplot(seasonal, aes(x=pCO2))+
  geom_density(aes(fill="All samplings"),alpha=0.5) +
  geom_text(aes(label="N=152",x =12000, y=0.4))+
  geom_density(data=seasonal_2, aes(fill="TA > 1 meq/L"),alpha=0.5)+
  geom_text(aes(label="N=102",x = 2000,y=0.4))+
  scale_fill_grey()+
  scale_x_log10()+
  ggtitle("D")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.15,.90),
        legend.text = element_text(size =12),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16))+
  labs(x= expression('pCO'[2]*' ('*mu*'atm)'), y = "Density",fill = element_blank())

grid.arrange(S2A,S2B,S2C,S2D,ncol=2)

#Figure S3.: See PCA results section in this script

#Figure S4:Frequency of pCO2 values based on sampling day period (morning or afternoon).
## This figure was made with Graph Prism

#Figure S5: Box-plots comparing pCO2 values between (A) reservoirs and lakes and between (B) humid and semiarid systems 
##5SA
t.test(spatial_2$pCO2~spatial_2$type) #student t-test
A5S<-ggplot(spatial_2, aes(type,pCO2))+
  geom_boxplot()+
  geom_jitter()+
  scale_y_log10()+
  ggtitle("A)")+
  theme_bw()+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.13,.13),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = expression('pCO'[2]*' ('*mu*'atm)'),shape = element_blank())

##5SB
t.test(spatial_2$pCO2~spatial_2$region) #student t-test
B5S<-ggplot(spatial_2, aes(region,pCO2))+
  geom_boxplot()+
  geom_jitter()+
  scale_y_log10()+
  ggtitle("B)")+
  theme_bw()+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.13,.13),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = "",shape = element_blank())

grid.arrange(A5S, B5S,ncol=2)

#Figure S6: Scatter plots: regressions with significant internal variables (spatial and seasonal)
## CO2 x chla relationship: eutrophication effect on pCO2
#S6A: spatial
m <-  lm (log10(pCO2) ~ log10(chla), data=dochla_data)
summary(m)
S6A<-ggplot(dochla_data, aes(chla, pCO2))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter(aes(colour=dayperiod),size=2,alpha=0.8)+
  scale_colour_manual(values = c("#000000", "#009E73", "#e79f00"))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10,10^5))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(1,10^2.75))+
  geom_hline(yintercept = 390,color="black",linetype="dashed")+
  ggtitle("A")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.18,.18),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = expression('pCO'[2]*' ('*mu*'atm)'),colour = element_blank())+
  annotate("text",x=10^2,y=10^5,label=lm_eqn(m),parse=TRUE)

#S6B: Seasonal
m <-  lm (log10(pCO2) ~ log10(chla), data=seasonal_2)#Bonfim is overestimated, thus it was excluded from the dataset for this analysis
summary(m)
S6B<-ggplot(seasonal_2, aes(chla,pCO2))+
  stat_smooth(method = "lm",fill="lightgrey", alpha = 0.5,size = 1,color="black") +
  geom_jitter()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10,10^5))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(1,10^2.75))+
  geom_hline(yintercept = 390,color="black",linetype="dashed")+
  ggtitle("B")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = c(.13,.15),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= expression("Chlorophyll-a ("*mu*"g/L)"), y = expression('pCO'[2]*' ('*mu*'atm)'))+
  annotate("text",x=10^2,y=10^5,label=lm_eqn(m),parse=TRUE)

#S6C) pCO2 vs Vol
m <-lm (log10(pCO2) ~ log10(vol), data=spatial_2)
summary(m)
S6C<-ggplot(spatial_2, aes(vol, pCO2))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter(size=2,alpha=0.7)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10,10^5))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10^2.8,10^7.7))+
  geom_hline(yintercept = 390,color="black",linetype="dashed")+
  ggtitle("C")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.13,.13),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = expression('pCO'[2]*' ('*mu*'atm)'),shape = element_blank())+
  annotate("text",x=10^6.3,y=10^5,label=lm_eqn(m),parse=TRUE)

#S6D) Chla vs vol

m <-  lm (log10(chla) ~ log10(vol), data=dochla_data)
summary(m)
S6D<-ggplot(dochla_data, aes(vol, chla))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter(size=2,alpha=0.7)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(3,10^3))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10^2.8,10^7.7))+
  ggtitle("D")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.13,.13),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "Water volume (m³)", y = expression("Chlorophyll-a ("*mu*"g/L)"),shape = element_blank())+
  annotate("text",x=10^6.2,y=10^3,label=lm_eqn(m),parse=TRUE)

plots1 <- list(S6A,S6B)
plots2 <- list(S6C,S6D)
grobs1 = lapply(plots1, ggplotGrob)
grobs2 = lapply(plots2, ggplotGrob)
g1 = do.call(cbind, c(grobs1, size="first"))
g2 = do.call(cbind, c(grobs2, size="first"))
g3 = do.call(rbind, c(list(g1,g2), size="first"))
grid.newpage()
grid.draw(g3)

#Figure S7: pCO2 x DO
#S7A)
m <- lm(log10(pCO2) ~ log10(Dosat), data=spatial_2[-20,]) #r2=-25.4; p=0.0003
summary(m)
AS7<-ggplot(spatial_2[-20,], aes(Dosat, pCO2))+ #outlier excluded(lakeID=61;DO=3.2mg/L)
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(14,10^5))+
  geom_hline(yintercept = 390,color="black",linetype="dashed")+
  ggtitle("A) Spatial")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "Dissolved Oxygen (%)", y = expression('pCO'[2]*' ('*mu*'atm)'))+
  annotate("text",x=160,y=10^5,label=lm_eqn(m),parse=TRUE)

#S7B)
m<- lm(log10(pCO2) ~ log10(DO_sat), data=seasonal_2) 
summary(m)
BS7<-ggplot(seasonal_2, aes(DO_sat,pCO2))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(14,10^5))+
  geom_hline(yintercept = 390,color="black",linetype="dashed")+
  ggtitle("B) Seasonal")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16))+
  labs(x= "Dissolved Oxygen (%)", y = expression('pCO'[2]*' ('*mu*'atm)'))+
  annotate("text",x=280,y=10^5,label=lm_eqn(m),parse=TRUE)

#S7C
m <- lm(log10(pCO2) ~ log10(DO_sat), data=garg) 
summary(m)
S7C<-ggplot(garg, aes(DO_sat,pCO2))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10,10^5))+
  geom_hline(yintercept = 390,color="black",linetype="dashed")+
  ggtitle("C) Gargalheiras")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16))+
  labs(x= "Dissolved Oxygen (%)", y = expression('pCO'[2]*' ('*mu*'atm)'))+
  annotate("text",x=240,y=10^5,label=lm_eqn(m),parse=TRUE)

#S7D
m <- lm(log10(pCO2) ~ log10(DO_sat), data=cruz) 
summary(m)
S7D<-ggplot(cruz, aes(DO_sat,pCO2))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10,10^5))+
  geom_hline(yintercept = 390,color="black",linetype="dashed")+
  ggtitle("D) Cruzeta")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16))+
  labs(x= "Dissolved Oxygen (%) ", y = expression('pCO'[2]*' ('*mu*'atm)'))+
  annotate("text",x=280,y=10^5,label=lm_eqn(m),parse=TRUE)

ggarrange(S7A,S7B,S7C,S7D,ncol=2,nrow=2,widths = 1,heights = 1) #fine adjustments were made using Adobe Illustrator

#Figure S8: Volume X Lake Area x pCO2
#S58: Lake area vs Volume
m <-  lm(log10(vol) ~ log10(lake_area_m2), data=spatial_2)
summary(m)
S8A<-ggplot(spatial_2, aes(lake_area_m2, vol))+
  stat_smooth(method = "lm",alpha = 0.5, fill="lightgrey", size = 1,color="black") +
  geom_jitter()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10^3.5,10^7.3))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10^2.8,10^7.8))+
  ggtitle("A")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = "Water volume (m³)")+
  annotate("text",x=10^4.5,y=10^7.7,label=lm_eqn(m),parse=TRUE)

#S8B)Lake area
m <-  lm(log10(pCO2) ~ log10(lake_area_m2), data=spatial_2)
summary(m)

S8B<-ggplot(spatial_2, aes(lake_area_m2, pCO2))+
  stat_smooth(method = "lm",alpha = 0.5, fill="lightgrey", size = 1,color="black") +
  geom_jitter()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(14,10^5))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10^3.5,10^7.3))+
  ggtitle("B")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "Lake area (m²)", y = expression('pCO'[2]*' ('*mu*'atm)'))+
  annotate("text",x=10^4.5,y=10^5,label=lm_eqn(m),parse=TRUE)

grid.arrange(S8A, S8B, ncol=1)

#FigS9: Internal processes 2: pCO2 vs DOC, pCO2 vs DOC/chla, DOC vs Chla, a365 vs DOC
#S9A) pCO2 x DOC
m <-  lm (log10(pCO2) ~ log10(DOC), data=spatial_2)
summary(m)
S9A<-ggplot(spatial_2, aes(DOC, pCO2))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10,10^5))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10,100))+
  geom_hline(yintercept = 390,color="black",linetype="dashed")+
  ggtitle("A")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "DOC (mg/L) ", y = expression('pCO'[2]*' ('*mu*'atm)'))+
  annotate("text",x=10^1.28,y=10^5,label=lm_eqn(m),parse=TRUE)

#S9B)
m <-  lm (log10(DOC) ~ log10(chla), data=spatial_2)
summary(m)
S9B<-ggplot(spatial_2, aes(chla,DOC))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(2,10^2.75))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(8,110))+
  ggtitle("B")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16))+
  labs(y= "DOC (mg/L)", x=expression("Chlorophyll-a ("*mu*"g/L)"))+
  annotate("text",x=10^1,y=10^2.03,label=lm_eqn(m),parse=TRUE)

#S9C) pCO2 vs BR:Chla
m <-  lm(log10(pCO2) ~ log10(BR_chla), data=PR_data)
summary(m)
S9C<-ggplot(PR_data, aes(BR_chla,pCO2))+
  stat_smooth(method = "lm",alpha = 0.5, fill="lightgrey", size = 1,color="black") +
  geom_jitter()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10,10^5))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10^-3.2,10^-0.3))+
  ggtitle("C")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= expression("BR ("*mu*"g C h⁻¹): Chl-a ("*mu*"g C)"), y = expression('pCO'[2]*' ('*mu*'atm)'))+
  annotate("text",x=10^-2.3,y=10^5,label=lm_eqn(m),parse=TRUE)

#S(D)) DOC x a365 (CDOM):strong positive relationship indicates high aromatic character of DOM
m<-  lm(log10(a365) ~ log10(DOC), data=spatial_2)
summary(m)
S9D<-ggplot(spatial_2, aes(DOC,a365))+
  stat_smooth(method = "lm",alpha = 0.5, fill="lightgrey", size = 1,color="black") +
  geom_jitter()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10^0.5,10^2))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(8,110))+
  ggtitle("D")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16))+
  labs(y= expression("a"[365]*~(~m^-1)), x="DOC (mg/L)")+
  annotate("text",x=10^1.25,y=10^2,label=lm_eqn(m),parse=TRUE)

plots1 <- list(S9A,S9B)
plots2 <- list(S9C,S9D)
grobs1 = lapply(plots1, ggplotGrob)
grobs2 = lapply(plots2, ggplotGrob)
g1 = do.call(cbind, c(grobs1, size="first"))
g2 = do.call(cbind, c(grobs2, size="first"))
g3 = do.call(rbind, c(list(g1,g2), size="first"))
grid.newpage()
grid.draw(g3)

#FigS10: Sediment influence: Surface pCO2 x Bottom pCO2 in Gargalheiras and Cruzeta
# Surface vs Bottom 
m <- lm(log10(pCO2_s) ~ log10(pCO2_b), sur_bot) #R2=40.5; p<0.0005
summary(m)

ggplot(sur_bot, aes(pCO2_b,pCO2_s))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter(size=2,alpha=0.7, aes(shape=system))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(800,10^4.6))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10,14000))+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(y= expression('Surface pCO'[2]*' ('*mu*'atm)'), x= expression('Bottom pCO'[2]*' ('*mu*'atm)'),shape = "System")+
  annotate("text",x=10^3.3,y=10^4,label=lm_eqn(m),parse=TRUE)

#Figure S11: diel variation in Gargalheiras
## This figure was made with STATISTICA 9.0

#Figure S12: Land use impact on variables of interest
##S12A) Chla x Forest
m <- lm(log10(chla) ~ fo_logit, chla_data)
summary(m)
AS12<-ggplot(chla_data, aes(fo_logit,chla))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter(size=2,alpha=0.7)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(3,10^3))+
  ggtitle("A")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.13,.13),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = expression("Chlorophyll-a ("*mu*"g/L)"),shape = element_blank())+
  annotate("text",x=0,y=10^3,label=lm_eqn(m),parse=TRUE)

##S12B) Chla x Anthropogenic
m <- lm(log10(chla) ~ anth_logit, chla_data) #R2=26.9; p<0.001
summary(m)

BS12<-ggplot(chla_data, aes(anth_logit,chla))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter(size=2,alpha=0.7)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(3,10^3))+
  ggtitle("B")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.13,.13),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = "",shape = element_blank())+
  annotate("text",x=2,y=10^3,label=lm_eqn(m),parse=TRUE)

##S12C) DOC:Chla x Forest
m <- lm(log10(DOC_chla) ~ fo_logit, chla_data) #r2=19.8, p<0.005
summary(m)
C8<-ggplot(chla_data, aes(fo_logit,DOC_chla))+ #relação negativa, extranho, contradiz o anterior
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter(size=2,alpha=0.7)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(0.0025,0.5))+
  ggtitle("C")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.13,.13),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = expression("DOC:Chl-a (mgC/L)"),shape = element_blank())+
  annotate("text",x=0,y=10^-0.32,label=lm_eqn(m),parse=TRUE)

##S12D) DOC:Chla x Anthropogenic
m <- lm(log10(DOC_chla) ~ anth_logit, chla_data)
summary(m)
DS12<-ggplot(chla_data, aes(anth_logit,DOC_chla))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter(size=2,alpha=0.7)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(0.0025,0.5))+
  ggtitle("D")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.13,.13),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = "",shape = element_blank())+
  annotate("text",x=1.8,y=10^-0.4,label=lm_eqn(m),parse=TRUE)

##S12E) SOC:Chla x Forest
m <- lm(log10(csoil_chla) ~ fo_logit, chla_data)
summary(m)
ES12<-ggplot(chla_data, aes(fo_logit,csoil_chla))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter(size=2,alpha=0.7)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(0.1,10^1.5))+
  ggtitle("E")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.13,.13),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = "SOC:Chl-a",shape = element_blank())+
  annotate("text",x=0,y=10^1.5,label=lm_eqn(m),parse=TRUE)

##S12F) SOC:Chla x Anthropogenic
m <- lm(log10(csoil_chla) ~ anth_logit, chla_data)
summary(m)

FS12<-ggplot(chla_data, aes(anth_logit,csoil_chla))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter(size=2,alpha=0.7)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(0.1,10^1.5))+
  ggtitle("F")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.13,.13),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "", y = "",shape = element_blank())+
  annotate("text",x=2,y=10^1.5,label=lm_eqn(m),parse=TRUE)

##S12G) BR:Chla x Forest
#Forest
m <- lm(log10(BR_chla) ~ fo_logit, chla_data) #+R2= 12.31, p<0.05
summary(m)
GS12<-ggplot(chla_data, aes(fo_logit,BR_chla))+ 
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter(size=2,alpha=0.7)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(0.0008,0.3))+
  ggtitle("G")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.13,.13),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "Logit (Proportion of Forested area)", y = expression("BR ("*mu*"g C h⁻¹):Chl-a ("*mu*"g C)"),shape = element_blank())+
  annotate("text",x=0,y=10^-0.55,label=lm_eqn(m),parse=TRUE)

##S12H) BR:Chla x Anthropogenic
m <- lm(log10(BR_chla) ~ anth_logit, chla_data) #R2=-12.5, p<0.05
summary(m)
HS12<-ggplot(chla_data, aes(anth_logit,BR_chla))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter(size=2,alpha=0.7)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(0.0008,0.3))+
  ggtitle("H")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position=c(.13,.13),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "Logit (Proportion of Anthropogenic area)", y = "",shape = element_blank())+
  annotate("text",x=1.6,y=10^-0.55,label=lm_eqn(m),parse=TRUE)

grid.arrange(AS12,BS12,CS12,DS12,ES12,FS12,GS12,HS12, ncol=2)