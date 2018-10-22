##Manucript: CO2 Supersaturation in Tropical Eutrophic Waters
###Script used to run the statistical analysis presented in the results and the supplementary material throughout the manuscript
### Authors: Pedro C. Junger et al.
###### 1.1 - loading required libraries #####
library(dplyr)
library(tidyr)
library(ggpubr)
library(pastecs)
library(gtools)
library(visreg)
library(AICcmodavg)
library(car)
library(Rmisc)
library(lme4)
library(lmerTest)

###### 2.1 - Data input ######

spatial <- read.csv("spatial_dataset.csv",header = TRUE) #full spatial dataset from 102 ecosystems used in this study
seasonal <- read.csv("seasonal_dataset.csv", sep = ",",header = TRUE) #full temporal dataset used in this study

##### 2.1.1 - Creating new columns with useful variables #####
spatial$BR_chla <- spatial$BR / (spatial$chla*40) #creating column with BR:chla ratio
spatial$chla_m3_C <- (spatial$chla*40)/1000
spatial$csoil_chla <- spatial$SOC / spatial$chla_m3_C

#creating column with anthropogenic land-use (agriculture+pasture+urban areas)
spatial["anth_ca"]<- spatial$agrpast_ca + spatial$urban_ca
#creating column w/ proportion of land-use
spatial["ca_fo"]<- spatial$fo_ca/100
spatial["ca_agr"]<- spatial$agrpast_ca/100
spatial["ca_urb"]<- spatial$urban_ca/100
spatial["ca_anth"]<- spatial$ca_agr + spatial$ca_urb

#replacing zeros by the lowest value by the minimum value in each category
spatial$ca_fo <- ifelse(spatial$ca_fo==0, minpositive(spatial$ca_fo), spatial$ca_fo)
spatial$ca_anth <- ifelse(spatial$ca_anth==0, minpositive(spatial$ca_anth), spatial$ca_anth)

#replacing 1 by 0.99999 before logit transformation
spatial$ca_fo <- ifelse(spatial$ca_fo==1, 0.99999, spatial$ca_fo)
spatial$ca_anth <- ifelse(spatial$ca_anth==1, 0.99999, spatial$ca_anth)

#logit transformation of land-use data
spatial["fo_logit"]<- logit(spatial$ca_fo)
spatial["anth_logit"]<- logit(spatial$ca_anth)

###### 2.1.2 - Filtering the spatial dataset by excluding observations with Total Alkalinity < 1meq/L ######
#Obs.: See detailed reasons in the methods (section 2.5) of the manuscript
spatial<- spatial %>%
  filter(lakeID!=67) #excluding outlier "Riacho da Cruz"
spatial_2 <- spatial %>%
  filter(TA>1) #excluding superestimated pCO2 values from the full spatial dataset
seasonal_2 <- seasonal %>%
  filter(TA>1) #excluding superestimated pCO2 values from the full seasonal dataset

###### 2.1.3 - Separating seasonal dataset by system ####
cruz <-subset(seasonal, system=="Cruzeta")
garg <-subset(seasonal, system=="Gargalheiras")
extr <-subset(seasonal, system=="Extremoz")
bon <-subset(seasonal, system=="Bonfim")

#### 3 - Statistical Analysis ####
##the codes for linear regressions presented in the figures are in the "script_figures.R" file available in the same github repository
##3.1. Best-fitting models selected with AICc (Table 2)
# 3.1.1. Spatial dataset - significant models:
m1 <- lm(log10(pCO2) ~ 1, data=spatial_2) #null model
m2 <- lm(log10(pCO2) ~ log10(depth), data=spatial_2)
m3 <- lm(log10(pCO2) ~ log10(vol), data=spatial_2)
m4 <- lm(log10(pCO2) ~ log10(chla), data=spatial_2)
m5 <- lm(log10(pCO2) ~ log10(DN), data=spatial_2)
m6 <- lm(log10(pCO2) ~ log10(DOC), data=spatial_2)
m7 <- lm(log10(pCO2) ~ log10(lake_area_m2), data=spatial_2)
m8 <- lm(log10(pCO2) ~ log10(real_SW_area_m2), data=spatial_2)
m9 <- lm(log10(pCO2) ~ log10(C_soil), data=spatial_2)
m10 <- lm(log10(pCO2) ~ log10(csoil_chla), data=spatial_2)
m11 <- lm(log10(pCO2) ~ log10(BR_chla), data=spatial_2)
m12 <- lm(log10(pCO2) ~ log10(chla)+log10(DOC), data=spatial_2)
m13 <- lm(log10(pCO2) ~ log10(chla)+log10(DN), data=spatial_2)
m14 <- lm(log10(pCO2) ~ log10(chla)+log10(vol), data=spatial_2)
m15 <- lm(log10(pCO2) ~ log10(chla)+log10(C_soil), data=spatial_2)
m16 <- lm(log10(pCO2) ~ log10(C_soil)+log10(real_SW_area_m2), data=spatial_2)
m17 <- lm(log10(pCO2) ~ log10(DOC)+log10(C_soil), data=spatial_2)
m18 <- lm(log10(pCO2) ~ log10(csoil_chla)+log10(DOC), data=spatial_2)
m19 <- lm(log10(pCO2) ~ log10(csoil_chla)+log10(DN), data=spatial_2)
m20 <- lm(log10(pCO2) ~ log10(csoil_chla)+log10(DN), data=spatial_2)
m21 <- lm(log10(pCO2) ~ log10(csoil_chla)+log10(vol), data=spatial_2)
m22 <- lm(log10(pCO2) ~ log10(csoil_chla)+log10(real_SW_area_m2), data=spatial_2)
m23 <- lm(log10(pCO2) ~ log10(csoil_chla)+log10(lake_area_m2), data=spatial_2)
m24 <- lm(log10(pCO2) ~ log10(BR_chla)+log10(DOC), data=spatial_2)
m25 <- lm(log10(pCO2) ~ log10(BR_chla)+log10(DN), data=spatial_2)
m26 <- lm(log10(pCO2) ~ log10(BR_chla)+log10(DOC), data=spatial_2)
m27 <- lm(log10(pCO2) ~ log10(BR_chla)+log10(lake_area_m2), data=spatial_2)
m28 <- lm(log10(pCO2) ~ log10(BR_chla)+log10(real_SW_area_m2), data=spatial_2)
m29 <- lm(log10(pCO2) ~ log10(BR_chla)+log10(vol), data=spatial_2)
m30 <- lm(log10(pCO2) ~ log10(BR_chla)+log10(C_soil), data=spatial_2)
m31 <- lm(log10(pCO2) ~ log10(BR_chla)+log10(vol), data=spatial_2)
m32 <- lm(log10(pCO2) ~ log10(DOC)+log10(C_soil)+log10(chla), data=spatial_2)
m33 <- lm(log10(pCO2) ~ log10(DOC)+log10(C_soil)+log10(chla)+log10(vol), data=spatial_2)
m34 <- lm(log10(pCO2) ~ log10(real_SW_area_m2)+log10(C_soil)+log10(chla)+log10(vol), data=spatial_2)
m35 <- lm(log10(pCO2) ~ log10(real_SW_area_m2)+log10(C_soil)+log10(DOC)+log10(vol), data=spatial_2)
m36 <- lm(log10(pCO2) ~ log10(real_SW_area_m2)+log10(C_soil)+log10(DOC)+log(chla)+log10(vol), data=spatial_2)
m37 <- lm(log10(pCO2) ~ log10(real_SW_area_m2)+log10(csoil_chla)+log10(DOC)+log(chla)+log10(vol), data=spatial_2)
m38 <- lm(log10(pCO2) ~ log10(real_SW_area_m2)+log10(BR_chla)+log10(DOC)+log(chla)+log10(vol), data=spatial_2)

mod_spatial <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, 
                    m15, m16, m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,
                    m30,m31,m32,m33,m34,m35,m36,m37,m38)

aictab(mod_spatial) #selecting best-fitting models

#Selected models:
summary(m17)
summary(m32)
summary(m15)
summary(m33)

#3.1.2. Seasonal dataset - significant models:
##3.1.2.1. Gargalheiras
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

##3.1.2.2. Cruzetas:
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

##3.1.2.3. Extremoz:
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

#3.1.2.4. Bonfim:
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

#3.2. General linear mixed-models (GLMM) of the analysis of variance for the effects of seasonality,reservoirs and their interaction on pCO2 (Table 3)
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

#3.3. Independent and interactive effects (mean ±1 SE) of seasonality and sampling stations on pCO2 for each ecosystem (Figure 9)
##3.3.1. Gargalheiras
m_g <- lm (log10(pCO2) ~ season + station, data = garg)
summary(m_g)#r2=35.16;p<0.0005
m_g <- lm (log10(pCO2) ~ season, data = garg)
summary(m_g)#r2=16.65;p<0.004
m_g <- lm (log10(pCO2) ~ station, data = garg)
summary(m_g)#r2=17.19;p<0.01

##3.3.2. Cruzeta
m_c <- lm (log10(pCO2) ~ season + station, data = cruz)
summary(m_c)
m_c <- lm (log10(pCO2) ~ season, data = cruz)
summary(m_c)
m_c <- lm (log10(pCO2) ~ station, data = cruz)
summary(m_c)

#3.3.3. Extremoz
m_e <- lm (log10(pCO2) ~ season + station, data = extr)
summary(m_e)#r2=16.1;p=0.05
m_e <- lm (log10(pCO2) ~ season, data = extr)
summary(m_e)#r2=14.73;p<0.05
m_e <- lm (log10(pCO2) ~ station, data = extr)
summary(m_e)

#3.3.4. Bonfim
m_b <- lm (log10(pCO2) ~ season + station, data = bon)
summary(m_b)
m_b <- aov(log10(pCO2) ~ season, data = bon)
summary(m_b)
m_b <- lm (log10(pCO2) ~ station, data = bon)
summary(m_b)

#3.4. Water volume variation in response to rainfall (Figure 9B)
mdvolprecip_garg<-  lm (log10(d_vol+1-min(d_vol)) ~ log10(precip+1-min(precip)), data=garg) 
summary(mdvolprecip_garg)
mdvolprecip_cruz <-lm(log10(d_vol+1-min(d_vol)) ~ log10(precip+1-min(precip)), data=cruz)
summary(mdvolprecip_cruz)
mdvolprecip_extr <- lm(log10(d_vol+1-min(d_vol)) ~ log10(precip+1-min(precip)), data=extr)
summary(mdvolprecip_extr)
mdvolprecip_bon<-  lm (log10(d_vol+1-min(d_vol)) ~ log10(precip+1-min(precip)), data=bon)
summary(mdvolprecip_bon)
