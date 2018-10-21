##Manucript: CO2 Supersaturation in Tropical Eutrophic Waters
###Script used to build up the figures used in the main text and the supplementary material throughout the manuscript
### Authors: Pedro C. Junger et al.

###### 1.1 - loading required libraries #####
library(ggplot2)
library(gtable)
library(scales)
library(gridExtra)
library(dplyr)
library(tidyr)
library(ggrepel)
library(grid)
library(ggpubr)
library(gtools)

###### 1.2 - External Functions #####
##plain: scientific notation for graphics
plain <- function(x,...) {
  format(x, ..., scientific = FALSE, trim = TRUE)
}

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

###### 2.1 - Data input ######

spatial <- read.csv("spatial_dataset_2.csv",header = TRUE) #full spatial dataset from 102 ecosystems used in this study
seasonal <- read.csv("seasonal_dataset.csv", sep = ",",header = TRUE) #full temporal dataset used in this study
mean <- read.csv("mean_data.csv", sep=",",header = TRUE) #dataset with mean values used in figures in supplementary material

spatial$DOC_chla<- spatial$DOC / spatial$chla*0.04 #creating column with doc:chla ratio
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

###### 2.2 - Filtering the spatial dataset by excluding observations with Total Alkalinity < 1meq/L ######
#Obs.: See detailed reasons in the methods (section 2.5) of the manuscript
spatial<- spatial %>%
  filter(lakeID!=67) #excluding outlier "Riacho da Cruz"
spatial_2 <- spatial %>%
  filter(TA>1) #excluding superestimated pCO2 values from the full spatial dataset
seasonal_2 <- seasonal %>%
  filter(TA>1) #excluding superestimated pCO2 values from the full seasonal dataset

dochla_data<-subset(spatial_2, DOC_chla<10) #excluding DOC:chla outlier
chla_data<-subset(spatial_2, chla>3.08) #excluding Chl-a outlier from the regression
PR_data<-subset(spatial_2, PR>7) #excluding PR outlier from the regression

###### 2.4 - Organizing dataset ######
seasonal$system <- factor(seasonal$system, levels = c("Gargalheiras", "Cruzeta","Extremoz","Bonfim"))
seasonal_2$system <- factor(seasonal_2$system, levels = c("Gargalheiras","Extremoz","Cruzeta","Bonfim"))
seasonal$date <- factor(seasonal$date, levels =c ("10-Dec","11-Jan","11-Feb","11-Mar","11-Apr","11-May",
                                                  "11-Jun","11-Jul","11-Aug","11-Sep","11-Oct","11-Nov",
                                                  "11-Dec","12-Jan","12-Jul","12-Aug","12-Sep","12-Oct","12-Nov","12-Dec",
                                                  "13-Jan","13-Feb","13-Mar","13-Apr","13-May","13-Jun",
                                                  "13-Jul","14-Sep","14-Oct","14-Nov","14-Dec","15-Jan","15-Feb",
                                                  "15-Mar","15-Apr","15-May","15-Jun","15-Jul","15-Aug",
                                                  "15-Sep","15-Oct"))

mean$system <- factor(mean$system, levels = c("Gargalheiras", "Cruzeta","Extremoz","Bonfim"))
###### Separating seasonal dataset by system ####
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

##### 3.1- Figures #####
#Figure 1: Map made with ArcGis

#Figure 2: Violin plots showing spatial and seasonal variability of the main variables
#Fig 2A
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

#Fig 2B
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

#Figure 3: Bargraph: prevalance of heterotrophy over space and time
spatial$trophic_state <- factor(spatial$trophic_state, levels=c("Hypereutrophic","Eutrophic","Mesotrophic","Oligotrophic"))
mycolors<-c("green4", "darkolivegreen3", "lightsteelblue3","darkblue")
micro <- expression(bold("390 "*mu*"atm"))

ggplot(spatial,aes(co2_rel, fill=trophic_state))+
  geom_histogram(binwidth = 4,alpha = 0.5,color="black")+
  geom_vline(xintercept = 0, linetype="dashed")+
  scale_y_continuous(expand=c(0,0), limits = c(0,18), breaks = seq(0,15,by=3))+
  scale_x_continuous(breaks = seq(-20,120,by=20))+
  scale_fill_manual(values=mycolors)+
  annotate("text",x=20,y=17,label=as.character(micro),size=5,color="red", parse=TRUE)+
  geom_label(aes(x = 0, y = 17, label = as.character(micro)), fill = "white",parse=TRUE)+
  annotate("segment", x = 0, xend = 10, y = 17, yend = 17, colour = "red", size=1, arrow=arrow())+
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

#Figure 4: Scatter plots: regressions with significant internal variables (spatial and seasonal)
## CO2 x chla relationship: eutrophication effect on pCO2
#4A: spatial
m <-  lm (log10(pCO2) ~ log10(chla), data=dochla_data)
summary(m)
A4<-ggplot(dochla_data, aes(chla, pCO2))+
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

#4B: Seasonal
m <-  lm (log10(pCO2) ~ log10(chla), data=seasonal_2)#Bonfim is overestimated, thus it was excluded from the dataset for this analysis
summary(m)
B4<-ggplot(seasonal_2, aes(chla,pCO2))+
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

#4C) pCO2 vs Vol
m <-lm (log10(pCO2) ~ log10(vol), data=spatial_2)
summary(m)
C4<-ggplot(spatial_2, aes(vol, pCO2))+
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

#4D) Chla vs vol

m <-  lm (log10(chla) ~ log10(vol), data=dochla_data)
summary(m)
D4<-ggplot(dochla_data, aes(vol, chla))+
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

plots1 <- list(A4,C4)
plots2 <- list(B4,D4)
grobs1 = lapply(plots1, ggplotGrob)
grobs2 = lapply(plots2, ggplotGrob)
g1 = do.call(cbind, c(grobs1, size="first"))
g2 = do.call(cbind, c(grobs2, size="first"))
g3 = do.call(rbind, c(list(g1,g2), size="first"))
grid.newpage()
grid.draw(g3) #fine adjustments were made using Adobe Illustrator

#Figure 5: Internal processes 2: pCO2 vs DOC, pCO2 vs DOC/chla, DOC vs Chla, a365 vs DOC
#5A) pCO2 x DOC
m <-  lm (log10(pCO2) ~ log10(DOC), data=spatial_2)
summary(m)
A5<-ggplot(spatial_2, aes(DOC, pCO2))+
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

#5B
m <-  lm (log10(DOC) ~ log10(chla), data=spatial_2)
summary(m)
B5<-ggplot(spatial_2, aes(chla,DOC))+
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

#5C) pCO2 vs BR:Chla
m <-  lm(log10(pCO2) ~ log10(BR_chla), data=PR_data)
summary(m)
C5<-ggplot(PR_data, aes(BR_chla,pCO2))+
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

#5D) DOC x a365 (CDOM):strong positive relationship indicates high aromatic character of DOM
m<-  lm(log10(a365) ~ log10(DOC), data=spatial_2)
summary(m)
D5<-ggplot(spatial_2, aes(DOC,a365))+
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

plots1 <- list(A5,B5)
plots2 <- list(C5,D5)
grobs1 = lapply(plots1, ggplotGrob)
grobs2 = lapply(plots2, ggplotGrob)
g1 = do.call(cbind, c(grobs1, size="first"))
g2 = do.call(cbind, c(grobs2, size="first"))
g3 = do.call(rbind, c(list(g1,g2), size="first"))
grid.newpage()
grid.draw(g3) #fine adjustments were made using Adobe Illustrator

#Figure 6: pCO2 x DO
m <- lm(log10(pCO2) ~ log10(Dosat), data=spatial_2[-20,]) #r2=-25.4; p=0.0003
summary(m)
A6<-ggplot(spatial_2[-20,], aes(Dosat, pCO2))+ #outlier excluded(lakeID=61;DO=3.2mg/L)
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

#6B)
m<- lm(log10(pCO2) ~ log10(DO_sat), data=seasonal_2) 
summary(m)
B6<-ggplot(seasonal_2, aes(DO_sat,pCO2))+
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

#6C
m <- lm(log10(pCO2) ~ log10(DO_sat), data=garg) 
summary(m)
C6<-ggplot(garg, aes(DO_sat,pCO2))+
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

#6D
m <- lm(log10(pCO2) ~ log10(DO_sat), data=cruz) 
summary(m)
D6<-ggplot(cruz, aes(DO_sat,pCO2))+
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

ggarrange(A6,B6,C6,D6,ncol=2,nrow=2,widths = 1,heights = 1) #fine adjustments were made using Adobe Illustrator

#Figure 7: External processes: scatter plots of regressions
##### Morphometric measurements 
#7A) pCO2 vs CA
m<-  lm(log10(pCO2) ~ log10(CA_m2), data=spatial_2)
summary(m)
A7<-ggplot(spatial_2, aes(CA_m2, pCO2))+
  stat_smooth(method = "lm",alpha = 0.5, fill="lightgrey", size = 1,color="black") +
  geom_jitter()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10,10^5))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10^5,10^9.5))+
  ggtitle("A")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= expression("Catchment area (m²)"), y = expression('pCO'[2]*' ('*mu*'atm)'))+
  annotate("text",x=10^8.5,y=10^5,label=lm_eqn(m),parse=TRUE)

#7B)
m <-  lm(log10(pCO2) ~ log10(csoil_chla), data=spatial_2)
summary(m)
B7<-ggplot(spatial_2, aes(csoil_chla, pCO2))+
  stat_smooth(method = "lm",alpha = 0.5, fill="lightgrey", size = 1,color="black") +
  geom_jitter()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(10,10^5))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(0.005/0.04,1.25/0.04))+
  ggtitle("B")+
  theme(plot.title = element_text(size=18,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black", size=0.2),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title=element_text(size=16)) +
  labs(x= "Soil OC : Chl-a ratio", y = "")+
  annotate("text",x=10^1,y=10^5,label=lm_eqn(m),parse=TRUE)

g1 <- ggplotGrob(A7)
g2 <- ggplotGrob(B7)
g <-cbind(g1, g2, size = "first")
g$heights <- unit.pmax(g1$heights,g2$heights)
grid.newpage()
grid.draw(g)

#Figure 9: Sampling station and Seasonality effect on pCO2 using seasonal dataset
#9A) Boxplots with season and stations
A9.1<-ggplot(garg, aes(station,pCO2))+
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

B9.1<-ggplot(cruz, aes(station,pCO2))+
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

C9.1<-ggplot(extr, aes(station,pCO2))+
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

D9.1<-ggplot(bon, aes(station,pCO2))+
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

A8<-ggarrange(A8.1,C8.1,B8.1,D8.1,ncol=2,nrow=2,widths = 1,heights = 1,common.legend=TRUE,legend="right")

#9B) volume responds to rainfall in the reservoirs, but not in the coastal lagoon
#Delta volume x Precipitation
B9<-ggplot(seasonal, aes(log10(precip+1-min(precip)),log10(d_vol+1-min(d_vol))))+
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

#Obs.:The figures 9A and 9B were gathered using Adobe Illustrator

#Figure 8: Land use impact on variables of interest
##8A) Chla x Forest
m <- lm(log10(chla) ~ fo_logit, chla_data)
summary(m)
A8<-ggplot(chla_data, aes(fo_logit,chla))+
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

##8B) Chla x Anthropogenic
m <- lm(log10(chla) ~ anth_logit, chla_data) #R2=26.9; p<0.001
summary(m)

B8<-ggplot(chla_data, aes(anth_logit,chla))+
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

##8C) DOC:Chla x Forest
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

##8D) DOC:Chla x Anthropogenic
m <- lm(log10(DOC_chla) ~ anth_logit, chla_data)
summary(m)
D8<-ggplot(chla_data, aes(anth_logit,DOC_chla))+
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

##8E) SOC:Chla x Forest
m <- lm(log10(csoil_chla) ~ fo_logit, chla_data)
summary(m)
E8<-ggplot(chla_data, aes(fo_logit,csoil_chla))+
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

##8F) SOC:Chla x Anthropogenic
m <- lm(log10(csoil_chla) ~ anth_logit, chla_data)
summary(m)

F8<-ggplot(chla_data, aes(anth_logit,csoil_chla))+
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

##8G) BR:Chla x Forest
#Forest
m <- lm(log10(BR_chla) ~ fo_logit, chla_data) #+R2= 12.31, p<0.05
summary(m)
G8<-ggplot(chla_data, aes(fo_logit,BR_chla))+ 
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

##8H) BR:Chla x Anthropogenic
m <- lm(log10(BR_chla) ~ anth_logit, chla_data) #R2=-12.5, p<0.05
summary(m)
H8<-ggplot(chla_data, aes(anth_logit,BR_chla))+
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

grid.arrange(A9,B9,C9,D9,E9,F9,G9,H9, ncol=2) # fine adjustments were made using Adobe Illustrator

####3.2 - Figures in Suplementary Material####

#Fig S1: Seasonal trends in pCO2, volume, precipitation
## pCO2, chla, rainfall and volume over time in all systems
#precipitação
A1s <- ggplot(mean, aes(date,precip))+
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
B1s <- ggplot(mean, aes(date,vol/1000))+
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

C1s <- ggplot(mean, aes(date,pCO2))+
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

D1s <- ggplot(mean, aes(date,chla))+
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

grid.arrange(A3s,B3s,C3s,D3s, nrow=2)

#Figure S2: density plots with two datasets for comparison
A2s <- ggplot(spatial, aes(TA,pCO2))+
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

B2s <- ggplot(spatial, aes(x=pCO2))+
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

C2s <- ggplot(seasonal, aes(TA,pCO2))+
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

D2s <- ggplot(seasonal, aes(x=pCO2))+
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

grid.arrange(A2s,B2s,C2s,D2s,ncol=2)

#Figure S3:Frequency of pCO2 values based on sampling day period (morning or afternoon).
## This figure was made with Graph Prism

#Figure S4: Box-plots comparing pCO2 values between (A) reservoirs and lakes and between (B) humid and semiarid systems 
##4Sa
t.test(spatial_2$pCO2~spatial_2$type) #student t-test
A4S<-ggplot(spatial_2, aes(type,pCO2))+
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

##4Sb
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

grid.arrange(s8a, s8b,ncol=2)

#Figure S5: Volume X Lake Area x pCO2
#S5A: Lake area vs Volume
m <-  lm(log10(vol) ~ log10(lake_area_m2), data=spatial_2)
summary(m)
S5A<-ggplot(spatial_2, aes(lake_area_m2, vol))+
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

#S5B)Lake area
m <-  lm(log10(pCO2) ~ log10(lake_area_m2), data=spatial_2)
summary(m)

S5B<-ggplot(spatial_2, aes(lake_area_m2, pCO2))+
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

grid.arrange(S5A, S5B, ncol=1)

#Fig6S: pCO2 vs land use
#A) Forested
m<- lm(log10(pCO2) ~ fo_logit, spatial_2) #tendencia, mas nada significativo
summary(m)
S6A<-ggplot(spatial_2, aes(fo_logit,pCO2))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter(size=2,alpha=0.7)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(14,10^5))+
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
  labs(x= "Logit (Proportion of Forested area)", y = expression('pCO'[2]*' ('*mu*'atm)'))

#B)Anthropogenic
m<- lm(log10(pCO2) ~ anth_logit, spatial_2)
summary(m)
S6B<-ggplot(spatial_2, aes(anth_logit,pCO2))+
  stat_smooth(method = "lm",alpha = 0.5, fill= "lightgrey",size = 1,color="black") +
  geom_jitter(size=2,alpha=0.7)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits =c(14,10^5))+
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
  labs(x= "Logit (Proportion of Anthropogenic area)", y = "")

grid.arrange(S6A, S6B, ncol=2)
#Fig7S: Sediment influence: Surface pCO2 x Bottom pCO2 in Gargalheiras and Cruzeta
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

#Figure S8: diel variation in Gargalheiras
## This figure was made with STATISTICA 9.0

