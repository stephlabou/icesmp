
#clear all variables in workspace and close all open windows
rm(list = ls()) 
graphics.off()

library(reshape2);library(ggplot2)
library(zoo)
library(lme4)
library(plyr)
library(gridExtra)
library(reshape2)

###################################
#set directories
base_dir<-paste(c("C:/Users/Steve/Desktop/Sandboxes/2015Aug27"), collapse="")
#data_dir<-paste(c(base_dir,"/Data"),collapse="")
data_dir<-base_dir
setwd(base_dir)

###############################
#read in data
setwd(data_dir)
full_under_ice.df<-read.csv("full_under_ice_data_20150613.csv")
full_under_ice.df
full_under_ice.df$iceonTF<-0
full_under_ice.df$iceonTF[which(full_under_ice.df$season=="iceon")]<-1
full_under_ice.df$WG<--1
full_under_ice.df$WG[which(full_under_ice.df$season=="iceoff")]<-1
names(full_under_ice.df)
###########################3
aggr<-aggregate(full_under_ice.df[,1],by=list(full_under_ice.df$lakename),FUN="length")
aggr
lakenames_nlim<-aggr[,1][which(aggr[,2]>=5)]
lakenames_nlim<-aggr[,1][which(aggr[,2]>=6)]

#drop Lake Muddus which only has one row of data point (ice off, with no corresponding ice on)
full_under_ice.df<-subset(full_under_ice.df,full_under_ice.df$lakename!="Lake Muddus")
nlim_under_ice.df<-subset(full_under_ice.df,full_under_ice.df$lakename %in% as.character(lakenames_nlim))

#################################
#lag function
bw1<-function(x)c(NA,x[1:(length(x)-1)])
bw2<-function(x)c(NA,NA,x[1:(length(x)-2)])
fw1<-function(x)c(x[2:(length(x))],NA)
fw2<-function(x)c(x[3:(length(x))],NA,NA)

######################################
lakenames<-unique(full_under_ice.df$lakename)
#vars<-c("avechla","maxchla")
vars<-names(full_under_ice.df)[grep("ave",names(full_under_ice.df))]
vars2<-names(full_under_ice.df)[grep("max",names(full_under_ice.df))]
vars3<-names(full_under_ice.df)[grep("prop",names(full_under_ice.df))]
vars4<-names(full_under_ice.df)[grep("cv",names(full_under_ice.df))]
vars<-unique(c(vars,vars2,vars3,vars4))
vars<-vars[-which(vars %in% c("avebenalgalmass","maxbenalgalmass","cvbenalgalmass","avebenchla","maxbenchla","cvbenchla"))]
#vars<-vars[-which(vars %in% c("avebenchla"))]
###########################################

id_cols<-1:38
long.df<-melt(full_under_ice.df,id.vars=c(names(full_under_ice.df)[id_cols]))
#long.means.df<-ddply(long.df,.(lakename,season,lakemaxdepth,variable,value),summarise,mean(value,na.rm=TRUE))
long.df<-long.df[which(is.na(long.df$value)==FALSE),]
long.df<-long.df[-grep("cv",long.df$variable),]
long.df<-long.df[-grep("max",long.df$variable),]

varlake.iceon<-unique(long.df[which(long.df$season=="iceon"),which(names(long.df) %in% c("variable","lakename"))])
varlake.iceoff<-unique(long.df[which(long.df$season=="iceoff"),which(names(long.df) %in% c("variable","lakename"))])
varlake<-rbind(varlake.iceon,varlake.iceoff)
varlake.WG<-varlake[which(duplicated(do.call(paste,varlake))),]

long.df<-subset(long.df,do.call(paste,data.frame(long.df$lakename,long.df$variable)) %in% do.call(paste,varlake.WG))

unique(long.df$lakename)

chla.Super<-subset(long.df,long.df$variable=="avechla" & long.df$lakename=="Lake Superior" & long.df$stationname=="ThunderBay")
#chla.Trout<-subset(long.df,long.df$variable=="avechla" & long.df$lakename=="Trout Lake")
#chla.Crystal<-subset(long.df,long.df$variable=="avechla" & long.df$lakename=="Crystal Lake")
#chla.CrystalBog<-subset(long.df,long.df$variable=="avechla" & long.df$lakename=="Crystal Bog")
data<-chla.Super
#data<-chla.Trout
#data<-chla.Crystal
#data<-chla.CrystalBog

data.paired<-reshape(data, v.names="value",timevar = "season", idvar = c("lakename","year"),direction = "wide")
data2<-data.frame(year=data.paired$year,iceon=data.paired$value.iceon,iceoff=data.paired$value.iceoff)

#Use the Granger causality test to see what follows what (iceoff follows iceon vs. iceon follows iceoff)
#see http://www.r-bloggers.com/chicken-or-the-egg-granger-causality-for-the-masses/

#get the differenced data
diceon1<-diff(as.numeric(data2$iceon))
diceoff2<-diff(as.numeric(data2$iceoff))

#test the hypothesis that ice off follows iceon
grangertest(diceoff2~diceon1,order=2)

#SIGNIFICANT for Lake Superior Chl a (H: ICEOFF CONDITION FOLLOWS ICEON CONDITION)
#Granger causality test

#Model 1: diceoff2 ~ Lags(diceoff2, 1:2) + Lags(diceon1, 1:2)
#Model 2: diceoff2 ~ Lags(diceoff2, 1:2)
#  Res.Df Df      F  Pr(>F)  
#1     21                    
#2     23 -2 4.5685 0.02253 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#now do the reverse test, that iceon follows iceoff...
#requires shifting iceoff column by one year
data3<-data.frame(year=data.paired$year,iceon=data.paired$value.iceon,iceoff=c(NA,data.paired$value.iceoff[1:(-1+length(data.paired[,2]))]))

diceoff1<-diff(as.numeric(data3$iceoff))
diceon2<-diff(as.numeric(data3$iceon))

grangertest(diceon2~diceoff1,order=2)

#OTHER WAY AROUND (H:ICEON CONDITION FOLLOWS ICEOFF CONDITION) IS NOT SIGNIF for Lake Superior Chl a
#Granger causality test
#
#Model 1: diceon2 ~ Lags(diceon2, 1:2) + Lags(diceoff1, 1:2)
#Model 2: diceon2 ~ Lags(diceon2, 1:2)
#  Res.Df Df      F Pr(>F)
#1     20                 
#2     22 -2 0.4132 0.6671



