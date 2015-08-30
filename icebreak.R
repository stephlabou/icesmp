
#clear all variables in workspace and close all open windows
rm(list = ls()) 
graphics.off()

library(reshape);library(ggplot2)
library(zoo)
#library(nlme)
library(lme4)
#library(dplyr)
library(plyr)
#library(data.table)
library(gridExtra)

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

varlake.iceon<-long.df[which(long.df$season=="iceon"),which(names(long.df) %in% c("variable","lakename"))]
varlake.iceoff<-long.df[which(long.df$season=="iceoff"),which(names(long.df) %in% c("variable","lakename"))]
varlake.merge<-merge(varlake.iceoff,varlake.iceon,by=c("variable","lakename"))


whichnames.chem<-c(grep(paste(c("phos","nitro","doc","suva"),collapse="|"),long.df$variable))
whichnames.bio<-c(grep(paste(c("chla","phyt","zoop","prop","ben"),collapse="|"),long.df$variable))
whichnames.phys<-c(grep(paste(c("secchi","temp","radiation"),collapse="|"),long.df$variable))
whichnames.cladocera<-c(grep(paste(c("propothercladoc"),collapse="|"),long.df$variable))

long.df$varname<-as.character(long.df$variable)
long.df$varname[whichnames.chem]<-paste("2chem.",long.df$varname[whichnames.chem],sep="")
long.df$varname[whichnames.bio]<-paste("3bio.",long.df$varname[whichnames.bio],sep="")
long.df$varname[whichnames.phys]<-paste("1phys.",long.df$varname[whichnames.phys],sep="")
long.df$varname[whichnames.cladocera]<-paste("3bio.","propothercladoc",sep="")

long.means.df<-aggregate(as.numeric(value)~lakename+season+lakeregloc+lakecountry+lakearea+lakemaxdepth+stationlat+variable+varname,data=long.df,FUN="mean")
names(long.means.df)[which(names(long.means.df) =="as.numeric(value)")]<-"value"

nlakes<-aggregate(long.means.df$value,by=list(long.means.df$varname,long.means.df$season),FUN="length")
which<-which(nlakes[,3]>10 & !as.character(nlakes[,1]) %in% c("WG","iceonTF"))
use_names<-nlakes[which,1]
use_names

#use_names<-"avechla"
#test<-subset(long.means.df,as.character(long.means.df$variable) %in% use_names)
#test

windows()
bplot <- ggplot(subset(long.means.df,as.character(long.means.df$varname) %in% use_names), aes(factor(season), value,colour=abs(stationlat)))
bplot<-bplot + geom_boxplot()+geom_point()+geom_jitter(position = position_jitter(height=0,width = 0.1))+ylab("value")
bplot<-bplot+theme(strip.text.x=element_text(size=9))+xlab("")+scale_colour_continuous(name="Lat(abs)")
#+theme(legend.position="none")
bplot<-bplot+facet_wrap(~varname,scales="free",ncol=7)
bplot

########################################
np.iceon.df<- cast(subset(long.means.df,long.means.df$season=="iceon"),lakename~variable,fun.aggregate="mean")
np.iceon.df$np<-(np.iceon.df$avetotdissnitro/14)/(np.iceon.df$avetotdissphos/31)
np.iceon.df<-data.frame(lakename=np.iceon.df$lakename,np=np.iceon.df$np)

np.iceon.df<-merge(subset(long.means.df,long.means.df$season=="iceon"),np.iceon.df,by=c("lakename"))

nlakes<-aggregate(long.means.df$value,by=list(long.means.df$varname,long.means.df$season),FUN="length")
which<-which(nlakes[,3]>6 & !as.character(nlakes[,1]) %in% c("WG","iceonTF"))
use_names<-nlakes[which,1]
#use_names<-use_names[grep("bio.",use_names)]

windows()
plot <- ggplot(subset(np.iceon.df,as.character(np.iceon.df$varname) %in% use_names & np.iceon.df$np>0), aes(x=np, value,colour=abs(stationlat)))
plot<-plot +geom_point()+ylab("value")
plot<-plot+theme(strip.text.x=element_text(size=9))+xlab("")+scale_colour_continuous(name="Lat(abs)")
#+theme(legend.position="none")
plot<-plot+facet_wrap(~varname,scales="free_y",ncol=5)
plot

windows()
plot <- ggplot(subset(np.iceon.df,as.character(np.iceon.df$varname) %in% use_names), aes(x=lakemaxdepth, value,colour=abs(stationlat)))
plot<-plot +geom_point()+ylab("value")
plot<-plot+theme(strip.text.x=element_text(size=9))+xlab("")+scale_colour_continuous(name="Lat(abs)")
#+theme(legend.position="none")
plot<-plot+facet_wrap(~varname,scales="free",ncol=7)
plot


############################################

vars_use<-c()
for(i in 1:length(vars)){
vari<-vars[i]
groupi<-full_under_ice.df$season
yi<-full_under_ice.df[,which(names(full_under_ice.df) %in% vari)]
t.testi<-t.test(yi~groupi)
p.vali<-t.testi$p.value
if(p.vali<=0.05){vars_use<-c(vars_use,vari)}
}

sigvars<-c()
for(i in 1:length(vars)){
varsi<-vars[i]
coli<-which(names(full_under_ice.df)==varsi)
which_na<-which(is.na(full_under_ice.df[,coli])==TRUE)
dataWGi<-full_under_ice.df[-which_na,]
dataWGi$val<-dataWGi[,coli]
#dataWGi$value<-dataWGi[,coli]
#dataWGi$groups<-dataWGi$lakename

if(length(unique(dataWGi$lakename))<=10){next()}

#order data frames
dataWGi<-dataWGi[order(dataWGi$lakename,dataWGi$year,dataWGi$iceonTF),]

dataWGnlimi<-subset(dataWGi,dataWGi$lakename %in% lakenames_nlim)

#dataWGnlimi$bw1<-unlist(tapply(dataWGnlimi[,coli],dataWGnlimi$lakename,bw1))
#dataWGnlimi$bw2<-unlist(tapply(dataWGnlimi[,coli],dataWGnlimi$lakename,bw2))
#dataWGnlimi$fw1<-unlist(tapply(dataWGnlimi[,coli],dataWGnlimi$lakename,fw1))
#dataWGnlimi$fw2<-unlist(tapply(dataWGnlimi[,coli],dataWGnlimi$lakename,fw2))


bivari <- reshape(dataWGi, v.names = "val",idvar = c("lakename","year"),
                timevar = "season", direction = "wide")

longi<-reshape(bivari, idvar = c("lakename","year"), timevar = "season",
        v.names = "val", direction = "long")

bivariplot<-bivari
bivariplot$val.iceon[which(bivariplot$val.iceon==0)]<-
0.5*min(na.omit(bivariplot$val.iceon[which(bivariplot$val.iceon!=0)]))
bivariplot$val.iceoff[which(bivariplot$val.iceoff==0)]<-
0.5*min(na.omit(bivariplot$val.iceoff[which(bivariplot$val.iceoff!=0)]))

longiplot<-longi
longiplot$val[which(longiplot$val==0)]<-
0.5*min(na.omit(longiplot$val[which(longiplot$val!=0)]))


#windows()
bplot <- ggplot(longiplot, aes(factor(season), log10(val),colour=log10(lakemaxdepth)))
bplot<-bplot + geom_boxplot()+geom_point()+geom_jitter(position = position_jitter(height=0,width = 0.1))+ylab(paste("log10",varsi))+theme(legend.position="none")

xylim<-log10(range(c(na.omit(bivariplot$val.iceoff),na.omit(bivariplot$val.iceon))))
plot <- ggplot(bivariplot, aes(x=log10(val.iceoff), y=log10(val.iceon),colour=log10(lakemaxdepth)))
#plot <- ggplot(bivari, aes(x=log10(val.iceoff), y=log10(val.iceon),colour=log10(lakearea)))
plot<-plot + geom_point()+geom_abline(slope=1)+xlim(xylim)+ylim(xylim)#ylim(0,max(val))#+coord_fixed(ratio=5)
#print(bplot)
#windows(width=11,height=4)
tiff(file=paste(varsi,".bplot.tif",sep=""),height=4,width=9,units="in",res=180)
grid.arrange(bplot,plot,ncol=2)
dev.off()

#}

#windows()
#plot <- ggplot(dataWGnlimi, aes(x=bw1, val,colour=season))
#plot<-plot + geom_point()+ylab(varsi)+facet_wrap(~season)
#print(plot)

#windows()
#plot <- ggplot(dataWGnlimi, aes(x=fw1, val,colour=season))
#plot<-plot + geom_point()+ylab(varsi)+facet_wrap(~season)
#print(plot)

#windows()

tiff(file=paste(varsi,".acf.tif",sep=""),height=6,width=4,units="in",res=180)
par(mfrow=c(3,2))
acf(dataWGnlimi$val[which(dataWGnlimi$iceonTF==1)],main="iceon, acf", xlab="Lag (years)")
pacf(dataWGnlimi$val[which(dataWGnlimi$iceonTF==1)],main="iceon, pacf", xlab="Lag (years)")
acf(dataWGnlimi$val[which(dataWGnlimi$iceonTF==0)],main="iceoff, acf", xlab="Lag (years)")
pacf(dataWGnlimi$val[which(dataWGnlimi$iceonTF==0)],main="iceoff, pacf", xlab="Lag (years)")
acf(dataWGnlimi$val,main="iceon+iceoff, acf", xlab="Lag (iceon/iceoff seasons)")
pacf(dataWGnlimi$val,main="iceon+iceoff, pacf", xlab="Lag (iceon/iceoff seasons)")
mtext(side=3,outer=TRUE,text=varsi,padj=1.5)
dev.off()

modeli<- lmer(log10(full_under_ice.df[,coli]+.001) ~ iceonTF + (sampledepth | lakename), full_under_ice.df)

print(varsi)
print(summary(modeli))
print(anova(modeli))

#

}




