#Script for comparative analysis of iceon/iceoff data, by SMP who holds no responsibility
#Notes...
#The script aggregates all data by variable, lake, and season, meaning one ice on value per variable per lake, and once ice off
#-->Lakes and variables lacking iceon data are discarded
#-->Variables with data for fewer than x lakes are discarded
#The graphical outputs are... 
#1) paired iceoff vs. iceon boxplots, by variable (bp.tif)
#2) scatterplots of iceon measurements vs. dissolved n:p ratio, by variable (np.covar.tif)
#Unfinished items as of 2015 Aug 31... 
#-->mixed models (lmer?) using lakename and/or region as random variables, and also year if the full under ice data set is used

#clear all variables in workspace and close all open windows
rm(list = ls()) 
graphics.off()

###################
#load packages
library(reshape);library(ggplot2)
library(zoo)
#library(nlme)
library(lme4)
library(plyr)
library(gridExtra)
###################################
#set directories (adjust to machine)
base_dir<-paste(c("C:/Users/Steve/Desktop/Sandboxes/2015Aug31"), collapse="")
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
#full_under_ice.df$WG<--1
#full_under_ice.df$WG[which(full_under_ice.df$season=="iceoff")]<-1
names(full_under_ice.df)
##############################
#drop Lake Muddus which only has one row of data point (ice off, with no corresponding ice on)
full_under_ice.df<-subset(full_under_ice.df,full_under_ice.df$lakename!="Lake Muddus")
#################################
#lag functions, backward and forward
bw1<-function(x)c(NA,x[1:(length(x)-1)])
bw2<-function(x)c(NA,NA,x[1:(length(x)-2)])
fw1<-function(x)c(x[2:(length(x))],NA)
fw2<-function(x)c(x[3:(length(x))],NA,NA)
######################################
#choose variables
vars<-names(full_under_ice.df)[grep("ave",names(full_under_ice.df))]
#vars2<-names(full_under_ice.df)[grep("max",names(full_under_ice.df))]
vars3<-names(full_under_ice.df)[grep("prop",names(full_under_ice.df))]
#vars4<-names(full_under_ice.df)[grep("cv",names(full_under_ice.df))]
#vars<-unique(c(vars,vars2,vars3,vars4))
vars<-unique(c(vars,vars3))
vars<-vars[-which(vars %in% c("avebenalgalmass","maxbenalgalmass","cvbenalgalmass","avebenchla","maxbenchla","cvbenchla"))]
#vars<-vars[-which(vars %in% c("avebenchla"))]
###########################################
#convert data to long format and restrict to key variables
id_cols<-1:38
long.df<-melt(full_under_ice.df,id.vars=c(names(full_under_ice.df)[id_cols]))
long.df<-long.df[which(is.na(long.df$value)==FALSE),]
long.df<-long.df[-grep("cv",long.df$variable),]
long.df<-long.df[-grep("max",long.df$variable),]

#restrict to variables that were measured during both ice on and ice off at a given lake
varlake.iceon<-unique(long.df[which(long.df$season=="iceon"),which(names(long.df) %in% c("variable","lakename"))])
varlake.iceoff<-unique(long.df[which(long.df$season=="iceoff"),which(names(long.df) %in% c("variable","lakename"))])
varlake<-rbind(varlake.iceon,varlake.iceoff)
varlake.WG<-varlake[which(duplicated(do.call(paste,varlake))),]
long.df<-subset(long.df,do.call(paste,data.frame(long.df$lakename,long.df$variable)) %in% do.call(paste,varlake.WG))

#add variable prefixes for easier sorting in multi-panel plots
whichnames.chem<-c(grep(paste(c("phos","nitro","doc","suva"),collapse="|"),long.df$variable))
whichnames.bio<-c(grep(paste(c("chla","phyt","zoop","prop","ben"),collapse="|"),long.df$variable))
whichnames.phys<-c(grep(paste(c("secchi","temp","radiation"),collapse="|"),long.df$variable))
whichnames.cladocera<-c(grep(paste(c("propothercladoc"),collapse="|"),long.df$variable))
long.df$varname<-as.character(long.df$variable)
long.df$varname[whichnames.chem]<-paste("2chem.",long.df$varname[whichnames.chem],sep="")
long.df$varname[whichnames.bio]<-paste("3bio.",long.df$varname[whichnames.bio],sep="")
long.df$varname[whichnames.phys]<-paste("1phys.",long.df$varname[whichnames.phys],sep="")
long.df$varname[whichnames.cladocera]<-paste("3bio.","propothercladoc",sep="")

#aggregate mean data by variable, lake, season
long.means.df<-aggregate(as.numeric(value)~lakename+season+lakeregloc+lakecountry+lakearea+lakemaxdepth+stationlat+variable+varname,data=long.df,FUN="mean")
names(long.means.df)[which(names(long.means.df) =="as.numeric(value)")]<-"value"

#choose variables that have ice on/ice off measurements for more than 10 lakes 
nlakes<-aggregate(long.means.df$value,by=list(long.means.df$varname,long.means.df$season),FUN="length")
which<-which(nlakes[,3]>10 & !as.character(nlakes[,1]) %in% c("WG","iceonTF"))
use_names<-nlakes[which,1]
use_names

#plot "grand" boxplots
tiff(file=paste("bp.tif",sep=""),height=8.5,width=13,units="in",res=360)
bplot <- ggplot(subset(long.means.df,as.character(long.means.df$varname) %in% use_names), aes(factor(season), value,colour=abs(stationlat)))
bplot<-bplot + geom_boxplot()+geom_point()+geom_jitter(position = position_jitter(height=0,width = 0.1))+ylab("value")
bplot<-bplot+theme(strip.text.x=element_text(size=8.5))+xlab("")+scale_colour_continuous(name="abs(Lat.)")
bplot<-bplot+facet_wrap(~varname,scales="free",ncol=7)
bplot
dev.off()
########################################
#reformat long data so that n:p ratio is a covariate in every row
np.iceon.df<- cast(subset(long.means.df,long.means.df$season=="iceon"),lakename~variable,fun.aggregate="mean")
np.iceon.df$np<-(np.iceon.df$avetotdissnitro/14)/(np.iceon.df$avetotdissphos/31)
np.iceon.df<-data.frame(lakename=np.iceon.df$lakename,np=np.iceon.df$np)
np.iceon.df<-merge(subset(long.means.df,long.means.df$season=="iceon"),np.iceon.df,by=c("lakename"))

#limit iceon data to variables having more than x lakes
nlakes <-aggregate(long.means.df$value,by=list(long.means.df$varname,long.means.df$season),FUN="length")
which<-which(nlakes[,3]>6 & !as.character(nlakes[,1]) %in% c("WG","iceonTF"))
use_names<-nlakes[which,1]

#plot iceon phys, chem, bio data against n:p ratio
tiff(file=paste("np.covar.tif",sep=""),height=8.5,width=11,units="in",res=360)
plot <- ggplot(subset(np.iceon.df,as.character(np.iceon.df$varname) %in% use_names & np.iceon.df$np>0), aes(x=np, value,colour=abs(stationlat)))
plot<-plot +geom_point()+ylab("value")
plot<-plot+theme(strip.text.x=element_text(size=10))+xlab("Dissolved n:p (molar ratio)")+scale_colour_continuous(name="abs(Lat.)")
plot<-plot+facet_wrap(~varname,scales="free_y",ncol=5)
plot
dev.off()
############################################
####################################
#do variable by variable stats
#warning: code below here is very in finished and contains relect sections

#not run if do=0
do<-0
if(do==1){
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
}
