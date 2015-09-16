##################################################################
####### Aggregate Zoop data sent by email from Corinna Gries on 2015 Sep 15 ##################
##################################################################

setwd("/Users/steve.powers/Desktop/Sandboxes/2015Sep15")

##################################################################
#################### SMP started code on 2015 Sep 15 ###########################
##################################################################

library("ncdf")
library("plyr")
library("dplyr")
library("zoo")

## Load data
zoop0 <- read.csv("zoop_density.csv", stringsAsFactors = FALSE)
taxanames0 <- read.csv("zoop_taxa.csv", stringsAsFactors = FALSE)
seasons0 <- read.csv("seasons.csv", stringsAsFactors = FALSE)
zoop<-zoop0
taxanames<-taxanames0
seasons<-seasons0

head(zoop)
head(taxanames)
names(zoop)
names(taxanames)

#assign taxon groups ("zoopgroup")
taxanames$zoopgroup<-"other"
taxanames$zoopgroup[which(taxanames$phylum=="Rotifera")]<-"rotifer"
taxanames$zoopgroup[which(taxanames$family=="Daphniidae")]<-"daphnia"
taxanames$zoopgroup[which(taxanames$order_=="Cyclopoida")]<-"cyclopoid"
taxanames$zoopgroup[which(taxanames$order_=="Calanoida")]<-"calanoid"

unique(data.frame(taxanames$order_,taxanames$family))

#merge zoop data with taxanames
zoop<-merge(zoop,taxanames,by.x="species_code",by.y="ntl_taxon_code")#,all=TRUE)

lakedates<-unique(data.frame(zoop[,c("lakeid","sample_date")]))
groups<-unique(taxanames$zoopgroup)
for(i in 1:length(groups)){
datai<-lakedates
datai$zoopgroup<-groups[i]  
if(i==1){lakedategroups<-datai}
if(i>1){
lakedategroups<-rbind(lakedategroups,datai)
}
}

zeroes<-lakedategroups
zeroes_merge<-merge(zeroes,zoop,by=c("lakeid","sample_date","zoopgroup"),all=TRUE)
zeroes_merge$number_per_liter<-zeroes_merge$number_per_liter
zeroes_merge$number_per_liter[which(is.na(zeroes_merge$number_per_liter)==TRUE)]<-0
#zeroes_merge<-zeroes_merge[,-which(names(zeroes_merge) %in% c("number_per_liter.x","number_per_liter.y"))]
zoop<-zeroes_merge
names(zoop)
head(zoop)

#calculate zoop total count (for each lake and date)
#zoopsfortot<-subset(zoop,zoop$zoopgroup!="calanoid")
zoopsfortot<-zoop
tots <- zoopsfortot %>% 
  group_by(lakeid,sample_date) %>% 
  dplyr::summarize(tot.number_per_liter = sum(number_per_liter, na.rm = TRUE))

head(tots)
head(zoop)

#merge zoop data and total count data
merge.df<-merge(zoop,tots,by=c("lakeid","sample_date"),all=TRUE)

#calculate zoop proportions (prop of count)
#merge.df$propofcount<-merge.df$number_per_liter/merge.df$sum.number_per_liter

#sum proportions by zoopgroup
#prop.df<-merge.df %>% 
#  group_by(lakeid,sample_date,zoopgroup) %>% 
#  dplyr::summarize(sum.propofcount = sum(propofcount, na.rm = TRUE))

sum.byzoopgroup.df<-merge.df %>% 
  group_by(lakeid,sample_date,zoopgroup,tot.number_per_liter) %>% 
  dplyr::summarize(sum.number_per_liter = sum(number_per_liter, na.rm = TRUE))

#data.frame(merge.df[,c("lakeid","sample_date","zoopgroup","tot.number_per_liter")])

#head(prop.df)

#merge zoop proportions with season (iceon/iceoff dates)
#merge.all<-merge(prop.df,seasons,by="lakeid")
merge.all<-merge(sum.byzoopgroup.df,seasons,by="lakeid",all=TRUE)
#reformat date columns
merge.all$sample_date <- format(as.POSIXct(merge.all$sample_date,format='%m/%d/%Y'),format='%Y-%m-%d')
merge.all$ice_on<- format(as.POSIXct(merge.all$ice_on,format='%m/%d/%Y'),format='%Y-%m-%d')
merge.all$ice_off<- format(as.POSIXct(merge.all$ice_off,format='%m/%d/%Y'),format='%Y-%m-%d')

#add year column
merge.all$year <- substr(merge.all$sample_date,1,4)

#assign seasons to merged zoop data 
merge.all$season<-""
merge.all$season[which(merge.all$sample_date> merge.all$ice_on & merge.all$sample_date<merge.all$ice_off)]<-"ice_on"
merge.all$season[which(merge.all$sample_date> merge.all$ice_off & merge.all$year==substr(merge.all$ice_off,1,4))]<-"ice_off" 

#discard unneeded rows
merge.real<-merge.all[which(merge.all$season %in% c("ice_on","ice_off")),]
#reassign year based on ice off date
merge.real$year<-substr(merge.real$ice_off,1,4)
zoop_season.df<-merge.real

#calculate seasonal averages
avg_zoop_season.df<-zoop_season.df %>% 
  group_by(lakeid,year,season,zoopgroup) %>% 
  dplyr::summarize(prop.seasonalavg= mean(sum.number_per_liter/tot.number_per_liter),
  tot.seasonalavg= mean(tot.number_per_liter))

#print to file for checking
write.csv(avg_zoop_season.df,"check.csv")
