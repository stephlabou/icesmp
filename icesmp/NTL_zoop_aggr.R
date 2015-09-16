##################################################################
####### Aggregate Zoop data sent by email from Corinna Gries on 2015 Sep 15 ##################
##################################################################

setwd("/Users/steve.powers/Desktop/Sandboxes/2015Sep15")

##################################################################
#################### SMP started code on 2015 Sep 15 ###########################
##################################################################

library("reshape2")
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

#assign taxon groups ("zoopgroup") - adds new column with zoopgroup
taxanames$zoopgroup<-"other"
taxanames$zoopgroup[which(taxanames$phylum=="Rotifera")]<-"rotifer"
taxanames$zoopgroup[which(taxanames$sub_class=="Cladocera"|taxanames$order %in% c("Ctenopoda", "Onychopoda", "Haplopoda", "Anomopoda"))]<-"otherclad"
taxanames$zoopgroup[which(taxanames$family=="Daphniidae")]<-"daphnia"
taxanames$zoopgroup[which(taxanames$order_=="Cyclopoida")]<-"cyclopoid"
taxanames$zoopgroup[which(taxanames$order_=="Calanoida")]<-"calanoid"


#unique(data.frame(taxanames$order_,taxanames$f))

#keep only needed columns in taxanames
taxanames <- select(taxanames, ntl_taxon_code, zoopgroup)

#merge zoop data with taxanames
zoop_taxa <- merge(zoop,taxanames,by.x="species_code",by.y="ntl_taxon_code")


#generate lake-date combos for each zoopgroup, used to add proportion=0 when not counted in sample
lakedates<-unique(data.frame(zoop_taxa[,c("lakeid","sample_date")]))
groups<-unique(taxanames$zoopgroup)
for(i in 1:length(groups)){
datai<-lakedates
datai$zoopgroup<-groups[i]  
if(i==1){lakedategroups<-datai}
if(i>1){
lakedategroups<-rbind(lakedategroups,datai)
}
}

#more steps to add zero data
zeroes<-lakedategroups

#create df with all lake/sample/zoopgroup combos - have nubmer per liter if in orig data
zeroes_merge<-merge(zeroes,zoop_taxa,by=c("lakeid","sample_date","zoopgroup"),all=TRUE)


#make NA values for number per liter 0 (these are samples that do not exist)
zeroes_merge$number_per_liter[which(is.na(zeroes_merge$number_per_liter)==TRUE)]<-0

#zeroes_merge<-zeroes_merge[,-which(names(zeroes_merge) %in% c("number_per_liter.x","number_per_liter.y"))]

zoop_all<-zeroes_merge
names(zoop_all)
head(zoop_all)

#remove depth columns and species code

zoop_all <- select(zoop_all, lakeid, sample_date, zoopgroup, number_per_liter)


#calculate zoop total count (for each lake and date)

#zoopsfortot<-subset(zoop,zoop$zoopgroup!="calanoid")

zoopsfortot<-zoop_all

#for a given sample date in a lake, find total zoop number per liter
tots <- zoopsfortot %>% 
  group_by(lakeid,sample_date) %>% 
  dplyr::summarize(tot.number_per_liter = sum(number_per_liter, na.rm = TRUE))

head(tots)
#head(zoop)

tots <- as.data.frame(tots)

#merge zoop data and total count data
merge.df<-merge(zoop_all,tots,by=c("lakeid","sample_date"),all=TRUE)

#calculate zoop proportions (prop of count)
#merge.df$propofcount<-merge.df$number_per_liter/merge.df$sum.number_per_liter

#sum proportions by zoopgroup
#prop.df<-merge.df %>% 
#  group_by(lakeid,sample_date,zoopgroup) %>% 
#  dplyr::summarize(sum.propofcount = sum(propofcount, na.rm = TRUE))

#find total number per liter for each zoop group (within lake, sample date combo)
sum.byzoopgroup.df<-merge.df %>% 
  group_by(lakeid,sample_date,tot.number_per_liter, zoopgroup) %>% 
  dplyr::summarize(sum.number_per_liter = sum(number_per_liter, na.rm = TRUE))

#data.frame(merge.df[,c("lakeid","sample_date","zoopgroup","tot.number_per_liter")])

#head(prop.df)

#merge zoop proportions with season (iceon/iceoff dates)
#merge.all<-merge(prop.df,seasons,by="lakeid")

#each lake will have ALL possible seasons
merge.all<-merge(sum.byzoopgroup.df,seasons,by="lakeid",all=TRUE)

#rename season to year

merge.all <- rename(merge.all, year = season)

#reformat date columns
merge.all$sample_date <- format(as.POSIXct(merge.all$sample_date,format='%m/%d/%Y'),format='%Y-%m-%d')
merge.all$ice_on_start<- format(as.POSIXct(merge.all$ice_on,format='%m/%d/%Y'),format='%Y-%m-%d')
merge.all$ice_on_end<- format(as.POSIXct(merge.all$ice_off,format='%m/%d/%Y'),format='%Y-%m-%d')
merge.all$ice_off_start <- format(as.POSIXct(merge.all$stratifon,format='%m/%d/%Y'),format='%Y-%m-%d')
merge.all$ice_off_end <- format(as.POSIXct(merge.all$stratifoff,format='%m/%d/%Y'),format='%Y-%m-%d')

#clean up df

merge.all <- select(merge.all, lakeid, sample_date, tot.number_per_liter, zoopgroup,
                    sum.number_per_liter, year, ice_on_start, ice_on_end,
                    ice_off_start, ice_off_end)

#add year column
#merge.all$year <- substr(merge.all$sample_date,1,4)


#assign seasons to merged zoop data 

merge.full <- merge.all

merge.full$season<-""

#merge.full$season[which(merge.full$sample_date >= merge.full$ice_on_start & merge.full$sample_date <= merge.full$ice_on_end)]<-"ice_on"

#merge.full$season[which(merge.full$sample_date >= merge.full$ice_off_start & merge.full$year <= merge.full$ice_off_end)]<-"ice_off" 

merge.full <- merge.full %>% 
              mutate(season = ifelse(sample_date >= ice_on_start & sample_date <= ice_on_end, "ice_on", 
                              ifelse(sample_date >= ice_off_start & sample_date <= ice_off_end, "ice_off", NA)))


#discard unneeded rows
merge.final <- filter(merge.full, !is.na(season))

#merge.real<-merge.all[which(merge.all$season %in% c("ice_on","ice_off")),]


#reassign year based on ice off date
#merge.real$year<-substr(merge.real$ice_off,1,4)

zoop_season.df<-merge.final

#keep only rows of interest

zoop_season.df <- select(zoop_season.df, lakeid, year, season, zoopgroup, tot.number_per_liter, sum.number_per_liter)


#calculate average, cv, and max
zoop_nums.df<-zoop_season.df %>% 
  group_by(lakeid, year, season) %>% 
  #avezoopcount
  mutate(avezoopcount = mean(tot.number_per_liter),
  #cvzoopcount
        cvzoopcount = sd(tot.number_per_liter)/avezoopcount,
  #maxzoopcount
        maxzoopcount = max(tot.number_per_liter))

#make as dataframe
zoop_nums.df <- as.data.frame(zoop_nums.df)

#calculate averages for zoopgroups

zoop_props <- zoop_nums.df %>% 
  group_by(lakeid, year, season, zoopgroup) %>% 
  mutate(prop.seasonalavg= sum(sum.number_per_liter)/sum(tot.number_per_liter))
    
#reorganize so props are with zoopgroups

zoop_props <- select(zoop_props, lakeid, year, season, avezoopcount, cvzoopcount, maxzoopcount, prop.seasonalavg, zoopgroup)

zoop_props <- as.data.frame(zoop_props)

#take zoop and prop columns from long format to wide format

zoop_props$ID <- 1:nrow(zoop_props)

zoop_reshape <- select(zoop_props, prop.seasonalavg, zoopgroup, ID)

zcast <- dcast(zoop_reshape, ID ~ zoopgroup, value.var = "prop.seasonalavg")

#merge back with full zoop set

zoop_correct <- merge(zcast, zoop_props, by = "ID")

#get rid of excess columns and consolidate rows

zoop_correct <- select(zoop_correct, lakeid, year, season, avezoopcount, cvzoopcount, maxzoopcount, 
                       calanoid, cyclopoid, daphnia, rotifer, other, otherclad)


zoop_final <- zoop_correct %>% 
              group_by(lakeid, year, season) %>% 
              mutate(propcalanoid = sum(calanoid, na.rm = TRUE),
                     propcyclopoid = sum(cyclopoid, na.rm = TRUE),
                     propdaphnia = sum(daphnia, na.rm = TRUE),
                     propothercladoc = sum(otherclad, na.rm = TRUE),
                     proprotifer = sum(rotifer, na.rm = TRUE),
                     propotherzoop = sum(other, na.rm = TRUE))
              
  
zoop_final <- as.data.frame(zoop_final)
  
zoop_formerge <- select(zoop_final, lakeid, year, season, avezoopcount, cvzoopcount, maxzoopcount,
                       propdaphnia, propothercladoc, propcyclopoid, propcalanoid, proprotifer, propotherzoop)
  
  
zoop_formerge <- unique(zoop_towrite)

zoop_formerge <- zoop_formerge[order(zoop_formerge$lakeid, zoop_formerge$year, rev(zoop_formerge$season)),]


#print to file for checking
write.csv(zoop_formerge, "NTL_correct_zoops.csv", row.names = FALSE)
