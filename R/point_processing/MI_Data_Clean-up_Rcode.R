
#The following is code Jody has used for southern Michigan and northern Michigan to do
#quality control checks to find things that need to be corrected.

#########################################################
###### Northern Michigan - check for duplicates #########
#########################################################

setwd("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected")
rm(list=ls())

options(digits = 12)
MI = read.csv("C:/Users/jmurray7/Dropbox/GIS PalEON/UMW_Point Data from Goring 2016 PlosOne/raw UMW posted to wiki/northernmi_pls_projected_v1.2.csv", header = TRUE, stringsAsFactors = FALSE)
head(MI)

#create uniqueIDs by concatenating all tree 1 and 2 information into one vector
uniquetrees = paste0(MI$SPP1,MI$DBH1,MI$AZIMUTH,MI$DIST1,MI$SPP2,MI$DBH2,MI$AZIMUTH2,MI$DIST2)
MI$uniquetrees = uniquetrees


######################## FIND DUPLICATES ##############################

# find the duplicate uniqueIDs
dups <-MI[duplicated(MI$uniquetrees)|duplicated(MI$uniquetrees, fromLast=TRUE),]
dups
head(dups)
write.csv(dups, file = "northernMI_Duplicates.csv", row.names = FALSE)



# Selecting just one of the duplicates from the list of Duplicate Corners
# that were checked in GIS
duplicates = read.csv("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/large tree corrections from Charlie 8-18-18/northernMichigan_duplicates/northernMI_Duplicates_checked.csv", header=T, stringsAsFactors = F)

dup_sub <- subset(duplicates, keep_1 %in% c(0,2))
head(dup_sub)
unique(dup_sub$keep_1)

unique_id <- unique(duplicates$uniquetree)
library(plyr)
# Use which() to figure out which ones aren't 2
duplicate.counts = ddply(dup_sub, "uniquetree", summarize, n = length(uniquetree))
write.csv(duplicate.counts, "C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/large tree corrections from Charlie 8-18-18/northernMichigan_duplicates/Duplicates_Counts.csv", row.names = F)
which(duplicate.counts$n >2)

# Create storage
mins <- rep(NA, length(unique_id))
maxs <- rep(NA, length(unique_id))
sub_data <- rep(NA,length(unique_id))

for (i in 1:length(unique_id)){
  sub_data <- subset(dup_sub, uniquetree == unique_id[i])
  mins[i] <- min(sub_data$FID_, na.rm = T)
  maxs[i] <- max(sub_data$FID_, na.rm = T)
}

todelete = as.data.frame(cbind(maxs,mins))
colnames(todelete) = c("Delete", "Keep")
write.csv(todelete, "C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/large tree corrections from Charlie 8-18-18/northernMichigan_duplicates/1_2_ToDeleteList.csv", row.names = F)



########################################################################################
### Combining files with northern Michigan Northern and Lower Peninsula labeled FIDs ###
########################################################################################

#the UP and LP files were created by selecting all the UP corners and all the LP corners in GIS
#and saving each selection as their own csvs.

#the code below subsets the data to remove the corners that were removed in the clean v1.2 file (had been 
#left in the v1.2_interim file) and then combines the FIDs labeled by their LP and UP locations

setwd("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected")
UP = read.csv("northernMI_UP_FIDs.csv", header = TRUE, stringsAsFactors = FALSE)
LP = read.csv("northernMI_LP_FIDs.csv", header = TRUE, stringsAsFactors = FALSE)
library(dplyr)

UPsub = subset(UP,Delete_in %in% c("Leave in v1.2_interim. Remove in v1.2 clean. Perhaps add back to v1.3"))
LPsub = subset(LP,Delete_in %in% c("Leave in v1.2_interim. Remove in v1.2 clean. Perhaps add back to v1.3"))

UPtrue = subset(UP,Delete_in %in% c("#N/A", ""))
LPtrue = subset(LP,Delete_in %in% c("#N/A"))

UP_LP_true = rbind(UPtrue, LPtrue)
write.csv(UP_LP_true, file = "northernMI_UP_LP_FIDs.csv", row.names = FALSE)



######################################################
### add recnums to the northern mi database ##########
######################################################

#the recnum values continually get truncated when saving the csvs or when importing into GIS
#this is code to add the recnum values back in using the FIDs
library(dplyr)
nmi_v1.2 = read.csv("./UMW_Point Data from Goring 2016 PlosOne/raw UMW posted to wiki/northernmi_pls_projected_v1.2.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(nmi)

#get recnum values from northern michigan v1.1 data that is on the wiki
nmi1 = read.csv("C:/Users/jmurray7/Desktop/northernmi_pls_projected_v1.1.csv", header = TRUE, stringsAsFactors = FALSE)
options(scipen = 999) #use this to deal with the scientific notation issues.
head(nmi1)

#turn recnum and recnum_c into point values
nmi1 <- nmi1 %>% mutate(newrecnum = format(pmax(RECNUM,RECNUM_C), scientific = FALSE))  ## format prevents NAs when 9-11 digits are 0s
nmi1 <- nmi1 %>% mutate(point = as.numeric(substr(newrecnum, 9, 11)))

point.summary = as.data.frame(table(nmi1$point))
point.summary

#save the new point values in a csv. Go through that using the point.summary  to look at weird entries 
#for example: 0, 129, 131, 193, 202, 270, 379, 685, 756 
#point values 87, 89, 91, 93, 95, 96, 122, 123, 124, 125, 127 are all north/west borders
myvars <- c("FID", "RECNUM", "RECNUM_C", "newrecnum","point")
nmi1.points <- nmi1[myvars]
write.csv(nmi1.points, "recnum_point_summary.csv", row.names = FALSE)

#Jody went through the recnum_points_summary.csv and checked the weird point entries
#she added a notes column to desribe why the entries were weird and what was done about it.
#For example, there were 23 entries that were added to the northernMI_toDelete.csv because they
#were extra entries that were not at the 1/4 section or section corner, but there were corners nearby 
#that were on the corner and should be used instead or the entries were on the west border.
#There were also 7 entries that had weird point values, that when Jody looked them up on the GIS map, 
#they just needed their point number corrected for example the point had been 685 but was really point 85,
#or had been 202 and is really point 20, etc.


#################################################################
### Northern Michigan Level 0 to Level 3a conversion check ######
#################################################################

library(readr)
library(dplyr)
library(fields)

mi = read.csv("northernmi_pls_projected_v1.2.csv", header = TRUE, stringsAsFactors = FALSE)

spec_codes <- read_csv("level0_to_level3a_v0.4-9.csv", guess_max = 1000) %>% 
  filter(domain == "Upper Midwest") %>%
  select(level1, level3a) 

mi <- mi %>% 
  left_join(spec_codes, by = c('SPP1' = 'level1')) %>% rename(L1_tree1 = SPP1, L3_tree1 = level3a) %>%
  left_join(spec_codes, by = c('SPP2' = 'level1')) %>% rename(L1_tree2 = SPP2, L3_tree2 = level3a) 


unique(mi$L3_tree1)

unknownL3 = mi[which(is.na(mi$L3_tree1)),]
unknownL3
write.csv(unknownL3, "unknown_tree1_check.csv", row.names = FALSE)
#ran it once, and checked/corrected the unknown tree 1s to match the codes in the conversion file
#then ran the code again and saved it as the file below
write.csv(unknownL3, "unknown_tree1_check2.csv", row.names = FALSE)


unique(mi$L3_tree2)

unknown.tree2.L3 = mi[which(is.na(mi$L3_tree2)),]
write.csv(unknown.tree2.L3, "unknown_tree2_check.csv", row.names = FALSE)
#ran it once, and checked/corrected the unknown tree 2s to match the codes in the conversion file
#then ran the code again and saved it as the file below
write.csv(unknown.tree2.L3, "unknown_tree2_check2.csv", row.names = FALSE)



###################################################################################
## Southern Michigan checking Huron, Tuscola, Sanilac, etc Townships ##############
## to make sure Charlie's Correction Factors will be applied to the right corners##
###################################################################################

smi = read.csv("southernmi_projected_v1.5.csv", header = TRUE, stringsAsFactors = FALSE)

#Subset data from townships >=11N in Tuscola, Sanilac and Huron, and >=18N in Bay, Arenac, Gladwin and Clare

#Townships to pull the corners from

towns = read.csv("NETwns_EMeridian_southernMI", header = TRUE, stringsAsFactors = FALSE)

colnames(towns)

#select corners in smi by twnrng
smi.sub = subset(smi, twnrng %in% towns$TWNRNG)
#save as csv to map in GIS and check that all the corners fall where they are supposed to
write.csv(smi.sub, "NEcorners.csv", row.names = FALSE)

#########################################################################################################
###### Southern Michigan - find nearest 1/4 section/section corners to southern_MI_v0.3 corners #########
#########################################################################################################
setwd("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected")
rm(list=ls())

MI_corners = read.csv("MIcorners.csv", header = TRUE, stringsAsFactors = FALSE)
MI_data = read.csv("southernMI_projected_v0.3.csv", header = TRUE, stringsAsFactors = FALSE)
MI_data = read.csv("southernMI_projected_v0.8.csv", header = TRUE, stringsAsFactors = FALSE)


#MI_corners = read.csv("MIcorners_test.csv", header = TRUE, stringsAsFactors = FALSE)
#MI_data = read.csv("southernMI_test.csv", header = TRUE, stringsAsFactors = FALSE)


#add cornertype column to MI_data
cornertype = 1:length(MI_data$FID)
MI_data$cornertype = cornertype


pos = MI_corners[,c("POINT_X", "POINT_Y")]
new.pos = MI_data[,c("POINT_X", "POINT_Y")]


for (i in 1:length(MI_data$cornertype)) {
  new.pos2 = c(new.pos[i,1],new.pos[i,2])
  nearest.idx = which.min(colSums((t(pos)-new.pos2)^2))
  nearest.idx 
  MI_data[i,]$cornertype = MI_corners[nearest.idx,13]
}


write.csv(MI_data, "southernMI_projected_v0.4.csv", row.names = FALSE)


#test code for finding a single x,y in the list of MI_corners
new.pos = c(957402.5, 655507.1)
new.pos2 = c(new.pos[i,1],new.pos[i,2])
nearest.idx = which.min(colSums((t(pos)-new.pos2)^2))
nearest.idx
MI_corners[nearest.idx,35]

new.pos[5,]

#########################################################################
###### Southern Michigan - check for duplicates, diameters, etc #########
#########################################################################
setwd("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected")
rm(list=ls())

options(digits = 12)
MI = read.csv("southernMI.csv", stringsAsFactors = FALSE)
#southernMI.csv came from the southernMI_projected_v0.3.shp shapefile - the attribute table was exported
#as a csv
head(MI)

##### CHECK FOR DUPLICATES IN UNIQUEIDS #######
#check to see how many unique IDs you have in the uniqueID column and the FID column if this equals the length
#of the uniqueID column. If it doesn't you have a duplicate uniqueID.
uniqueIDcount = unique(MI$uniqueID)
FIDcount = unique(MI$FID)
length(uniqueIDcount)
length(uniqueIDcount) == length(MI$uniqueID)
#### IF THIS EQUALS FALSE YOU HAVE A UNIQUEID THAT IS DUPLICATED!!!######

######################## FIX THIS PROBLEM! Use Code Below ##############################

#create object with just FID and uniqueIDs
temporary = MI[,c("FID","uniqueID")]
head(temporary)
# find the duplicate uniqueIDs
dups <-temporary[duplicated(temporary$uniqueID)|duplicated(temporary$uniqueID, fromLast=TRUE),]
dups
head(dups)
write.csv(dups, file = "SouthernMI_Duplicates.csv", row.names = FALSE)

dupsV=cbind(c(duplicated(sort(MI$uniqueID))),MI$FID)
hmmV=cbind(c(which(apply(dupsV=="TRUE",1,any))))
hmmV
tail(MI$FID)


##################################################################
#Southern Michigan duplicate tree info, but different quadrats
##################################################################

#read in the southernMI attribute table
options(digits = 12)
MI = read.csv("southernMI.csv", stringsAsFactors = FALSE)
#southernMI.csv came from the southernMI_projected_v0.3.shp shapefile - the attribute table was exported
#as a csv. I think added the columns uniquetree1, uniquetree2, uniquetree3, and uniquetree4.
#these uniquetreeX were the concatenation of speciesX,distX,azX,diamZ. The QX was not included in the concatenation.

head(MI)

uniquetree1 = paste0(MI$species1,MI$dist1,MI$az1,MI$diam1)
uniquetree2 = paste0(MI$species2,MI$dist2,MI$az2,MI$diam2)

MI$uniquetree1 = uniquetree1
MI$uniquetree2 = uniquetree2

temporary = MI[,c("FID","uniquetree1","uniquetree2")]

out = 1:length(temporary$uniquetree1)
for(i in 1:length(temporary$uniquetree1)){
out[i] = match(temporary[i,]$uniquetree1,temporary$uniquetree2)
}

y = temporary[!is.na(out),]

head(y)


table(out)

write.csv(matchtemp, file = "SouthernMI_Tree1_2_Match.csv", row.names = T)


for(i in 1:length(temporary$uniquetree1)){
  out[i] = temporary$uniquetree1[i]==temporary$uniquetree2}

tree1 = MI[,c("FID","uniquetree1")]
tree2 = MI[,c("FID","uniquetree2")]

##################################################################################
##Code to look at diameters and distances in southernMI_projected_v0.7 shapefile##
##################################################################################

#The next 3 lines of code were used to look at diameters and distances for version 0.7 from the GIS dbf
#library(foreign)
#dbf = read.dbf("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/southernMI_projected_v0.3.dbf", as.is = TRUE)
#head(dbf)


MI_data = read.csv("southernMI_projected_v0.7.csv", header = TRUE, stringsAsFactors = FALSE)

##################
### Diam1 ########
##################

#remove 9999 from diam1
diam1 = MI_data$diam1
sort(diam1)
remove = 9999
values = sort(unique(diam1))
table(diam1)
diam1_no9999 = diam1[!diam1%in%remove]
#diam1_9999 = diam1[diam1%in%remove]

#create historgram of diam1 with 9999 removed
hist(diam1_no9999,breaks = 200, xlim=range(1:1500))
diam1.summary = summary(diam1_no9999)
diam1.summary
diam1.count = length(diam1_no9999)

#table of diam1 values with the 9999 entries removed
diam1_table = as.data.frame(table(diam1_no9999))
diam1_table

#select only the diam1 that are <= 200 and create a histogram
newdiam1 = dbf[which(dbf$diam1<=200),]
summary(newdiam1$diam1)
hist(newdiam1$diam1)

png(paste("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/RCode_Results/diam1hist.png", sep = " "),   height = 768, width=1024)
hist(newdiam1$diam1, xlab = "Diam1", main = "Histogram of Diam1 Values 0-200")
dev.off()


##################
### Diam2 ########
##################

#remove 9999 from diam2
diam2 = MI_data$diam2
sort(diam2)
remove = 9999
values = sort(unique(diam2))
values
table(diam2)
diam2_no9999 = diam2[!diam2%in%remove]

#create historgram of diam2 with 9999 removed
hist(diam2_no9999,breaks = 200, xlim=range(1:1500))
diam2.summary = summary(diam2_no9999)
diam2.summary
diam2.count = length(diam2_no9999)


#table of diam2 values with the 9999 entries removed
diam2_table = as.data.frame(table(diam2_no9999))
diam2_table

#select only the diam2 that are <= 200 and create a histogram
newdiam2 = dbf[which(dbf$diam2<=200),]
summary(newdiam2$diam2)
hist(newdiam2$diam2)

png(paste("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/RCode_Results/diam2hist.png", sep = " "),   height = 768, width=1024)
hist(newdiam2$diam2, xlab = "Diam2", main = "Histogram of Diam2 Values 0-200")
dev.off()

newdiam2.a = dbf[which(dbf$diam2>0 & dbf$diam2<9999),]
summary(newdiam2.a$diam2)
hist(newdiam2.a$diam2, breaks = 200, xlim=range(1:1200))

png(paste("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/RCode_Results/diam2hist_1-200.png", sep = " "),   height = 768, width=1024)
hist(newdiam2.a$diam2, breaks = 200, xlim=range(1:1200), xlab = "Diam2", main = "Histogram of Diam2 Values 1-1200")
dev.off()

newdiam1.a = dbf[which(dbf$diam1>0 & dbf$diam1<9999),]
summary(newdiam1.a$diam1)
hist(newdiam1.a$diam1, breaks = 200, xlim=range(1:1200))

png(paste("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/RCode_Results/diam1hist_1-200.png", sep = " "),   height = 768, width=1024)
hist(newdiam1.a$diam1, breaks = 200, xlim=range(1:1200), xlab = "Diam2", main = "Histogram of Diam1 Values 1-1200")
dev.off()



##################
### Diam3 ########
##################

#remove 9999 from diam3
diam3 = MI_data$diam3
remove = 9999
values = sort(unique(diam3))
values
table(diam3)
diam3_no9999 = diam3[!diam3%in%remove]

#create historgram of diam3 with 9999 removed
hist(diam3_no9999,breaks = 200, xlim=range(1:300))
diam3.summary = summary(diam3_no9999)
diam3.summary
diam3.count = length(diam3_no9999)

png(paste("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/RCode_Results/diam3_no9999hist.png", sep = " "),   height = 768, width=1024)
hist(diam3_no9999,breaks = 200, xlim=range(1:300), xlab = "Diam3", main = "Histogram of Diam3 Values 0-300")
dev.off()


#table of diam3 values with the 9999 entries removed
diam3_table = as.data.frame(table(diam3_no9999))
diam3_table

#select only the diam3 that are <= 200 and create a histogram
newdiam3 = dbf[which(dbf$diam3>0 & dbf$diam3<9999),]
summary(newdiam3$diam3)
hist(newdiam3$diam3, breaks = 200, xlim=range(1:1200))

png(paste("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/RCode_Results/diam3hist.png", sep = " "),   height = 768, width=1024)
hist(newdiam3$diam3, breaks = 200, xlim=range(1:1200), xlab = "Diam3", main = "Histogram of Diam3 Values 1-1200")
dev.off()

##################
### Diam4 ########
##################

#remove 9999 from diam4
diam4 = MI_data$diam4
remove = 9999
values = sort(unique(diam4))
values
table(diam4)
diam4_no9999 = diam4[!diam4%in%remove]

#create historgram of diam4 with 9999 removed
hist(diam4_no9999,breaks = 200, xlim=range(1:300))
diam4.summary = summary(diam4_no9999)
diam4.summary
diam4.count = length(diam4_no9999)

png(paste("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/RCode_Results/diam4_no9999hist.png", sep = " "),   height = 768, width=1024)
hist(diam4_no9999,breaks = 200, xlim=range(1:300), xlab = "Diam4", main = "Histogram of Diam4 Values 0-300")
dev.off()


#table of diam4 values with the 9999 entries removed
diam4_table = as.data.frame(table(diam4_no9999))
diam4_table

#select only the diam4 that are <= 200 and create a histogram
newdiam4 = dbf[which(dbf$diam4>0 & dbf$diam4<9999),]
summary(newdiam4$diam4)
hist(newdiam4$diam4, breaks = 200, xlim=range(1:1200))

png(paste("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/RCode_Results/diam4hist.png", sep = " "),   height = 768, width=1024)
hist(newdiam4$diam4, breaks = 200, xlim=range(1:1200), xlab = "Diam4", main = "Histogram of Diam4 Values 1-1200")
dev.off()
length(newdiam4$diam4)



png(paste("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/RCode_Results/diam_histograms.png", sep = " "),   height = 768, width=1024)
par(mfrow = c(2, 2))
hist(newdiam1.a$diam1, breaks = 200, xlim=range(1:1200), xlab = "Diam2", main = "Histogram of Diam1 Values 1-1200")
hist(newdiam2.a$diam2, breaks = 200, xlim=range(1:1200), xlab = "Diam2", main = "Histogram of Diam2 Values 1-1200")
hist(newdiam3$diam3, breaks = 200, xlim=range(1:1200), xlab = "Diam3", main = "Histogram of Diam3 Values 1-1200")
hist(newdiam4$diam4, breaks = 200, xlim=range(1:1200), xlab = "Diam4", main = "Histogram of Diam4 Values 1-1200")
dev.off()


##################
### Dist1 ########
##################

#remove 9999 from dist1
dist1 = MI_data$dist1
sort(diam1)
remove = 9999
values = sort(unique(dist1))
values
table(dist1)
dist1_no9999 = dist1[!dist1%in%remove]

#create historgram of dist1 with 9999 removed
hist(dist1_no9999,breaks = 200, xlim=range(1:500))
table = as.data.frame(order(table(dist1_no9999)))
table
dist1.summary = summary(dist1_no9999)
dist1.summary


png(paste("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/RCode_Results/diam1hist.png", sep = " "),   height = 768, width=1024)
hist(newdiam1$diam1, xlab = "Diam1", main = "Histogram of Diam1 Values 0-200")
dev.off()


##################
### Dist2 ########
##################

#remove 9999 from dist2
dist2 = MI_data$dist2
remove = 9999
values = sort(unique(dist2))
values
table(dist2)
dist2_no9999 = dist2[!dist2%in%remove]

#create historgram of diam2 with 9999 removed
hist(dist2_no9999,breaks = 200, xlim=range(1:250))
diam2.summary = summary(diam2_no9999)
diam2.summary
diam2.count = length(diam2_no9999)


#table of diam2 values with the 9999 entries removed
diam2_table = as.data.frame(table(diam2_no9999))
diam2_table

#select only the diam2 that are <= 200 and create a histogram
newdiam2 = dbf[which(dbf$diam2<=200),]
summary(newdiam2$diam2)
hist(newdiam2$diam2)

png(paste("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/RCode_Results/diam2hist.png", sep = " "),   height = 768, width=1024)
hist(newdiam2$diam2, xlab = "Diam2", main = "Histogram of Diam2 Values 0-200")
dev.off()


##################
### Dist3 ########
##################

#remove 9999 from dist3
dist3 = MI_data$dist3
remove = 9999
values = sort(unique(dist3))
values
table(dist3)
dist3_no9999 = dist3[!dist3%in%remove]

#create historgram of diam3 with 9999 removed
hist(dist3_no9999,breaks = 200, xlim=range(1:150))
diam3.summary = summary(diam3_no9999)
diam3.summary
diam3.count = length(diam3_no9999)

png(paste("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/RCode_Results/diam3_no9999hist.png", sep = " "),   height = 768, width=1024)
hist(diam3_no9999,breaks = 200, xlim=range(1:300), xlab = "Diam3", main = "Histogram of Diam3 Values 0-300")
dev.off()


#table of diam3 values with the 9999 entries removed
diam3_table = as.data.frame(table(diam3_no9999))
diam3_table

#select only the diam2 that are <= 200 and create a histogram
newdiam3 = dbf[which(dbf$diam3>0 & dbf$diam3<9999),]
summary(newdiam3$diam3)
hist(newdiam3$diam3, breaks = 200, xlim=range(1:1200))

png(paste("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/RCode_Results/diam3hist.png", sep = " "),   height = 768, width=1024)
hist(newdiam3$diam3, breaks = 200, xlim=range(1:1200), xlab = "Diam3", main = "Histogram of Diam3 Values 1-1200")
dev.off()

##################
### Dist4 ########
##################

#remove 9999 from diam4
dist4 = MI_data$dist4
remove = 9999
values = sort(unique(dist4))
values
table(dist4)
dist4_no9999 = dist4[!dist4%in%remove]

#create historgram of diam3 with 9999 removed
hist(dist4_no9999,breaks = 200, xlim=range(1:200))
diam4.summary = summary(diam4_no9999)
diam4.summary
diam4.count = length(diam4_no9999)

png(paste("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/RCode_Results/diam4_no9999hist.png", sep = " "),   height = 768, width=1024)
hist(diam4_no9999,breaks = 200, xlim=range(1:300), xlab = "Diam4", main = "Histogram of Diam4 Values 0-300")
dev.off()



#table of diam3 values with the 9999 entries removed
diam4_table = as.data.frame(table(diam4_no9999))
diam4_table

#select only the diam2 that are <= 200 and create a histogram
newdiam4 = dbf[which(dbf$diam4>0 & dbf$diam4<9999),]
summary(newdiam4$diam4)
hist(newdiam4$diam4, breaks = 200, xlim=range(1:1200))

png(paste("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/RCode_Results/diam4hist.png", sep = " "),   height = 768, width=1024)
hist(newdiam4$diam4, breaks = 200, xlim=range(1:1200), xlab = "Diam4", main = "Histogram of Diam4 Values 1-1200")
dev.off()
length(newdiam4$diam4)




output <- list(diam1.summary = diam1.summary,
               diam1.count = diam1.count,
               diam2.summary = diam1.summary,
               diam2.count = diam1.count,
               diam3.summary = diam1.summary,
               diam3.count = diam1.count,
               diam4.summary = diam1.summary,
               diam4.count = diam1.count
)



# Print output 

sink(file = paste("./RCode_Results/Diam_Dist_Output.txt", sep=""))
output
sink() # stops sinking


####################################################################################
##Code to create csv with large diameters from southernMI_projected_v0.3 shapefile##
####################################################################################
rm(list=ls())

library(foreign)
dbf = read.dbf("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/southernMI_projected_v0.3.dbf", as.is = TRUE)
head(dbf)
length(dbf$fid_1)

dbf["bigdiam"] = NA

#created dataframes for each diameter
diam4 = dbf[which(dbf$diam4 >254),] #95
diam3 = dbf[which(dbf$diam3 >254),] #106
diam2 = dbf[which(dbf$diam2 >254),] #958
diam1 = dbf[which(dbf$diam1 >254),] #1214

bigdiam = rbind(diam4,diam3,diam2,diam1)

bigdiam = bigdiam[,-36]

#Find the number of uniqueIDs to see if you have duplicate entries
uniqueIDcount = unique(bigdiam$uniqueID)
length(uniqueIDcount)
length(uniqueIDcount) == length(bigdiam$uniqueID) #if FALSE you have duplicate entries

#Find duplicate entries
dups <-bigdiam[duplicated(bigdiam$uniqueID)|duplicated(bigdiam$uniqueID, fromLast=TRUE),]
dups

#this returns a table with the uniqueID label and how many times it is duplicated.
table(dups$uniqueID)







for(i in 1:length(dbf$bigdiam)){
  if(dbf$diam1[i] > 254){
    dbf$bigdiam[i] <- 'big'
  }
  else if(dbf$diam2[i] > 254){
    dbf$bigdiam[i] <- 'big'
  }
  else if(dbf$diam3[i] >254){
    dbf$bigdiam[i] <- 'big'
  }
}



#check that there are no duplicates in the x,y coordinates
cornersPointXcount = unique(MI_data$point_x)
length(cornersPointXcount)
length(cornersPointXcount) == length(MI_data$point_x) #needs to say TRUE

cornersPointYcount = unique(MI_data$point_y)
length(cornersPointYcount)
length(cornersPointYcount) == length(MI_data$point_y) #needs to say TRUE


xdups <-MI_data[duplicated(MI_data$point_x)|duplicated(MI_data$point_x, fromLast=TRUE),]
xdups
xdups2 = xdups[order(xdups$point_x),]


ydups <-MI_data[duplicated(MI_data$point_y)|duplicated(MI_data$point_y, fromLast=TRUE),]
ydups
ydups2 = ydups[order(ydups$point_y),]

xdups <-MI_data[duplicated(MI_data$uniqueID)|duplicated(MI_data$unique, fromLast=TRUE),]
xdups
xdups2 = xdups[order(xdups$point_x),]




################################################################################################
## Create a Level0 to Level3a conversion table that can be used to update the conversion file###
################################################################################################

#this was done above when the L1 & L3 trees were checked. But can use the code again here.
rm(list=ls())
setwd("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected")
state = read.csv("./southernMI_projected_v1.2.csv", header = TRUE, stringsAsFactors = FALSE)

#combine all the L1 and L3 trees, to get a table with counts of all 4 trees L1/L3 labels
#select L1 & L3 of just tree1
L1.L3tree1 = state[,c("species1","level3a_tree1")]
colnames(L1.L3tree1) = c("L1_tree", "L3_tree")


#select L1 & L3 of just tree2 and remove NAs
L1.L3tree2 = state[,c("species2","level3a_tree2")]
colnames(L1.L3tree2) = c("L1_tree", "L3_tree")
L1.L3tree2 = L1.L3tree2[complete.cases(L1.L3tree2),]


#select L1 & L3 of just tree3 and remove NAs
L1.L3tree3 = state[,c("species3","level3a_tree3")]
L1.L3tree3 = L1.L3tree3[complete.cases(L1.L3tree3),]
colnames(L1.L3tree3) = c("L1_tree", "L3_tree")


#select L1 & L3 of just tree4 and remove NAs
L1.L3tree4= state[,c("species4","level3a_tree4")]
L1.L3tree4 = L1.L3tree4[complete.cases(L1.L3tree4),]
colnames(L1.L3tree4) = c("L1_tree", "L3_tree")

#combine all 4 L1&L3 trees
combined = rbind(L1.L3tree1,L1.L3tree2,L1.L3tree3,L1.L3tree4)

#create a table of the counts of L1 trees in the L3 categories
library(dplyr)
L1.L3combined = combined %>% group_by(L1_tree) %>% tally()

write.csv(L1.L3combined, file = "./southernMI_projected_v1.2_L1tree_summary.csv", row.names = FALSE)



#join Level2 and Comments from the L0 to L3 conversion file.
#read in the conversion file
conversion = read.csv("C:/Users/jmurray7/Dropbox/PalEON2/Conversion Tables - Allometry, PEcAn, PLS, FIA/conversion files uploaded to wiki/level0_to_level3a_v0.4-7_NEED sMI udpate.csv",header = TRUE, stringsAsFactors = FALSE)
#subset just the domain you want
#use unique(conversion$domain) to find the name of the Domains
sMIconversion = conversion[which(conversion$domain == "Southern MI"),]


#merge the L1.L3 count file with the conversion file - keeping all the L1.L3 entries (left outer join)
sMIconversion_merge = merge(x=L1.L3combined, y=sMIconversion, by.x = "L1_tree", by.y = "level1", all.x = TRUE)

#now the column headings are not in the same order as the IL conversion so we won't be able to seamlessly combine
#the updated conversion data with the old conversione file. SO select each column and then join into a database
#that is in the same order as the conversion file
level0 = sMIconversion_merge$level0
level1 = sMIconversion_merge$L1_tree
level2 = sMIconversion_merge$level2
check = sMIconversion_merge$check
level3a = sMIconversion_merge$level3a
count = sMIconversion_merge$n
domain = sMIconversion_merge$domain
comments = sMIconversion_merge$comments

conversiontable = data.frame(cbind(level0,level1, level2, check, level3a, count, domain, comments), stringsAsFactors = FALSE)
colnames(conversiontable) = c("level0","level1","level2","check","level3a","count","domain","comments")

write.csv(conversiontable, file = "./southernMI_projected_v1.2_conversion.csv", row.names = FALSE)

#################################################################################
## Dealing with corners with 4 trees and 3 trees. Make them be 2 tree corners. ##
#################################################################################

#10-2-18
#In southern Michigan we have a number of corners with 4 trees or 3 trees. 
#Charlie says these are really all 2 tree corners. He checked all the interior corners (n=68). 
#We need to also fix exterior township corners (n>4000). Because there were so many Charlie came up with 
#the following rule to help us remove the additional 2 trees.
#For S border corners we keep trees in quadrats 1 and 4, for E border corners keep trees in quadrats 3 and 4.

#For the exterior corners, we needed to label them as South or East border which JP did in GIS.
#Shapefiles of just the corner entries from southernmi_pls_projected_v1.5_interim.csv (the working version that had these 4 and 3 tree corners cleaned up - you can find them in v1.4 if need be)
#with 4 trees or 3 trees were created in GIS. Then JP went through and labeled them as "s" or "e" Border
#corresponding with their location on the township borders. 
#Then the attribute tables from the shapefiles were saved as 4tree_corners_tocorrect.csv and 3tree_corners_tocorrect.csv


#The following code determines what quadrats trees are in at each corner.
library(dplyr)

setwd("C:/Users/jmurray7/Dropbox/GIS PalEON/Michigan PLSS/Michigan Projected/southernMI_4tree,3tree corners")

#Selecting S and E borders for corners with 4 trees
tree.4 = read.csv("4tree_corners_tocorrect.csv", header = TRUE, stringsAsFactors = FALSE)

tree.4.s = subset(tree.4,Border %in% c("s"))
tree.4.e = subset(tree.4,Border %in% c("e"))

quadrat.4 = paste0(tree.4.s$q1,tree.4.s$q2,tree.4.s$q3,tree.4.s$q4)
quadrat.4.e = paste0(tree.4.e$q1,tree.4.e$q2,tree.4.e$q3,tree.4.e$q4)

tree.4.s$quadrat = quadrat.4
tree.4.e$quadrat = quadrat.4.e


#selecting 4 tree corners with south border quadrats to keep
tree4.quadrats = as.data.frame(sort(table(tree.4.s$quadrat)))
write.csv(tree4.quadrats,file = "tree4_quadrat_summary.csv", row.names = FALSE)

tree.4.s.correctquadrats = subset(tree.4.s,quadrat %in% c("1004", "1024", "1034","1240","1340","1403","1432","2314","2412","3412","4123",
                                                          "4321","1433","4312","1204","1342","2413","2134","2341","1243",
                                                          "1324","1224","1423","1234","1334"))

tree.4.s.not_correctquadrats = subset(tree.4.s,!quadrat %in% c("1004", "1024", "1034","1240","1340","1403","1432","2314","2412","3412","4123",
                                                          "4321","1433","4312","1204","1342","2413","2134","2341","1243",
                                                          "1324","1224","1423","1234","1334"))

trees.4.s.correctquadrats.sort = tree.4.s.correctquadrats[order(tree.4.s.correctquadrats$quadrat),]

write.csv(tree.4.s.correctquadrats,file="tree4SBorder_CORRECT_to_2trees.csv",row.names = FALSE)

#selecting 4 tree corners with east border quadrats to keep
tree4.e.quadrats = as.data.frame(sort(table(tree.4.e$quadrat)))
write.csv(tree4.e.quadrats,file = "tree4_quadrat_Eborder_summary.csv", row.names = FALSE)

tree.4.e.correctquadrats = subset(tree.4.e,quadrat %in% c("1314", "3412", "3421","4322","2314","2340",
                                                          "1342","4312"," 234","1423","2413","1034","1134","1324",
                                                          "2134","2341","2234","4123","1243","1234"))

tree.4.e.not_correctquadrats = subset(tree.4.e,!quadrat %in% c("1314", "3412", "3421","4322", "2314","2340",
                                                          "1342","4312","234","1423","2413","1034","1134","1324",
                                                          "2134","2341","2234","4123","1243","1234"))

trees.4.e.correctquadrats.sort = tree.4.e.correctquadrats[order(tree.4.e.correctquadrats$quadrat),]

write.csv(tree.4.e.correctquadrats,file="tree4EBorder_CORRECT_to_2trees.csv",row.names = FALSE)

#Selecting S and E borders for corners with 3 trees
tree.3 = read.csv("3tree_corners_tocorrect.csv", header = TRUE, stringsAsFactors = FALSE)

tree.3.s = subset(tree.3,Border %in% c("s")) 
tree.3.e = subset(tree.3,Border %in% c("e"))

quadrat.3.s = paste0(tree.3.s$q1,tree.3.s$q2,tree.3.s$q3,tree.3.s$q4)
quadrat.3.e = paste0(tree.3.e$q1,tree.3.e$q2,tree.3.e$q3,tree.3.e$q4)

tree.3.s$quadrat = quadrat.3.s
tree.3.e$quadrat = quadrat.3.e

#selecting 3 tree corners with south border quadrats to keep
tree3.quadrats = as.data.frame(sort(table(tree.3.s$quadrat)))
write.csv(tree3.quadrats,file = "tree3_quadrat_SBorder_summary.csv", row.names = FALSE)

tree.3.s.correctquadrats = subset(tree.3.s,quadrat %in% c("2140", "4120","4130","4310","1040","1430","1240","1340"))

tree.3.s.not_correctquadrats = subset(tree.3.s,!quadrat %in% c("2140", "4120","4130","4310","1040","1430","1240","1340"))

trees.3.s.correctquadrats.sort = tree.3.s.correctquadrats[order(tree.3.s.correctquadrats$quadrat),]

write.csv(tree.3.s.correctquadrats,file="tree3SBorder_CORRECT_to_2trees.csv",row.names = FALSE)

#selecting 3 tree corners with east border quadrats to keep
tree3.e.quadrats = as.data.frame(sort(table(tree.3.e$quadrat)))
write.csv(tree3.e.quadrats,file = "tree3_quadrat_Eborder_summary.csv", row.names = FALSE)

tree.3.e.correctquadrats = subset(tree.3.e,quadrat %in% c("3410","4320","4230","4310","1430","1340","2340"))

tree.3.e.not_correctquadrats = subset(tree.3.e,!quadrat %in% c("3410","4320","4230","4310","1430","1340","2340"))

trees.3.e.correctquadrats.sort = tree.3.e.correctquadrats[order(tree.3.e.correctquadrats$quadrat),]

write.csv(tree.3.e.correctquadrats,file="tree3EBorder_CORRECT_to_2trees.csv",row.names = FALSE)


#combine all the corners that have trees in multiple quadrats 
#(e.g., two trees in quad 1 or quad 4 (i.e.,1124,1434, etc) for S border)
#these cannot be easily be fixed and need someone to go back to the original notes to deterrmine 
#which two trees to keep

#remove the extra column from tree 3 file in order to rbind tree 3 and tree 4 data
tree.3.e.not_correctquadrats = tree.3.e.not_correctquadrats[,-2]
tree.3.s.not_correctquadrats = tree.3.s.not_correctquadrats[,-2]

tree.4.s.not_correctquadrats = tree.4.s.not_correctquadrats[,-26]
tree.4.e.not_correctquadrats = tree.4.e.not_correctquadrats[,-26]
colnames(tree.4.s.not_correctquadrats)[44] <- "Border"
colnames(tree.4.e.not_correctquadrats)[44] <- "Border"

trees.to.fix = rbind(tree.4.e.not_correctquadrats, tree.4.s.not_correctquadrats,tree.3.s.not_correctquadrats,tree.3.e.not_correctquadrats)
write.csv(trees.to.fix, file="Corners_to_CORRECT_for_v1.6.csv", row.names = FALSE)


