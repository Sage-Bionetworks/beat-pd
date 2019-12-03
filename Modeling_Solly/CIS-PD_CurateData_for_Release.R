setwd("~/Documents/PDDB Challenge 2/CIS-PD")

library(synapser)
library(gdata)
library(plyr)


#Load labels
scores_query <- synTableQuery("select * from syn20489608")
scores <- scores_query$asDataFrame()
scores$on_off[scores$on_off==-1]<-NA

# Load features
#features<-read.csv("Analyses/Features/CIS-PD_Watch_Features_Full_NShawen.csv", header=T, as.is=T)


# Load Nick's sensor QC summary data
synfile<-synGet("syn21071756", downloadLocation="Analyses/Features/")
msmtqc<-read.csv(synfile$path, header=T, as.is=T)

winthresh<-4
scores<-scores[scores$measurement_id%in%msmtqc$measurement_id[msmtqc$Highest...of.Consecutive.Above.Threshold.Windows>=winthresh],]

#scores<-scores[!(is.na(scores$on_off)&is.na(scores$dyskinesia)&is.na(scores$tremor)),]


indidvec<-unique(scores$subject_id)


# Load Random Splits from 
synfile<-synGet("syn21095189", downloadLocation="Analyses/")
splits<-read.csv(synfile$path, header=T, as.is=T)

# Curate Training, Test and Ancillary Data
splitid<-7

training<-splits[splits[,paste("training", splitid, sep="")], 1:5]
test<-splits[!splits[,paste("training", splitid, sep="")], 1:2]

testinternal<-splits[!splits[,paste("training", splitid, sep="")], 1:5]

ancillary<-scores[!scores$subject_id%in%splits$subject_id,]
ancillary<-ancillary[,c("measurement_id", "subject_id", "on_off", "dyskinesia", "tremor")]

#Versions to Share
write.csv(training, "../Data Release/CIS-PD_Training_Data_IDs_Labels.csv", row.names=F, quote=F)
write.csv(test, "../Data Release/CIS-PD_Test_Data_IDs.csv", row.names=F, quote=F)
write.csv(ancillary, "../Data Release/CIS-PD_Ancillary_Data_IDs_Labels.csv", row.names=F, quote=F)

#Internal ONLY
write.csv(testinternal, "../Data Release/CIS-PD_Test_Data_IDs_Labels_INTERNALDoNOTShare.csv", row.names=F, quote=F)


storageloc<-"syn21291569"

synfile<-File("../Data Release/CIS-PD_Training_Data_IDs_Labels.csv", parentId=storageloc)
synfile<-synStore(synfile)
synfile<-File("../Data Release/CIS-PD_Test_Data_IDs.csv", parentId=storageloc)
synfile<-synStore(synfile)
synfile<-File("../Data Release/CIS-PD_Test_Data_IDs_Labels_INTERNALDoNOTShare.csv", parentId=storageloc)
synfile<-synStore(synfile)
synfile<-File("../Data Release/CIS-PD_Ancillary_Data_IDs_Labels.csv", parentId=storageloc)
synfile<-synStore(synfile)

idsforrelease<-sort(unique(c(training$subject_id, ancillary$subject_id)))


# Download and curate Demographic data and UPDRS Part 3

updrs_query <- synTableQuery("select * from syn18435297")
updrs <- updrs_query$asDataFrame()
updrs<-updrs[,-c(1:2)]
names(updrs)[1]<-"subject_id"
updrs<-updrs[updrs$subject_id%in%idsforrelease,]
updrs<-updrs[updrs$Visit%in%c("Baseline", "2 Weeks: Time 0", "2 Weeks: Time 60"),]
names(updrs)[-c(1:3)]<-paste("UPDRS_", names(updrs)[-c(1:3)], sep="")
updrs<-updrs[,1:(dim(updrs)[2]-1)]
# Keep baseline only for those that don't have week 2
tmptbl<-table(updrs$subject_id)
tmpids<-names(tmptbl[tmptbl==1])
updrs<-rbind(updrs[updrs$Visit!="Baseline",],updrs[updrs$subject_id%in%tmpids&updrs$Visit=="Baseline",])
updrs<-updrs[order(updrs$subject_id),]
# Harmonize Hoen & Yahr to numeric
updrs$UPDRS_3.20[updrs$UPDRS_3.20=="Unilateral involvement only"]<-1
updrs$UPDRS_3.20[updrs$UPDRS_3.20=="Bilateral involvement without impairment of balance"]<-2
updrs$UPDRS_3.20[updrs$UPDRS_3.20=="Mild to moderate involvement; some postural instability but physically independent; needs assistance to recover from pull test."]<-3
updrs$UPDRS_3.20<-as.numeric(updrs$UPDRS_3.20)

write.csv(updrs, "../Data Release/CIS-PD_UPDRS_Part3.csv", row.names=F, quote=F)
synfile<-File("../Data Release/CIS-PD_UPDRS_Part3.csv", parentId=storageloc)
synfile<-synStore(synfile)
write.csv(names(updrs), "../Data Release/CIS-PD_updrsIII_names.csv", row.names = F, quote=F)

demos_query <- synTableQuery("select * from syn18418057")
demos <- demos_query$asDataFrame()
demos<-demos[,-c(1:2)]
names(demos)[1]<-"subject_id"
demos<-demos[,c("subject_id", "Age", "Gender", "Race", "Ethnicity")]
demos<-demos[demos$subject_id%in%idsforrelease,]
write.csv(demos, "../Data Release/CIS-PD_Demographics.csv", row.names=F, quote=F)
synfile<-File("../Data Release/CIS-PD_Demographics.csv", parentId=storageloc)
synfile<-synStore(synfile)


updrs124_query <- synTableQuery("select * from syn18418786")
updrs124 <- updrs124_query$asDataFrame()
updrs124<-updrs124[,-c(1:2)]
names(updrs124)[1]<-"subject_id"
updrs124<-updrs124[updrs124$subject_id%in%idsforrelease,]
updrs124<-updrs124[updrs124$Visit=="Baseline",]
updrs124$UPDRS_PartI_Total<-apply(updrs124[,substr(names(updrs124), 1,1)=="1"],1,sum)
updrs124$UPDRS_PartII_Total<-apply(updrs124[,substr(names(updrs124), 1,1)=="2"],1,sum)
updrs124<-updrs124[,!substr(names(updrs124), 1,1)%in%c(1,2)]
updrs124<-updrs124[,-c(3:5)]
updrs124<-updrs124[,c(1:2, 9:10, 3:8)]
names(updrs124)[substr(names(updrs124),1,1)=="4"]<-paste("UPDRS_", names(updrs124)[substr(names(updrs124),1,1)=="4"], sep="")
write.csv(updrs124, "../Data Release/CIS-PD_UPDRS_Part1_2_4.csv", row.names=F, quote=F)
synfile<-File("../Data Release/CIS-PD_UPDRS_Part1_2_4.csv", parentId=storageloc)
synfile<-synStore(synfile)
