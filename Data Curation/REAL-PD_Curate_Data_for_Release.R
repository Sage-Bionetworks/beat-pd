setwd("~/Documents/PDDB Challenge 2/REAL-PD")

library(synapser)
library(gdata)
library(plyr)
library(dplyr)
library(lubridate)


alltraits<-c("tremor", "on_off", "dyskinesia")

#Load labels
scores_query <- synTableQuery("select * from syn20822276")
scores <- scores_query$asDataFrame()


#scores<-scores[!apply(is.na(scores[,alltraits]), 1, all),]


# Load Nick's sensor QC summary data
synfile<-synGet("syn21123211", downloadLocation="Analyses/Features/")
msmtqcphone<-read.csv(synfile$path, header=T, as.is=T)

synfile<-synGet("syn21123212", downloadLocation="Analyses/Features/")
msmtqcwatch<-read.csv(synfile$path, header=T, as.is=T)

msmtqcmerge<-full_join(msmtqcwatch, msmtqcphone, by=c("measurement_id","subject_id","medication_state","slowness_walking","tremor","on_off","dyskinesia"),
                       suffix=c(".watch", ".phone"))


# Threshold based on enough watch or phone data
winthresh<-4
scores<-scores[scores$measurement_id%in%msmtqcmerge$measurement_id[(msmtqcmerge$Highest...of.Consecutive.Above.Threshold.Windows.watch >= winthresh&!is.na(msmtqcmerge$Highest...of.Consecutive.Above.Threshold.Windows.watch))|(msmtqcmerge$Highest...of.Consecutive.Above.Threshold.Windows.phone >= winthresh&!is.na(msmtqcmerge$Highest...of.Consecutive.Above.Threshold.Windows.phone))],]

#scores<-scores[!(is.na(scores$on_off)&is.na(scores$dyskinesia)&is.na(scores$tremor)),]


# Load Random Splits from 
synfile<-synGet("syn21141640", downloadLocation="Analyses/")
splits<-read.csv(synfile$path, header=T, as.is=T)
# Transform Tremor to tremor-1 to put it on 0-4 scale
splits$tremor<-splits$tremor-1
#scores$tremor<-scores$tremor-1

# Transform on/off by 1-on_off to be consistent with CIS-PD
# where 0=ON and 4=OFF
splits$on_off<-1-splits$on_off
scores$on_off<-1-scores$on_off

# Curate Training, Test and Ancillary Data
splitid<-1

training<-splits[splits[,paste("training", splitid, sep="")], 1:5]
test<-splits[!splits[,paste("training", splitid, sep="")], 1:2]

testinternal<-splits[!splits[,paste("training", splitid, sep="")], 1:5]

ancillary<-scores[!scores$subject_id%in%splits$subject_id,]
ancillary<-ancillary[,c("measurement_id", "subject_id", "on_off", "dyskinesia", "tremor")]

#Versions to Share
write.csv(training, "../Data Release/REAL-PD_Training_Data_IDs_Labels.csv", row.names=F, quote=F)
write.csv(test, "../Data Release/REAL-PD_Test_Data_IDs.csv", row.names=F, quote=F)
write.csv(ancillary, "../Data Release/REAL-PD_Ancillary_Data_IDs_Labels.csv", row.names=F, quote=F)

#Internal ONLY
write.csv(testinternal, "../Data Release/REAL-PD_Test_Data_IDs_Labels_INTERNALDoNOTShare.csv", row.names=F, quote=F)


storageloc<-"syn21291570"

synfile<-File("../Data Release/REAL-PD_Training_Data_IDs_Labels.csv", parentId=storageloc)
synfile<-synStore(synfile)
synfile<-File("../Data Release/REAL-PD_Test_Data_IDs.csv", parentId=storageloc)
synfile<-synStore(synfile)
synfile<-File("../Data Release/REAL-PD_Test_Data_IDs_Labels_INTERNALDoNOTShare.csv", parentId=storageloc)
synfile<-synStore(synfile)
synfile<-File("../Data Release/REAL-PD_Ancillary_Data_IDs_Labels.csv", parentId=storageloc)
synfile<-synStore(synfile)

idsforrelease<-sort(unique(c(training$subject_id, ancillary$subject_id)))


#
# Read in Clinical and Survey data
#

survey<-read.delim("Other/Home-based_validation_csv_export_20181129160617/Home-based_validation_Uw_ervaringen_export_20181129.csv", sep=";", header=T, as.is=T)
survey<-survey[survey$Castor.Record.ID%in%idsforrelease,]
survey<-survey[,c("Castor.Record.ID", "User_experience_23", "Smartphon_position_other")]
names(survey)<-c("subject_id", "Most_Common_Smartphone_Location", "Smartphone_Location_Other")
survey$Most_Common_Smartphone_Location[survey$Most_Common_Smartphone_Location==1]<-"Front pocket"
survey$Most_Common_Smartphone_Location[survey$Most_Common_Smartphone_Location==2]<-"Back pocket"
survey$Most_Common_Smartphone_Location[survey$Most_Common_Smartphone_Location==3]<-"Jacket"
survey$Most_Common_Smartphone_Location[survey$Most_Common_Smartphone_Location==4]<-"Shoulder bag"
survey$Most_Common_Smartphone_Location[survey$Most_Common_Smartphone_Location==5]<-"Backpack"
survey$Most_Common_Smartphone_Location[survey$Most_Common_Smartphone_Location==6]<-"Handbag"
survey$Most_Common_Smartphone_Location[survey$Most_Common_Smartphone_Location==7]<-"Other"
survey<-survey[order(survey$subject_id),]

clin<-read.delim("Other/Home-based_validation_csv_export_20181129160617/Home-based_validation_export_20181129.csv", sep=";", header=T, as.is=T)
names(clin)[1]<-"subject_id"
clin<-clin[clin$subject_id%in%idsforrelease,]
clin<-clin[,-c(2:3)]

demos<-clin[,c("subject_id", "date_screening", "gender", "year_of_birth")]
demos$Age<-year(dmy(demos$date_screening))-demos$year_of_birth
demos<-demos[,c("subject_id", "Age", "gender")]
names(demos)[3]<-"Gender"
demos$Gender[demos$Gender==0]<-"Male"
demos$Gender[demos$Gender==1]<-"Female"
write.csv(demos, "../Data Release/REAL-PD_Demographics.csv", row.names=F, quote=F)
synfile<-File("../Data Release/REAL-PD_Demographics.csv", parentId=storageloc)
synfile<-synStore(synfile)


metadat<-clin[,c("subject_id", "smartphone_brand", "smartphone_model", "android_version")]
metadat<-merge.data.frame(metadat, survey, by="subject_id", all=T)
write.csv(metadat, "../Data Release/REAL-PD_Smartphone_Metadata.csv", row.names=F, quote=F)
synfile<-File("../Data Release/REAL-PD_Smartphone_Metadata.csv", parentId=storageloc)
synfile<-synStore(synfile)



names(clin)<-gsub("UPDRS_3_22", "UPDRS_3_18", names(clin))

tmp<-names(clin)
offupdrsnames<-66:101
onupdrsnames<-198:233

tmp2<-data.frame(clin[,c(1, offupdrsnames)], ParticipantState="Off", stringsAsFactors = F)
names(tmp2)<-gsub("OFF_", "", names(tmp2))
tmp3<-data.frame(clin[,c(1, onupdrsnames)], ParticipantState="On", stringsAsFactors = F)
names(tmp3)<-gsub("ON_", "", names(tmp3))

updrs<-rbind.fill(tmp2,tmp3)
updrs<-updrs[,c(1,dim(updrs)[2], 2:(dim(updrs)[2]-1))]

cisupdrsnames<-unlist(read.csv("../Data Release/CIS-PD_updrsIII_names.csv", header=T, as.is=T))
cisupdrsnames<-cisupdrsnames[-2]
cisupdrsnames<-cisupdrsnames[-c(3:6)]

names(updrs)<-cisupdrsnames
updrs[updrs==-99]<-NA

write.csv(updrs, "../Data Release/REAL-PD_UPDRS_Part3.csv", row.names=F, quote=F)
synfile<-File("../Data Release/REAL-PD_UPDRS_Part3.csv", parentId=storageloc)
synfile<-synStore(synfile)

write.csv(metadat, "../Data Release/REAL-PD_Smartphone_Metadata.csv", row.names=F, quote=F)
synfile<-File("../Data Release/REAL-PD_Smartphone_Metadata.csv", parentId=storageloc)
synfile<-synStore(synfile)


updrs1.<-read.delim("Other/Home-based_validation_csv_export_20181129160617/Home-based_validation_Klachten_in_het_dagelijks_leven_export_20181129.csv", sep=";", header=T, as.is=T)
names(updrs1.)[2]<-"subject_id"
updrs1.<-updrs1.[updrs1.$subject_id%in%idsforrelease,]
#table(unlist(updrs1.[,grep("MDS_UPDRS", names(updrs1.))]))
updrs1.$UPDRS_PartI_Total<-apply(updrs1.[,grep("MDS_UPDRS_1", names(updrs1.))],1,sum)
updrs1.$UPDRS_PartII_Total<-apply(updrs1.[,grep("MDS_UPDRS_2", names(updrs1.))],1,sum)
updrs1.<-updrs1.[,c("subject_id", "UPDRS_PartI_Total", "UPDRS_PartII_Total")]

u1names<-grep("UPDRS_1", names(clin))
u4names<-grep("UPDRS_4", names(clin))
updrs4.<-clin[,c(1, u1names, u4names)]
updrs4.<-updrs4.[,-grep("UPDRS_1A", names(updrs4.))]
updrs4.$U1.1tot<-apply(updrs4.[,grep("UPDRS_1", names(updrs4.))],1,sum)
updrs4.<-updrs4.[,-grep("UPDRS_1", names(updrs4.))]

updrs1.<-merge.data.frame(updrs1., updrs4., by=1, all=T)
updrs1.$UPDRS_PartI_Total<-updrs1.$UPDRS_PartI_Total+updrs1.$U1.1tot
updrs1.<-updrs1.[,-dim(updrs1.)]
substr(names(updrs1.)[grep("4", names(updrs1.))],8,8)<-"."

write.csv(updrs1., "../Data Release/REAL-PD_UPDRS_Part1_2_4.csv", row.names=F, quote=F)
synfile<-File("../Data Release/REAL-PD_UPDRS_Part1_2_4.csv", parentId=storageloc)
synfile<-synStore(synfile)


