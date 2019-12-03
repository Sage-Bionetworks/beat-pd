setwd("~/Documents/PDDB Challenge 2/")

library(synapser)

realsynid<-"syn21292051"
cissynid<-"syn21291582"

realsyn<-synGet(realsynid)
real<-read.csv(realsyn$path, header=T, as.is=T)

cissyn<-synGet(cissynid)
cis<-read.csv(cissyn$path, header=T, as.is=T)

onoff<-rbind(cis[,c("measurement_id", "on_off")], real[,c("measurement_id", "on_off")])
onoff<-onoff[!is.na(onoff$on_off),]
onoff$prediction<-NA
onoff<-onoff[,c("measurement_id", "prediction")]

dyskinesia<-rbind(cis[,c("measurement_id", "dyskinesia")], real[,c("measurement_id", "dyskinesia")])
dyskinesia<-dyskinesia[!is.na(dyskinesia$dyskinesia),]
dyskinesia$prediction<-NA
dyskinesia<-dyskinesia[,c("measurement_id", "prediction")]

tremor<-rbind(cis[,c("measurement_id", "tremor")], real[,c("measurement_id", "tremor")])
tremor<-tremor[!is.na(tremor$tremor),]
tremor$prediction<-NA
tremor<-tremor[,c("measurement_id", "prediction")]

synloc<-"syn21344932"

write.csv(onoff, "Data Release/BEAT-PD_SC1_OnOff_Submission_Template.csv", row.names=F, quote=F)
synfile<-File("Data Release/BEAT-PD_SC1_OnOff_Submission_Template.csv", parentId=synloc)
synfile<-synStore(synfile)

write.csv(dyskinesia, "Data Release/BEAT-PD_SC2_Dyskinesia_Submission_Template.csv", row.names=F, quote=F)
synfile<-File("Data Release/BEAT-PD_SC2_Dyskinesia_Submission_Template.csv", parentId=synloc)
synfile<-synStore(synfile)

write.csv(tremor, "Data Release/BEAT-PD_SC3_Tremor_Submission_Template.csv", row.names=F, quote=F)
synfile<-File("Data Release/BEAT-PD_SC3_Tremor_Submission_Template.csv", parentId=synloc)
synfile<-synStore(synfile)
