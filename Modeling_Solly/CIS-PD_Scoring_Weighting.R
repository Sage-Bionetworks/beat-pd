setwd("~/Documents/PDDB Challenge 2/CIS-PD")

library(dplyr)
library(tidyr)
library(ggplot2)

files<-list.files("Analyses/Features/solly-splits/", full.names = TRUE)
files<-files[grep("csv",files)]

#On/Off

on_off_res<-NULL
for(i in files[grep("on_off", files)]){
  data<-read.csv(i, header = TRUE, as.is=TRUE)
  res <- data %>% mutate(nwtdmse=n_test*mse, nwtdmse_null=n_test*null_mse, nwtdmae=n_test*mae, nwtdmae_null=n_test*null_mae,
                  sqrtnwtdmse=sqrt(n_test)*mse, sqrtnwtdmse_null=sqrt(n_test)*null_mse, sqrtnwtdmae=sqrt(n_test)*mae, sqrtnwtdmae_null=sqrt(n_test)*null_mae,
                  lognwtdmse=log(n_test)*mse, lognwtdmse_null=log(n_test)*null_mse, lognwtdmae=log(n_test)*mae, lognwtdmae_null=log(n_test)*null_mae) %>%
    group_by(split_id) %>% summarize(nmse=sum(nwtdmse)/sum(n_test), nmse_null=sum(nwtdmse_null)/sum(n_test), nmae=sum(nwtdmae)/sum(n_test), nmae_null=sum(nwtdmae_null)/sum(n_test),
                                     sqrtnmse=sum(sqrtnwtdmse)/sum(sqrt(n_test)), sqrtnmse_null=sum(sqrtnwtdmse_null)/sum(sqrt(n_test)), sqrtnmae=sum(sqrtnwtdmae)/sum(sqrt(n_test)), sqrtnmae_null=sum(sqrtnwtdmae_null)/sum(sqrt(n_test)),
                                     lognmse=sum(lognwtdmse)/sum(log(n_test)), lognmse_null=sum(lognwtdmse_null)/sum(log(n_test)), lognmae=sum(lognwtdmae)/sum(log(n_test)), lognmae_null=sum(lognwtdmae_null)/sum(log(n_test)))
  res<-data.frame(model=gsub("_on_off.csv", "", strsplit(i, "//")[[1]][2]), res, stringsAsFactors = FALSE)
  on_off_res<-rbind(on_off_res, res)
}

on_off_res<-on_off_res %>% transmute(Model=gsub("classif-", "", model), split_id=split_id+1, n.MSE=nmse_null-nmse,
                                  n.MAE=nmae_null-nmae, sqrt.MSE=sqrtnmse_null-sqrtnmse,
                                  sqrt.MAE=sqrtnmae_null-sqrtnmae, log.MSE=lognmse_null-lognmse,
                                  log.MAE=lognmae_null-lognmae)

on_off_res <- gather(on_off_res, stat, lift, n.MSE:log.MAE) %>% separate(stat, c("Weight", "Stat"))

on_off_res$Weight<-factor(on_off_res$Weight, levels = c("log", "sqrt", "n"), ordered = TRUE)

ggplot(on_off_res)+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight)
ggsave("Analyses/Plots/CIS-PD_on_off_weighting_Alexs_models.png", height=5.5, width=10)

on_off_summary<-aggregate(on_off_res$lift, by=list(Stat=on_off_res$Stat, Weight=on_off_res$Weight, Model=on_off_res$Model), FUN="summary")
on_off_summary<-on_off_summary[order(on_off_summary$Stat, on_off_summary$Weight, on_off_summary$Model),]
names(on_off_summary)<-gsub("x.", "", names(on_off_summary))

write.csv(on_off_summary, "Analyses/CIS-PD_on_off_weighting_Alexs_models.csv", row.names=F, quote=F)



#Dyskinesia

dyskinesia_res<-NULL
for(i in files[grep("dyskinesia", files)]){
  data<-read.csv(i, header = TRUE, as.is=TRUE)
  res <- data %>% mutate(nwtdmse=n_test*mse, nwtdmse_null=n_test*null_mse, nwtdmae=n_test*mae, nwtdmae_null=n_test*null_mae,
                         sqrtnwtdmse=sqrt(n_test)*mse, sqrtnwtdmse_null=sqrt(n_test)*null_mse, sqrtnwtdmae=sqrt(n_test)*mae, sqrtnwtdmae_null=sqrt(n_test)*null_mae,
                         lognwtdmse=log(n_test)*mse, lognwtdmse_null=log(n_test)*null_mse, lognwtdmae=log(n_test)*mae, lognwtdmae_null=log(n_test)*null_mae) %>%
    group_by(split_id) %>% summarize(nmse=sum(nwtdmse)/sum(n_test), nmse_null=sum(nwtdmse_null)/sum(n_test), nmae=sum(nwtdmae)/sum(n_test), nmae_null=sum(nwtdmae_null)/sum(n_test),
                                     sqrtnmse=sum(sqrtnwtdmse)/sum(sqrt(n_test)), sqrtnmse_null=sum(sqrtnwtdmse_null)/sum(sqrt(n_test)), sqrtnmae=sum(sqrtnwtdmae)/sum(sqrt(n_test)), sqrtnmae_null=sum(sqrtnwtdmae_null)/sum(sqrt(n_test)),
                                     lognmse=sum(lognwtdmse)/sum(log(n_test)), lognmse_null=sum(lognwtdmse_null)/sum(log(n_test)), lognmae=sum(lognwtdmae)/sum(log(n_test)), lognmae_null=sum(lognwtdmae_null)/sum(log(n_test)))
  res<-data.frame(model=gsub("_dyskinesia.csv", "", strsplit(i, "//")[[1]][2]), res, stringsAsFactors = FALSE)
  dyskinesia_res<-rbind(dyskinesia_res, res)
}

dyskinesia_res<-dyskinesia_res %>% transmute(Model=gsub("classif-", "", model), split_id=split_id+1, n.MSE=nmse_null-nmse,
                                     n.MAE=nmae_null-nmae, sqrt.MSE=sqrtnmse_null-sqrtnmse,
                                     sqrt.MAE=sqrtnmae_null-sqrtnmae, log.MSE=lognmse_null-lognmse,
                                     log.MAE=lognmae_null-lognmae)

dyskinesia_res <- gather(dyskinesia_res, stat, lift, n.MSE:log.MAE) %>% separate(stat, c("Weight", "Stat"))

dyskinesia_res$Weight<-factor(dyskinesia_res$Weight, levels = c("log", "sqrt", "n"), ordered = TRUE)

ggplot(dyskinesia_res)+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight)
ggsave("Analyses/Plots/CIS-PD_dyskinesia_weighting_Alexs_models.png", height=5.5, width=10)

dyskinesia_summary<-aggregate(dyskinesia_res$lift, by=list(Stat=dyskinesia_res$Stat, Weight=dyskinesia_res$Weight, Model=dyskinesia_res$Model), FUN="summary")
dyskinesia_summary<-dyskinesia_summary[order(dyskinesia_summary$Stat, dyskinesia_summary$Weight, dyskinesia_summary$Model),]
names(dyskinesia_summary)<-gsub("x.", "", names(dyskinesia_summary))

write.csv(dyskinesia_summary, "Analyses/CIS-PD_dyskinesia_weighting_Alexs_models.csv", row.names=F, quote=F)



#Tremor

tremor_res<-NULL
for(i in files[grep("tremor", files)]){
  data<-read.csv(i, header = TRUE, as.is=TRUE)
  res <- data %>% mutate(nwtdmse=n_test*mse, nwtdmse_null=n_test*null_mse, nwtdmae=n_test*mae, nwtdmae_null=n_test*null_mae,
                         sqrtnwtdmse=sqrt(n_test)*mse, sqrtnwtdmse_null=sqrt(n_test)*null_mse, sqrtnwtdmae=sqrt(n_test)*mae, sqrtnwtdmae_null=sqrt(n_test)*null_mae,
                         lognwtdmse=log(n_test)*mse, lognwtdmse_null=log(n_test)*null_mse, lognwtdmae=log(n_test)*mae, lognwtdmae_null=log(n_test)*null_mae) %>%
    group_by(split_id) %>% summarize(nmse=sum(nwtdmse)/sum(n_test), nmse_null=sum(nwtdmse_null)/sum(n_test), nmae=sum(nwtdmae)/sum(n_test), nmae_null=sum(nwtdmae_null)/sum(n_test),
                                     sqrtnmse=sum(sqrtnwtdmse)/sum(sqrt(n_test)), sqrtnmse_null=sum(sqrtnwtdmse_null)/sum(sqrt(n_test)), sqrtnmae=sum(sqrtnwtdmae)/sum(sqrt(n_test)), sqrtnmae_null=sum(sqrtnwtdmae_null)/sum(sqrt(n_test)),
                                     lognmse=sum(lognwtdmse)/sum(log(n_test)), lognmse_null=sum(lognwtdmse_null)/sum(log(n_test)), lognmae=sum(lognwtdmae)/sum(log(n_test)), lognmae_null=sum(lognwtdmae_null)/sum(log(n_test)))
  res<-data.frame(model=gsub("_tremor.csv", "", strsplit(i, "//")[[1]][2]), res, stringsAsFactors = FALSE)
  tremor_res<-rbind(tremor_res, res)
}

tremor_res<-tremor_res %>% transmute(Model=gsub("classif-", "", model), split_id=split_id+1, n.MSE=nmse_null-nmse,
                                     n.MAE=nmae_null-nmae, sqrt.MSE=sqrtnmse_null-sqrtnmse,
                                     sqrt.MAE=sqrtnmae_null-sqrtnmae, log.MSE=lognmse_null-lognmse,
                                     log.MAE=lognmae_null-lognmae)

tremor_res <- gather(tremor_res, stat, lift, n.MSE:log.MAE) %>% separate(stat, c("Weight", "Stat"))

tremor_res$Weight<-factor(tremor_res$Weight, levels = c("log", "sqrt", "n"), ordered = TRUE)

ggplot(tremor_res)+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight)
ggsave("Analyses/Plots/CIS-PD_tremor_weighting_Alexs_models.png", height=5.5, width=10)

tremor_summary<-aggregate(tremor_res$lift, by=list(Stat=tremor_res$Stat, Weight=tremor_res$Weight, Model=tremor_res$Model), FUN="summary")
tremor_summary<-tremor_summary[order(tremor_summary$Stat, tremor_summary$Weight, tremor_summary$Model),]
names(tremor_summary)<-gsub("x.", "", names(tremor_summary))

write.csv(tremor_summary, "Analyses/CIS-PD_tremor_weighting_Alexs_models.csv", row.names=F, quote=F)

save.image("../Code/CIS-PD Modeling/Modeling_Solly/CIS-PD_Scoring_Weighting.Rdata")


#
# Try Diskinesia w/o 1044
#

ftr<-c(1044)

dyskinesia_res_filter<-NULL
for(i in files[grep("dyskinesia", files)]){
  data<-read.csv(i, header = TRUE, as.is=TRUE)
  res <- data %>% filter(!subject_id %in% ftr) %>% mutate(nwtdmse=n_test*mse, nwtdmse_null=n_test*null_mse, nwtdmae=n_test*mae, nwtdmae_null=n_test*null_mae,
                         sqrtnwtdmse=sqrt(n_test)*mse, sqrtnwtdmse_null=sqrt(n_test)*null_mse, sqrtnwtdmae=sqrt(n_test)*mae, sqrtnwtdmae_null=sqrt(n_test)*null_mae,
                         lognwtdmse=log(n_test)*mse, lognwtdmse_null=log(n_test)*null_mse, lognwtdmae=log(n_test)*mae, lognwtdmae_null=log(n_test)*null_mae) %>%
    group_by(split_id) %>% summarize(nmse=sum(nwtdmse)/sum(n_test), nmse_null=sum(nwtdmse_null)/sum(n_test), nmae=sum(nwtdmae)/sum(n_test), nmae_null=sum(nwtdmae_null)/sum(n_test),
                                     sqrtnmse=sum(sqrtnwtdmse)/sum(sqrt(n_test)), sqrtnmse_null=sum(sqrtnwtdmse_null)/sum(sqrt(n_test)), sqrtnmae=sum(sqrtnwtdmae)/sum(sqrt(n_test)), sqrtnmae_null=sum(sqrtnwtdmae_null)/sum(sqrt(n_test)),
                                     lognmse=sum(lognwtdmse)/sum(log(n_test)), lognmse_null=sum(lognwtdmse_null)/sum(log(n_test)), lognmae=sum(lognwtdmae)/sum(log(n_test)), lognmae_null=sum(lognwtdmae_null)/sum(log(n_test)))
  res<-data.frame(model=gsub("_dyskinesia.csv", "", strsplit(i, "//")[[1]][2]), res, stringsAsFactors = FALSE)
  dyskinesia_res_filter<-rbind(dyskinesia_res_filter, res)
}


dyskinesia_res_filter<-dyskinesia_res_filter %>% 
  transmute(Model=gsub("classif-", "", model), split_id=split_id+1, n.MSE=nmse_null-nmse,
            n.MAE=nmae_null-nmae, sqrt.MSE=sqrtnmse_null-sqrtnmse,
            sqrt.MAE=sqrtnmae_null-sqrtnmae, log.MSE=lognmse_null-lognmse,
            log.MAE=lognmae_null-lognmae)

dyskinesia_res_filter <- gather(dyskinesia_res_filter, stat, lift, n.MSE:log.MAE) %>% separate(stat, c("Weight", "Stat"))

dyskinesia_res_filter$Weight<-factor(dyskinesia_res$Weight, levels = c("log", "sqrt", "n"), ordered = TRUE)

dyskinesia_summary_filter<-aggregate(dyskinesia_res_filter$lift, by=list(Stat=dyskinesia_res_filter$Stat, Weight=dyskinesia_res_filter$Weight, Model=dyskinesia_res_filter$Model), FUN="summary")
dyskinesia_summary_filter<-dyskinesia_summary_filter[order(dyskinesia_summary_filter$Stat, dyskinesia_summary_filter$Weight, dyskinesia_summary_filter$Model),]
names(dyskinesia_summary_filter)<-gsub("x.", "", names(dyskinesia_summary_filter))



#
# Plot point for Specific Model
#
mod<-1

ggplot(on_off_res)+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight) + geom_point(data=filter(on_off_res, split_id==mod), aes(x=Model, y=lift), position=position_dodge(width=0.75), col=2)
ggsave("Analyses/Plots/CIS-PD_on_off_weighting_Alexs_models_1.png", height=5.5, width=10)
ggplot(dyskinesia_res)+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight) + geom_point(data=filter(dyskinesia_res, split_id==mod), aes(x=Model, y=lift), position=position_dodge(width=0.75), col=2)
ggsave("Analyses/Plots/CIS-PD_dyskinesia_weighting_Alexs_models_1.png", height=5.5, width=10)
ggplot(tremor_res)+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight) + geom_point(data=filter(tremor_res, split_id==mod), aes(x=Model, y=lift), position=position_dodge(width=0.75), col=2)
ggsave("Analyses/Plots/CIS-PD_tremor_weighting_Alexs_models_1.png", height=5.5, width=10)



#
# Percent Lift
#

#On/Off

on_off_res_pct<-NULL
for(i in files[grep("on_off", files)]){
  data<-read.csv(i, header = TRUE, as.is=TRUE)
  res <- data %>% mutate(nwtdmse=n_test*(null_mse-mse)/null_mse, nwtdmae=n_test*(null_mae-mae)/null_mae, 
                         sqrtnwtdmse=sqrt(n_test)*(null_mse-mse)/null_mse, sqrtnwtdmae=sqrt(n_test)*(null_mae-mae)/null_mae,
                         lognwtdmse=log(n_test)*(null_mse-mse)/null_mse, lognwtdmae=log(n_test)*(null_mae-mae)/null_mae) %>%
    group_by(split_id) %>% summarize(nmse=sum(nwtdmse)/sum(n_test), nmae=sum(nwtdmae)/sum(n_test), 
                                     sqrtnmse=sum(sqrtnwtdmse)/sum(sqrt(n_test)), sqrtnmae=sum(sqrtnwtdmae)/sum(sqrt(n_test)),
                                     lognmse=sum(lognwtdmse)/sum(log(n_test)), lognmae=sum(lognwtdmae)/sum(log(n_test)))
  res<-data.frame(model=gsub("_on_off.csv", "", strsplit(i, "//")[[1]][2]), res, stringsAsFactors = FALSE)
  on_off_res_pct<-rbind(on_off_res_pct, res)
}

on_off_res_pct<-on_off_res_pct %>% transmute(Model=gsub("classif-", "", model), split_id=split_id+1, n.MSE=nmse,
                                     n.MAE=nmae, sqrt.MSE=sqrtnmse,
                                     sqrt.MAE=sqrtnmae, log.MSE=lognmse,
                                     log.MAE=lognmae)

on_off_res_pct <- gather(on_off_res_pct, stat, lift, n.MSE:log.MAE) %>% separate(stat, c("Weight", "Stat"))

on_off_res_pct$Weight<-factor(on_off_res_pct$Weight, levels = c("log", "sqrt", "n"), ordered = TRUE)

ggplot(on_off_res_pct)+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight)+ylab("Percent Lift")+ geom_point(data=filter(on_off_res_pct, split_id==mod), aes(x=Model, y=lift), position=position_dodge(width=0.75), col=2)
ggsave("Analyses/Plots/CIS-PD_on_off_weighting_pct_Alexs_models.png", height=5.5, width=10)

on_off_summary_pct<-aggregate(on_off_res_pct$lift, by=list(Stat=on_off_res_pct$Stat, Weight=on_off_res_pct$Weight, Model=on_off_res_pct$Model), FUN="summary")
on_off_summary_pct<-on_off_summary_pct[order(on_off_summary_pct$Stat, on_off_summary_pct$Weight, on_off_summary_pct$Model),]
names(on_off_summary_pct)<-gsub("x.", "", names(on_off_summary_pct))

write.csv(on_off_summary, "Analyses/CIS-PD_on_off_weighting_pct_Alexs_models.csv", row.names=F, quote=F)

#Dyskinesia

dyskinesia_res_pct<-NULL
for(i in files[grep("dyskinesia", files)]){
  data<-read.csv(i, header = TRUE, as.is=TRUE)
  res <- data %>% mutate(nwtdmse=n_test*(null_mse-mse)/null_mse, nwtdmae=n_test*(null_mae-mae)/null_mae, 
                         sqrtnwtdmse=sqrt(n_test)*(null_mse-mse)/null_mse, sqrtnwtdmae=sqrt(n_test)*(null_mae-mae)/null_mae,
                         lognwtdmse=log(n_test)*(null_mse-mse)/null_mse, lognwtdmae=log(n_test)*(null_mae-mae)/null_mae) %>%
    group_by(split_id) %>% summarize(nmse=sum(nwtdmse)/sum(n_test), nmae=sum(nwtdmae)/sum(n_test), 
                                     sqrtnmse=sum(sqrtnwtdmse)/sum(sqrt(n_test)), sqrtnmae=sum(sqrtnwtdmae)/sum(sqrt(n_test)),
                                     lognmse=sum(lognwtdmse)/sum(log(n_test)), lognmae=sum(lognwtdmae)/sum(log(n_test)))
  res<-data.frame(model=gsub("_dyskinesia.csv", "", strsplit(i, "//")[[1]][2]), res, stringsAsFactors = FALSE)
  dyskinesia_res_pct<-rbind(dyskinesia_res_pct, res)
}

dyskinesia_res_pct<-dyskinesia_res_pct %>% transmute(Model=gsub("classif-", "", model), split_id=split_id+1, n.MSE=nmse,
                                             n.MAE=nmae, sqrt.MSE=sqrtnmse,
                                             sqrt.MAE=sqrtnmae, log.MSE=lognmse,
                                             log.MAE=lognmae)

dyskinesia_res_pct <- gather(dyskinesia_res_pct, stat, lift, n.MSE:log.MAE) %>% separate(stat, c("Weight", "Stat"))

dyskinesia_res_pct$Weight<-factor(dyskinesia_res_pct$Weight, levels = c("log", "sqrt", "n"), ordered = TRUE)

ggplot(dyskinesia_res_pct)+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight)+ylab("Percent Lift")+ geom_point(data=filter(dyskinesia_res_pct, split_id==mod), aes(x=Model, y=lift), position=position_dodge(width=0.75), col=2)
ggsave("Analyses/Plots/CIS-PD_dyskinesia_weighting_pct_Alexs_models.png", height=5.5, width=10)

dyskinesia_summary_pct<-aggregate(dyskinesia_res_pct$lift, by=list(Stat=dyskinesia_res_pct$Stat, Weight=dyskinesia_res_pct$Weight, Model=dyskinesia_res_pct$Model), FUN="summary")
dyskinesia_summary_pct<-dyskinesia_summary_pct[order(dyskinesia_summary_pct$Stat, dyskinesia_summary_pct$Weight, dyskinesia_summary_pct$Model),]
names(dyskinesia_summary_pct)<-gsub("x.", "", names(dyskinesia_summary_pct))

write.csv(dyskinesia_summary, "Analyses/CIS-PD_dyskinesia_weighting_pct_Alexs_models.csv", row.names=F, quote=F)


#Tremor

tremor_res_pct<-NULL
for(i in files[grep("tremor", files)]){
  data<-read.csv(i, header = TRUE, as.is=TRUE)
  res <- data %>% mutate(nwtdmse=n_test*(null_mse-mse)/null_mse, nwtdmae=n_test*(null_mae-mae)/null_mae, 
                         sqrtnwtdmse=sqrt(n_test)*(null_mse-mse)/null_mse, sqrtnwtdmae=sqrt(n_test)*(null_mae-mae)/null_mae,
                         lognwtdmse=log(n_test)*(null_mse-mse)/null_mse, lognwtdmae=log(n_test)*(null_mae-mae)/null_mae) %>%
    group_by(split_id) %>% summarize(nmse=sum(nwtdmse)/sum(n_test), nmae=sum(nwtdmae)/sum(n_test), 
                                     sqrtnmse=sum(sqrtnwtdmse)/sum(sqrt(n_test)), sqrtnmae=sum(sqrtnwtdmae)/sum(sqrt(n_test)),
                                     lognmse=sum(lognwtdmse)/sum(log(n_test)), lognmae=sum(lognwtdmae)/sum(log(n_test)))
  res<-data.frame(model=gsub("_tremor.csv", "", strsplit(i, "//")[[1]][2]), res, stringsAsFactors = FALSE)
  tremor_res_pct<-rbind(tremor_res_pct, res)
}

tremor_res_pct<-tremor_res_pct %>% transmute(Model=gsub("classif-", "", model), split_id=split_id+1, n.MSE=nmse,
                                             n.MAE=nmae, sqrt.MSE=sqrtnmse,
                                             sqrt.MAE=sqrtnmae, log.MSE=lognmse,
                                             log.MAE=lognmae)

tremor_res_pct <- gather(tremor_res_pct, stat, lift, n.MSE:log.MAE) %>% separate(stat, c("Weight", "Stat"))

tremor_res_pct$Weight<-factor(tremor_res_pct$Weight, levels = c("log", "sqrt", "n"), ordered = TRUE)

ggplot(tremor_res_pct)+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight)+ylab("Percent Lift")+ geom_point(data=filter(tremor_res_pct, split_id==mod), aes(x=Model, y=lift), position=position_dodge(width=0.75), col=2)
ggsave("Analyses/Plots/CIS-PD_tremor_weighting_pct_Alexs_models.png", height=5.5, width=10)

tremor_summary_pct<-aggregate(tremor_res_pct$lift, by=list(Stat=tremor_res_pct$Stat, Weight=tremor_res_pct$Weight, Model=tremor_res_pct$Model), FUN="summary")
tremor_summary_pct<-tremor_summary_pct[order(tremor_summary_pct$Stat, tremor_summary_pct$Weight, tremor_summary_pct$Model),]
names(tremor_summary_pct)<-gsub("x.", "", names(tremor_summary_pct))

write.csv(tremor_summary, "Analyses/CIS-PD_tremor_weighting_pct_Alexs_models.csv", row.names=F, quote=F)



#
# Plot by Individual
#


on_off_res_ind<-NULL
for(i in files[grep("on_off", files)]){
  data<-read.csv(i, header = TRUE, as.is=TRUE)
  res <- data %>% transmute(subject_id=subject_id, split_id=split_id+1, n_total=n_total, n_train=n_train, n_test=n_test,
                            MSE=null_mse-mse, MAE=null_mae-mae) 
   res<-data.frame(model=gsub("classif-", "", gsub("_on_off.csv", "", strsplit(i, "//")[[1]][2])), res, stringsAsFactors = FALSE)
  on_off_res_ind<-rbind(on_off_res_ind, res)
}

on_off_res_ind<-pivot_longer(on_off_res_ind, cols=MSE:MAE, names_to="Statistic", values_to = "Lift")


dyskinesia_res_ind<-NULL
for(i in files[grep("dyskinesia", files)]){
  data<-read.csv(i, header = TRUE, as.is=TRUE)
  res <- data %>% transmute(subject_id=subject_id, split_id=split_id+1, n_total=n_total, n_train=n_train, n_test=n_test,
                            MSE=null_mse-mse, MAE=null_mae-mae) 
  res<-data.frame(model=gsub("classif-", "", gsub("_dyskinesia.csv", "", strsplit(i, "//")[[1]][2])), res, stringsAsFactors = FALSE)
  dyskinesia_res_ind<-rbind(dyskinesia_res_ind, res)
}

dyskinesia_res_ind<-pivot_longer(dyskinesia_res_ind, cols=MSE:MAE, names_to="Statistic", values_to = "Lift")


tremor_res_ind<-NULL
for(i in files[grep("tremor", files)]){
  data<-read.csv(i, header = TRUE, as.is=TRUE)
  res <- data %>% transmute(subject_id=subject_id, split_id=split_id+1, n_total=n_total, n_train=n_train, n_test=n_test,
                            MSE=null_mse-mse, MAE=null_mae-mae) 
  res<-data.frame(model=gsub("classif-", "", gsub("_tremor.csv", "", strsplit(i, "//")[[1]][2])), res, stringsAsFactors = FALSE)
  tremor_res_ind<-rbind(tremor_res_ind, res)
}

tremor_res_ind<-pivot_longer(tremor_res_ind, cols=MSE:MAE, names_to="Statistic", values_to = "Lift")

res_ind<-rbind(data.frame(Phenotype="On/off", on_off_res_ind, stringsAsFactors = F),
               data.frame(Phenotype="Dyskinesia", dyskinesia_res_ind, stringsAsFactors = F),
               data.frame(Phenotype="Tremor", tremor_res_ind, stringsAsFactors = F))



 stat_box_data <- function(x, lower_limit = NA){
  return( 
    data.frame(
      y = lower_limit*0.9,
      label = paste(" (", round(mean(x), 1), ")", sep="")
    )
  )
}

# Plot
res_ind$subject_id<-factor(res_ind$subject_id)
res_ind$Phenotype<-factor(res_ind$Phenotype, levels=c("On/off", "Dyskinesia", "Tremor"), ordered = TRUE)
ggplot(subset(res_ind, model=="rf"), aes(x=subject_id, y=Lift)) + geom_boxplot() + 
  facet_grid(Statistic~Phenotype) + theme(axis.text.x = element_text(angle = 90)) +
  xlab("Subject ID") + ylab("Lift") + 
  stat_summary(aes(x=subject_id, y=n_total), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90, position = position_dodge(width = 0.8),  fun.args = list(lower_limit=min(res_ind$Lift[res_ind$model=="rf"])*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/CIS-PD_AlexRF_mse_mae.png", width = 8, height = 6.5)

