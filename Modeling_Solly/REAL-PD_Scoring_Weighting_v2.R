setwd("~/Documents/PDDB Challenge 2/REAL-PD")

library(dplyr)
library(tidyr)
library(ggplot2)

files<-list.files("Analyses/Features/real-solly-nick/", full.names = TRUE)
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

on_off_res<-on_off_res %>% transmute(Model=model, split_id=split_id+1, n.MSE=nmse_null-nmse,
                                     n.MAE=nmae_null-nmae, sqrt.MSE=sqrtnmse_null-sqrtnmse,
                                     sqrt.MAE=sqrtnmae_null-sqrtnmae, log.MSE=lognmse_null-lognmse,
                                     log.MAE=lognmae_null-lognmae)

on_off_res <- gather(on_off_res, stat, lift, n.MSE:log.MAE) %>% separate(stat, c("Weight", "Stat"))

on_off_res$Weight<-factor(on_off_res$Weight, levels = c("log", "sqrt", "n"), ordered = TRUE)

ggplot(on_off_res)+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight)+ theme(axis.text.x = element_text(angle = 90))
ggsave("Analyses/Plots/REAL-PD_on_off_weighting_Alexs_models.png", height=5.5, width=10)


ggplot(on_off_res[on_off_res$Model!="regress-mlp",])+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight)+ theme(axis.text.x = element_text(angle = 90))
ggsave("Analyses/Plots/REAL-PD_on_off_weighting_Alexs_models_v2.png", height=5.5, width=10)




on_off_summary<-aggregate(on_off_res$lift, by=list(Stat=on_off_res$Stat, Weight=on_off_res$Weight, Model=on_off_res$Model), FUN="summary")
on_off_summary<-on_off_summary[order(on_off_summary$Stat, on_off_summary$Weight, on_off_summary$Model),]
names(on_off_summary)<-gsub("x.", "", names(on_off_summary))

write.csv(on_off_summary, "Analyses/REAL-PD_on_off_weighting_Alexs_models_v2.csv", row.names=F, quote=F)



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

dyskinesia_res<-dyskinesia_res %>% transmute(Model=model, split_id=split_id+1, n.MSE=nmse_null-nmse,
                                             n.MAE=nmae_null-nmae, sqrt.MSE=sqrtnmse_null-sqrtnmse,
                                             sqrt.MAE=sqrtnmae_null-sqrtnmae, log.MSE=lognmse_null-lognmse,
                                             log.MAE=lognmae_null-lognmae)

dyskinesia_res <- gather(dyskinesia_res, stat, lift, n.MSE:log.MAE) %>% separate(stat, c("Weight", "Stat"))

dyskinesia_res$Weight<-factor(dyskinesia_res$Weight, levels = c("log", "sqrt", "n"), ordered = TRUE)

ggplot(dyskinesia_res)+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight)+ theme(axis.text.x = element_text(angle = 90))
ggsave("Analyses/Plots/REAL-PD_dyskinesia_weighting_Alexs_models.png", height=5.5, width=10)

ggplot(dyskinesia_res[dyskinesia_res$Model!="regress-mlp",])+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight)+ theme(axis.text.x = element_text(angle = 90))
ggsave("Analyses/Plots/REAL-PD_dyskinesia_weighting_Alexs_models_v2.png", height=5.5, width=10)



dyskinesia_summary<-aggregate(dyskinesia_res$lift, by=list(Stat=dyskinesia_res$Stat, Weight=dyskinesia_res$Weight, Model=dyskinesia_res$Model), FUN="summary")
dyskinesia_summary<-dyskinesia_summary[order(dyskinesia_summary$Stat, dyskinesia_summary$Weight, dyskinesia_summary$Model),]
names(dyskinesia_summary)<-gsub("x.", "", names(dyskinesia_summary))

write.csv(dyskinesia_summary, "Analyses/REAL-PD_dyskinesia_weighting_Alexs_models_v2.csv", row.names=F, quote=F)


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

tremor_res<-tremor_res %>% transmute(Model=model, split_id=split_id+1, n.MSE=nmse_null-nmse,
                                     n.MAE=nmae_null-nmae, sqrt.MSE=sqrtnmse_null-sqrtnmse,
                                     sqrt.MAE=sqrtnmae_null-sqrtnmae, log.MSE=lognmse_null-lognmse,
                                     log.MAE=lognmae_null-lognmae)

tremor_res <- gather(tremor_res, stat, lift, n.MSE:log.MAE) %>% separate(stat, c("Weight", "Stat"))

tremor_res$Weight<-factor(tremor_res$Weight, levels = c("log", "sqrt", "n"), ordered = TRUE)

ggplot(tremor_res)+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight)+ theme(axis.text.x = element_text(angle = 90))
ggsave("Analyses/Plots/REAL-PD_tremor_weighting_Alexs_models.png", height=5.5, width=10)

ggplot(tremor_res[tremor_res$Model!="regress-mlp",])+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight)+ theme(axis.text.x = element_text(angle = 90))
ggsave("Analyses/Plots/REAL-PD_tremor_weighting_Alexs_models_v2.png", height=5.5, width=10)




tremor_summary<-aggregate(tremor_res$lift, by=list(Stat=tremor_res$Stat, Weight=tremor_res$Weight, Model=tremor_res$Model), FUN="summary")
tremor_summary<-tremor_summary[order(tremor_summary$Stat, tremor_summary$Weight, tremor_summary$Model),]
names(tremor_summary)<-gsub("x.", "", names(tremor_summary))

write.csv(tremor_summary, "Analyses/REAL-PD_tremor_weighting_Alexs_models_v2.csv", row.names=F, quote=F)

save.image("../Code/REAL-PD Modeling/Modeling_Solly/REAL-PD_Scoring_Weighting_v2.Rdata")



mod<-1

ggplot(on_off_res[on_off_res$Model!="regress-mlp",])+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight) + geom_point(data=filter(on_off_res[on_off_res$Model!="regress-mlp",], split_id==mod), aes(x=Model, y=lift), position=position_dodge(width=0.75), col=2)+ theme(axis.text.x = element_text(angle = 90))+ geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/REAL-PD_on_off_weighting_Alexs_models_1_v2.png", height=5.5, width=10)
ggplot(dyskinesia_res[dyskinesia_res$Model!="regress-mlp",])+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight) + geom_point(data=filter(dyskinesia_res[dyskinesia_res$Model!="regress-mlp",], split_id==mod), aes(x=Model, y=lift), position=position_dodge(width=0.75), col=2)+ theme(axis.text.x = element_text(angle = 90))+ geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/REAL-PD_dyskinesia_weighting_Alexs_models_1_v2.png", height=5.5, width=10)
ggplot(tremor_res[tremor_res$Model!="regress-mlp",])+geom_boxplot(aes(x=Model, y=lift))+facet_grid(Stat~Weight) + geom_point(data=filter(tremor_res[tremor_res$Model!="regress-mlp",], split_id==mod), aes(x=Model, y=lift), position=position_dodge(width=0.75), col=2)+ theme(axis.text.x = element_text(angle = 90))+ geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/REAL-PD_tremor_weighting_Alexs_models_1_v2.png", height=5.5, width=10)



