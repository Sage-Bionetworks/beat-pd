setwd("~/Documents/PDDB Challenge 2/REAL-PD/")

library(caret)
library(synapser)
library(gdata)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)




alltraits<-c("slowness_walking", "tremor", "on_off", "dyskinesia")

#Load labels
scores_query <- synTableQuery("select * from syn20822276")
scores <- scores_query$asDataFrame()
scores<-scores[!apply(is.na(scores[,alltraits]), 1, all),]

# Load features
featphone<-synGet("syn20928641", downloadLocation="Analyses/Features")
featwatchaccel<-synGet("syn20928640", downloadLocation="Analyses/Features")
featwatchgyro<-synGet("syn20928636", downloadLocation="Analyses/Features")

features_phone<-read.csv(featphone$path, header=T, as.is=T)
features_watchaccel<-read.csv(featwatchaccel$path, header=T, as.is=T)
features_watchgyro<-read.csv(featwatchgyro$path, header=T, as.is=T)


# Load Nick's sensor QC summary data
synfile<-synGet("syn21123211", downloadLocation="Analyses/Features/")
msmtqcphone<-read.csv(synfile$path, header=T, as.is=T)

synfile<-synGet("syn21123212", downloadLocation="Analyses/Features/")
msmtqcwatch<-read.csv(synfile$path, header=T, as.is=T)

scores_filtered<-scores[scores$measurement_id%in%features_phone$ID|scores$measurement_id%in%features_watchaccel$ID|scores$measurement_id%in%features_watchgyro$ID,]

rm(features_phone)
rm(features_watchaccel)
rm(features_watchgyro)

# Look for correlation between missingnesss and traits
returncortest<-function(X, Y, method="kendall"){
  if(sum(!is.na(X)&!is.na(Y))>3&var(X, na.rm=T)>0&var(Y, na.rm=T)>0){
    res<-cor.test(X, Y, use="pairwise.complete", method = method)$p.value
  } else {
    res<-NA
  }
  return(res)
}

QCcorphone<-msmtqcphone %>% group_by(subject_id) %>% summarize(subth.on_off=returncortest(Sub.Threshold.Windows, on_off),
                                                     subth.dyskinesia=returncortest(Sub.Threshold.Windows, dyskinesia),
                                                     subth.tremor=returncortest(Sub.Threshold.Windows, tremor),
                                                     consec.on_off=returncortest(Highest...of.Consecutive.Above.Threshold.Windows, on_off),
                                                     consec.dyskinesia=returncortest(Highest...of.Consecutive.Above.Threshold.Windows, dyskinesia),
                                                     consec.tremor=returncortest(Highest...of.Consecutive.Above.Threshold.Windows, tremor))


QCcorwatch<-msmtqcwatch %>% group_by(subject_id) %>% summarize(subth.on_off=returncortest(Sub.Threshold.Windows, on_off),
                                                               subth.dyskinesia=returncortest(Sub.Threshold.Windows, dyskinesia),
                                                               subth.tremor=returncortest(Sub.Threshold.Windows, tremor),
                                                               consec.on_off=returncortest(Highest...of.Consecutive.Above.Threshold.Windows, on_off),
                                                               consec.dyskinesia=returncortest(Highest...of.Consecutive.Above.Threshold.Windows, dyskinesia),
                                                               consec.tremor=returncortest(Highest...of.Consecutive.Above.Threshold.Windows, tremor))


plot_qq<-function(trait, traitname=""){
  trait<-trait[!is.na(trait)]
  n <- length(trait)
  p <- (1 : n) / n - 0.5 / n
  ggplot() + geom_point(aes(x = -log10(p), y = -log10(sort(trait)))) + ylim(0, NA) +
    xlab("-log10(Expected)") + ylab("-log10(Observed)") + geom_abline(intercept = 0, slope = 1, color = "red") +
    geom_hline(yintercept = -log10(0.05/n), color = "black", linetype = "dashed") + ggtitle(traitname)
}

plot_qq(QCcorphone$subth.on_off, "Num Sub-threshold/On-Off")
ggsave("Analyses/Plots/REAL-PD_phone_onoff_vs_numSubThresh_qq.png")
plot_qq(QCcorphone$subth.dyskinesia, "Num Sub-threshold/Dyskinesia")
ggsave("Analyses/Plots/REAL-PD_phone_dyskinesia_vs_numSubThresh_qq.png")
plot_qq(QCcorphone$subth.tremor, "Num Sub-threshold/Tremor")
ggsave("Analyses/Plots/REAL-PD_phone_tremor_vs_numSubThresh_qq.png")
plot_qq(QCcorphone$consec.on_off, "Best Consecutive/On-Off")
ggsave("Analyses/Plots/REAL-PD_phone_onoff_vs_bestConsecutiveWindow_qq.png")
plot_qq(QCcorphone$consec.dyskinesia, "Best Consecutive/Dyskinesia")
ggsave("Analyses/Plots/REAL-PD_phone_dyskinesia_vs_bestConsecutiveWindow_qq.png")
plot_qq(QCcorphone$consec.tremor, "Best Consecutive/Tremor")
ggsave("Analyses/Plots/REAL-PD_phone_tremor_vs_bestConsecutiveWindow_qq.png")


plot_qq(QCcorwatch$subth.on_off, "Num Sub-threshold/On-Off")
ggsave("Analyses/Plots/REAL-PD_watch_onoff_vs_numSubThresh_qq.png")
plot_qq(QCcorwatch$subth.dyskinesia, "Num Sub-threshold/Dyskinesia")
ggsave("Analyses/Plots/REAL-PD_watch_dyskinesia_vs_numSubThresh_qq.png")
plot_qq(QCcorwatch$subth.tremor, "Num Sub-threshold/Tremor")
ggsave("Analyses/Plots/REAL-PD_watch_tremor_vs_numSubThresh_qq.png")
plot_qq(QCcorwatch$consec.on_off, "Best Consecutive/On-Off")
ggsave("Analyses/Plots/REAL-PD_watch_onoff_vs_bestConsecutiveWindow_qq.png")
plot_qq(QCcorwatch$consec.dyskinesia, "Best Consecutive/Dyskinesia")
ggsave("Analyses/Plots/REAL-PD_watch_dyskinesia_vs_bestConsecutiveWindow_qq.png")
plot_qq(QCcorwatch$consec.tremor, "Best Consecutive/Tremor")
ggsave("Analyses/Plots/REAL-PD_watch_tremor_vs_bestConsecutiveWindow_qq.png")





qplot(msmtqcwatch$Sub.Threshold.Windows, geom="histogram", xlab = "# Subthreshold Windows", main = "Watch Measurements")
ggsave("Analyses/Plots/REAL-PD_watch_NumSubThresh_hist.png")
qplot(msmtqcwatch$Highest...of.Consecutive.Above.Threshold.Windows, geom="histogram", xlab = "Longest Above Threshold Window", main = "Watch Measurements")
ggsave("Analyses/Plots/REAL-PD_watch_BestCohsecutiveWindow_hist.png")

qplot(msmtqcphone$Sub.Threshold.Windows, geom="histogram", xlab = "# Subthreshold Windows", main = "Phone Measurements")
ggsave("Analyses/Plots/REAL-PD_phone_NumSubThresh_hist.png")
qplot(msmtqcphone$Highest...of.Consecutive.Above.Threshold.Windows, geom="histogram", xlab = "Longest Above Threshold Window", main = "Phone Measurements")
ggsave("Analyses/Plots/REAL-PD_phone_BestCohsecutiveWindow_hist.png")


# Try various thresholds for data filtering

thrsh<-c(0, 4, 6, 8, 10)

msmtqcmerge<-full_join(msmtqcwatch, msmtqcphone, by=c("measurement_id","subject_id","medication_state","slowness_walking","tremor","on_off","dyskinesia"),
                       suffix=c(".watch", ".phone"))


is_enough_data<-function(x){
  if(sum(!is.na(x))>40 & (sum(table(x)>=10)>1 | sum(table(x)>=5)>2)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


res<-lapply(thrsh, function(x, data){
  res<-data %>% filter((Highest...of.Consecutive.Above.Threshold.Windows.watch >= x&!is.na(Highest...of.Consecutive.Above.Threshold.Windows.watch))|(Highest...of.Consecutive.Above.Threshold.Windows.phone >= x&!is.na(Highest...of.Consecutive.Above.Threshold.Windows.phone)))  %>% group_by(subject_id) %>% summarize(on_off=is_enough_data(on_off), dyskinesia=is_enough_data(dyskinesia), tremor=is_enough_data(tremor))
  res<-data.frame(Threshold=x, res, stringsAsFactors = F)
  return(res)
}, data=msmtqcmerge)
res<-do.call('rbind', res)

res<-pivot_wider(res, id_cols = subject_id, names_from = Threshold, values_from = c(on_off, dyskinesia, tremor))
res<-res[!apply(!res[,-1], 1, all),]


res_watch<-lapply(thrsh, function(x, data){
  res<-data %>% filter(Highest...of.Consecutive.Above.Threshold.Windows >= x) %>% group_by(subject_id) %>% summarize(on_off=is_enough_data(on_off), dyskinesia=is_enough_data(dyskinesia), tremor=is_enough_data(tremor))
  res<-data.frame(Threshold=x, res, stringsAsFactors = F)
  return(res)
}, data=msmtqcwatch)
res_watch<-do.call('rbind', res_watch)

res_watch<-pivot_wider(res_watch, id_cols = subject_id, names_from = Threshold, values_from = c(on_off, dyskinesia, tremor))
res_watch<-res_watch[!apply(!res_watch[,-1], 1, all),]


res_numtotal<-lapply(thrsh, function(x, data){
  res<-data %>% group_by(subject_id) %>% 
    summarize(on_off=sum(((Highest...of.Consecutive.Above.Threshold.Windows.watch >= x&!is.na(Highest...of.Consecutive.Above.Threshold.Windows.watch))|(Highest...of.Consecutive.Above.Threshold.Windows.phone >= x&!is.na(Highest...of.Consecutive.Above.Threshold.Windows.phone)))&!is.na(on_off)), dyskinesia = sum(((Highest...of.Consecutive.Above.Threshold.Windows.watch >= x&!is.na(Highest...of.Consecutive.Above.Threshold.Windows.watch))|(Highest...of.Consecutive.Above.Threshold.Windows.phone >= x&!is.na(Highest...of.Consecutive.Above.Threshold.Windows.phone)))&!is.na(dyskinesia)), tremor=sum(((Highest...of.Consecutive.Above.Threshold.Windows.watch >= x&!is.na(Highest...of.Consecutive.Above.Threshold.Windows.watch))|(Highest...of.Consecutive.Above.Threshold.Windows.phone >= x&!is.na(Highest...of.Consecutive.Above.Threshold.Windows.phone)))&!is.na(tremor)))
  res<-data.frame(Threshold=x, res, stringsAsFactors = F)
  return(res)
}, data=msmtqcmerge)
res_numtotal<-do.call('rbind', res_numtotal)

res_numtotal<-pivot_wider(res_numtotal, id_cols = subject_id, names_from = Threshold, values_from = c(on_off, dyskinesia, tremor))
res_numtotal<-res_numtotal[res_numtotal$subject_id%in%res$subject_id,]


onofftotal<-colSums(res_numtotal[, grep("on_off", names(res_numtotal))]*res$on_off_0)
onofftotal[1]-onofftotal

dyskinesiatotal<-colSums(res_numtotal[, grep("dyskinesia", names(res_numtotal))]*res$dyskinesia_0)
dyskinesiatotal[1]-dyskinesiatotal

tremortotal<-colSums(res_numtotal[, grep("tremor", names(res_numtotal))]*res$tremor_0)
tremortotal[1]-tremortotal


plot_qq(QCcorphone$subth.on_off[QCcorphone$subject_id%in%res$subject_id[res$on_off_4]], "Num Sub-threshold/On-Off")
ggsave("Analyses/Plots/REAL-PD_phone_onoff_vs_numSubThresh_qq_filtered.png")
plot_qq(QCcorphone$subth.dyskinesia[QCcorphone$subject_id%in%res$subject_id[res$dyskinesia_4]], "Num Sub-threshold/Dyskinesia")
ggsave("Analyses/Plots/REAL-PD_phone_dyskinesia_vs_numSubThresh_qq_filtered.png")
plot_qq(QCcorphone$subth.tremor[QCcorphone$subject_id%in%res$subject_id[res$tremor_4]], "Num Sub-threshold/Tremor")
ggsave("Analyses/Plots/REAL-PD_phone_tremor_vs_numSubThresh_qq_filtered.png")
plot_qq(QCcorphone$consec.on_off[QCcorphone$subject_id%in%res$subject_id[res$on_off_4]], "Best Consecutive/On-Off")
ggsave("Analyses/Plots/REAL-PD_phone_onoff_vs_bestConsecutiveWindow_qq_filtered.png")
plot_qq(QCcorphone$consec.dyskinesia[QCcorphone$subject_id%in%res$subject_id[res$dyskinesia_4]], "Best Consecutive/Dyskinesia")
ggsave("Analyses/Plots/REAL-PD_phone_dyskinesia_vs_bestConsecutiveWindow_qq_filtered.png")
plot_qq(QCcorphone$consec.tremor[QCcorphone$subject_id%in%res$subject_id[res$tremor_4]], "Best Consecutive/Tremor")
ggsave("Analyses/Plots/REAL-PD_phone_tremor_vs_bestConsecutiveWindow_qq_filtered.png")


plot_qq(QCcorwatch$subth.on_off[QCcorwatch$subject_id%in%res$subject_id[res$on_off_4]], "Num Sub-threshold/On-Off")
ggsave("Analyses/Plots/REAL-PD_watch_onoff_vs_numSubThresh_qq_filtered.png")
plot_qq(QCcorwatch$subth.dyskinesia[QCcorwatch$subject_id%in%res$subject_id[res$dyskinesia_4]], "Num Sub-threshold/Dyskinesia")
ggsave("Analyses/Plots/REAL-PD_watch_dyskinesia_vs_numSubThresh_qq_filtered.png")
plot_qq(QCcorwatch$subth.tremor[QCcorwatch$subject_id%in%res$subject_id[res$tremor_4]], "Num Sub-threshold/Tremor")
ggsave("Analyses/Plots/REAL-PD_watch_tremor_vs_numSubThresh_qq_filtered.png")
plot_qq(QCcorwatch$consec.on_off[QCcorwatch$subject_id%in%res$subject_id[res$on_off_4]], "Best Consecutive/On-Off")
ggsave("Analyses/Plots/REAL-PD_watch_onoff_vs_bestConsecutiveWindow_qq_filtered.png")
plot_qq(QCcorwatch$consec.dyskinesia[QCcorwatch$subject_id%in%res$subject_id[res$dyskinesia_4]], "Best Consecutive/Dyskinesia")
ggsave("Analyses/Plots/REAL-PD_watch_dyskinesia_vs_bestConsecutiveWindow_qq_filtered.png")
plot_qq(QCcorwatch$consec.tremor[QCcorwatch$subject_id%in%res$subject_id[res$tremor_4]], "Best Consecutive/Tremor")
ggsave("Analyses/Plots/REAL-PD_watch_tremor_vs_bestConsecutiveWindow_qq_filtered.png")



save.image("Analyses/REAL-PD_SampleCounts_with_Thresholding.Rdata")

