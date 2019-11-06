
setwd("~/Documents/PDDB Challenge 2/CIS-PD")

library(synapser)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load Nick's sensor QC summary data
synfile<-synGet("syn21071756", downloadLocation="Analyses/Features/")
msmtqc<-read.csv(synfile$path, header=T, as.is=T)

#Load labels
scores_query <- synTableQuery("select * from syn20489608")
scores <- scores_query$asDataFrame()
scores$on_off[scores$on_off==-1]<-NA

# Load features
features<-read.csv("Analyses/Features/CIS-PD_Watch_Features_Full_NShawen.csv", header=T, as.is=T)
scores<-scores[scores$measurement_id%in%features$ID,]

scores<-scores[!(is.na(scores$on_off)&is.na(scores$dyskinesia)&is.na(scores$tremor)),]
# Check for scores not present in Nick' sensor QC file
#write.csv(scores[!scores$measurement_id %in% msmtqc$measurement_id,], "Analyses/Features/MeasurementIds_missing_from_VarianceThresholdResults_CIS-PD.csv", row.names=F, quote=F)

rm(features)

# Look for correlation between missingnesss and traits
QCcor<-msmtqc %>% group_by(subject_id) %>% summarize(subth.on_off=cor.test(Sub.Threshold.Windows, on_off, use="pairwise.complete", method = "kendall")$p.value,
                                              subth.dyskinesia=cor.test(Sub.Threshold.Windows, dyskinesia, use="pairwise.complete", method = "kendall")$p.value,
                                              subth.tremor=cor.test(Sub.Threshold.Windows, tremor, use="pairwise.complete", method = "kendall")$p.value,
                                              consec.on_off=cor.test(Highest...of.Consecutive.Above.Threshold.Windows, on_off, use="pairwise.complete", method = "kendall")$p.value,
                                              consec.dyskinesia=cor.test(Highest...of.Consecutive.Above.Threshold.Windows, dyskinesia, use="pairwise.complete", method = "kendall")$p.value,
                                              consec.tremor=cor.test(Highest...of.Consecutive.Above.Threshold.Windows, tremor, use="pairwise.complete", method = "kendall")$p.value)


plot_qq<-function(trait, traitname=""){
  trait<-trait[!is.na(trait)]
  n <- length(trait)
  p <- (1 : n) / n - 0.5 / n
  ggplot() + geom_point(aes(x = -log10(p), y = -log10(sort(trait)))) +
    xlab("-log10(Expected)") + ylab("-log10(Observed)") + geom_abline(intercept = 0, slope = 1, color = "red") +
    geom_hline(yintercept = -log10(0.05/n), color = "black", linetype = "dashed") + ggtitle(traitname)
}
plot_qq(QCcor$subth.on_off, "Num Sub-threshold/On-Off")
ggsave("Analyses/Plots/CIS-PD_onoff_vs_numSubThresh_qq.png")
plot_qq(QCcor$subth.dyskinesia, "Num Sub-threshold/Dyskinesia")
ggsave("Analyses/Plots/CIS-PD_dyskinesia_vs_numSubThresh_qq.png")
plot_qq(QCcor$subth.tremor, "Num Sub-threshold/Tremor")
ggsave("Analyses/Plots/CIS-PD_tremor_vs_numSubThresh_qq.png")
plot_qq(QCcor$consec.on_off, "Best Consecutive/On-Off")
ggsave("Analyses/Plots/CIS-PD_onoff_vs_bestConsecutiveWindow_qq.png")
plot_qq(QCcor$consec.dyskinesia, "Best Consecutive/Dyskinesia")
ggsave("Analyses/Plots/CIS-PD_dyskinesia_vs_bestConsecutiveWindow_qq.png")
plot_qq(QCcor$consec.tremor, "Best Consecutive/Tremor")
ggsave("Analyses/Plots/CIS-PD_tremor_vs_bestConsecutiveWindow_qq.png")

qplot(msmtqc$Sub.Threshold.Windows, geom="histogram", xlab = "# Subthreshold Windows")
ggsave("Analyses/Plots/CIS-PD_NumSubThresh_hist.png")
qplot(msmtqc$Highest...of.Consecutive.Above.Threshold.Windows, geom="histogram", xlab = "Longest Above Threshold Window")
ggsave("Analyses/Plots/CIS-PD_BestCohsecutiveWindow_hist.png")



# Try various thresholds for data filtering

thrsh<-c(0, 4, 6, 8, 10)


is_enough_data<-function(x){
  if(sum(!is.na(x))>40 & (sum(table(x)>=10)>1 | sum(table(x)>=5)>2)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

res<-lapply(thrsh, function(x, data){
  res<-data %>% filter(Highest...of.Consecutive.Above.Threshold.Windows >= x) %>% group_by(subject_id) %>% summarize(on_off=is_enough_data(on_off), dyskinesia=is_enough_data(dyskinesia), tremor=is_enough_data(tremor))
  res<-data.frame(Threshold=x, res, stringsAsFactors = F)
  return(res)
}, data=msmtqc)
res<-do.call('rbind', res)

res<-pivot_wider(res, id_cols = subject_id, names_from = Threshold, values_from = c(on_off, dyskinesia, tremor))
res<-res[!apply(!res[,-1], 1, all),]



res_numtotal<-lapply(thrsh, function(x, data){
  res<-data %>% group_by(subject_id) %>% summarize(on_off=sum(Highest...of.Consecutive.Above.Threshold.Windows >= x&!is.na(on_off)), dyskinesia_tremor=sum(Highest...of.Consecutive.Above.Threshold.Windows >= x))
  res<-data.frame(Threshold=x, res, stringsAsFactors = F)
  return(res)
}, data=msmtqc)
res_numtotal<-do.call('rbind', res_numtotal)

res_numtotal<-pivot_wider(res_numtotal, id_cols = subject_id, names_from = Threshold, values_from = c(on_off, dyskinesia_tremor))
res_numtotal<-res_numtotal[res_numtotal$subject_id%in%res$subject_id,]


onofftotal<-colSums(res_numtotal[, grep("on_off", names(res_numtotal))]*res$on_off_0)
onofftotal[1]-onofftotal

dyskinesiatotal<-colSums(res_numtotal[, grep("dyskinesia", names(res_numtotal))]*res$dyskinesia_0)
dyskinesiatotal[1]-dyskinesiatotal

tremortotal<-colSums(res_numtotal[, grep("tremor", names(res_numtotal))]*res$tremor_0)
tremortotal[1]-tremortotal

