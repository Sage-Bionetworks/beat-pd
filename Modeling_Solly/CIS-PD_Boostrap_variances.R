
library(tidyr)
library(e1071)

setwd("~/Documents/PDDB Challenge 2/CIS-PD")
scores_withsplits_optimized<-read.csv("Analyses/CISPD_Labels_training_splits.csv", header=T, as.is=T)


boots<-10000
phenotypes<-c("on_off", "dyskinesia", "tremor")
trainingvars<-names(scores_withsplits_optimized)[grep("training",names(scores_withsplits_optimized))]

bootstrapvar<-function(test, training, boots){
  meanest<-mean(training)
  msereps<-replicate(boots, mean((sample(test, length(test), replace = TRUE)-meanest)^2))
  medest<-median(training)
  maereps<-replicate(boots, mean(abs(sample(test, length(test), replace = TRUE)-medest)))
  msereps<-var(msereps)*(boots-1)/boots
  maereps<-var(maereps)*(boots-1)/boots
  return(data.frame(MSE=msereps, MAE=maereps, stringsAsFactors = F))
}


wrapper_bootstrapvar<-function(trainingvar, pheno, scoresi, boots){
#  print(trainingvar)
  test<-scoresi[!scoresi[,trainingvar],pheno]
  test<-test[!is.na(test)]
  training <- scoresi[scoresi[,trainingvar],pheno]
  training <- training[!is.na(training)]
  res<-bootstrapvar(test=test, training = training, boots)
  res<-data.frame(Training_split=trainingvar, res, stringsAsFactors = F)
  return(res)
}

compute_bootstrap_var<- function(indid, scoresi, phenos, trainingvars, boots){
#  print(indid)
  scoresi<-scoresi[scoresi$subject_id==indid,]
  phenologic<-apply(is.na(scoresi[,phenos]), 2, all)
  phenos<-phenos[!phenologic]
  res<-NULL
  for(i in phenos){
#    print(i)
    temp<-do.call('rbind', lapply(trainingvars, wrapper_bootstrapvar, pheno=i, scoresi=scoresi, boots=boots))
    temp<-data.frame(Trait=i, temp, stringsAsFactors = F)
    res<-rbind(res, temp)
#    print(paste(i, "done"))
  }
  res<-data.frame(subject_id=indid, res, stringsAsFactors = F)
#  print(paste(indid, "done"))
  return(res)
}

set.seed(29)

res.var<-lapply(unique(scores_withsplits_optimized$subject_id), compute_bootstrap_var, scoresi=scores_withsplits_optimized, phenos=phenotypes, trainingvars=trainingvars, boots=boots)
res.var<-do.call('rbind', res.var)

mean_res_var<-res.var %>% group_by(subject_id, Trait) %>% summarize(MSE=mean(MSE), MAE=mean(MAE))
n_test<-scores_withsplits_optimized %>% group_by(subject_id) %>% summarize(ninverse=1/sum(!training1), n2inverse=1/sum(!training1)^2)

mean_res_var<-merge.data.frame(mean_res_var, n_test, by="subject_id")

data_var<- scores_withsplits_optimized %>% group_by(subject_id) %>% summarize(on_off=var(on_off, na.rm = T), dyskinesia=var(dyskinesia), tremor=var(tremor))
data_var <- gather(data_var, Trait, Var, on_off:tremor)

mean_res_var<-full_join(mean_res_var, data_var)

data_kurt<- scores_withsplits_optimized %>% group_by(subject_id) %>% summarize(on_off=kurtosis(on_off, na.rm = T), dyskinesia=kurtosis(dyskinesia), tremor=kurtosis(tremor))
data_kurt <- gather(data_kurt, Trait, Kurtosis, on_off:tremor)

mean_res_var<-full_join(mean_res_var, data_kurt)

mean_res_var$varovern<-mean_res_var$ninverse*mean_res_var$Var

data_skew<- scores_withsplits_optimized %>% group_by(subject_id) %>% summarize(on_off=skewness(on_off, na.rm = T), dyskinesia=skewness(dyskinesia), tremor=skewness(tremor))
data_skew <- gather(data_skew, Trait, Skew, on_off:tremor)

mean_res_var<-full_join(mean_res_var, data_skew)

par(mfrow=c(1,3))
plot(mean_res_var$Var[mean_res_var$Trait=="on_off"], mean_res_var$MSE[mean_res_var$Trait=="on_off"], xlab="Data Variance", ylab="Bootstrap MSE Variance", main ="On/Off")
plot(mean_res_var$Var[mean_res_var$Trait=="dyskinesia"], mean_res_var$MSE[mean_res_var$Trait=="dyskinesia"], xlab="Data Variance", ylab="Bootstrap MSE Variance", main ="Dyskinesia")
plot(mean_res_var$Var[mean_res_var$Trait=="tremor"], mean_res_var$MSE[mean_res_var$Trait=="tremor"], xlab="Data Variance", ylab="Bootstrap MSE Variance", main ="Tremor")

par(mfrow=c(1,3))
plot(mean_res_var$Skew[mean_res_var$Trait=="on_off"], mean_res_var$MSE[mean_res_var$Trait=="on_off"], xlab="Data Skewness", ylab="Bootstrap MSE Variance", main ="On/Off")
plot(mean_res_var$Skew[mean_res_var$Trait=="dyskinesia"], mean_res_var$MSE[mean_res_var$Trait=="dyskinesia"], xlab="Data Skewness", ylab="Bootstrap MSE Variance", main ="Dyskinesia")
plot(mean_res_var$Skew[mean_res_var$Trait=="tremor"], mean_res_var$MSE[mean_res_var$Trait=="tremor"], xlab="Data Skewness", ylab="Bootstrap MSE Variance", main ="Tremor")

par(mfrow=c(1,3))
plot(mean_res_var$Kurtosis[mean_res_var$Trait=="on_off"], mean_res_var$MSE[mean_res_var$Trait=="on_off"], xlab="Data Kurtosis", ylab="Bootstrap MSE Variance", main ="On/Off")
plot(mean_res_var$Kurtosis[mean_res_var$Trait=="dyskinesia"], mean_res_var$MSE[mean_res_var$Trait=="dyskinesia"], xlab="Data Kurtosis", ylab="Bootstrap MSE Variance", main ="Dyskinesia")
plot(mean_res_var$Kurtosis[mean_res_var$Trait=="tremor"], mean_res_var$MSE[mean_res_var$Trait=="tremor"], xlab="Data Kurtosis", ylab="Bootstrap MSE Variance", main ="Tremor")

par(mfrow=c(1,3))
plot(mean_res_var$ninverse[mean_res_var$Trait=="on_off"], mean_res_var$MSE[mean_res_var$Trait=="on_off"], xlab=expression(1/n), ylab="Bootstrap MSE Variance", main ="On/Off")
plot(mean_res_var$ninverse[mean_res_var$Trait=="dyskinesia"], mean_res_var$MSE[mean_res_var$Trait=="dyskinesia"], xlab=expression(1/n), ylab="Bootstrap MSE Variance", main ="Dyskinesia")
plot(mean_res_var$ninverse[mean_res_var$Trait=="tremor"], mean_res_var$MSE[mean_res_var$Trait=="tremor"], xlab=expression(1/n), ylab="Bootstrap MSE Variance", main ="Tremor")



par(mfrow=c(1,3))
plot(mean_res_var$varovern[mean_res_var$Trait=="on_off"], mean_res_var$MSE[mean_res_var$Trait=="on_off"], xlab=expression(hat(sigma)^2/n), ylab="Bootstrap MSE Variance", main ="On/Off")
abline(a=0, b=1, lty=2, col=2)
plot(mean_res_var$varovern[mean_res_var$Trait=="dyskinesia"], mean_res_var$MSE[mean_res_var$Trait=="dyskinesia"], xlab=expression(hat(sigma)^2/n), ylab="Bootstrap MSE Variance", main ="Dyskinesia")
abline(a=0, b=1, lty=2, col=2)
plot(mean_res_var$varovern[mean_res_var$Trait=="tremor"], mean_res_var$MSE[mean_res_var$Trait=="tremor"], xlab=expression(hat(sigma)^2/n), ylab="Bootstrap MSE Variance", main ="Tremor")
abline(a=0, b=1, lty=2, col=2)


par(mfrow=c(1,3))
plot(1/mean_res_var$varovern[mean_res_var$Trait=="on_off"], 1/mean_res_var$MSE[mean_res_var$Trait=="on_off"], xlab=expression(n/hat(sigma)^2), ylab="1/(Bootstrap MSE Variance)", main ="On/Off")
abline(a=0, b=1, lty=2, col=2)
plot(1/mean_res_var$varovern[mean_res_var$Trait=="dyskinesia"], 1/mean_res_var$MSE[mean_res_var$Trait=="dyskinesia"], xlab=expression(n/hat(sigma)^2), ylab="1/(Bootstrap MSE Variance)", main ="Dyskinesia")
abline(a=0, b=1, lty=2, col=2)
plot(1/mean_res_var$varovern[mean_res_var$Trait=="tremor"], 1/mean_res_var$MSE[mean_res_var$Trait=="tremor"], xlab=expression(n/hat(sigma)^2), ylab="1/(Bootstrap MSE Variance)", main ="Tremor")
abline(a=0, b=1, lty=2, col=2)


par(mfrow=c(1,3))
plot(mean_res_var$varovern[mean_res_var$Trait=="on_off"], mean_res_var$MAE[mean_res_var$Trait=="on_off"], xlab=expression(hat(sigma)^2/n), ylab="Bootstrap MAE Variance", main ="On/Off")
plot(mean_res_var$varovern[mean_res_var$Trait=="dyskinesia"], mean_res_var$MAE[mean_res_var$Trait=="dyskinesia"], xlab=expression(hat(sigma)^2/n), ylab="Bootstrap MAE Variance", main ="Dyskinesia")
plot(mean_res_var$varovern[mean_res_var$Trait=="tremor"], mean_res_var$MAE[mean_res_var$Trait=="tremor"], xlab=expression(hat(sigma)^2/n), ylab="Bootstrap MAE Variance", main ="Tremor")

par(mfrow=c(1,3))
plot(1/mean_res_var$varovern[mean_res_var$Trait=="on_off"], 1/mean_res_var$MAE[mean_res_var$Trait=="on_off"], xlab=expression(n/hat(sigma)^2), ylab="1/(Bootstrap MAE Variance)", main ="On/Off")
plot(1/mean_res_var$varovern[mean_res_var$Trait=="dyskinesia"], 1/mean_res_var$MAE[mean_res_var$Trait=="dyskinesia"], xlab=expression(n/hat(sigma)^2), ylab="1/(Bootstrap MAE Variance)", main ="Dyskinesia")
plot(1/mean_res_var$varovern[mean_res_var$Trait=="tremor"], 1/mean_res_var$MAE[mean_res_var$Trait=="tremor"], xlab=expression(n/hat(sigma)^2), ylab="1/(Bootstrap MAE Variance)", main ="Tremor")




par(mfrow=c(1,3))
hist(1/mean_res_var$ninverse[mean_res_var$Trait=="on_off"], xlab="n",  main ="On/Off")
hist(1/mean_res_var$ninverse[mean_res_var$Trait=="dyskinesia"], xlab="n", main ="Dyskinesia")
hist(1/mean_res_var$ninverse[mean_res_var$Trait=="tremor"], xlab="n",  main ="Tremor")


par(mfrow=c(1,3))
hist(1/sqrt(mean_res_var$ninverse[mean_res_var$Trait=="on_off"]), xlab=expression(sqrt(n)),  main ="On/Off")
hist(1/sqrt(mean_res_var$ninverse[mean_res_var$Trait=="dyskinesia"]), xlab=expression(sqrt(n)), main ="Dyskinesia")
hist(1/sqrt(mean_res_var$ninverse[mean_res_var$Trait=="tremor"]), xlab=expression(sqrt(n)),  main ="Tremor")

# library(sm)
# 
# trait.f <- factor(mean_res_var$Trait, levels= c("on_off","dyskinesia","tremor"),
#                 labels = c("On/Off", "Dyskinesia", "Tremor"))
# 
# tmp<-1/mean_res_var$ninverse
# sm.density.compare(tmp, mean_res_var$Trait, xlab="n")
# title(main="N Distribution by Trait")

par(mfrow=c(1,3))
plot(density(1/mean_res_var$ninverse[mean_res_var$Trait=="on_off"], na.rm=T), xlab="n",  main ="On/Off")
plot(density(1/mean_res_var$ninverse[mean_res_var$Trait=="dyskinesia"], na.rm=T), xlab="n",  main ="Dyskinesia")
plot(density(1/mean_res_var$ninverse[mean_res_var$Trait=="tremor"], na.rm=T), xlab="n",  main ="Tremor")

write.csv(mean_res_var %>% group_by(Trait) %>% summarize(min.mse=min(1/MSE, na.rm=T), max.mse=max(1/MSE, na.rm=T), min.var=min(1/varovern, na.rm=T), max.var=max(1/varovern, na.rm=T), min.n=min(1/ninverse, na.rm=T), max.n=max(1/ninverse, na.rm=T), min.sqrtn=min(1/sqrt(ninverse), na.rm=T), max.sqrtn=max(1/sqrt(ninverse), na.rm=T)), "CIS-PD_variance_weighting.csv", row.names=F, quote=F)


