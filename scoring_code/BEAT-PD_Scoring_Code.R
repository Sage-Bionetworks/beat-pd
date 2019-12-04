#
# BEAT-PD Scoring Function:
# Metric is MSE, scored within individual, 
# and weighted by either log or sqrt
#
# Input: filename = csv filepath,
#         trait = "on_off", "dyskinesia" or "tremor"
# Output: Vector contianing sqrt-weighted 
#         score and log-weighted score
#

library(synapser)
library(dplyr)

weightedMSE<-function(filename, trait){
  # Get predictions and truth
  pred<-read.csv(filename, header=T, as.is=T)
  dat<-getTruth(trait)

  dat$pred<-pred$prediction[match(dat$measurement_id, pred$measurement_id)]
  
  mseframe<- dat %>% group_by(subject_id) %>% summarize(n=length(truth), mse=getMSE(truth, pred))
  res<-c(sqrt_weighted_mse=weighted.mean(mseframe$mse, sqrt(mseframe$n)), log_weighted_mse=weighted.mean(mseframe$mse, log(mseframe$n)))
  return(res)
}

getTruth<-function(trait){
  realsynid<-"syn21292051"
  cissynid<-"syn21291582"
  
  realsyn<-synGet(realsynid)
  real<-read.csv(realsyn$path, header=T, as.is=T)
  cissyn<-synGet(cissynid)
  cis<-read.csv(cissyn$path, header=T, as.is=T)
  
  truth<-rbind(cis, real)
  truth<-truth[, c("measurement_id", "subject_id", trait)]
  names(truth)[3]<-"truth"
  truth<-truth[!is.na(truth$truth),]
  return(truth)
}

getMSE<-function(x, y){
  return(mean((x-y)^2))
}
