###########################
# Build CART and ENET models ignoring any classes of size <6
# then score for Macro and Traditionally weighted MSE and MAE,
# as well as MSE and MAE within label class.
#
# This is compared to the Null estimate where the Null predictor
# is estimated from the Training split.
###########################

setwd("~/Documents/PDDB Challenge 2/CIS-PD")

library(caret)
library(synapser)
library(gdata)
library(pROC)
library(tidyr)
library(dplyr)
library(ggplot2)
#library(plyr)

#Load labels
scores_query <- synTableQuery("select * from syn20489608")
scores <- scores_query$asDataFrame()
scores$on_off[scores$on_off==-1]<-NA

# Load features
features<-read.csv("Analyses/Features/CIS-PD_Watch_Features_Full_NShawen.csv", header=T, as.is=T)
individual_summaries<-read.csv("Analyses/CIS-PD_phenotype_counts_by_individual.csv", header=T, as.is=T)
scores_filtered<-scores[scores$measurement_id%in%features$ID,]

indidvec<-unique(scores_filtered$subject_id)

#Generate statistics on reporting time by individual
on_off_time_summary<- scores %>% filter(!is.na(on_off)) %>% 
  mutate(time_diff = on_off_reported_timestamp-timestamp) %>% 
  group_by(subject_id) %>% summarize(timediff_mean= mean(time_diff), timediff_median= median(time_diff), timediff_sd= sd(time_diff), prop_gt_1hr=sum(time_diff>60*60)/sum(!is.na(time_diff)), prop_gt_2hr=sum(time_diff>2*60*60)/sum(!is.na(time_diff)))

dyskinesia_time_summary<- scores %>% filter(!is.na(dyskinesia)) %>% 
  mutate(time_diff = dyskinesia_reported_timestamp-timestamp) %>% 
  group_by(subject_id) %>% summarize(timediff_mean= mean(time_diff), timediff_median= median(time_diff), timediff_sd= sd(time_diff), prop_gt_1hr=sum(time_diff>60*60)/sum(!is.na(time_diff)), prop_gt_2hr=sum(time_diff>2*60*60)/sum(!is.na(time_diff)))

tremor_time_summary<- scores %>% filter(!is.na(tremor)) %>% 
  mutate(time_diff = tremor_reported_timestamp-timestamp) %>% 
  group_by(subject_id) %>% summarize(timediff_mean= mean(time_diff), timediff_median= median(time_diff), timediff_sd= sd(time_diff), prop_gt_1hr=sum(time_diff>60*60)/sum(!is.na(time_diff)), prop_gt_2hr=sum(time_diff>2*60*60)/sum(!is.na(time_diff)))


# Fit rpart return prob vector
fit_individual_model_return_prob<-function(indid, trait, scores, features){
  print(indid)
  scoresi<-scores[scores$subject_id==indid,]
  scoresi<-scoresi[!is.na(scoresi[,trait]),]
  if(dim(scoresi)[1]> 40){
    tbscores<-table(scoresi[,trait])
    if(any(tbscores)<6){
      scoresi[scoresi[,trait]%in%as.numeric(names(which(tbscores<6))),trait]<-NA
    }
  }
  scoresi<-scoresi[!is.na(scoresi[,trait]),]
  tbscores<-table(scoresi[,trait])
  print(tbscores)
  
  if(sum(tbscores>=6)>1 & length(scoresi[,trait]) > 40){
    scoresi[,trait]<-factor(scoresi[,trait], ordered = TRUE)
    
    intrain <- createDataPartition(y = scoresi[,trait], p = 0.8, list = FALSE)
    
    features_train<-features[features$ID%in%scoresi$measurement_id[intrain],]
    scores_train<-scoresi[match(features_train$ID, scoresi$measurement_id),]
    
    features_test<-features[features$ID%in%scoresi$measurement_id[-intrain],]
    scores_test<-scoresi[match(features_test$ID, scoresi$measurement_id),]
    
    myfolds<-groupKFold(scores_train$measurement_id, k = 5)
    
    preProcessInTrain <-c("center", "scale")
    
    
    rtmodel <- train(features_train[,-dim(features_train)[2]], scores_train[,trait], 
                     method = "rpart", tuneLength = 10, preProc = preProcessInTrain,
                     trControl = trainControl(verboseIter = TRUE, savePredictions = F, index=myfolds),
                     metric="Kappa")
    
    predicts2 <- predict(rtmodel, features_test[,-dim(features_test)[2]], type = "prob")
    
    predicts3<-aggregate(predicts2, list(scores_test$measurement_id), mean)    
    scores_test3<-scores_test[match(predicts3[,1], scores_test$measurement_id),trait]
    predicts3$mean.train<-mean(as.numeric(as.character(scoresi[intrain, trait])))
    predicts3$median.train<-median(as.numeric(as.character(scoresi[intrain, trait])))
    predicts3$mean.class.train<-mean(as.numeric(as.character(unique(scoresi[intrain, trait]))))
    predicts3$median.class.train<-median(as.numeric(as.character(unique(scoresi[intrain, trait]))))
    predicts3$Truth<-scores_test3
    
    return(predicts3)
    
  } else {
    return(NULL)
  }
}

nreps<-100
# Run model fitting for 100 training/test bootstraps of rpart model
on_off_results_predictions<-lapply(indidvec, function(x, trait, scores, features, nreps){
  replicate(nreps, fit_individual_model_return_prob(x, trait, scores, features), simplify =F)
}, "on_off", scores_filtered, features, nreps)

dyskinesia_results_predictions<-lapply(indidvec, function(x, trait, scores, features, nreps){
  replicate(nreps, fit_individual_model_return_prob(x, trait, scores, features), simplify =F)
}, "dyskinesia", scores_filtered, features, nreps)

tremor_results_predictions<-lapply(indidvec, function(x, trait, scores, features, nreps){
  replicate(nreps, fit_individual_model_return_prob(x, trait, scores, features), simplify =F)
}, "tremor", scores_filtered, features, nreps)



# Code to estrapolate numeric estimate from probability vectors
# 2 different statistics: max probability class, and weighted average
estimatecategory<-function(pred){
  if(!is.null(pred)){
    tmp<-pred[,-1]
    tmp<-tmp[,!names(tmp)%in%c("mean.train", "median.train", "mean.class.train", "median.class.train", "Truth")]
    nms<-as.numeric(names(tmp))
    pred$Predict.cat<-nms[apply(tmp,1, which.max)]
    pred$Predict.weightedave<-rowSums(sweep(tmp, MARGIN=2, nms, `*`))
  }
  return(pred)
}

#Functions to compute mse and mae
mse<-function(true, pred){
  if(!is.null(true)&!is.null(pred)){
    wh<-which(!is.na(true)&!is.na(pred))
    true<-as.numeric(as.character(true[wh]))
    pred<-as.numeric(as.character(pred[wh]))
    mse<-sum((true-pred)^2)/length(true)
    return(mse)
  }
}

mae<-function(true, pred){
  if(!is.null(true)&!is.null(pred)){
    wh<-which(!is.na(true)&!is.na(pred))
    true<-as.numeric(as.character(true[wh]))
    pred<-as.numeric(as.character(pred[wh]))
    mse<-sum(abs(true-pred))/length(true)
    return(mse)
  }
}

# Functions to compute MSE/MAE by label, macroaveraged, and traditional
compute_macroaveragemse<-function(true, pred){
  tmp<-data.frame(true, pred)
  tmp2<-tmp %>% group_by(true) %>% summarize(N=sum(!is.na(true)&!is.na(pred)), Value=mse(true, pred))
  tmp2<-as.data.frame(tmp2)
  tmp2<-rbind(tmp2, data.frame(true="Macro", N=sum(tmp2$N), Value=mean(tmp2$Value), stringsAsFactors = F))
  tmp2<-rbind(tmp2, data.frame(true="Traditional", N=tmp2$N[dim(tmp2)[1]], Value=mse(true, pred), stringsAsFactors = F))
  names(tmp2)[1]<-"Grouping.or.Version"
  tmp2<-data.frame(Statistic="MSE", tmp2, stringsAsFactors = F)
  return(tmp2)
}

compute_macroaveragemae<-function(true, pred){
  tmp<-data.frame(true, pred)
  tmp2<-tmp %>% group_by(true) %>% summarize(N=sum(!is.na(true)&!is.na(pred)), Value=mae(true, pred))
  tmp2<-as.data.frame(tmp2)
  tmp2<-rbind(tmp2, data.frame(true="Macro", N=sum(tmp2$N), Value=mean(tmp2$Value), stringsAsFactors = F))
  tmp2<-rbind(tmp2, data.frame(true="Traditional", N=tmp2$N[dim(tmp2)[1]], Value=mae(true, pred), stringsAsFactors = F))
  names(tmp2)[1]<-"Grouping.or.Version"
  tmp2<-data.frame(Statistic="MAE", tmp2, stringsAsFactors = F)
  return(tmp2)
}

compute_mse_mae_results<-function(res){
  if(!is.null(res)){
    mse.weighted=data.frame(Estimate="Weighted Average", compute_macroaveragemse(res$Truth, res$Predict.weightedave), stringsAsFactors = F)
    mse.categorical=data.frame(Estimate="Max Prob", compute_macroaveragemse(res$Truth, res$Predict.cat), stringsAsFactors = F)
    mse.null=compute_macroaveragemse(res$Truth,res$mean.train)
    mse.null<-data.frame(Estimate="NULL", mse.null[mse.null$Grouping.or.Version=="Traditional",], stringsAsFactors = F)
    tmp<-data.frame(Estimate="NULL", compute_macroaveragemse(res$Truth,res$mean.class.train), stringsAsFactors = F)
    tmp<-tmp[tmp$Grouping.or.Version!="Traditional",]
    mse.null<-rbind(tmp, mse.null)
    
    mae.weighted=data.frame(Estimate="Weighted Average", compute_macroaveragemae(res$Truth, res$Predict.weightedave), stringsAsFactors = F)
    mae.categorical=data.frame(Estimate="Max Prob", compute_macroaveragemae(res$Truth, res$Predict.cat), stringsAsFactors = F)
    mae.null=compute_macroaveragemae(res$Truth,res$mean.train)
    mae.null<-data.frame(Estimate="NULL", mae.null[mae.null$Grouping.or.Version=="Traditional",], stringsAsFactors = F)
    tmp<-data.frame(Estimate="NULL", compute_macroaveragemae(res$Truth,res$median.class.train), stringsAsFactors = F)
    tmp<-tmp[tmp$Grouping.or.Version!="Traditional",]
    mae.null<-rbind(tmp, mae.null)
    return(rbind(mse.categorical, mse.weighted,mse.null,
                 mae.categorical,mae.weighted,
                 mae.null))
  } 
  else {
    return(data.frame(Estimate=character(), Statistic=character(), Grouping.or.Version=character(),
                      N=integer(), Value=double(), stringsAsFactors = F))
  }
}

#Function to add subject_id column to dataframe elements in a list
annotate_id<-function(x, dat, columnname) { 
  dat<-dat[[x]]
  if(dim(dat)[1]>0){
    dat[,columnname]<-as.numeric(x)
  } else {
    dat[,columnname]<-integer()
  }
  return(dat)
}

# Estimate categories for all 3 phenotypes
on_off_results_predictions<-lapply(on_off_results_predictions, function(x) { lapply(x, estimatecategory) })
dyskinesia_results_predictions<-lapply(dyskinesia_results_predictions, function(x) { lapply(x, estimatecategory) })
tremor_results_predictions<-lapply(tremor_results_predictions, function(x) { lapply(x, estimatecategory) })

# USE THIS IF YOU NEED TO REMOVE ESTIMATED CATEGORIES FROM LISTS
# on_off_results_predictions<-lapply(on_off_results_predictions, function(x) { lapply(x, function(y){
#   if(!is.null(y)){
#     y<-y[,-dim(y)[2]]
#     y<-y[,-dim(y)[2]]
#   }
#   return(y)
# }) })

# Compute MSE/MAE statistics return large dataframe with all results
on_off_results_mse_mae<-lapply(on_off_results_predictions, function(x){ return(lapply(x, compute_mse_mae_results)) })
on_off_results_mse_mae<-lapply(on_off_results_mse_mae, function(X){ return( do.call('rbind', lapply(1:nreps, annotate_id, dat=X, columnname="Rep")) ) })
names(on_off_results_mse_mae)<-indidvec
on_off_results_mse_mae<-do.call('rbind', lapply(as.character(indidvec), annotate_id, dat=on_off_results_mse_mae, columnname="subject_id"))

dyskinesia_results_mse_mae<-lapply(dyskinesia_results_predictions, function(x){ return(lapply(x, compute_mse_mae_results)) })
dyskinesia_results_mse_mae<-lapply(dyskinesia_results_mse_mae, function(X){ return( do.call('rbind', lapply(1:nreps, annotate_id, dat=X, columnname="Rep")) ) })
names(dyskinesia_results_mse_mae)<-indidvec
dyskinesia_results_mse_mae<-do.call('rbind', lapply(as.character(indidvec), annotate_id, dat=dyskinesia_results_mse_mae, columnname="subject_id"))

tremor_results_mse_mae<-lapply(tremor_results_predictions, function(x){ return(lapply(x, compute_mse_mae_results)) })
tremor_results_mse_mae<-lapply(tremor_results_mse_mae, function(X){ return( do.call('rbind', lapply(1:nreps, annotate_id, dat=X, columnname="Rep")) ) })
names(tremor_results_mse_mae)<-indidvec
tremor_results_mse_mae<-do.call('rbind', lapply(as.character(indidvec), annotate_id, dat=tremor_results_mse_mae, columnname="subject_id"))


#Compute difference from NULL model
on_off_results_mse_mae<-spread(on_off_results_mse_mae, Estimate, Value)
on_off_results_mse_mae$MaxProb.Null.diff<-on_off_results_mse_mae$`Max Prob` - on_off_results_mse_mae$`NULL`
on_off_results_mse_mae$WeightedAve.Null.diff<-on_off_results_mse_mae$`Weighted Average` - on_off_results_mse_mae$`NULL`

dyskinesia_results_mse_mae<-spread(dyskinesia_results_mse_mae, Estimate, Value)
dyskinesia_results_mse_mae$MaxProb.Null.diff<-dyskinesia_results_mse_mae$`Max Prob` - dyskinesia_results_mse_mae$`NULL`
dyskinesia_results_mse_mae$WeightedAve.Null.diff<-dyskinesia_results_mse_mae$`Weighted Average` - dyskinesia_results_mse_mae$`NULL`

tremor_results_mse_mae<-spread(tremor_results_mse_mae, Estimate, Value)
tremor_results_mse_mae$MaxProb.Null.diff<-tremor_results_mse_mae$`Max Prob` - tremor_results_mse_mae$`NULL`
tremor_results_mse_mae$WeightedAve.Null.diff<-tremor_results_mse_mae$`Weighted Average` - tremor_results_mse_mae$`NULL`

###############
# Plot results CART Models
###############

#function to display N stats
stat_box_data <- function(x, lower_limit = NA){
  return( 
    data.frame(
      y = lower_limit*0.95,
      label = paste("n =", round(mean(x), 1), '\n')
    )
  )
}

# On/Off
on_off_results_mse_mae$Grouping.or.Version<-factor(on_off_results_mse_mae$Grouping.or.Version)
on_off_results_mse_mae$subject_id<-factor(on_off_results_mse_mae$subject_id)
tmp<-on_off_results_mse_mae[on_off_results_mse_mae$Statistic=="MSE",]
tmp$plotgroup<-0
tmp$plotgroup[!tmp$Grouping.or.Version%in%c("Macro", "Traditional")]<-1
tmp1<-tmp[tmp$plotgroup==1,]
ggplot(tmp1, aes(x=subject_id, y=WeightedAve.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MSE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90, position = position_dodge(width = 0.8),  fun.args = list(lower_limit=min(tmp1$WeightedAve.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/On_Off_CART_mse_components.png", width = 11, height = 6.5)
tmp0<-tmp[tmp$plotgroup==0,]
ggplot(tmp0, aes(x=subject_id, y=WeightedAve.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MSE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90,  position = position_dodge(width = 0.9),  fun.args = list(lower_limit=min(tmp0$WeightedAve.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/On_Off_CART_mse_summary.png", width = 9.5, height = 5.75)
#MAE
tmp<-on_off_results_mse_mae[on_off_results_mse_mae$Statistic=="MAE",]
tmp$plotgroup<-0
tmp$plotgroup[!tmp$Grouping.or.Version%in%c("Macro", "Traditional")]<-1
tmp1<-tmp[tmp$plotgroup==1,]
ggplot(tmp1, aes(x=subject_id, y=WeightedAve.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MAE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90, position = position_dodge(width = 0.8),  fun.args = list(lower_limit=min(tmp1$WeightedAve.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/On_Off_CART_mae_components.png", width = 11, height = 6.5)
tmp0<-tmp[tmp$plotgroup==0,]
ggplot(tmp0, aes(x=subject_id, y=WeightedAve.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MAE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90,  position = position_dodge(width = 0.9),  fun.args = list(lower_limit=min(tmp0$WeightedAve.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/On_Off_CART_mae_summary.png", width = 9.5, height = 5.75)

tmp<-on_off_results_mse_mae[on_off_results_mse_mae$Statistic=="MAE",]
tmp$plotgroup<-0
tmp$plotgroup[!tmp$Grouping.or.Version%in%c("Macro", "Traditional")]<-1
tmp0<-tmp[tmp$plotgroup==0,]
ggplot(tmp0, aes(x=subject_id, y=MaxProb.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MAE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90,  position = position_dodge(width = 0.9),  fun.args = list(lower_limit=min(tmp0$MaxProb.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")



# Dyskinesia
dyskinesia_results_mse_mae$Grouping.or.Version<-factor(dyskinesia_results_mse_mae$Grouping.or.Version)
dyskinesia_results_mse_mae$subject_id<-factor(dyskinesia_results_mse_mae$subject_id)
tmp<-dyskinesia_results_mse_mae[dyskinesia_results_mse_mae$Statistic=="MSE",]
tmp$plotgroup<-0
tmp$plotgroup[!tmp$Grouping.or.Version%in%c("Macro", "Traditional")]<-1
tmp1<-tmp[tmp$plotgroup==1,]
ggplot(tmp1, aes(x=subject_id, y=WeightedAve.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MSE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90, position = position_dodge(width = 0.8),  fun.args = list(lower_limit=min(tmp1$WeightedAve.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Dyskinesia_CART_mse_components.png", width = 11, height = 6.5)
tmp0<-tmp[tmp$plotgroup==0,]
ggplot(tmp0, aes(x=subject_id, y=WeightedAve.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MSE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90,  position = position_dodge(width = 0.9),  fun.args = list(lower_limit=min(tmp0$WeightedAve.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Dyskinesia_CART_mse_summary.png", width = 9.5, height = 5.75)
#MAE
tmp<-dyskinesia_results_mse_mae[dyskinesia_results_mse_mae$Statistic=="MAE",]
tmp$plotgroup<-0
tmp$plotgroup[!tmp$Grouping.or.Version%in%c("Macro", "Traditional")]<-1
tmp1<-tmp[tmp$plotgroup==1,]
ggplot(tmp1, aes(x=subject_id, y=WeightedAve.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MAE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90, position = position_dodge(width = 0.8),  fun.args = list(lower_limit=min(tmp1$WeightedAve.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Dyskinesia_CART_mae_components.png", width = 11, height = 6.5)
tmp0<-tmp[tmp$plotgroup==0,]
ggplot(tmp0, aes(x=subject_id, y=WeightedAve.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MAE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90,  position = position_dodge(width = 0.9),  fun.args = list(lower_limit=min(tmp0$WeightedAve.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Dyskinesia_CART_mae_summary.png", width = 9.5, height = 5.75)


# Tremor
tremor_results_mse_mae$Grouping.or.Version<-factor(tremor_results_mse_mae$Grouping.or.Version)
tremor_results_mse_mae$subject_id<-factor(tremor_results_mse_mae$subject_id)
tmp<-tremor_results_mse_mae[tremor_results_mse_mae$Statistic=="MSE",]
tmp$plotgroup<-0
tmp$plotgroup[!tmp$Grouping.or.Version%in%c("Macro", "Traditional")]<-1
tmp1<-tmp[tmp$plotgroup==1,]
ggplot(tmp1, aes(x=subject_id, y=WeightedAve.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MSE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90, position = position_dodge(width = 0.8),  fun.args = list(lower_limit=min(tmp1$WeightedAve.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Tremor_CART_mse_components.png", width = 11, height = 6.5)
tmp0<-tmp[tmp$plotgroup==0,]
ggplot(tmp0, aes(x=subject_id, y=WeightedAve.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MSE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90,  position = position_dodge(width = 0.9),  fun.args = list(lower_limit=min(tmp0$WeightedAve.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Tremor_CART_mse_summary.png", width = 9.5, height = 5.75)
#MAE
tmp<-tremor_results_mse_mae[tremor_results_mse_mae$Statistic=="MAE",]
tmp$plotgroup<-0
tmp$plotgroup[!tmp$Grouping.or.Version%in%c("Macro", "Traditional")]<-1
tmp1<-tmp[tmp$plotgroup==1,]
ggplot(tmp1, aes(x=subject_id, y=WeightedAve.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MAE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90, position = position_dodge(width = 0.8),  fun.args = list(lower_limit=min(tmp1$WeightedAve.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Tremor_CART_mae_components.png", width = 11, height = 6.5)
tmp0<-tmp[tmp$plotgroup==0,]
ggplot(tmp0, aes(x=subject_id, y=WeightedAve.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MAE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90,  position = position_dodge(width = 0.9),  fun.args = list(lower_limit=min(tmp0$WeightedAve.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Tremor_CART_mae_summary.png", width = 9.5, height = 5.75)





###############
# Fit ENET Models
###############

fit_enet_model_return_prob<-function(indid, trait, scores, features){
  print(indid)
  scoresi<-scores[scores$subject_id==indid,]
  scoresi<-scoresi[!is.na(scoresi[,trait]),]
  
  
  if(dim(scoresi)[1]> 40){
    tbscores<-table(scoresi[,trait])
    if(any(tbscores)<6){
      scoresi[scoresi[,trait]%in%as.numeric(names(which(tbscores<6))),trait]<-NA
    }
  }
  scoresi<-scoresi[!is.na(scoresi[,trait]),]
  tbscores<-table(scoresi[,trait])
  print(tbscores)
  
  if(sum(tbscores>=6)>1 & length(scoresi[,trait]) > 40){
    scoresi[,trait]<-factor(scoresi[,trait])
    intrain <- createDataPartition(y = scoresi[,trait], p = 0.8, list = FALSE)
    scoresi[,trait]<-as.numeric(scoresi[,trait])
    
    features_train<-features[features$ID%in%scoresi$measurement_id[intrain],]
    scores_train<-scoresi[match(features_train$ID, scoresi$measurement_id),]
    
    features_test<-features[features$ID%in%scoresi$measurement_id[-intrain],]
    scores_test<-scoresi[match(features_test$ID, scoresi$measurement_id),]
    
    myfolds<-groupKFold(scores_train$measurement_id, k = 5)
    
    preProcessInTrain <-c("center", "scale")
    
    enetmodel <- train(features_train[,-dim(features_train)[2]], scores_train[,trait], 
                       method = "glmnet", family = "gaussian", preProc = preProcessInTrain,
                       trControl = trainControl(verboseIter = TRUE, savePredictions = T, index=myfolds),
                       metric="RMSE")
    
    
    predicts2 <- predict(enetmodel, features_test[,-dim(features_test)[2]])
    #    rmse2<-RMSE(predicts2,scores_test[,trait])
    
    predicts3<-aggregate(predicts2, list(scores_test$measurement_id), mean)
    predicts3med<-aggregate(predicts2, list(scores_test$measurement_id), median)    
    scores_test3<-scores_test[match(predicts3[,1], scores_test$measurement_id),trait]
    #    auroc3<-as.numeric(multiclass.roc(scores_test3, predicts3[,-1])$auc)
    predicts3<-data.frame(Predict.mean=predicts3, Predict.median=predicts3med, stringsAsFactors = F)
    
    predicts3$mean.train<-mean(as.numeric(as.character(scoresi[intrain, trait])))
    predicts3$median.train<-median(as.numeric(as.character(scoresi[intrain, trait])))
    predicts3$mean.class.train<-mean(as.numeric(as.character(unique(scoresi[intrain, trait]))))
    predicts3$median.class.train<-median(as.numeric(as.character(unique(scoresi[intrain, trait]))))
    predicts3$Truth<-scores_test3
    
    
    return(predicts3)
  } 
}


set.seed(5000)
on_off_enet_results_predictions<-lapply(indidvec, function(x, trait, scores, features, nreps){
  replicate(nreps, fit_enet_model_return_prob(x, trait, scores, features), simplify =F)
}, "on_off", scores_filtered, features, nreps)

dyskinesia_enet_results_predictions<-lapply(indidvec, function(x, trait, scores, features, nreps){
  replicate(nreps, fit_enet_model_return_prob(x, trait, scores, features), simplify =F)
}, "dyskinesia", scores_filtered, features, nreps)


tremor_enet_results_predictions<-lapply(indidvec, function(x, trait, scores, features, nreps){
  replicate(nreps, fit_enet_model_return_prob(x, trait, scores, features), simplify =F)
}, "tremor", scores_filtered, features, nreps)


# Function to compute MSE/MAE scores for this enet results
compute_mse_mae_enet_results<-function(res){
  if(!is.null(res)){
    mse.weighted=data.frame(Estimate="Mean Summarized", compute_macroaveragemse(res$Truth, res$Predict.mean.x), stringsAsFactors = F)
    mse.categorical=data.frame(Estimate="Median Summarized", compute_macroaveragemse(res$Truth, res$Predict.median.x), stringsAsFactors = F)
    mse.null=compute_macroaveragemse(res$Truth,res$mean.train)
    mse.null<-data.frame(Estimate="NULL", mse.null[mse.null$Grouping.or.Version=="Traditional",], stringsAsFactors = F)
    tmp<-data.frame(Estimate="NULL", compute_macroaveragemse(res$Truth,res$mean.class.train), stringsAsFactors = F)
    tmp<-tmp[tmp$Grouping.or.Version!="Traditional",]
    mse.null<-rbind(tmp, mse.null)
    
    mae.weighted=data.frame(Estimate="Mean Summarized", compute_macroaveragemae(res$Truth, res$Predict.mean.x), stringsAsFactors = F)
    mae.categorical=data.frame(Estimate="Median Summarized", compute_macroaveragemae(res$Truth, res$Predict.median.x), stringsAsFactors = F)
    mae.null=compute_macroaveragemae(res$Truth,res$mean.train)
    mae.null<-data.frame(Estimate="NULL", mae.null[mae.null$Grouping.or.Version=="Traditional",], stringsAsFactors = F)
    tmp<-data.frame(Estimate="NULL", compute_macroaveragemae(res$Truth,res$median.class.train), stringsAsFactors = F)
    tmp<-tmp[tmp$Grouping.or.Version!="Traditional",]
    mae.null<-rbind(tmp, mae.null)
    return(rbind(mse.categorical, mse.weighted,mse.null,
                 mae.categorical,mae.weighted,
                 mae.null))
  } 
  else {
    return(data.frame(Estimate=character(), Statistic=character(), Grouping.or.Version=character(),
                      N=integer(), Value=double(), stringsAsFactors = F))
  }
}

# Compute MSE/MAE statistics for Enet return large dataframe with all results
on_off_enet_results_mse_mae<-lapply(on_off_enet_results_predictions, function(x){ return(lapply(x, compute_mse_mae_enet_results)) })
on_off_enet_results_mse_mae<-lapply(on_off_enet_results_mse_mae, function(X){ return( do.call('rbind', lapply(1:nreps, annotate_id, dat=X, columnname="Rep")) ) })
names(on_off_enet_results_mse_mae)<-indidvec
on_off_enet_results_mse_mae<-do.call('rbind', lapply(as.character(indidvec), annotate_id, dat=on_off_enet_results_mse_mae, columnname="subject_id"))

dyskinesia_enet_results_mse_mae<-lapply(dyskinesia_enet_results_predictions, function(x){ return(lapply(x, compute_mse_mae_enet_results)) })
dyskinesia_enet_results_mse_mae<-lapply(dyskinesia_enet_results_mse_mae, function(X){ return( do.call('rbind', lapply(1:nreps, annotate_id, dat=X, columnname="Rep")) ) })
names(dyskinesia_enet_results_mse_mae)<-indidvec
dyskinesia_enet_results_mse_mae<-do.call('rbind', lapply(as.character(indidvec), annotate_id, dat=dyskinesia_enet_results_mse_mae, columnname="subject_id"))

tremor_enet_results_mse_mae<-lapply(tremor_enet_results_predictions, function(x){ return(lapply(x, compute_mse_mae_enet_results)) })
tremor_enet_results_mse_mae<-lapply(tremor_enet_results_mse_mae, function(X){ return( do.call('rbind', lapply(1:nreps, annotate_id, dat=X, columnname="Rep")) ) })
names(tremor_enet_results_mse_mae)<-indidvec
tremor_enet_results_mse_mae<-do.call('rbind', lapply(as.character(indidvec), annotate_id, dat=tremor_enet_results_mse_mae, columnname="subject_id"))

#Compute difference from NULL model
on_off_enet_results_mse_mae<-spread(on_off_enet_results_mse_mae, Estimate, Value)
on_off_enet_results_mse_mae$MeanSummarized.Null.diff<-on_off_enet_results_mse_mae$`Mean Summarized` - on_off_enet_results_mse_mae$`NULL`
on_off_enet_results_mse_mae$MedianSummarized.Null.diff<-on_off_enet_results_mse_mae$`Median Summarized` - on_off_enet_results_mse_mae$`NULL`

dyskinesia_enet_results_mse_mae<-spread(dyskinesia_enet_results_mse_mae, Estimate, Value)
dyskinesia_enet_results_mse_mae$MeanSummarized.Null.diff<-dyskinesia_enet_results_mse_mae$`Mean Summarized` - dyskinesia_enet_results_mse_mae$`NULL`
dyskinesia_enet_results_mse_mae$MedianSummarized.Null.diff<-dyskinesia_enet_results_mse_mae$`Median Summarized` - dyskinesia_enet_results_mse_mae$`NULL`

tremor_enet_results_mse_mae<-spread(tremor_enet_results_mse_mae, Estimate, Value)
tremor_enet_results_mse_mae$MeanSummarized.Null.diff<-tremor_enet_results_mse_mae$`Mean Summarized` - tremor_enet_results_mse_mae$`NULL`
tremor_enet_results_mse_mae$MedianSummarized.Null.diff<-tremor_enet_results_mse_mae$`Median Summarized` - tremor_enet_results_mse_mae$`NULL`


###############
# Plot results ENET Models
###############


# On/Off
on_off_enet_results_mse_mae$Grouping.or.Version<-factor(on_off_enet_results_mse_mae$Grouping.or.Version)
on_off_enet_results_mse_mae$subject_id<-factor(on_off_enet_results_mse_mae$subject_id)
tmp<-on_off_enet_results_mse_mae[on_off_enet_results_mse_mae$Statistic=="MSE",]
tmp$plotgroup<-0
tmp$plotgroup[!tmp$Grouping.or.Version%in%c("Macro", "Traditional")]<-1
tmp1<-tmp[tmp$plotgroup==1,]
ggplot(tmp1, aes(x=subject_id, y=MeanSummarized.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MSE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90, position = position_dodge(width = 0.8),  fun.args = list(lower_limit=min(tmp1$MeanSummarized.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/On_Off_enet_mse_components.png", width = 11, height = 6.5)
tmp0<-tmp[tmp$plotgroup==0,]
ggplot(tmp0, aes(x=subject_id, y=MeanSummarized.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MSE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90,  position = position_dodge(width = 0.9),  fun.args = list(lower_limit=min(tmp0$MeanSummarized.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/On_Off_enet_mse_summary.png", width = 9.5, height = 5.75)
#MAE
tmp<-on_off_enet_results_mse_mae[on_off_enet_results_mse_mae$Statistic=="MAE",]
tmp$plotgroup<-0
tmp$plotgroup[!tmp$Grouping.or.Version%in%c("Macro", "Traditional")]<-1
tmp1<-tmp[tmp$plotgroup==1,]
ggplot(tmp1, aes(x=subject_id, y=MeanSummarized.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MAE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90, position = position_dodge(width = 0.8),  fun.args = list(lower_limit=min(tmp1$MeanSummarized.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/On_Off_enet_mae_components.png", width = 11, height = 6.5)
tmp0<-tmp[tmp$plotgroup==0,]
ggplot(tmp0, aes(x=subject_id, y=MeanSummarized.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MAE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90,  position = position_dodge(width = 0.9),  fun.args = list(lower_limit=min(tmp0$MeanSummarized.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/On_Off_enet_mae_summary.png", width = 9.5, height = 5.75)


# Dyskinesia
dyskinesia_enet_results_mse_mae$Grouping.or.Version<-factor(dyskinesia_enet_results_mse_mae$Grouping.or.Version)
dyskinesia_enet_results_mse_mae$subject_id<-factor(dyskinesia_enet_results_mse_mae$subject_id)
tmp<-dyskinesia_enet_results_mse_mae[dyskinesia_enet_results_mse_mae$Statistic=="MSE",]
tmp$plotgroup<-0
tmp$plotgroup[!tmp$Grouping.or.Version%in%c("Macro", "Traditional")]<-1
tmp1<-tmp[tmp$plotgroup==1,]
ggplot(tmp1, aes(x=subject_id, y=MeanSummarized.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MSE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90, position = position_dodge(width = 0.8),  fun.args = list(lower_limit=min(tmp1$MeanSummarized.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Dyskinesia_enet_mse_components.png", width = 11, height = 6.5)
tmp0<-tmp[tmp$plotgroup==0,]
ggplot(tmp0, aes(x=subject_id, y=MeanSummarized.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MSE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90,  position = position_dodge(width = 0.9),  fun.args = list(lower_limit=min(tmp0$MeanSummarized.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Dyskinesia_enet_mse_summary.png", width = 9.5, height = 5.75)
#MAE
tmp<-dyskinesia_enet_results_mse_mae[dyskinesia_enet_results_mse_mae$Statistic=="MAE",]
tmp$plotgroup<-0
tmp$plotgroup[!tmp$Grouping.or.Version%in%c("Macro", "Traditional")]<-1
tmp1<-tmp[tmp$plotgroup==1,]
ggplot(tmp1, aes(x=subject_id, y=MeanSummarized.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MAE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90, position = position_dodge(width = 0.8),  fun.args = list(lower_limit=min(tmp1$MeanSummarized.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Dyskinesia_enet_mae_components.png", width = 11, height = 6.5)
tmp0<-tmp[tmp$plotgroup==0,]
ggplot(tmp0, aes(x=subject_id, y=MeanSummarized.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MAE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90,  position = position_dodge(width = 0.9),  fun.args = list(lower_limit=min(tmp0$MeanSummarized.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Dyskinesia_enet_mae_summary.png", width = 9.5, height = 5.75)


# Tremor
tremor_enet_results_mse_mae$Grouping.or.Version<-factor(tremor_enet_results_mse_mae$Grouping.or.Version)
tremor_enet_results_mse_mae$subject_id<-factor(tremor_enet_results_mse_mae$subject_id)
tmp<-tremor_enet_results_mse_mae[tremor_enet_results_mse_mae$Statistic=="MSE",]
tmp$plotgroup<-0
tmp$plotgroup[!tmp$Grouping.or.Version%in%c("Macro", "Traditional")]<-1
tmp1<-tmp[tmp$plotgroup==1,]
ggplot(tmp1, aes(x=subject_id, y=MeanSummarized.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MSE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90, position = position_dodge(width = 0.8),  fun.args = list(lower_limit=min(tmp1$MeanSummarized.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Tremor_enet_mse_components.png", width = 11, height = 6.5)
tmp0<-tmp[tmp$plotgroup==0,]
ggplot(tmp0, aes(x=subject_id, y=MeanSummarized.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MSE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90,  position = position_dodge(width = 0.9),  fun.args = list(lower_limit=min(tmp0$MeanSummarized.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Tremor_enet_mse_summary.png", width = 9.5, height = 5.75)
#MAE
tmp<-tremor_enet_results_mse_mae[tremor_enet_results_mse_mae$Statistic=="MAE",]
tmp$plotgroup<-0
tmp$plotgroup[!tmp$Grouping.or.Version%in%c("Macro", "Traditional")]<-1
tmp1<-tmp[tmp$plotgroup==1,]
ggplot(tmp1, aes(x=subject_id, y=MeanSummarized.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MAE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90, position = position_dodge(width = 0.8),  fun.args = list(lower_limit=min(tmp1$MeanSummarized.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Tremor_enet_mae_components.png", width = 11, height = 6.5)
tmp0<-tmp[tmp$plotgroup==0,]
ggplot(tmp0, aes(x=subject_id, y=MeanSummarized.Null.diff, fill=Grouping.or.Version)) + geom_boxplot() + 
  xlab("Subject ID") + ylab("MAE Difference (Model-Null)")+
  stat_summary(aes(x=subject_id, y=N), fun.data = stat_box_data, geom = "text", size = 3, hjust = 0.5, vjust=0.9, angle=90,  position = position_dodge(width = 0.9),  fun.args = list(lower_limit=min(tmp0$MeanSummarized.Null.diff)*1.15)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
ggsave("Analyses/Plots/Tremor_enet_mae_summary.png", width = 9.5, height = 5.75)



# Write Out and Share Results

CART_results<-rbind(data.frame(Phenotype="On/Off", Model="CART", on_off_results_mse_mae, stringsAsFactors = F),
                    data.frame(Phenotype="Dyskinesia", Model="CART", dyskinesia_results_mse_mae, stringsAsFactors = F),
                    data.frame(Phenotype="Tremor", Model="CART", tremor_results_mse_mae, stringsAsFactors = F))

write.csv(CART_results, "Analyses/CISPD_CART_model_scores_total_and_byLabel.csv", row.names = F, quote=F)

ENET_results<-rbind(data.frame(Phenotype="On/Off", Model="ENET", on_off_enet_results_mse_mae, stringsAsFactors = F),
                    data.frame(Phenotype="Dyskinesia", Model="ENET", dyskinesia_enet_results_mse_mae, stringsAsFactors = F),
                    data.frame(Phenotype="Tremor", Model="ENET", tremor_enet_results_mse_mae, stringsAsFactors = F))

write.csv(ENET_results, "Analyses/CISPD_ENET_model_scores_total_and_byLabel.csv", row.names = F, quote=F)

rm(list=c("tmp", "tmp0", "tmp1"))

save.image("Analyses/Features/CIS-PD_Build_Models_Try2.Rdata")


#Store to Synapse
act<-Activity(name='CIS-PD Modeling/Scoring', description='Individual model building bootstrapping and scoring via (macro and traditional) MSE/MAE and within label class')
act$used(c('syn20712268', 'syn20489608'))
act$executed('https://raw.githubusercontent.com/sieberts/pddb2/master/Modeling_Solly/CIS-PD_Build_Models_v2.R')

syncart<-File("Analyses/CISPD_CART_model_scores_total_and_byLabel.csv", parentId='syn20833961')
syncart<-synStore(syncart, activity=act)

synenet<-File("Analyses/CISPD_ENET_model_scores_total_and_byLabel.csv", parentId='syn20833961')
synenet<-synStore(synenet, activity=act)
