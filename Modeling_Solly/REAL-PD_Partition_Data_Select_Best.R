
setwd("~/Documents/PDDB Challenge 2/REAL-PD")

library(synapser)
library(gdata)
library(plyr)
library(dplyr)

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

msmtqcmerge<-full_join(msmtqcwatch, msmtqcphone, by=c("measurement_id","subject_id","medication_state","slowness_walking","tremor","on_off","dyskinesia"),
                       suffix=c(".watch", ".phone"))


# Threshold based on enough watch or phone data
winthresh<-4
scores<-scores[scores$measurement_id%in%msmtqcmerge$measurement_id[(msmtqcmerge$Highest...of.Consecutive.Above.Threshold.Windows.watch >= winthresh&!is.na(msmtqcmerge$Highest...of.Consecutive.Above.Threshold.Windows.watch))|(msmtqcmerge$Highest...of.Consecutive.Above.Threshold.Windows.phone >= winthresh&!is.na(msmtqcmerge$Highest...of.Consecutive.Above.Threshold.Windows.phone))],]

scores<-scores[!(is.na(scores$on_off)&is.na(scores$dyskinesia)&is.na(scores$tremor)),]





indidvec<-unique(scores$subject_id)

is_enough_data<-function(x){
  if(sum(!is.na(x))>40 & (sum(table(x)>=10)>1 | sum(table(x)>=5)>2)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
subsample <- function(dat, p) {
  if (nrow(dat) == 1) {
    out <- dat$index
  }
  else {
    num <- ceiling(nrow(dat) * p)
    out <- sample(dat$index, size = num)
  }
  out
}

createsplits_v2<-function(indid, scores, p, times){
  print(indid)
  scoresi<-scores[scores$subject_id==indid,]
  traits<-c("on_off", "dyskinesia", "tremor")
  whichtraits<-sapply(traits, function(x){ return(is_enough_data(scoresi[,x])) }, simplify = T)
  scoresi[,traits[!whichtraits]]<-NA
  
  scoresi2<-scoresi[,c("measurement_id", "subject_id", traits)]
  traits<-traits[whichtraits]
  
  if(any(whichtraits)){
    for(i in 1:times){
      colname<-paste("training", i, sep="")
      scoresi2[,colname]<-F
      for(tr in traits){
        tbl<-table(scoresi2[,tr])
        if(any(tbl==1)){
          scoresi2[scoresi2[,tr]%in%as.numeric(names(tbl)[tbl==1]),colname]<-T
        }
      }
      fixed<-scoresi2$measurement_id[scoresi2[,colname]]
      
      num <- ceiling(nrow(scoresi2) * p) - length(fixed)
      out <- sample(scoresi2$measurement_id, size = num)
      
      scoresi2[scoresi2$measurement_id%in%out,colname]<-T
      
      for(tr in traits){
        tbl<-table(scoresi2[,tr], scoresi2[,colname])
        if(any(tbl[,2]==0)){      # Make sure there are no values in test but not training
          for(k in as.numeric(names(tbl[,2])[tbl[,2]==0])){
            chng<-scoresi2$measurement_id[scoresi2[,tr]==k&!is.na(scoresi2[,tr])]
            chng<-sample(chng)[1]
            chng2<-sample(scoresi2$measurement_id[scoresi2[,colname]==TRUE&!scoresi$measurement_id%in%fixed])[1]
            scoresi2[scoresi2$measurement_id==chng, colname]<-TRUE
            scoresi2[scoresi2$measurement_id==chng2, colname]<-FALSE
          }
        }
        
        tbl<-table(scoresi2[,tr], scoresi2[,colname])
        if(any(tbl[,1]==0&tbl[,2]>1)){      # If there are 2 or more observations, make sure at least one is in test
          for(k in as.numeric(names(tbl[,1])[tbl[,1]==0&tbl[,2]>1])){
            chng<-scoresi2$measurement_id[scoresi2[,tr]==k&!is.na(scoresi2[,tr])]
            chng<-sample(chng)[1]
            chng2<-sample(scoresi2$measurement_id[scoresi2[,colname]==FALSE&!scoresi$measurement_id%in%fixed])[1]
            scoresi2[scoresi2$measurement_id==chng, colname]<-FALSE
            scoresi2[scoresi2$measurement_id==chng2, colname]<-TRUE
          }
        }
      }
    }
    return(scoresi2)
  }
}
  
set.seed(1000)
times<-10000
p<-0.75
scores_withsplits_v2<-lapply(indidvec, createsplits_v2, scores=scores, p=p, times=times)
scores_withsplits_v2<-do.call('rbind', scores_withsplits_v2)


totalprop2<-ddply(scores_withsplits_v2, .(subject_id), summarize, prop.table(table(training1))[1])

tbls2<-lapply(unique(scores_withsplits_v2$subject_id), function(x, scores){
  lapply(c("on_off", "dyskinesia", "tremor"), function(y, x, scores){
    scoresi<-scores[scores$subject_id==x,]
    return(cbind(table(scoresi[,y], scoresi$training1), prop.table(table(scoresi[,y], scoresi$training1), margin = 1)))
  }, x=x, scores=scores)
}, scores=scores_withsplits_v2)



comp_dev<-function(X, Y, p){
  tt<-table(Y, X)
  tp<-prop.table(tt, margin = 1)
  tot<-rowSums(tt)
  va<-sum((tp[,2]-p)^2*tot)/sum(tot)
  dif<-abs(prop.table(table(X))["TRUE"]-p)
  return(Var=va)
}

comp_lowestvar_splits<-function(indid, scoresi, p){
  scoresi<-scoresi[scoresi$subject_id==indid,]
  res<-do.call('rbind', lapply(c("on_off", "dyskinesia", "tremor"), function(y, scoresi, p){
    Y<-scoresi[,y]
    res<-data.frame(Trait=y,apply(scoresi[,grep("training", names(scoresi))], MARGIN=2, FUN = comp_dev, Y=Y, p=p), stringsAsFactors = F)
    res<-data.frame(Rep=rownames(res), res, stringsAsFactors = F)
    return(res)
  }, scoresi=scoresi, p=p))
  names(res)[dim(res)[2]]<-"Var"
  res<-reshape(res, idvar = "Rep", timevar = "Trait", direction="wide")
  res<-res[,!apply(is.na(res), 2, all)]
  tmp<-as.data.frame(apply(res[-1], 2, rank))
  names(tmp)<-gsub("Var", "Rank", names(tmp))
  res<-cbind(res, tmp)
  if(length(grep("Rank", names(res)))>1){
   res$Mean.rank<-apply(res[,grep("Rank", names(res))], 1, mean)
  } else {
   res$Mean.rank<-res[,grep("Rank", names(res))]
  }
  res$Rank.mean.rank<-rank(res$Mean.rank)
  return(res)
}

split_ranks<-lapply(unique(scores_withsplits_v2$subject_id),comp_lowestvar_splits,
                    scoresi=scores_withsplits_v2, p=p)  



names(split_ranks)<-unique(scores_withsplits_v2$subject_id)
scores_withsplits_optimized<-lapply(names(split_ranks), function(indid, split_ranks, scoresi, nreps=10){
  reps <- split_ranks[[as.character(indid)]] %>% top_n(nreps, -Rank.mean.rank) %>% select(Rep, Rank.mean.rank)
  if(dim(reps)[1]>nreps){
    remv<-sample(reps$Rep[reps$Rank.mean.rank==max(reps$Rank.mean.rank)], dim(reps)[1]-nreps)
    reps<-reps[!reps$Rep%in%remv,]
  }
  cols<-c(names(scoresi)[-grep("training", names(scoresi))], unlist(reps$Rep))
  scoresi<-scoresi[scoresi$subject_id==indid,cols]
  names(scoresi)<-c(names(scoresi)[-grep("training", names(scoresi))], paste("training", 1:nreps, sep=""))
  return(scoresi)
}, split_ranks=split_ranks, scoresi=scores_withsplits_v2)
scores_withsplits_optimized<-do.call('rbind', scores_withsplits_optimized)


tbls3<-lapply(unique(scores_withsplits_optimized$subject_id), function(x, scores){
  lapply(c("on_off", "dyskinesia", "tremor"), function(y, x, scores){
    scoresi<-scores[scores$subject_id==x,]
    return(cbind(table(scoresi[,y], scoresi$training1), prop.table(table(scoresi[,y], scoresi$training1), margin = 1)))
  }, x=x, scores=scores)
}, scores=scores_withsplits_optimized)


write.csv(scores_withsplits_optimized, "Analyses/REALPD_Labels_training_splits_v2.csv", row.names = F, quote=F)

#Store to Synapse
act<-Activity(name='REAL-PD Partition Data', description='Split into 75/25 training/test splits across phenotypes, 10 reps, after filtering low variance segments')
act$used(c('syn20822276', 'syn21123211', 'syn21123212'))
act$executed('https://raw.githubusercontent.com/sieberts/pddb2/master/Modeling_Solly/REAL-PD_Partition_Data_Select_Best.R')

syncart<-File("Analyses/REALPD_Labels_training_splits_v2.csv", parentId='syn20836185')
syncart<-synStore(syncart, activity=act)





