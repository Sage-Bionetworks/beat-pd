setwd("~/Documents/PDDB Challenge 2/CIS-PD")

library(synapser)
library(gdata)
library(plyr)


#Load labels
scores_query <- synTableQuery("select * from syn20489608")
scores <- scores_query$asDataFrame()
scores$on_off[scores$on_off==-1]<-NA

# Load features
features<-read.csv("Analyses/Features/CIS-PD_Watch_Features_Full_NShawen.csv", header=T, as.is=T)
scores<-scores[scores$measurement_id%in%features$ID,]

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

createsplits<-function(indid, scores, p, times){
  scoresi<-scores[scores$subject_id==indid,]
  traits<-c("on_off", "dyskinesia", "tremor")
  whichtraits<-sapply(traits, function(x){ return(is_enough_data(scoresi[,x])) }, simplify = T)
  scoresi[,traits[!whichtraits]]<-NA

  if(any(whichtraits)){
    for(i in 1:times){
      tmp <- unlist(dlply(data.frame(scoresi[,c("measurement_id", "subject_id", names(whichtraits)[whichtraits])], index = seq(along = scoresi[,names(whichtraits)[1]])), 
             names(whichtraits)[whichtraits], subsample, p = p))
      scoresi[,paste("training", i, sep="")]<-F
      scoresi[tmp,paste("training", i, sep="")]<-T
    }
    return(scoresi)
  }
}

set.seed(1000)

p<-0.71
times<-20

scores_withsplits<-lapply(indidvec, createsplits, scores=scores, p=p, times=times)
scores_withsplits<-do.call('rbind', scores_withsplits)

write.csv(scores_withsplits, "Analyses/CISPD_Labels_training_splits.csv", row.names = F, quote=F)

#Store to Synapse
act<-Activity(name='CIS-PD Partition Data', description='Split into 75/25 training/test splits across phenotypes, 20 reps')
act$used(c('syn20712268', 'syn20489608'))
act$executed('https://raw.githubusercontent.com/sieberts/pddb2/master/Modeling_Solly/CIS-PD_Partition_Data.R')

syncart<-File("Analyses/CISPD_Labels_training_splits.csv", parentId='syn20552049')
syncart<-synStore(syncart, activity=act)





tbls<-lapply(unique(scores_withsplits$subject_id), function(x, scores){
  lapply(c("on_off", "dyskinesia", "tremor"), function(y, x, scores){
    scoresi<-scores[scores$subject_id==x,]
    return(cbind(table(scoresi[,y], scoresi$training1), prop.table(table(scoresi[,y], scoresi$training1), margin = 1)))
  }, x=x, scores=scores)
}, scores=scores_withsplits)
