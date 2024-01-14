#core code
####machine learning model####
library(gbm)
library(xgboost)
library(randomForest)
library(kernlab)
library(pROC)
library(tidyverse)
library(sampling)
library(reshape2)
library(caret)
###parameter tuning
#gbm
gbmGrid <-  expand.grid(
  interaction.depth = 3*(3:7),
  n.trees = 1000*(3:6),
  shrinkage = 0.005,
  n.minobsinnode = 3*(2:6)
)
set.seed(seed)
gbmFit <- train(
  x = x.train,
  y = y.train,
  method = 'gbm',
  trControl = fitControl, 
  verbose = FALSE,
  tuneGrid = gbmGrid,
  metric = "ROC"
)
bestTune <- gbmFit$bestTune
#xgb
xgbGrid <-  expand.grid(
  nrounds=2000*(1:3), 
  max_depth=5*(1:3), 
  eta=0.005, 
  gamma=0, 
  subsample=1, 
  colsample_bytree=0.1*(8:10), 
  min_child_weight=2^(0:4)
)
set.seed(seed)
xgbFit <- train(
  x = x.train,
  y = y.train,
  method = "xgbTree",
  trControl = fitControl, 
  verbose = FALSE,
  tuneGrid = xgbGrid,
  metric = "ROC"
)
bestTune <- xgbFit$bestTune
#svm
svmGrid <-  expand.grid(
  sigma = 2^(-10:-1),
  C = 1:10
)
set.seed(seed)
svmFit <- train(
  x = x.train,
  y = y.train,
  method = "svmRadial",
  trControl = fitControl, 
  verbose = FALSE,
  tuneGrid = svmGrid,
  metric = "ROC"
)
bestTune <- svmFit$bestTune
#rf
rfGrid <-  expand.grid(
  mtry = 2:20
)
set.seed(seed)
rfFit <- train(
  x = x.train,
  y = y.train,
  method = "rf",
  trControl = fitControl, 
  verbose = FALSE,
  tuneGrid = rfGrid,
  metric = "ROC"
)
bestTune <- rfFit$bestTune
###sampling
#train:test(4:1)
n <- count(train_data_all,train_data_all$dis==1) %>% as.data.frame()
set.seed(seed)  
train_test_id <- strata(train_data_all,stratanames = 'dis',size = c(0.2*n[2,2],0.2*n[1,2]),method = 'srswor')
train_test <- train_data_all[train_test_id$ID_unit,]
train_train <- train_data_all[-train_test_id$ID_unit,]
#different situations in spatial modelling
if (length(positivepoint$DDD)<376 & length(positivepoint$DDD)>124)  {
  set.seed(seed)
  sample0 <- sample(negativepoint$DDD,500-length(positivepoint$DDD),replace = F,prob = negativepoint$weightindex)%>%as.data.frame()
  names(sample0) <- 'DDD'
  sample1 <- left_join(sample0,allpoint,by='DDD') 
  sample1 <- rbind(sample1,positivepoint)
}
if (length(positivepoint$DDD)<125) {
  set.seed(seed)
  sample0 <- sample(negativepoint$DDD,length(positivepoint$DDD)*3,replace = F,prob = negativepoint$weightindex)%>%as.data.frame()
  names(sample0) <- 'DDD'
  sample1 <- left_join(sample0,allpoint,by='DDD') 
  sample1 <- rbind(sample1,positivepoint)
}
if (length(positivepoint$DDD)>375) {
  set.seed(seed)
  sample0=strata(allpoint,stratanames="dis",size=c(375,125),method="systematic",pik=allpoint$weightindex)
  sample1 <- getdata(allpoint,sample0)
}
###modelling
#train
set.seed(seed+i)
ML.train <- train(
  x = x.train,
  y = y.train,
  method = 'gbm',#"xgbTree","svmRadial","rf"
  trControl = fitControl,
  verbose = FALSE,
  tuneGrid = bestTune,
  metric = "ROC"
)
#prediction
predict(ML.train, newdata = pred_data,type='prob')[2]
#contribution
varImp(ML.train,scale = F)$importance[1]
#assessment
rocobj1 <- roc(dataset1$case,dataset1$pred,auc=T)
assessment<- coords(rocobj1, "best",
                    best.method ="youden",
                    ret=c("threshold", "specificity", "sensitivity", "accuracy",
                          "precision", "recall","youden","tn","tp","fn","fp"), transpose = FALSE)
assessment$auc <- rocobj1$auc
assessment$f1score <- 2*assessment$precision*assessment$recall/(assessment$precision+assessment$recall)
#cutoff
cutoff <- assessment$threshold
####entropy weight####
Rescale = function(x) {
  rng = range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1]) + 1/10000
}
Entropy_Weight = function(X) {
  X = lapply(X, Rescale)
  P = data.frame(lapply(X, function(x) x / sum(x)))
  e = sapply(P, function(x) sum(x * log(x)) *(-1/log(nrow(P))))
  d = 1 - e
  w = d / sum(d)
  return(w)
}
####durbin####
library(spdep)
library(spatialreg)
#data prepare
nbProvince <- read.gal('map/provincefinal.gal',override.id = T)
OBJECTID <- attr(nbProvince,'region.id')
nbProvince1 <- read.gwt2nb('map/provincefinal1.gwt',region.id = OBJECTID)
ndists <- attr(nbProvince1,'GeoDa')$dist
invdist.province <- lapply(ndists,function(x) 1/x)
lwProvince <- nb2listw(nbProvince1,glist = invdist.province)
#moran test
morancase <- moran.test(seq1$cases,lwProvince)
#GNS
gns <- sacsarlm(fm,data = seq1,lwProvince,Durbin = T)
gns_sum <- summary(gns)
#impact
imps <- impacts(gns,listw=lwProvince,R=2000)
imps_sum <- summary(imps,zstats=T,short=T)
####sequence analysis####
##distance
library(ape)
library(phangorn)
NS_fasta <- read.dna("NS.fasta",format = 'fasta')
NS_fasta_1 <- phyDat(NS_fasta)
NS_dist <- dist.ml(NS_fasta_1)
##similarity
library(Biostrings)
pair_wise <- pairwiseAlignment(sequence1, sequence2)
similarity <- pid(pair_wise)
