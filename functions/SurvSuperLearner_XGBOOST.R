
#----------------------------------------------------------------------
## Author: Emir S
## Created: Feb 28, 2023
## Version: 
## Last-Updated: July 11, 2023
##           By: 
##     Update #: 
#----------------------------------------------------------------------
## 
### Commentary: XGBOOST Library for SurvSuperLearner
## 
### Change Log:
#----------------------------------------------------------------------


survSL.xgboost<-function(time, event, X, newX, new.times, obsWeights, ...){
  
  require(xgboost)
  # wrapper for xgboost
  require(survXgboost)
  # xgb dependency
  require(Ckmeans.1d.dp)
  
  xbg_time<-ifelse(event == 1, time, -time)
  
  # model matrices
  train_x<-model.matrix(~.,data=X)[,-1]
  test_x<-model.matrix(~.,data=newX)[,-1]
  
  # scale each column if it is not binary
  train_x<-apply(train_x,2,function(x){
    noNA<-x[!is.na(x)]
    
    if(length(unique(noNA))>2){
      x<-scale(x)
    }
    return(x)
  })
  test_x<-apply(test_x,2,function(x){
    noNA<-x[!is.na(x)]
    
    if(length(unique(noNA))>2){
      x<-scale(x)
    }
    return(x)
  })
  
  # "labels", in this case the status for melanoma
  labels<-xbg_time
  
  # data matrices that the xgb functions would accept
  xgb_matrix <- xgb.DMatrix(data = train_x, label=labels)
  
  watchlist<-list("train"=xgb_matrix)
  
  surv_xgboost_model <- xgb.train.surv(
    params = list(
      objective = "survival:cox",
      eval_metric = "cox-nloglik",
      eta = 0.05 # larger eta leads to algorithm not converging, resulting in NaN predictions
    ),watchlist=watchlist,early_stopping_rounds = 2, data = train_x, label = labels,
    nrounds = 20,verbose=0)
  
  pred <- predict(object = surv_xgboost_model, newdata = test_x, type = "surv",
                  times = new.times)
  
  fit <- list(object = surv_xgboost_model)
  class(fit) <- c("survSL.xgb")
  out <- list(pred = pred, fit = fit)
  
  return(out)
  
}

# predict functionality
predict.survSL.xgb<-function(object, newX, new.times, ...){
  test_x<-model.matrix(~.,data=newX)[,-1]
  pred<-predict(object = object$object, newdata = test_x, type = "surv",
                times = new.times)
  return(pred)
}
