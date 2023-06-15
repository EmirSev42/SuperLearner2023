#----------------------------------------------------------------------
## Author: Emir S
## Created: Feb 28, 2023
## Version: 
## Last-Updated: 
##           By: 
##     Update #: 
#----------------------------------------------------------------------
## 
### Commentary: XGBOOST Library for SurvSuperLearner

## notes: XGBoost requires some data processing (rescale continuous vars, and code discrete vars)
# and have a time variable set as POSITIVE if it is NOT censored, and NEGATIVE if it IS censored
# we could pre-process it and then run; the problem is the survSuperLearner function does not like that;
# it wants to see a single consistent time, event and X (covariate matrix) input for ALL libraries considered
# thus we have to any data changes for a specific library inside the library function itself, as we do here

# the caveat is we have to select the variables to be considered manually each time
#----------------------------------------------------------------------


survSL.xgboost<-function(time, event, X, newX, new.times, obsWeights, ...){
  
  require(xgboost)
  # surv wrapper for xgboost
  require(survXgboost)
  # xgb dependency
  require(Ckmeans.1d.dp)
  require(dplyr)

  # time (event) positive: not censored, negative: censored
  xbg_time<-ifelse(event == 1, time, -time)

  # SELECT VARIABLES HERE
  train_cont_predictors<-X%>%select(c("age","gfr21","pACRc"))
  train_disc_predictors<-X%>%select("male","dm","cvd")
  
  test_cont_predictors<-newX%>%select(c("age","gfr21","pACRc"))
  test_disc_predictors<-newX%>%select("male","dm","cvd")
  
  # create the dummy variables object using caret
  train_dummy <- caret::dummyVars(" ~ .", data=train_disc_predictors)
  test_dummy <- caret::dummyVars(" ~ .", data=test_disc_predictors)
  
  train_disc_predictors_encoded <- data.frame(predict(train_dummy, newdata = train_disc_predictors))
  test_disc_predictors_encoded <- data.frame(predict(test_dummy, newdata = test_disc_predictors))
  
  train_cont_predictors_rescaled<-data.frame(lapply(train_cont_predictors,scale))
  test_cont_predictors_rescaled<-data.frame(lapply(test_cont_predictors,scale))
  
  
  # final data matrices for train and test
  train_x<-data.frame(train_disc_predictors_encoded,
                      train_cont_predictors_rescaled)%>%as.matrix()
  test_x<-data.frame(test_disc_predictors_encoded,
                     test_cont_predictors_rescaled)%>%as.matrix()
  
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
