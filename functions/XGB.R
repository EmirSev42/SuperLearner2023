
# XGB.R
#----------------------------------------------------------------------
## Author: Emir S
## Created: July 12, 2023
## Version: 
## Last-Updated: 
##           By: Emir S
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: Xgboost, predict risk and model fit based on 
# formula input and a predictRisk function for the appropriate class

## Xgboost by itself doesn't have formula input, nor does it keep the call
# so we wrap it such that it has both, so that it may be used in the Score fct
### Change Log:
#----------------------------------------------------------------------

Xgb <- function(formula,
                 data,
                 ...){
  require(xgboost)
  # wrapper for xgboost
  require(survXgboost)
  # xgb dependency
  require(Ckmeans.1d.dp)
  call <- match.call(expand.dots=TRUE)
  formula.names <- try(all.names(formula),silent=TRUE)
  if (!(formula.names[2]=="Hist")) stop("The left hand side of formula should look like this: Hist(time,event).")
  actual.terms <- terms(formula,data=data)
  formula <- eval(call$formula)
  response <- model.response(model.frame(formula,data))
  Time <- as.numeric(response[,"time"])
  Event <- as.numeric(response[,"status"])
  
  xbg_time<-ifelse(Event== 1, Time, -Time)
  
  ## remove intercept
  X <- model.matrix(actual.terms,data=data)[,-c(1),drop=FALSE]
  
  train_x<-apply(X,2,function(x){
    noNA<-x[!is.na(x)]
    
    if(length(unique(noNA))>2){
      x<-scale(x)
    }
    return(x)
  })
  
  labels<-xbg_time
  
  # data matrices that the xgb functions would accept
  xgb_matrix <- xgb.DMatrix(data = train_x, label=labels)
  
  # watchlist
  watchlist<-list("train"=xgb_matrix)
  
  # call xgb
  surv_xgboost_model <- xgb.train.surv(
    params = list(
      objective = "survival:cox",
      eval_metric = "cox-nloglik",
      eta = 0.05 # larger eta leads to algorithm not converging, resulting in NaN predictions
    ),watchlist=watchlist,early_stopping_rounds = 2, data = train_x, label = labels,
    nrounds = 50,verbose=0)
  
  output = list(xgb_fit = surv_xgboost_model,
                call = match.call(),
                times = Time)
  class(output) <- "Xgb"
  return(output)
  
}

predictRisk.Xgb <- function(object,newdata,times,...){
  formula <- as.formula(object$call$formula)
  actual.terms <- terms(formula,data=newdata)
  formula <- eval(as.formula(object$call$formula))
  
  ## remove intercept
  X <- model.matrix(actual.terms,data=newdata)[,-c(1),drop=FALSE]
  
  test_x<-apply(X,2,function(x){
    noNA<-x[!is.na(x)]
    
    if(length(unique(noNA))>2){
      x<-scale(x)
    }
    return(x)
  })
  
  
  preds<-1-predict(object = object$xgb_fit, newdata = test_x, type = "surv",
          times = times)
  return(preds)
}

# example:
library(riskRegression)

# generate sample data
set.seed(42)
d <- sampleData(1000,outcome = "survival")
head(d)

# fit xgboost
fit_xgb<-Xgb(formula = Hist(time,event)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data = d)

# predict risks
preds<-predictRisk.Xgb(object=fit_xgb,newdata=d,times=1:12)

# Score
x<-Score(list("Xgboost"=fit_xgb),
         data=data.frame(d),
         formula=Hist(time,event)~1,times=1:12,metrics=c("AUC","Brier"),
         plots = "Calibration")

# Calibration
plotCalibration(x,times=12,cens.method="local",
                xlim=c(0,1),ylim=c(0,1),brier.in.legend=T,auc.in.legend=T)


# try with loob split
x2<-Score(list("Xgboost"=fit_xgb),
          data=data.frame(d),
          formula=Hist(time,event)~1,times=1:3,null.model=F,se.fit=F,metrics="Brier",
          plots = "Calibration",B=50,split.method="loob",
          M=round(.632*NROW(d)),seed=42,conf.int=F)

# calibration with loob
plotCalibration(x2,times=3,cens.method="local",xlim=c(0,1),ylim=c(0,1))

