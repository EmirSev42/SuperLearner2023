
# Snet.R
#----------------------------------------------------------------------
## Author: Emir S
## Created: July 12, 2023
## Version: 
## Last-Updated: 
##           By: Emir S
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: Survnet, predict risk and model fit based on 
# formula input and a predictRisk function for the appropriate class

## Survnet by itself doesn't have formula input, nor does it keep the call
# so we wrap it such that it has both, so that it may be used in the Score fct
### Change Log:
#----------------------------------------------------------------------

Snet <- function(formula,
                data,
                ...){
  require("prodlim")
  require("survnet")
  require("survival")
  call <- match.call(expand.dots=TRUE)
  formula.names <- try(all.names(formula),silent=TRUE)
  if (!(formula.names[2]=="Hist")) stop("The left hand side of formula should look like this: Hist(time,event).")
  actual.terms <- terms(formula,data=data)
  formula <- eval(call$formula)
  response <- model.response(model.frame(formula,data))
  Time <- as.numeric(response[,"time"])
  Event <- as.numeric(response[,"status"])
  
  # response
  Y <- data.frame("time"=Time,"event"=Event)
  ## remove intercept
  X <- model.matrix(actual.terms,data=data)[,-c(1),drop=FALSE]
  
  train_x<-apply(X,2,function(x){
    noNA<-x[!is.na(x)]
    
    if(length(unique(noNA))>2){
      x<-scale(x)
    }
    return(x)
  })
  
  if(length(Time)<250){
    time.grid<-seq(from=1,to=max(Time),length.out=250)
  }

  if(length(Time>=250)){
    time.grid<-sort(Time)
  }
  
  
  ## calling survnet
  fit_nn <- (survnet::survnet(x=train_x,y=Y,
                              breaks=sort(time.grid),epochs=30,verbose=0))
  
  output = list(nnet_fit = fit_nn,
                call = match.call(),
                times = Time)
  class(output) <- "Snet"
  return(output)
  
  
}

predictRisk.Snet <- function(object,newdata,times,...){
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
  
  time.grid<-object$nnet_fit$breaks
  
  predict_nn<-predict(object$nnet_fit,newdata=test_x)
  
  pred<-matrix(NA,nrow=nrow(test_x),ncol=length(times))
  
  for(i in 1:nrow(test_x)){
    pred[i,]<-stats::approx(x=time.grid,predict_nn[i,], 
                            xout = times, method = 'constant', rule = 2)$y
  }
  return(pred)
}

# example:
library(riskRegression)

# generate sample data
set.seed(42)
d <- sampleData(1000,outcome = "survival")
head(d)

# fit survnet
fit_nn<-Snet(formula = Hist(time,event)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data = d)

fit_nn$nnet_fit$breaks%>%length()

# predict risks
preds<-predictRisk.Snet(object=fit_nn,newdata=d,times=1:12)
preds%>%head()


# Score
x<-Score(list("Survnet"=fit_nn),
         data=data.frame(d),
         formula=Hist(time,event)~1,times=1:12,metrics=c("AUC","Brier"),
         plots = "Calibration")


# Calibration
plotCalibration(x,times=5,cens.method="local",
                xlim=c(0,1),ylim=c(0,1),brier.in.legend=T,auc.in.legend=T)


# try with loob split
x2<-Score(list("Survnet"=fit_nn),
          data=data.frame(d),
          formula=Hist(time,event)~1,times=1:3,null.model=F,se.fit=F,metrics="Brier",
          plots = "Calibration",B=50,split.method="loob",
          M=round(.632*NROW(d)),seed=42,conf.int=F)

# calibration with loob
plotCalibration(x2,times=3,cens.method="local",xlim=c(0,1),ylim=c(0,1))

