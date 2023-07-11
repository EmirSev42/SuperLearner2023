### SSL.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 10 2023 (15:53) 
## Version: 
## Last-Updated: Jul 11 2023
##           By: Emir S
##     Update #: 9
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

SSL <- function(formula,
                data,
                event.SL.library,
                cens.SL.library,
                times, # prediction time horizons
                ...){
  require("survSuperLearner")
  require("prodlim")
  require("survival")
  call <- match.call(expand.dots=TRUE)
  formula.names <- try(all.names(formula),silent=TRUE)
  if (!(formula.names[2]=="Hist")) stop("The left hand side of formula should look like this: Hist(time,event).")
  actual.terms <- terms(formula,data=data)
  formula <- eval(call$formula)
  response <- model.response(model.frame(formula,data))
  Time <- as.numeric(response[,"time"])
  Event <- as.numeric(response[,"status"])
  ## remove intercept
  X <- model.matrix(actual.terms,data=data)[,-c(1),drop=FALSE]
  ## calling Westling 
  ssl = survSuperLearner(time = Time,
                         event = Event,
                         event.SL.library = event.SL.library,
                         cens.SL.library = cens.SL.library,
                         new.times = times,X=data.frame(X),
                         ...)
  output = list(ssl_fit = ssl,
                call = match.call(),
                times = times)
  class(output) <- "SSL"
  output
 
}

predictRisk.SSL <- function(object,newdata,times,...){
  # process newdata
  formula <- as.formula(object$call$formula)
  actual.terms <- terms(formula,data=newdata)
  formula <- eval(as.formula(object$call$formula))
  ## remove intercept
  newX <- model.matrix(actual.terms,data=newdata)[,-c(1),drop=FALSE]
  surv<-survSuperLearner:::predict.survSuperLearner(object$ssl_fit,
                                                     newdata = data.frame(newX),
                                              new.times = object$times,...)$event.SL.predict
  # risk
  risk<-1-surv
  return(risk)
}

# example:
library(riskRegression)

# generate sample data
set.seed(42)
d <- sampleData(203,outcome = "survival")

# SSL fit
set.seed(42)
ssl_fit<-SSL(formula = Hist(time,event)~X1+X8,data = d,
             event.SL.library=c("survSL.km","survSL.weibreg","survSL.coxph",
                                "survSL.rfsrc"),
             cens.SL.library=c("survSL.km","survSL.weibreg","survSL.coxph",
                               "survSL.rfsrc"),times=c(1,2))

# event coefficients
ssl_fit$ssl_fit$event.coef

# Predict
preds_ssl<-predictRisk.SSL(object=ssl_fit,newdata=d,times=c(1,2))

# Score
x<-Score(list(ssl_fit),
          data=data.frame(d),
          formula=Hist(time,event)~1,times=1:2,null.model=F,se.fit=F,metrics="Brier",
          plots = "Calibration")

# Calibration
plotCalibration(x,times=c(1,2),cens.method="local",xlim=c(0,0.5),ylim=c(0,0.5))

# try with loob split
x2<-Score(list(ssl_fit),
          data=data.frame(d),
          formula=Hist(time,event)~1,times=1:2,null.model=F,se.fit=F,metrics="Brier",
          plots = "Calibration",B=50,split.method="loob",
          M=round(.632*NROW(d)),seed=42,conf.int=F)

# calibration with loob
plotCalibration(x2,times=c(1,2),cens.method="local",xlim=c(0,0.5),ylim=c(0,0.5))

######################################################################
### SSL.R ends here
