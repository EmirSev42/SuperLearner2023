### SSL.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 10 2023 (15:53) 
## Version: 
## Last-Updated: Jul 10 2023 (16:22) 
##           By: Thomas Alexander Gerds
##     Update #: 7
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(riskRegression)
library(prodlim)
d <- sampleData(203,outcome = "survival")
SSL(formula = Hist(time,event)~X1+X8,data = d)

with(d,survSuperLearner(time = time,
                        event = event,
                        X = d[,.(X1,X8)],
                        event.SL.library = c("survSL.coxph", "survSL.rfsrc"),
                        cens.SL.library = c("survSL.coxph", "survSL.rfsrc")))


SSL <- function(formula,
                data,
                event.SL.library,
                cens.SL.library,
                times, # prediction time horizons
                ...){
    ## requireNamesSpace("survSuperLearner")
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
                           new.times = times,
                           ...)
    output = list(ssl_fit = ssl,
                  call = match.call(),
                  times = times)
    class(output) <- "SSL"
    output
}

predictRisk.SSL <- function(object,newdata,times,...){
    # process newdata
    formula = object$call$formula
    actual.terms <- terms(formula,data=newdata)
    formula <- eval(call$formula)
    ## remove intercept
    newX <- model.matrix(actual.terms,data=newdata)[,-c(1),drop=FALSE]
    predict(ssl_fit,newdata = newX,newtimes = object$times,...)
    # risk
    1-surv
}


######################################################################
### SSL.R ends here
