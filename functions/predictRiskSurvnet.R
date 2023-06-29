### predictRiskSurvnet.R --- 
#----------------------------------------------------------------------
## Author: Emir S
## Created: Dec 8 2021 (16:29) 
## Version: 
## Last-Updated: Dec 5 2022
##           By: Emir S
##     Update #: 5
#----------------------------------------------------------------------
## 
### Commentary: Attempting update the placeholder approach for Survnet
## 
### Change Log:
#----------------------------------------------------------------------
## 


### placeholder:
pk <- function(data){
    placeholder <- list()
    class(placeholder) <- "pk"
    placeholder$call <- match.call()
    placeholder
}


# regular predict risk function for survnet
# if we try to run this with score (with validation), then we have a warning:

# Error in FUN(X[[i]], ...) : model Survnet does not have a call argument.

predictRisk.survnet<-function(object,newdata,times,cause){
  cat("prepating data...\n")
  disc_predictors<-newdata[,c("male","dm","cvd")]
  dummy <- caret::dummyVars(" ~ .", data=disc_predictors)
  # use the predict function with the dummy object for the data
  disc_predictors_encoded <- data.frame(predict(dummy, newdata = disc_predictors))
  cont_predictors<-newdata[,c("age","gfr21","pACRc")]
  cont_predictors_rescaled<-data.frame(lapply(cont_predictors,scale))
  # data matrix for survnet
  x<-data.frame(disc_predictors_encoded,
                cont_predictors_rescaled)
  # acquire predictions
  cat("acquiring predictions...\n")
  predictions_matrix <- survnet:::predict.survnet(object=object,newdata=x,
                                                  cause=cause)
  
  predictions<-t(sapply(1:nrow(predictions_matrix), function(i) {
    stats::approx(c(1:max_horizon), c(predictions_matrix[i,]),
                  method = "constant", xout = times, 
                  rule = 2)$y }))
  return(predictions)
}

# predictRisk function for placeholder
# CATCH: the fit.object must be defaulted to the survnet model of choice
# for some reason, Score doesn't work otherwise

# to call: set object as the place holder ph, using ph<-pk(NULL) and set default for
# fit.object as the survnet object of your choice, in our case "fit_survnet"

predictRisk.pk <- function(object,fit.object=fit_survnet,newdata,times,cause){
  cat("\n prepating data...\n")
  disc_predictors<-newdata[,c("male","dm","cvd")]
  dummy <- caret::dummyVars(" ~ .", data=disc_predictors)
  # use the predict function with the dummy object for the data
  disc_predictors_encoded <- data.frame(predict(dummy, newdata = disc_predictors))
  cont_predictors<-newdata[,c("age","gfr21","pACRc")]
  cont_predictors_rescaled<-data.frame(lapply(cont_predictors,scale))
  
  # data matrix for survnet
  x<-data.frame(disc_predictors_encoded,
                cont_predictors_rescaled)
  # acquire predictions
  cat("\n acquiring predictions...\n")
  predictions_matrix <- survnet:::predict.survnet(fit.object,x)
  predictions<-t(sapply(1:nrow(predictions_matrix), function(i) {
    stats::approx(c(1:max_horizon), c(predictions_matrix[i,]),
                  method = "constant", xout = times, 
                  rule = 2)$y }))
  out<-predictions
  return(out)
}
