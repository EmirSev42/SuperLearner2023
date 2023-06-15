 
#----------------------------------------------------------------------
## Author: Emir S
## Created: Feb 28, 2023
## Version: 
## Last-Updated: 
##           By: 
##     Update #: 
#----------------------------------------------------------------------
## 
### Commentary: Survnet Library for SurvSuperLearner

## notes: Survnet requires some data processing (rescale continuous vars, and code discrete vars)
# have a single response y that contains the time and event variables, 
# we could pre-process it and then run; the problem is the survSuperLearner function does not like that;
# it wants to see a single consistent time, event and X (covariate matrix) input for ALL libraries considered
# thus we have to any data changes for a specific library inside the library function itself, as we do here

# the caveat is we have to select the variables to be considered manually each time
#----------------------------------------------------------------------


survSL.survnet<-function(time, event, X, newX, new.times, obsWeights, ...){
  require(survnet)
  require(dplyr)
  
 
  y <- data.frame("time"=time,"event"=event)

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
                      train_cont_predictors_rescaled)
  test_x<-data.frame(test_disc_predictors_encoded,
                     test_cont_predictors_rescaled)
  
  if(length(new.times==250)){
    time.grid<-new.times
  }
  
  if(length(new.times<250)){
    time.grid<-seq(from=1,to=max(time),length.out=250)
  }
  
  #cat("grid length:",length(time.grid),"\n\n")
  
  #cat("fitting survnet... \n\n")
  
  fit.nn <- survnet::survnet(x=train_x,y=y,breaks=sort(time.grid),epochs=20,verbose=0)
  
  #cat("no. breaks:",fit.nn$breaks%>%length(),"\n\n")
  
  #cat("predicting survnet... \n\n")
  
  predict_nn<-(1-(predict(fit.nn,newdata=test_x)))
  
  #cat("dimension of grid predictions:",dim(predict_nn),"\n\n")
  
  #cat("compiling time predictions... \n\n")
  
  #cat("nrow(X):",nrow(test_x),"\n\n")
  
  pred<-matrix(NA,nrow=nrow(test_x),ncol=length(new.times))
  
  for(i in 1:nrow(test_x)){
    pred[i,]<-stats::approx(x=time.grid,predict_nn[i,], 
                            xout = new.times, method = 'constant', rule = 2)$y
  }
  
  #cat("dimension of time predictions:",dim(pred),"\n\n")
  
  #cat("concluding... \n\n")
  
  fit <- list(object = fit.nn)
  class(fit) <- c("survSL.nn")
  out <- list(pred = pred, fit = fit)
  
  return(out)
}

