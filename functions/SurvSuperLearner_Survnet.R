 
#----------------------------------------------------------------------
## Author: Emir S
## Created: Feb 28, 2023
## Version: 
## Last-Updated: July 11, 2023
##           By: Emir S
##     Update #: 
#----------------------------------------------------------------------
## 
### Commentary: Survnet Library for SurvSuperLearner
## 
### Change Log:
#----------------------------------------------------------------------

# survnet funciton for SL library
survSL.survnet<-function(time, event, X, newX, new.times, obsWeights, ...){
  require(survnet)
  require(dplyr)
  
  y <- data.frame("time"=time,"event"=event)
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
  
  # time grid
  # survnet only performs well if the time input is relatively long, since it is recurrent
  # if the inputted time is of length less than 250, we create an appropriate time grid
  if(length(new.times==250)){
    time.grid<-new.times
  }
  
  if(length(new.times<250)){
    time.grid<-seq(from=1,to=max(time),length.out=250)
  }
  
  # the survnet call
  fit_nn <- survnet::survnet(x=train_x,y=y,breaks=sort(time.grid),epochs=20,verbose=0)

  # predictcall for survnet
  predict_nn<-(1-(predict(fit_nn,newdata=test_x)))
  
  # this empty matrix will contain the survival predictions
  pred<-matrix(NA,nrow=nrow(test_x),ncol=length(new.times))

  # fill up the prediction matrix, using interpolation at out prediction horizon points
  for(i in 1:nrow(test_x)){
    pred[i,]<-stats::approx(x=time.grid,predict_nn[i,], 
                            xout = new.times, method = 'constant', rule = 2)$y
  }

  # create output
  fit_nn$call<-match.call()
  fit <- list(object = fit_nn)
  class(fit) <- c("survSL.nn")
  out <- list(pred = pred, fit = fit)
  
  return(out)
}

# predict functionality
# this is the predict function for the class "SurvSL.nn" which we created above
predict.survSL.nn<-function(object, newX, new.times, ...){
  test_x<-model.matrix(~.,data=newX)[,-1]
  predict_nn<-(1-(predict(object$object,newdata=test_x)))
  time.grid<-object$object$breaks
  pred<-matrix(NA,nrow=nrow(test_x),ncol=length(new.times))
  
  for(i in 1:nrow(test_x)){
    pred[i,]<-stats::approx(x=time.grid,predict_nn[i,], 
                            xout = new.times, method = 'constant', rule = 2)$y
  }
  return(pred)
}
