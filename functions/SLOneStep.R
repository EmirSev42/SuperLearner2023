### SLOneStep.R ###
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
## Author: Emir S
## Created: Feb 6, 2023 
## Version: 
## Last-Updated: 
##           By:
##     Update #:
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
## 
### Commentary: 
# performs one step of the super learner algorithm (that is for the indicated i'th sample)

# inputs:

# samples_list: list; list of samples, as outputted by "LoobSample.R". See corresponding file.
# B: integer; the sample to consider. i.e the B'th subsample from the samples_list input will be considered.
# times: vector; time points to predict for.
# event_library: string vector; the event libraries to consider for the survSuperLearner function. 
# See survSuperLearner:::survlistWrappers() for available  default wrappers.
# cens_library: string vector; the censoring libraries to consider for the survSuperLearner function. 
# See survSuperLearner:::survlistWrappers() for available  default wrappers.

# output: list of 3. contents:

# preds: named list of predictions. Each member of the list is a (i x times) matrix of predictions 
# with a corresponding name (i.e preds$survSL.coxph for a cox model) where the unit i,j of the matrix is the risk prediction for the
# i'th person in the subsample for the j'th time. The meta learner predictions are alwats named "sl".

# scores: list of length "times" (see inputs),  named accordingly (t1, t2, t3..) that contains the model name, time and corresponding Brier score
# for each library along with the meta learner. By default arranged from lowest to highest Brier scores.

# coefs: single named vector of length equal to length(event_library). Contains the alpha (weight) of each library for the meta learner prediction.
# the libraries with coefficients 0 are dropped from the algorithm. Does not vary across time.

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

one_step<-function(samples_list,B,times=1:5,event_library,cens_library){
  # note run time
  st<-Sys.time()
  
  if(is.numeric(B)==F){
    stop("Invalid B")
  }
  
  if(!(B %in% 1:length(samples_list$cv_train))){
    stop("B out of bounds")
  }
  
  cat("preparing data... \n \n")

  # select the matrices for the train and test batches for the B'th subsample
  train_matrix<-samples_list$cv_train[[B]]%>%select(-c(years,death))
  test_matrix<-samples_list$cv_test[[B]]%>%select(-c(years,death))
  train_data<-samples_list$cv_train[[B]]
  test_data<-samples_list$cv_test[[B]]
  # select the corresponding time and event
  train_years<-samples_list$cv_train[[B]]$years
  train_event<-samples_list$cv_train[[B]]$death
  
  cat("fitting super learner... \n \n")

  # call survSuperLearner with the appropriate inputs
  sl_fit<-survSuperLearner(time=train_years,event=train_event,X=train_matrix,
                           new.times=times,event.SL.library=event_library,verbose=T,
                           cens.SL.library=cens_library,newX = test_matrix)
  
  
  cat("acquiring SL coefficients... \n \n")
  
  coefs <- sl_fit$event.coef
  
  cat("acquiring predictions... \n \n")
  
  # number of libraries used
  nlibraries<-sl_fit$event.SL.library$library%>%nrow()
  # risk predictions of the super learner (meta)
  sl_risk<-(1-(sl_fit$event.SL.predict))
  # list that will contain all the risks
  all_risks<-NULL
  for(i in 1:nlibraries){
    # extract library predictions
    # the event.library.predict from the fit object contains the predictions
    all_risks[[i]]<-1-((sl_fit$event.library.predict[,,i]))
  }
  # fill the last spot with the super learner risks
  all_risks[[nlibraries+1]]<-sl_risk
  # name the risk predictions correspondingly
  names(all_risks)<-c(event_library,"sl")
  
  # IMPORTANT STEP: IF THERE IS AN NA PREDICTION IN ANY OF THE LIBRARIES, REMOVE THAT LIBRARY FROM THE LIST OF PREDICTIONS!
  # otherwise the score function can fail
  for(i in 1:nlibraries){
    
    if(NA %in% all_risks[[i]]){
      cat("library", i, "has prodocued NA predictions, and is removed \n \n")
      all_risks<-all_risks[-i]
    }
  }
  
  preds<-all_risks
  
  cat("calculating Score... \n \n")
  
  x1<-Score(preds,
            data=test_data,
            formula=Hist(years,death)~1,times=times,null.model=F,se.fit=F,metrics="Brier")
  
  cat("compiling Brier Scores... \n \n")
  
  scores<-NULL
  for(i in times){
    scores[[i]]<-x1$Brier$score%>%filter(times==i)
  }
  names(scores)<-paste("t",times,sep="")
  
  cat("preparing output... \n \n")
  
  out<-list("preds"=preds,"scores"=scores,"coefs"=coefs)
  
  # save time
  et<-Sys.time()
  
  cat("Algorithm Complete. Run time: \n \n")
  print(et-st)
  
  return(out)
}


