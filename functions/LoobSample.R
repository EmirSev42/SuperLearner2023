
### LoobSample.R ###
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
# collect leave-one-out bootstrap samples
# input: data, B for the number of samples desired

# output: lists of cv_train, cv_test each of length B corresponding to the B'th sample
# matrix "oob" logical, showing if each subject was OOB in that sample or not

### Change Log:
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

loob<-function(data,B,rep=F,id=NULL){
  require(dplyr)
  
  # initialize outputs;
  cv_train<-NULL
  cv_test<-NULL
  OOB<-data.frame(matrix(NA,nrow=nrow(data),ncol=B))
  colnames(OOB)<-paste("sample",1:B,sep="")
  
  # create ID variable if id is left empty
  if(is.null(id)){
    obs_id<-1:nrow(data)
  }
  # rename ID if id is specified
  else{
    obs_id<-id
  }
  # assign ID
  data$obs_id<-obs_id
  
  # main loop
  for(i in 1:B){
    
    cat("now collecting sample", i, "...\n\n")
    
    # collect train batch for the i'th subsample
    train_id<-sample(1:(nrow(data)),size=(nrow(data)*0.632),replace=rep)
    # sort by ID
    cv_train[[i]]<-data[train_id,]%>%arrange(obs_id)
    # samples not in the train batch are allocated to the test batch
    cv_test[[i]]<-data[-unique(train_id),]%>%arrange(obs_id)
    # check if the person with that ID has been OOB
    OOB[,i]<-data$obs_id %in% cv_test[[i]]$obs_id
  }
  
  # check if some subjects are never OOB
  is_OOB<-NULL
  for(i in 1:nrow(OOB)){
    is_OOB[i]<-TRUE %in% OOB[i,]
  }
  # warn if some subjects are never OOB
  if(FALSE %in% is_OOB){
    cat("Some subjects are never OOB! Increase B")
  }
  # return content as list
  return(list("cv_train"=cv_train,"cv_test"=cv_test,
              "oob"=OOB))
}
