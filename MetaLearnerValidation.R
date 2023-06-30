#----------------------------------------------------------------------
## Author: Emir S
## Created: March 9, 2023
## Version: 
## Last-Updated: Jun 30 2023 (09:11) 
##           By: Thomas Alexander Gerds
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: Meta-learner approach for super learning, with leave-one-out bootstrap validation
###
###             EMIR: please state here in more detail what this script tries to achieve
###                   and perhaps where you would like to have my feedback/help/suggestions
## 
### Change Log:
#----------------------------------------------------------------------
rm(list=ls())
gc()
# libraries:

# meta-learner function
library(survSuperLearner)
# piping and others
library(tidyverse)
# Score function
library(riskRegression)
# survival
library(survival)
# surv neural network
library(survnet)
# for dummy variables
library(caret)
# xgboost functions
library(xgboost)
# wrapper for xgboost
library(survXgboost)
# xgb dependency
library(Ckmeans.1d.dp)

# optional set wd
try(setwd("V:/RESEARCH_EMPLOYEES/EMIR SEVINC/G3b4_2023"))
try(setwd("~/research/SuperVision/Emir/SuperLearner2023/"))

#----------------------------------------------------------------------
# STEP 1: LOAD DATA
#----------------------------------------------------------------------

df<-read.csv("./dataderived/SynthData.csv")
# select relevant columns
df<-df%>%select(c("years","death","age","gfr21","pACRc","male","dm","cvd"))

# random sample, for faster run, to check if the codes work
n<-10000
set.seed(42)
df<-df[sample(1:nrow(df),size=n,replace=F),]

# pick prediction horizon (max)
set.seed(42)
boot.max <- bootstrap::bootstrap(df$years,
                                 1000, # number of bootstrap samples
                                 theta=function(x){max(x)})

# this is the maximum prediction horizon in years
# extracted via bootstrap method
max_horizon<-boot.max$thetastar%>%min()%>%floor()
max_horizon

# ----------------------------------------------------------------------- #
# Step 2: COLLECT LOOB SAMPLE
# ----------------------------------------------------------------------- #

# self-written LOOB function
source("./functions/LoobSample.R")
# see corresponding file for ino

# set samples to collect
B<-50
set.seed(42)
# the output will indicate if some subjects are never OOB
sample<-loob(df,B=B)

# ----------------------------------------------------------------------- #
# Step 2.5: ONE STEP OF THE ALGORITHM
# ----------------------------------------------------------------------- #

# let's test one step of our algorithm to make sure it works...
# function that performs the super-learner algorithm for one subsample
source("./functions/SLOneStep.R")
# see corresponding file for info

# functions that allow using Survnet and XGBoost
source("./functions/SurvSuperLeaner_Survnet.R")
source("./functions/SurvSuperLeaner_XGBOOST.R")
# see corresponding files for info

# optional: list available wrappers
survSuperLearner:::survlistWrappers()

# before we run it, we must establish the library to be used
# let's keep it simple and use two basic libraries (cox,rf) plus one custom written one (xgb)
event.SL.library <- cens.SL.library <- c("survSL.coxph",
                                         "survSL.rfsrc","survSL.xgboost")

# start time
onestep_st<-Sys.time()
# call our one-step function for subsample 1
foo<-one_step(samples_list=sample,B=1,times=1:5,event_library=event.SL.library,
              cens_library=cens.SL.library)
# end time
onestep_et<-Sys.time()
# end time - start time
print(onestep_et-onestep_st)

# optional:  save for later use
# saveRDS(foo, file = "./output/SLOneStep.rds")

# load
# foo<-readRDS("./output/SLOneStep.rds")

# output of a single step:
foo%>%names()
# the preds outputs are based on sample (rows) and time (columns)
foo$preds%>%names()
foo$preds$survSL.survnet%>%head()
foo$preds$sl%>%head()
# note that this is exactly the nrow(OOB):
foo$preds$sl%>%dim()
nrow(df)*(1-0.632)

# scores are a list based on time
foo$scores%>%class()
foo$scores%>%length()
# scores for t=1
foo$scores$`t1`%>%arrange(Brier)%>% 
  mutate_if(is.numeric, round,digits=3)
# scores for t=5
foo$scores$`t5`%>%arrange(Brier)%>% 
  mutate_if(is.numeric, round,digits=3)
# coefs are simply the calculated event coefficients
# ie sl is the weighted sum of the library predictions, weighted by these coefs...
foo$coefs

# Let's see if they calibrate reasonably well...
calib_score<-Score(foo$preds,
          data=sample$cv_test[[1]],
          formula=Hist(years,death)~1,times=1:5,null.model=F,se.fit=F,metrics="Brier",
          plots = "Calibration")

# see them all
plotCalibration(calib_score,times=5)
# super learner only
plotCalibration(calib_score,times=5,models="sl")

# bar for SL
plotCalibration(calib_score,times=5,models = "sl",bars=T)

# ----------------------------------------------------------------------- #
# Step 3: MAIN LOOP
# ----------------------------------------------------------------------- #

# This would take a long time

# note time
main_st<-Sys.time()

# one simply iterates the one step function B times over the entire sample 
# to perform the full algorithm

# create empty lists that will contain the corresponding results for each subsample
all_results<-list()
all_preds<-list()
all_scores<-list()
all_coefs<-list()

# iterate from 1 to B: the total number of subsamples
k<-B
for(i in 1:k){
  # print
  cat("NOW COMPUTING SAMPLE", i, "...\n \n")
  all_results[[i]]<-one_step(samples_list=sample,B=i,times=1:5,
                             event_library=event.SL.library,cens_library=cens.SL.library)
  all_preds[[i]]<-all_results[[i]]$preds
  all_scores[[i]]<-all_results[[i]]$scores
  all_coefs[[i]]<-all_results[[i]]$coefs
}

names(all_coefs)<-paste("sample",1:k,sep="")
names(all_preds)<-paste("sample",1:k,sep="")
names(all_scores)<-paste("sample",1:k,sep="")

# end time
main_et<-Sys.time()
# print time
print(main_et-main_st)

# coefficients for subsamples...
all_coefs$sample1
all_coefs$sample2
all_coefs$sample3

# sl predictions for sample 1
all_preds$sample1$sl%>%head()

# scores for sample 1, time 5
all_scores$sample1$t5%>%arrange(Brier)

# Let's compute the average Brier scores of all methods for time 5
all_scores_matrix<-matrix(NA,nrow=0,ncol=length(event.SL.library)+1)

for(i in 1:k){
  all_scores_matrix<-rbind(all_scores_matrix,all_scores[[i]]$t5$Brier)
}

rownames(all_scores_matrix)<-paste("sample",1:k,sep="")
colnames(all_scores_matrix)<-all_preds$sample1%>%names()

# average Brier scores
all_scores_matrix%>%colMeans()%>%round(3)

# ----------------------------------------------------------------------- #
# FIXES (DO NOT RUN)
# ----------------------------------------------------------------------- #
