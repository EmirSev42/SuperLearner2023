### 02_RFS_Tuning.R ###
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
## Author: Pietro Ravani
## Created: Friday, Jan 27, 2023, (14:22) 
## Version: 
## Last-Updated: June 15, 2023
##           By: Emir S
##     Update #: 1
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
## 
### Commentary: range of hyperparameter values for RFS
### 
## Loading libraries
library(riskRegression)
library(survival)
library(tidyverse)
library(randomForestSRC)
library(data.table)

#----------------------------------------------------------------------
# STEP 1: LOAD DATA
#----------------------------------------------------------------------

df<-read.csv("./dataderived/SynthData.csv")
# select relevant columns
df06<-df%>%select(c("years","death","age","gfr21","pACRc","male","dm","cvd"))

# random sample, for faster run, to check if the codes work
n<-10000
set.seed(42)
df06<-df06[sample(1:nrow(df06),size=n,replace=F),]

# ----------------------------------------------------------------------- #
# STEP 2: TUNE MTRY AND NODESIZE
# ----------------------------------------------------------------------- #

# built in function that selects optimal mtry and nodesize
 system.time({
 o06 <- tune(Surv(years,death) ~ ., df06, doBest = TRUE,trace=T)
 }) # ? minutes
# save(o06, file="./models/o06.rda")
 
 # see selected parameters
 o06$optimal
 
 # put results of the search into dataframe for ggplot
 X<-as.data.frame(o06$results)
 # mtry used:
 X$mtry%>%unique()
 
# see performance for different mtry
X$mtry <- as.factor(X$mtry)
 ggplot(data=X, aes(x=nodesize,y=err, fill=mtry))+geom_point(aes(col=mtry))+
   geom_line(aes(col=mtry))
 
 # see them separately
 ggplot(data=X, aes(x=nodesize,y=err, fill=mtry))+geom_point(aes(col=mtry))+
   geom_line(aes(col=mtry))+facet_wrap(~mtry)
 
 # entries with the lowest errors
 X%>%arrange(err)%>%head()
 # mtry 2 and 3 are the winners

 # ----------------------------------------------------------------------- #
 # STEP 3: N OF TREES
 # ----------------------------------------------------------------------- #
 
 
 # 50 trees vs 100
 Wk050.2 <- rfsrc(formula = Surv(years,death) ~ age + male + gfr21 + pACRc + dm + cvd,
                  data = df06, mtry=2,
                  ntree = 50, nodesize = o06$optimal[["nodesize"]], nsplit = 10, seed = 132)
 Wk050.3 <- rfsrc(formula = Surv(years,death) ~ age + male + gfr21 + pACRc + dm + cvd,
                  data = df06, mtry=3,
                  ntree = 50, nodesize = o06$optimal[["nodesize"]], nsplit = 10, seed = 132)
 Wk100.2 <- rfsrc(formula = Surv(years,death) ~ age + male + gfr21 + pACRc + dm + cvd,
                  data = df06, mtry=2,
                  ntree = 100, nodesize = o06$optimal[["nodesize"]], nsplit = 10, seed = 132)
 Wk100.3 <- rfsrc(formula = Surv(years,death) ~ age + male + gfr21 + pACRc + dm + cvd,
                  data = df06, mtry=3,
                  ntree = 100, nodesize = o06$optimal[["nodesize"]], nsplit = 10, seed = 132)
 
 rsfmod <- list(Wk050.2=Wk050.2,Wk050.3=Wk050.3,
                Wk100.2=Wk100.2,Wk100.3=Wk100.3
 )
 
 # call score function with loob
 # NOTE: can be done in parallel
 x <- Score(rsfmod,
            formula=Surv(years,death)~1,
            data=df06,
            split.method="loob",
            B=100,
            M=round(.632*NROW(df06)),
            seed=132,
            times=5,
            se.fit=FALSE,
            conf.int=FALSE,
            metrics="brier",
            null.model = F)
            
 # arrange by Brier score to see winner
 x$Brier$score%>%arrange(Brier)
