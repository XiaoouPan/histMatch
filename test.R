rm(list = ls())
setwd("~/Dropbox/Mayo-intern/aim2/data")

library(dplyr)
library(BART)
library(readxl)
library(survival)
library(forestplot)
library(tikzDevice)
library(caret)

load("allTrials.RData") 

################################## propensity score matching #############################################
library(randomForest)
times <- all.trial$fu_mos
delta <- all.trial$status
study = all.trial$STUDY
all.trial$history = as.factor(all.trial$STUDY != 1203)
index = which(all.trial$arm == 0)
x.train.con <- all.trial[index, c("history", "SEX_ID", "ps", "WBC", "age_floor", "trant")]
times.con <- times[index]
delta.con <- delta[index]
## x.test.con for proximity score calculation
x.test.con = all.trial[index, c("SEX_ID", "ps", "WBC", "age_floor", "trant")]
study = study[index]

sapply(x.train.con, function(y) sum(is.na(y)))
miss.index = which(is.na(x.train.con$ps) | is.na(x.train.con$WBC))
x.train.con = x.train.con[-miss.index, ]
times.con = times.con[-miss.index]
delta.con = delta.con[-miss.index]
x.test.con = x.test.con[-miss.index, ]
study = study[-miss.index]

## A toy case with 20 samples
toy.index = c(1:10, 255:264)
x.train.con = x.train.con[toy.index, ]
times.con = times.con[toy.index]
delta.con = delta.con[toy.index]
x.test.con = x.test.con[toy.index, ]
study = study[toy.index]

rf = randomForest(
  history ~ .,
  data = x.train.con,
  ntree = 500,
  proximity = TRUE,
  na.action = na.omit
)
#rf
#plot(rf)
score = predict(rf, newdata = x.test.con)  
prox = rf$proximity ## proximity score: a 1457 * 1457 matrix

surv.proc = surv.pre.bart(times = times.con, delta = delta.con, x.train = x.train.con, x.test = x.train.con) 



