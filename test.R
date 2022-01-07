rm(list = ls())

library(dplyr)
library(BART)
library(readxl)
library(survival)
library(forestplot)
library(tikzDevice)
library(caret)
library(rjags)

load("~/Dropbox/Mayo-intern/aim2/data/allTrials.RData") 

################################## propensity score matching #############################################
library(randomForest)
times <- all.trial$fu_mos
delta <- all.trial$status
study = all.trial$STUDY
all.trial$history = as.factor(all.trial$STUDY != 1203)
index = which(all.trial$arm == 0)
#x.train.con <- all.trial[index, c("history", "SEX_ID", "ps", "WBC", "age_floor", "trant")]
x.train.con <- all.trial[index, c("history", "age_floor")]
times.con <- times[index]
delta.con <- delta[index]
## x.test.con for proximity score calculation
#x.test.con = all.trial[index, c("SEX_ID", "ps", "WBC", "age_floor", "trant")]
study = study[index]

sapply(x.train.con, function(y) sum(is.na(y)))
#miss.index = which(is.na(x.train.con$ps) | is.na(x.train.con$WBC))
x.train.con = x.train.con[-miss.index, ]
times.con = times.con[-miss.index]
delta.con = delta.con[-miss.index]
x.test.con = x.test.con[-miss.index, ]
study = study[-miss.index]

## A toy case with 20 samples
toy.index = c(1:10, 262:271)
x.train.con$history = as.numeric(x.train.con$history) - 1
x.train.con = data.matrix(x.train.con[toy.index, ])
times.con = times.con[toy.index]
delta.con = delta.con[toy.index]
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

## Survival data processing
surv.proc = surv.pre.bart(times = times.con, delta = delta.con, x.train = x.train.con, x.test = x.train.con, K = 3) 
y.train = surv.proc$y.train ## 40
x.train = surv.proc$tx.train ## 40 * 3
x.train = x.train[, c(1, 3)]

fit = hist_match(y = y.train, x = x.train, n1 = 19, n = 40, n.adapt = 1000, n.burn = 1000, n.iter = 1000)


