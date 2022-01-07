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
index.cur = which(study == 1203)
m = length(index.cur)
match = rep(0, m)
for (i in 1:m) {
  max.ind = which.max(prox[i, -index.cur])
  match[i] = m + max.ind
}
study = as.factor(study)
levels(study) = c("S0106", "S1203", "E1900", "C10201", "C10603")
table(study[match])


#################################### historical borrowing with 1 trial #####################################
trial.num = c(106, 10201, 10603, 1900)
trial.name = c("dummy_s0106", "dummy_c10201", "dummy_c10603", "dummy_e1900")
trial.id = c("S0106", "C10201", "C10603", "E1900")
num = 4
index = which(all.trial$STUDY == 1203  | all.trial$STUDY == trial.num[num])
index.test = which(all.trial$STUDY == 1203)
times <- all.trial$fu_mos[index]
delta <- all.trial$status[index]
x.train <- all.trial[index, c("arm", "SEX_ID", "ps", "WBC", "npm1", "flt3", "age_floor", "Cyto", "trant", trial.name[num])]
x.train.trt <- data.matrix(filter(x.train, arm == 1)[, -1])
x.train.con <- data.matrix(filter(x.train, arm == 0)[, -1])
x.test <- all.trial[index.test, c("arm", "SEX_ID", "ps", "WBC", "npm1", "flt3", "age_floor", "Cyto", "trant", trial.name[num])]
x.test.trt <- data.matrix(filter(x.test , arm == 1)[, -1])
x.test.con <- data.matrix(filter(x.test , arm == 0)[, -1])
times.trt <- times[x.train$arm == 1]
times.con <- times[x.train$arm == 0]
delta.trt <- delta[x.train$arm == 1]
delta.con <- delta[x.train$arm == 0]

bart.4.trt <- mc.surv.bart(x.train = x.train.trt, times = times.trt, 
                           delta = delta.trt, x.test = x.test.trt, 
                           K = 40, ntree = 100, mc.cores = 8, seed = 321) 
bart.4.con <- mc.surv.bart(x.train = x.train.con, times = times.con, 
                           delta = delta.con, x.test = x.test.con, 
                           K = 40, ntree = 100, mc.cores = 8, seed = 321) 

