rm(list = ls())

library(dplyr)
library(BART)
library(readxl)
library(survival)
library(forestplot)
library(tikzDevice)
library(caret)
library(rjags)
library(randomForest)

load("~/Dropbox/Mayo-intern/aim2/data/allTrials.RData") 

times <- all.trial$fu_mos
delta <- all.trial$status
study = all.trial$STUDY
all.trial$history = as.factor(all.trial$STUDY != 1203)
index = which(all.trial$arm == 0)
x.train.con <- all.trial[index, c("history", "SEX_ID", "ps", "WBC", "age_floor", "trant")]
#x.train.con <- all.trial[index, c("history", "age_floor")]
times.con <- times[index]
delta.con <- delta[index]
## x.test.con for proximity score calculation
#x.test.con = all.trial[index, c("SEX_ID", "ps", "WBC", "age_floor", "trant")]
study = study[index]

#sapply(x.train.con, function(y) sum(is.na(y)))
miss.index = which(is.na(x.train.con$ps) | is.na(x.train.con$WBC))
x.train.con = x.train.con[-miss.index, ]
times.con = times.con[-miss.index]
delta.con = delta.con[-miss.index]
#x.test.con = x.test.con[-miss.index, ]
study = study[-miss.index]

## A toy case with 20 samples
toy.index = 1:400
x.train.con$history = as.numeric(x.train.con$history) - 1
x.train.con$SEX_ID = as.numeric(x.train.con$SEX_ID) - 1
x.train.con$ps = as.numeric(x.train.con$ps) - 1
x.train.con$trant = as.numeric(x.train.con$trant) - 1
#x.train.con = data.matrix(x.train.con)
x.train.con = data.matrix(x.train.con[toy.index, ])
times.con = times.con[toy.index]
delta.con = delta.con[toy.index]
study = study[toy.index]

## Survival data processing
surv.proc = surv.pre.bart(times = times.con, delta = delta.con, x.train = x.train.con, x.test = x.train.con, K = 20) 
y.train = surv.proc$y.train ## 40
x.train = surv.proc$tx.train ## 40 * 7
n = length(y.train)
n1 = sum(x.train[, 2] == 0)
x.train = x.train[, -2]

fit = hist_match(y = y.train, x = x.train, n1 = n1, n = n, n.adapt = 1000, n.burn = 1000, n.iter = 1000)
w_post = fit$w_post
u0_post = fit$u0_post
u1_post = fit$u1_post

### survival probability at one year
test.index = which(x.train[, 1] == 11.6725)
w_est = c(rep(0, n1), rowMeans(w_post))
beta0 = rowMeans(u0_post)
beta1 = rowMeans(u1_post)
mu_pred = (1 - w_est[test.index]) * (x.train[test.index, ] %*% beta0) + w_est[test.index] * (x.train[test.index, ] %*% beta1)
prob_pred = pnorm(0, mu_pred, 1)
mean(prob_pred)
