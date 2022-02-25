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
library(psrwe)
library(xtable)

load("~/Dropbox/Mayo-intern/aim2/data/allTrials.RData") 

times <- all.trial$fu_mos
delta <- all.trial$status
## control
index = which(all.trial$arm == 0)
x.train.con <- all.trial[index, c("STUDY", "SEX_ID", "ps", "WBC", "age_floor", "trant")]
times.con <- times[index]
delta.con <- delta[index]
## treatment
index.trt = which(all.trial$arm == 1 & all.trial$STUDY == 1203)
times.trt <- times[index.trt]
delta.trt <- delta[index.trt]
event.trt = as.numeric(times.trt <= 12 & delta.trt == 1)

miss.index = which(is.na(x.train.con$ps) | is.na(x.train.con$WBC))
x.train.con = x.train.con[-miss.index, ]
times.con = times.con[-miss.index]
delta.con = delta.con[-miss.index]
event = as.numeric(times.con <= 12 & delta.con == 1)
history = as.numeric(x.train.con$STUDY == 1203)
x.train.con$history = history
x.train.con = as.data.frame(x.train.con)
prop.fit = glm(history ~ SEX_ID + ps + WBC + age_floor + trant, family = 'binomial', data = x.train.con)
prop.score = prop.fit$fitted.values
dist.1203 = prop.score[which(x.train.con$STUDY == 1203)]
dist.106 = prop.score[which(x.train.con$STUDY == 106)]
dist.10201 = prop.score[which(x.train.con$STUDY == 10201)]
dist.10603 = prop.score[which(x.train.con$STUDY == 10603)]
dist.1900 = prop.score[which(x.train.con$STUDY == 1900)]
dist = list(dist.1203, dist.106, dist.10201, dist.10603, dist.1900)
cor.dist = diag(5)
for (i in 1:4) {
  for (j in (i + 1):5) {
    cor.dist[i, j] = cor.dist[j, i] = get_distance(dist[[i]], dist[[j]], 'ovl')
  }
}
omega = solve(cor.dist) ## correlation matrix

size = c(length(dist.1203), length(dist.106), length(dist.10201), length(dist.10603), length(dist.1900))
y = c(sum(event[which(x.train.con$STUDY == 1203)]), 
      sum(event[which(x.train.con$STUDY == 106)]),
      sum(event[which(x.train.con$STUDY == 10201)]),
      sum(event[which(x.train.con$STUDY == 10603)]),
      sum(event[which(x.train.con$STUDY == 1900)]))
year = c(11, 2, 1, 6, 0)

## fit the mixed effects model for control group
fit = mix_effect(y = y, year = year, omega = omega, size = size, mu = rep(0, 5), n.adapt = 5000, n.burn = 5000, n.iter = 10000)
alpha_post = fit$alpha_post
beta_post = fit$beta_post
effect_post = fit$effect_post
tau_post = fit$tau_post

## posterior probability
temp = rowMeans(alpha_post) + mean(beta_post) * year + rowMeans(effect_post)
post_prob = 1 / (1 + exp(-temp))

## fit the naive model for the trreatment group
alpha.trt = naive(y = sum(event.trt), size = length(event.trt), n.adapt = 5000, n.burn = 5000, n.iter = 10000)
## posterior probability
temp = mean(alpha.trt)
post_prob_trt = 1 / (1 + exp(-temp))

post_prob_trt - post_prob[1] 







## previous codes for survival analysis

x.train.con$SEX_ID = as.numeric(x.train.con$SEX_ID) - 1
x.train.con$ps1 = as.numeric(x.train.con$ps == 1)
x.train.con$ps2 = as.numeric(x.train.con$ps == 2)
x.train.con$trant = as.numeric(x.train.con$trant) - 1
x.train.con = data.matrix(x.train.con[, -3])

## Survival data processing
surv.proc = surv.pre.bart(times = times.con, delta = delta.con, x.train = x.train.con, x.test = x.train.con, K = 20) 
y.train = surv.proc$y.train ## 15295
x.train = surv.proc$tx.train ## 15295 * 8
n = length(y.train)
size = rep(0, 5)
study = c(1203, 106, 10201, 10603, 1900)
for (i in 1:5) {
  size[i] = sum(x.train[, 2] == study[i])
}
for (i in 2:5) {
  size[i] = size[i - 1] + size[i]
}
x.train = x.train[, -2]

fit = hist_match(y = y.train, x = x.train, n1 = n1, n = n, n.adapt = 1000, n.burn = 1000, n.iter = 1000)
w_post = fit$w_post
u0_post = fit$u0_post
u1_post = fit$u1_post
tau_post = fit$tau_post

### survival probability
w_est = c(rep(0, n1), rowMeans(w_post))
beta0 = rowMeans(u0_post)
beta1 = rowMeans(u1_post)
tau = mean(tau_post)
time_int = unique(x.train[, 1])
prob_pred = rep(0, 20)
for (i in 1:20) {
  test.index = which(x.train[, 1] == time_int[i])
  mu_pred = (1 - w_est[test.index]) * (x.train[test.index, ] %*% beta0) + w_est[test.index] * (x.train[test.index, ] %*% beta1)
  prob_pred[i] = mean(pnorm(0, mu_pred, 1 / tau))
}
plot(time_int, prob_pred, type = "l")

