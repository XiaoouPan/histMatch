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

#################################### survival curves of control groups #####################################
getSurv = function(bart.con) {
  time.index <- seq(1, length(bart.con$surv.test.mean), by = length(bart.con$times)) #by K (number of time pts) 
  bart.con.surv <- rep(NA, length(bart.con$times)) 
  for (i in seq(length(bart.con$times)) - 1) { 
    bart.con.surv[i+1] <- mean(bart.con$surv.test.mean[time.index + i], na.rm = TRUE)
  }
  return (bart.con.surv)
}

trial.num = c(1203, 106, 10201, 10603, 1900)
trial.id = c("S1203", "S0106", "C10201", "C10603", "E1900")
num = 5
index = which(all.trial$STUDY == trial.num[num])
times <- all.trial$fu_mos[index]
delta <- all.trial$status[index]
x.train <- all.trial[index, c("arm", "SEX_ID", "ps", "WBC", "npm1", "flt3", "age_floor", "Cyto", "trant")]
x.train.con <- data.matrix(filter(x.train, arm == 0)[, -1])
times.con <- times[x.train$arm == 0]
delta.con <- delta[x.train$arm == 0]

bart.5.con <- mc.surv.bart(x.train = x.train.con, times = times.con, 
                           delta = delta.con, x.test = x.train.con, 
                           K = 40, ntree = 100, mc.cores = 8, seed = 321) 
bart.5.con.surv = getSurv(bart.5.con)




tikz("plot.tex", standAlone = TRUE, width = 7, height = 5)
plot(c(0, bart.1.con$times), c(1, bart.1.con.surv), type = "s", lwd = 5, cex = 1, col = "red", axes = FALSE, xlim = c(0, 121), 
     ylim = c(0, 1), xlab = "", ylab = "")
lines(c(0, bart.2.con$times), c(1, bart.2.con.surv), type = "s", lwd = 5, cex = 1, col = "blue") 
lines(c(0, bart.3.con$times), c(1, bart.3.con.surv), type = "s", lwd = 5, cex = 1, col = "darkorchid") 
lines(c(0, bart.4.con$times), c(1, bart.4.con.surv), type = "s", lwd = 5, cex = 1, col = "gold4") 
lines(c(0, bart.5.con$times), c(1, bart.5.con.surv), type = "s", lwd = 5, cex = 1, col = "seagreen") 
color = c("blue", "red", "gold4", "seagreen", "darkorchid")
labels = c("\\texttt{S0106}", "\\texttt{S1203}", "\\texttt{C10603}", "\\texttt{E1900}", "\\texttt{C10201}")
legend("topright", labels, col = color, lwd = 5, cex = 1.7, box.lwd = 1, bg = "white")
axis(1, line = 0, cex.axis = 1.5)
axis(2, line = 0, cex.axis = 1.5)
box()
title(xlab = "Time in Months", line = 3, cex.lab = 2)
title(ylab = "Survival Probability", line = 2.5, cex.lab = 1.6)
dev.off() 
tools::texi2dvi("plot.tex", pdf = T)



################################## propensity score matching #############################################
library(randomForest)
times <- all.trial$fu_mos
delta <- all.trial$status
study = all.trial$STUDY
all.trial$history = as.factor(all.trial$STUDY != 1203)
index = which(all.trial$arm == 0)
#x.train.con <- all.trial[index, c("history", "SEX_ID", "ps", "WBC", "npm1", "flt3", "age_floor", "Cyto", "trant")]
x.train.con <- all.trial[index, c("history", "SEX_ID", "ps", "WBC", "age_floor", "trant")]
times.con <- times[index]
delta.con <- delta[index]
#x.test.con = all.trial[index, c("SEX_ID", "ps", "WBC", "npm1", "flt3", "age_floor", "Cyto", "trant")]
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
rf
plot(rf)
score = predict(rf, newdata = x.test.con)
table(observed = x.train.con$history, predicted = score)


prox = rf$proximity
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
df.plot = data.frame(group = levels(study[match]), value = as.numeric(table(study[match])))
tikz("plot.tex", standAlone = TRUE, width = 7, height = 5)
ggplot(df.plot, aes(x = "", y = value, fill = group)) + geom_bar(stat = "identity") + coord_polar("y", start = 0) +
  xlab("Frequency") + ylab("") +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(0.5, "cm"),
        legend.background = element_rect(fill = alpha("white", 0)), axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20))
dev.off() 
tools::texi2dvi("plot.tex", pdf = T)


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

### borrowing with proximity score matching
trial.num = c(106, 10201, 10603, 1900)
trial.name = c("dummy_s0106", "dummy_c10201", "dummy_c10603", "dummy_e1900")
trial.id = c("S0106", "C10201", "C10603", "E1900")
index = union(which(all.trial$STUDY == 1203), match)
index.test = which(all.trial$STUDY == 1203)
times <- all.trial$fu_mos[index]
delta <- all.trial$status[index]
x.train <- all.trial[index, c("arm", "SEX_ID", "ps", "WBC", "age_floor", "trant")]
x.train.trt <- data.matrix(filter(x.train, arm == 1)[, -1])
x.train.con <- data.matrix(filter(x.train, arm == 0)[, -1])
x.test <- all.trial[index.test, c("arm", "SEX_ID", "ps", "WBC", "age_floor", "trant")]
x.test.trt <- data.matrix(filter(x.test , arm == 1)[, -1])
x.test.con <- data.matrix(filter(x.test , arm == 0)[, -1])
times.trt <- times[x.train$arm == 1]
times.con <- times[x.train$arm == 0]
delta.trt <- delta[x.train$arm == 1]
delta.con <- delta[x.train$arm == 0]

bart.trt <- mc.surv.bart(x.train = x.train.trt, times = times.trt, 
                           delta = delta.trt, x.test = x.test.trt, 
                           K = 40, ntree = 100, mc.cores = 8, seed = 321) 
bart.con <- mc.surv.bart(x.train = x.train.con, times = times.con, 
                           delta = delta.con, x.test = x.test.con, 
                           K = 40, ntree = 100, mc.cores = 8, seed = 321) 


#### all trials
index.test = which(all.trial$STUDY == 1203)
times <- all.trial$fu_mos
delta <- all.trial$status
x.train <- all.trial[, c("arm", "SEX_ID", "ps", "WBC", "npm1", "flt3", "age_floor", "Cyto", "trant", trial.name)]
x.train.trt <- data.matrix(filter(x.train, arm == 1)[, -1])
x.train.con <- data.matrix(filter(x.train, arm == 0)[, -1])
x.test <- all.trial[index.test, c("arm", "SEX_ID", "ps", "WBC", "npm1", "flt3", "age_floor", "Cyto", "trant", trial.name)]
x.test.trt <- data.matrix(filter(x.test , arm == 1)[, -1])
x.test.con <- data.matrix(filter(x.test , arm == 0)[, -1])
times.trt <- times[x.train$arm == 1]
times.con <- times[x.train$arm == 0]
delta.trt <- delta[x.train$arm == 1]
delta.con <- delta[x.train$arm == 0]

bart.all.trt <- mc.surv.bart(x.train = x.train.trt, times = times.trt, 
                           delta = delta.trt, x.test = x.test.trt, 
                           K = 40, ntree = 100, mc.cores = 8, seed = 321) 
bart.all.con <- mc.surv.bart(x.train = x.train.con, times = times.con, 
                           delta = delta.con, x.test = x.test.con, 
                           K = 40, ntree = 100, mc.cores = 8, seed = 321) 
