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

find_bart_surv_mat <- function(bart.obj) { # input is a BART object 
  surv.mc.mat <- bart.obj$surv.test 
  bart.surv.mat <- matrix(NA, nrow = dim(surv.mc.mat)[1], ncol = length(bart.obj$times)) 
  time.index <- seq(1, dim(surv.mc.mat)[2], 
                    by = length(bart.obj$times)) #by K (number of time pts) 
  for (i in seq(length(bart.obj$times)) - 1) {
    bart.surv.mat[, i+1] <- rowMeans(surv.mc.mat[, time.index + i, drop = FALSE])
  }
  bart.surv.mat
}

find_bart_stat <- function(bart.surv.vec, bart.obj) { 
  # input is a bart survival probability vector and a BART object 
  mc.mat.t <- find_bart_surv_mat(bart.obj)
  
  med.os <- bart.obj$times[bart.surv.vec <= 0.5][1] 
  med.os.lower <- bart.obj$times[apply(mc.mat.t, 2, 
                                       function(x) quantile(x, probs = 0.025)) <= 0.5][1] 
  med.os.upper <- bart.obj$times[apply(mc.mat.t, 2, 
                                       function(x) quantile(x, probs = 0.975)) <= 0.5][1] 
  # med.os.lower <- round(median(apply(mc.mat.t, 2, function(x) quantile(x, probs = 0.025))), 3)
  # med.os.upper <- round(median(apply(mc.mat.t, 2, function(x) quantile(x, probs = 0.975))), 3)
  med.os.combo <- paste0(med.os, " (", med.os.lower, ", ", med.os.upper, ")")
  
  one.yr.os <- round(bart.surv.vec[bart.obj$times >= 12][1], 3)
  one.yr.os.lower <- round(quantile(mc.mat.t[, bart.obj$times >= 12][, 1], probs = 0.025), 3)
  one.yr.os.upper <- round(quantile(mc.mat.t[, bart.obj$times >= 12][, 1], probs = 0.975), 3)
  one.yr.os.combo <- paste0(one.yr.os, " (", one.yr.os.lower, ", ", one.yr.os.upper, ")") 
  
  two.yr.os <- round(bart.surv.vec[bart.obj$times >= 24][1], 3)
  two.yr.os.lower <- round(quantile(mc.mat.t[, bart.obj$times >= 24][, 1], probs = 0.025), 3)
  two.yr.os.upper <- round(quantile(mc.mat.t[, bart.obj$times >= 24][, 1], probs = 0.975), 3)
  two.yr.os.combo <- paste0(two.yr.os, " (", two.yr.os.lower, ", ", two.yr.os.upper, ")") 
  
  c(Median_OS = med.os.combo, One_year_OS = one.yr.os.combo, Two_year_OS = two.yr.os.combo)
}

bart_confint_finder <- function(bart.obj.con, bart.obj.trt, year) {
  # input is two bart objects 
  # year is 1 or 2 for confidence intervals of the difference in survival at 1 year or 2 year 
  # Note: the length of the objects below are different for the control and the treatment survfit objects 

  mons = year * 12
  
  # find index 
  index.con <- length(bart.obj.con$times[bart.obj.con$times < mons]) + 1 
  index.trt <- length(bart.obj.trt$times[bart.obj.trt$times < mons]) + 1 
  
  bart.surv.mat.con <- find_bart_surv_mat(bart.obj.con)
  bart.surv.mat.trt <- find_bart_surv_mat(bart.obj.trt) 
  
  surv.diff.coln <- bart.surv.mat.trt[, index.trt] - bart.surv.mat.con[, index.con] 
  surv.diff <- mean(surv.diff.coln)
  surv.ci <- quantile(surv.diff.coln, probs = c(0.025, 0.975)) 
  
  c(surv.diff, surv.ci)
}

concor = function(bart.obj, times, delta, n, year = 1) {
  sVal = matrix(bart.obj$surv.test.mean, n, bart.obj$K, byrow = TRUE)
  index = length(bart.obj$times[bart.obj$times < year * 12]) + 1 
  score = sVal[, index]
  num = 0
  den = 0
  for (i in 1:n) {
    if (times[i] <  year * 12) {
      next
    }
    for (j in 1:n) {
      if (i == j || times[j] <  year * 12 || delta[j] == 0 || times[j] >= times[i]) {
        next
      }
      den = den + 1
      if (score[j] < score[i]) {
        num = num + 1
      }
    }
  }
  return (num / den)
}


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



## Some statistics, not necessary for forest plots
time.index <- seq(1, length(bart.1.con$surv.test.mean), 
                  by = length(bart.1.con$times)) #by K (number of time pts) 
bart.1.con.surv <- rep(NA, length(bart.1.con$times)) 
for (i in seq(length(bart.1.con$times)) - 1) { 
  bart.1.con.surv[i+1] <- mean(bart.1.con$surv.test.mean[time.index + i], 
                               na.rm = TRUE)
}
time.index <- seq(1, length(bart.1.trt$surv.test.mean), 
                  by = length(bart.1.trt$times)) #by K (number of time pts) 
bart.1.trt.surv <- rep(NA, length(bart.1.trt$times)) 
for (i in seq(length(bart.1.trt$times)) - 1) { 
  bart.1.trt.surv[i+1] <- mean(bart.1.trt$surv.test.mean[time.index + i], 
                               na.rm = TRUE)
}

bart.1.stat <- cbind(find_bart_stat(bart.1.con.surv, bart.1.con), 
                     find_bart_stat(bart.1.trt.surv, bart.1.trt))
colnames(bart.1.stat) <- c("Control", "Treatment")
bart.1.stat <- rbind(N = c(sum(filter(one.historical, arm == 0)$status), sum(filter(one.historical, arm == 1)$status)), 
                     bart.1.stat)

#################################### Forest plots ######################################

####### one year survival forest plot ########

confint.df <- readRDS(file = "conint_res_one_yr.rds")

year = 1.5
confint.df <- rbind(confint.df, 
                    "1203.BART.hist.10201" = bart_confint_finder(bart.2.con, bart.2.trt, year = year),
                    "1203.BART.hist.10603" = bart_confint_finder(bart.3.con, bart.3.trt, year = year),
                    "1203.BART.hist.0106" = bart_confint_finder(bart.1.con, bart.1.trt, year = year),
                    "1203.BART.hist.1900" = bart_confint_finder(bart.4.con, bart.4.trt, year = year),
                    "1203.BART.hist.all" = bart_confint_finder(bart.all.con, bart.all.trt, year = year))

tabletext <- cbind(c("Study", "C10201", "", "C10603", "", "S0106", "", "E1900", "", "S1203", "", "", "", "", "", ""), 
                   c("Method", rep(c("KM", "BART"), 5), paste("BART+", trial.id[2], sep = ""), paste("BART+", trial.id[3], sep = ""), 
                     paste("BART+", trial.id[1], sep = ""), paste("BART+", trial.id[4], sep = ""), "BART+ALL"), 
                   c("S(1)-S(0)", round(confint.df[, "Survival Difference"], 3)))

confint.forestplot <- data.frame(mean = c(NA, round(confint.df[, "Survival Difference"], 3)), 
                                 lower = c(NA, confint.df[, "Lower"]), 
                                 upper = c(NA, confint.df[, "Upper"]))

# box size is based on precision 
forestplot.bart1 <- confint.forestplot %>% 
  forestplot(labeltext = tabletext, 
             is.summary = c(TRUE, rep(FALSE, 15)), 
             xlab = "Difference in 1 yr survival (S(1)-S(0))")

tikz("plot.tex", standAlone = TRUE, width = 7, height = 5)
forestplot.bart1
dev.off() 
tools::texi2dvi("plot.tex", pdf = T)


####### two year survival forest plot ########

confint.df <- readRDS(file = "conint_res_two_yr.rds")

confint.df <- rbind(confint.df, 
                    "1203.BART.hist.10201" = bart_confint_finder(bart.2.con, bart.2.trt, year = 2),
                    "1203.BART.hist.10603" = bart_confint_finder(bart.3.con, bart.3.trt, year = 2),
                    "1203.BART.hist.0106" = bart_confint_finder(bart.1.con, bart.1.trt, year = 2),
                    "1203.BART.hist.1900" = bart_confint_finder(bart.4.con, bart.4.trt, year = 2),
                    "1203.BART.hist.all" = bart_confint_finder(bart.all.con, bart.all.trt, year = 2))

tabletext <- cbind(c("Study", "C10201", "", "C10603", "", "S0106", "", "E1900", "", "S1203", "", "", "", "", "", ""), 
                   c("Method", rep(c("KM", "BART"), 5), paste("BART+", trial.id[2], sep = ""), paste("BART+", trial.id[3], sep = ""), 
                     paste("BART+", trial.id[1], sep = ""), paste("BART+", trial.id[4], sep = ""), "BART+ALL"), 
                   c("S(1)-S(0)", round(confint.df[, "Survival Difference"], 3)))

confint.forestplot <- data.frame(mean = c(NA, round(confint.df[, "Survival Difference"], 3)), 
                                 lower = c(NA, confint.df[, "Lower"]), 
                                 upper = c(NA, confint.df[, "Upper"]))

# box size is based on precision 
two.yr.forestplot.bart1 <- confint.forestplot %>% 
  forestplot(labeltext = tabletext, 
             is.summary = c(TRUE, rep(FALSE, 15)), 
             xlab = "Difference in 2 yr survival (S(1)-S(0))")

tikz("plot.tex", standAlone = TRUE, width = 7, height = 5)
two.yr.forestplot.bart1
dev.off() 
tools::texi2dvi("plot.tex", pdf = T)



####### Forest plots with historical borrowing based on proximity score matching #####

confint.df <- readRDS(file = "conint_res_one_yr.rds")

year = 1
confint.df <- rbind(confint.df, 
                    "1203.BART.matching" = bart_confint_finder(bart.con, bart.trt, year = year))

tabletext <- cbind(c("Study", "C10201", "", "C10603", "", "S0106", "", "E1900", "", "S1203", "", ""), 
                   c("Method", rep(c("KM", "BART"), 5), "BART+Matching"), 
                   c("S(1)-S(0)", round(confint.df[, "Survival Difference"], 3)))

confint.forestplot <- data.frame(mean = c(NA, round(confint.df[, "Survival Difference"], 3)), 
                                 lower = c(NA, confint.df[, "Lower"]), 
                                 upper = c(NA, confint.df[, "Upper"]))

# box size is based on precision 
forestplot.bart1 <- confint.forestplot %>% 
  forestplot(labeltext = tabletext, 
             is.summary = c(TRUE, rep(FALSE, 15)), 
             xlab = "Difference in 1 yr survival (S(1)-S(0))")

tikz("plot.tex", standAlone = TRUE, width = 7, height = 5)
forestplot.bart1
dev.off() 
tools::texi2dvi("plot.tex", pdf = T)


####### two year survival forest plot ########

confint.df <- readRDS(file = "conint_res_two_yr.rds")

confint.df <- rbind(confint.df, 
                    "1203.BART.matching" = bart_confint_finder(bart.con, bart.trt, year = 2))

tabletext <- cbind(c("Study", "C10201", "", "C10603", "", "S0106", "", "E1900", "", "S1203", "", ""), 
                   c("Method", rep(c("KM", "BART"), 5), "BART+Matching"), 
                   c("S(1)-S(0)", round(confint.df[, "Survival Difference"], 3)))

confint.forestplot <- data.frame(mean = c(NA, round(confint.df[, "Survival Difference"], 3)), 
                                 lower = c(NA, confint.df[, "Lower"]), 
                                 upper = c(NA, confint.df[, "Upper"]))

# box size is based on precision 
two.yr.forestplot.bart1 <- confint.forestplot %>% 
  forestplot(labeltext = tabletext, 
             is.summary = c(TRUE, rep(FALSE, 15)), 
             xlab = "Difference in 2 yr survival (S(1)-S(0))")

tikz("plot.tex", standAlone = TRUE, width = 7, height = 5)
two.yr.forestplot.bart1
dev.off() 
tools::texi2dvi("plot.tex", pdf = T)






## Extract sample from s1203 for testing, and create a trial with insufficient individuals in control group
index = 1:471 
folds = createFolds(all.trial$arm[index], 5, FALSE)
test.index = which(folds == 1)

times <- all.trial$fu_mos[-test.index]
delta <- all.trial$status[-test.index]
x.train <- all.trial[-test.index, c("arm", "SEX_ID", "ps", "WBC", "npm1", "flt3", "age_floor", "Cyto", "trant", "dummy_s0106",
                         "dummy_c10201", "dummy_c10603", "dummy_e1900")]
x.train.trt <- data.matrix(filter(x.train, arm == 1)[, -1])
x.train.con <- data.matrix(filter(x.train, arm == 0)[, -1])
times.trt <- times[x.train$arm == 1]
times.con <- times[x.train$arm == 0]
delta.trt <- delta[x.train$arm == 1]
delta.con <- delta[x.train$arm == 0]
x.test <- all.trial[test.index, c("arm", "SEX_ID", "ps", "WBC", "npm1", "flt3", "age_floor", "Cyto", "trant", "dummy_s0106",
                                    "dummy_c10201", "dummy_c10603", "dummy_e1900")]
x.test.trt <- data.matrix(filter(x.test, arm == 1)[, -1])
x.test.con <- data.matrix(filter(x.test, arm == 0)[, -1])
times.test <- all.trial$fu_mos[test.index]
delta.test <- all.trial$status[test.index]
times.test.trt <- times.test[x.test$arm == 1]
times.test.con <- times.test[x.test$arm == 0]
delta.test.trt <- delta.test[x.test$arm == 1]
delta.test.con <- delta.test[x.test$arm == 0]

## Proportional hazard
df = data.frame(x.train.trt)
df$times = times.trt
df$delta = delta.trt
df$SEX_ID = as.numeric(df$SEX_ID) - 1
fit.coxph <- coxph(Surv(times, delta) ~ SEX_ID + ps + WBC + npm1 + flt3 + age_floor + Cyto + trant, data = df)


## Run BART
K = 40; ntree = 100
bart.1.trt <- mc.surv.bart(x.train = x.train.trt, times = times.trt, 
                           delta = delta.trt, x.test = x.test.trt, 
                           K = K, ntree = ntree, mc.cores = 8, seed = 321) 
concor(bart.1.trt, times.test.trt, delta.test.trt, length(times.test.trt), year = 0.5)

K = 40; ntree = 100
bart.1.con <- mc.surv.bart(x.train = x.train.con, times = times.con, 
                           delta = delta.con, x.test = x.test.con, K = K, ntree = ntree, mc.cores = 8, seed = 321) 
#pre = surv.pre.bart(x.train = x.train.con, times = times.con, delta = delta.con, K = K, x.test = x.test.con)
#bart.pred = predict(bart.1.con, newdata = pre$tx.test, mc.cores = 8)
concor(bart.1.con, times.test.con, delta.test.con, length(times.test.con), year = 0.5)


## Feature importance
M = 20
feature.imp = matrix(0, 4, M)
pb = txtProgressBar(style = 3)
for (m in 1:M) {
  pert = rbinom(nrow(x.test.con), 1, 0.5)
  ## S0106
  x.test.pert = x.test.con
  x.test.pert[, 9] = pert
  pre = surv.pre.bart(x.train = x.train.con, times = times.con, delta = delta.con, K = K, x.test = x.test.pert)
  bart.pred = predict(bart.1.con, newdata = pre$tx.test, mc.cores = 8)
  feature.imp[1, m] = concor(bart.pred, times.test.con, delta.test.con, length(times.test.con), year = 0.5)
  ## C10201
  x.test.pert = x.test.con
  x.test.pert[, 10] = pert
  pre = surv.pre.bart(x.train = x.train.con, times = times.con, delta = delta.con, K = K, x.test = x.test.pert)
  bart.pred = predict(bart.1.con, newdata = pre$tx.test, mc.cores = 8)
  feature.imp[2, m] = concor(bart.pred, times.test.con, delta.test.con, length(times.test.con), year = 0.5)
  ## C10603
  x.test.pert = x.test.con
  x.test.pert[, 11] = pert
  pre = surv.pre.bart(x.train = x.train.con, times = times.con, delta = delta.con, K = K, x.test = x.test.pert)
  bart.pred = predict(bart.1.con, newdata = pre$tx.test, mc.cores = 8)
  feature.imp[3, m] = concor(bart.pred, times.test.con, delta.test.con, length(times.test.con), year = 0.5)
  ## E1900
  x.test.pert = x.test.con
  x.test.pert[, 12] = pert
  pre = surv.pre.bart(x.train = x.train.con, times = times.con, delta = delta.con, K = K, x.test = x.test.pert)
  bart.pred = predict(bart.1.con, newdata = pre$tx.test, mc.cores = 8)
  feature.imp[4, m] = concor(bart.pred, times.test.con, delta.test.con, length(times.test.con), year = 0.5)
  
  setTxtProgressBar(pb, m / M)
}

rowMeans(feature.imp)


## S0106
x.test.pert = x.test.con
x.test.pert[, 9] = 1
pre = surv.pre.bart(x.train = x.train.con, times = times.con, delta = delta.con, K = K, x.test = x.test.pert)
bart.pred = predict(bart.1.con, newdata = pre$tx.test, mc.cores = 8)
concor(bart.pred, times.test.con, delta.test.con, length(times.test.con), year = 0.5)
## C10201
x.test.pert = x.test.con
x.test.pert[, 10] = 1
pre = surv.pre.bart(x.train = x.train.con, times = times.con, delta = delta.con, K = K, x.test = x.test.pert)
bart.pred = predict(bart.1.con, newdata = pre$tx.test, mc.cores = 8)
concor(bart.pred, times.test.con, delta.test.con, length(times.test.con), year = 0.5)
## C10603
x.test.pert = x.test.con
x.test.pert[, 11] = 1
pre = surv.pre.bart(x.train = x.train.con, times = times.con, delta = delta.con, K = K, x.test = x.test.pert)
bart.pred = predict(bart.1.con, newdata = pre$tx.test, mc.cores = 8)
concor(bart.pred, times.test.con, delta.test.con, length(times.test.con), year = 0.5)
## E1900
x.test.pert = x.test.con
x.test.pert[, 12] = 1
pre = surv.pre.bart(x.train = x.train.con, times = times.con, delta = delta.con, K = K, x.test = x.test.pert)
bart.pred = predict(bart.1.con, newdata = pre$tx.test, mc.cores = 8)
concor(bart.pred, times.test.con, delta.test.con, length(times.test.con), year = 0.5)
