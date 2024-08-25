library(SuperLearner)
library(dplyr)
library(scales)
library(mgcv)

icu <- read.csv("icu_pseudo_data.csv")
icu <- icu %>% mutate_at(vars(site, male, sepsis_dx:winter, v_cc1:v_cc_r5), as.factor)

### relevant covariates

covariates = c("site", "age", "male", "sepsis_dx",
               "periarrest", "weekend", "winter", "icnarc_score", "news_score",
               "sofa_score", "open_bin", "v_cc1", "v_cc2", "v_cc3", 
               "v_cc4", "v_cc5", "v_cc_r1", "v_cc_r2", "v_cc_r3", "v_cc_r4", 
               "v_cc_r5")
covariates_sub = c("site", "age", "male", "sepsis_dx", "weekend", "winter")

dat <- data.frame(A = icu$icu_bed, Y = 1 - icu$dead28) # recode to Y = survival
dat <- cbind.data.frame(dat, icu[, covariates])
rm(icu)

n <- nrow(dat)
sl.lib=c("SL.glm","SL.ranger","SL.rpart")

#### Procedure with cross-fitting ####

set.seed(538)
M <- 4 ## number of data splits for cross-fitting
test.indices <- list()
remaining <- 1:n
for (m in 1:(M-1)) {
  test.indices[[m]] <- sample(remaining, floor(n/M), replace = F)
  remaining <- remaining[! (remaining %in% test.indices[[m]])]
}
test.indices[[M]] <- remaining

dat$mu.hat.0 <- dat$mu.hat.1 <- dat$pi.hat <- NA
dat$tau.sub.hat <- dat$xi.sub.hat <- NA
dat$beta.subcont.hat <- NA

for(j in 1:M) {
  inds.test <- test.indices[[j]]
  inds.train <- (1:n)[-inds.test]
  
  ### fit primary nuisance models, mu and pi
  mu.mod <- SuperLearner(Y = dat[inds.train,]$Y,
                         X = dat[inds.train, c("A", covariates)],
                         family = "binomial",
                         SL.library = sl.lib)
  pi.mod <- SuperLearner(Y = dat[inds.train,]$A,
                         X = dat[inds.train, covariates],
                         family = "binomial",
                         SL.library = sl.lib)
  
  dat.1 <- dat.0 <- dat[inds.test, c("A", covariates)]
  dat.0$A <- 0; dat.1$A <- 1
  
  dat$mu.hat.0[inds.test] <- predict.SuperLearner(mu.mod, dat.0)$pred
  dat$mu.hat.1[inds.test] <- predict.SuperLearner(mu.mod, dat.1)$pred
  dat$pi.hat[inds.test] <- 
    predict.SuperLearner(pi.mod, 
                         newdata = dat[inds.test,covariates], 
                         type = "response")$pred
  
  ### fit nuisance models for subset approaches, tau_w, xi_w, check{beta}_w
  dat.1.train <- dat.0.train <- dat[inds.train, c("A", covariates)]
  dat.0.train$A <- 0; dat.1.train$A <- 1
  tau.hats <- predict.SuperLearner(mu.mod, dat.1.train)$pred -
    predict.SuperLearner(mu.mod, dat.0.train)$pred
  pi.hats <- predict.SuperLearner(pi.mod)$pred
  beta.hats <- tau.hats * (1 * (tau.hats > 0) - pi.hats)
  
  tau.sub.mod <- SuperLearner(Y = tau.hats,
                              X = dat[inds.train, covariates_sub],
                              SL.library = sl.lib)
  xi.sub.mod <- SuperLearner(Y = tau.hats * pi.hats,
                             X = dat[inds.train, covariates_sub],
                             SL.library = sl.lib)
  beta.subcont.mod <- SuperLearner(Y = beta.hats,
                                   X = dat[inds.train, covariates_sub],
                                   SL.library = sl.lib)
  
  dat$tau.sub.hat[inds.test] <- 
    predict.SuperLearner(tau.sub.mod, 
                         newdata = dat[inds.test,covariates_sub], 
                         type = "response")$pred
  dat$xi.sub.hat[inds.test] <- 
    predict.SuperLearner(xi.sub.mod, 
                         newdata = dat[inds.test,covariates_sub], 
                         type = "response")$pred
  dat$beta.subcont.hat[inds.test] <- 
    predict.SuperLearner(beta.subcont.mod, 
                         newdata = dat[inds.test,covariates_sub], 
                         type = "response")$pred
  
}

delta.width = 0.00005
deltas <- seq(0, 1, by = delta.width)

dat$tau.hat <- dat$mu.hat.1 - dat$mu.hat.0
dat$h.star.hat <- 1*(dat$tau.hat > 0)
dat$beta.hat <- dat$tau.hat * (dat$h.star.hat - dat$pi.hat)
dat$c.hat <- dat$h.star.hat * (1 - dat$pi.hat) + 
  (1 - dat$h.star.hat) * dat$pi.hat

q.hats <- quantile(dat$beta.hat, probs = 1 - deltas)

Delta.star.hats <- sapply(q.hats, function(q.hat) {
  1 * (dat$beta.hat > q.hat)
}, simplify = 0)
colnames(Delta.star.hats) <- paste("delta = ", deltas)

dat$phi.hat <- (dat$h.star.hat - dat$pi.hat) * 
  (dat$A / dat$pi.hat - (1 - dat$A) / (1 - dat$pi.hat)) *
  (dat$Y - ifelse(dat$A == 1, dat$mu.hat.1, dat$mu.hat.0)) +
  dat$tau.hat * (dat$h.star.hat - dat$A)

value.hats <- mean(dat$Y) +
  sapply(1:length(deltas), function(j) {
    mean(Delta.star.hats[,j] * dat$phi.hat)
  }, simplify = 0)

var.hats <- sapply(1:length(deltas), function(j) {
  var(dat$Y + Delta.star.hats[,j] * (dat$phi.hat - q.hats[j]))
}, simplify = 0)

### using subset of covariates for both Delta and h

dat$h.star.sub.hat <- 1*(dat$tau.sub.hat > 0)
dat$beta.sub.hat <- dat$tau.sub.hat * dat$h.star.sub.hat - dat$xi.sub.hat

r.hats <- quantile(dat$beta.sub.hat, probs = 1 - deltas)

Delta.star.sub.hats <- sapply(r.hats, function(r.hat) {
  1 * (dat$beta.sub.hat > r.hat)
}, simplify = 0)
colnames(Delta.star.sub.hats) <- paste("delta = ", deltas)

dat$phi.sub.hat <- (dat$h.star.sub.hat - dat$pi.hat) * 
  (dat$A / dat$pi.hat - (1 - dat$A) / (1 - dat$pi.hat)) *
  (dat$Y - ifelse(dat$A == 1, dat$mu.hat.1, dat$mu.hat.0)) +
  dat$tau.hat * (dat$h.star.sub.hat - dat$A)

value.sub.hats <- mean(dat$Y) +
  sapply(1:length(deltas), function(j) {
    mean(Delta.star.sub.hats[,j] * dat$phi.sub.hat)
  }, simplify = 0)

var.sub.hats <- sapply(1:length(deltas), function(j) {
  var(dat$Y + Delta.star.sub.hats[,j] * (dat$phi.sub.hat - r.hats[j]))
}, simplify = 0)


### using subset of covariates for Delta only

s.hats <- quantile(dat$beta.subcont.hat, probs = 1 - deltas)

Delta.star.subcont.hats <- sapply(s.hats, function(s.hat) {
  1 * (dat$beta.subcont.hat > s.hat)
}, simplify = 0)
colnames(Delta.star.subcont.hats) <- paste("delta = ", deltas)

value.subcont.hats <- mean(dat$Y) +
  sapply(1:length(deltas), function(j) {
    mean(Delta.star.subcont.hats[,j] * dat$phi.hat)
  }, simplify = 0)

var.subcont.hats <- sapply(1:length(deltas), function(j) {
  var(dat$Y + Delta.star.subcont.hats[,j] * (dat$phi.hat - s.hats[j]))
}, simplify = 0)


# save(dat, n, deltas, delta.width, M,
#      Delta.star.hats, q.hats,
#      value.hats, var.hats,
#      Delta.star.sub.hats, r.hats,
#      value.sub.hats, var.sub.hats,
#      Delta.star.subcont.hats, s.hats,
#      value.subcont.hats, var.subcont.hats,
#      file = "icu-4folds.Rdata")
# load("icu-4folds.Rdata")


### plotting results

par(mar = c(4,4.6,1,1))
sub.delts <- seq(1,length(deltas) - 1, length.out = 20)
## all covariates approach
reord <- order(value.hats)
plot(x = deltas, y = value.hats[reord], type = "l",
     xlab = expression("Probability of contact" ~ (delta)), 
     ylab = expression("Mean counterfactual"),
     ylim = c(0.74,0.86))
arrows(x0 = deltas[sub.delts], x1 = deltas[sub.delts],
       y0 = value.hats[reord][sub.delts] - 1.96 *
         sqrt(var.hats[reord][sub.delts] / n),
       y1 = value.hats[reord][sub.delts] + 1.96 *
         sqrt(var.hats[reord][sub.delts] / n),
       code = 3, angle = 90, length = 0.05,
       col = alpha("black", 0.35))
# segments(y0 = min(value.hats),y1 = max(value.hats),
#          x0=0,x1=1,lty="dashed", col='red')
# abline(h = mean(dat$Y), lty="dashed", col = "blue")
# abline(h = max(value.hats), lty="dashed", col = "blue")
## subset of covariates approach (both Delta and h restricted)
reord_sub <- order(value.sub.hats)
lines(x = deltas, y = value.sub.hats[reord_sub],
      col = 'red')
arrows(x0 = deltas[sub.delts], x1 = deltas[sub.delts],
       y0 = value.sub.hats[reord_sub][sub.delts] - 1.96 *
         sqrt(var.sub.hats[reord_sub][sub.delts] / n),
       y1 = value.sub.hats[reord_sub][sub.delts] + 1.96 *
         sqrt(var.sub.hats[reord_sub][sub.delts] / n),
       code = 3, angle = 90, length = 0.05,
       col = alpha('red', 0.35))
## subset of covariates approach (just Delta restricted)
reord_subcont <- order(value.subcont.hats)
lines(x = deltas, y = value.subcont.hats[reord_subcont],
      col = 'blue')
arrows(x0 = deltas[sub.delts], x1 = deltas[sub.delts],
       y0 = value.subcont.hats[reord_subcont][sub.delts] - 1.96 *
         sqrt(var.subcont.hats[reord_subcont][sub.delts] / n),
       y1 = value.subcont.hats[reord_subcont][sub.delts] + 1.96 *
         sqrt(var.subcont.hats[reord_subcont][sub.delts] / n),
       code = 3, angle = 90, length = 0.05,
       col = alpha('blue', 0.35))
legend("bottomright", lty = 1, col = c("black", "blue", "red"),
       bty = 'n',
       legend = c("Unrestricted treatment rule",
                  expression("Contact rule" ~ italic(W) ~ "dependent"),
                  expression("Contact rule and policy" ~ italic(W) ~ "dependent")))

# hist(dat$tau.hat)
# hist(dat$c.hat)
# plot(dat$c.hat, abs(dat$tau.hat))

### estimation and inference for AUC measures

## all covariates
AUPBC.hat <- sum(value.hats) * delta.width - 0.5 * value.hats[length(deltas)] - 0.5 * mean(dat$Y)
AUPBC.norm.hat <- 2 * (sum(value.hats) * delta.width - mean(dat$Y)) / 
  (value.hats[length(deltas)] - mean(dat$Y)) - 1

N <- sapply(1:n, function(i) {
  sum(Delta.star.hats[i,] * dat$phi.hat[i] - q.hats * 
        (Delta.star.hats[i,] - deltas)) * delta.width - 0.5 * dat$phi.hat[i]
},simplify = 0)
N.norm <- as.vector((2 * N - AUPBC.norm.hat * dat$phi.hat) / mean(dat$phi.hat))
AUPBC.hat + c(-1,1) * qnorm(0.975) * sqrt(var(N) / n)
AUPBC.norm.hat + c(-1,1) * qnorm(0.975) * sqrt(var(N.norm) / n)

## subset of covariates
AUPBC.sub.hat <- sum(value.subcont.hats) * delta.width - 
  0.5 * value.subcont.hats[length(deltas)] - 0.5 * mean(dat$Y)
AUPBC.norm.sub.hat <- 2 * (sum(value.subcont.hats) * delta.width -
                                 mean(dat$Y)) / 
  (value.subcont.hats[length(deltas)] - mean(dat$Y)) - 1

N.sub <- sapply(1:n, function(i) {
  sum(Delta.star.subcont.hats[i,] * dat$phi.hat[i] - s.hats * 
        (Delta.star.subcont.hats[i,] - deltas)) * delta.width - 
    0.5 * dat$phi.hat[i]
},simplify = 0)
N.sub.norm <- as.vector((2 * N - AUPBC.norm.sub.hat * dat$phi.hat) / mean(dat$phi.hat))
AUPBC.sub.hat + c(-1,1) * qnorm(0.975) * sqrt(var(N) / n)
AUPBC.norm.sub.hat + c(-1,1) * qnorm(0.975) * sqrt(var(N.norm) / n)

### DR-Learner for CPB estimation

dr.beta.mod <- gam(phi.hat ~ s(sofa_score), data = dat)
sofa.scores <- seq(0,14, by = 0.1)
cpb.hat <- predict(dr.beta.mod, 
                   newdata = data.frame(sofa_score = sofa.scores))
plot(sofa.scores, cpb.hat,type = "l",
     xlab = "SOFA Score", ylab = "Estimated CPB")
