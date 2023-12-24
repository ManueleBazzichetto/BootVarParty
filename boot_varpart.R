#Bootstrapping variance partitioning

library(MASS)
library(ggplot2)
library(ggpubr)

#the idea is to use bootstrap to compute a measure of uncertainty around estimated variance
#explained by predictors, which is usually derived through 'variance partitioning' (see Venn diagrams)

#any quantity computed from the data (i.e., the sample) estimates a true, fixed, unknown parameter (estimand) #through estimators (at least, this is the case in the frequentist world)
#however, being a function of the data, estimators' value changes across samples (estimators are random variables)
#therefore, it is usually the case that estimators are combined with measures of uncertainty (e.g. variance)
#which tells us how precisely the estimator estimates parameter(s)

#estimators with high variance will largely change from one sample to an another, which means they 
#have low precision (e.g., reflected by wide confidence intervals)

#here, my claim is that numbers reported in Venn diagrams to represent variance explained by individual predictors or groups of predictors are just an estimate of an unknown quantity, and, as such, 
#should be combined with a measure of variance which tells us how precisely we can estimate the 
#true proportion of explained variance

#using bootstrap I can get a (non-parametric) image of the sampling distribution of the estimator for the variance explained by individual or groups of predictors. From the estimated sampling distribution, I can get a measure of sampling variance (variance of the estimator)

#here, I simulate 9 scenarios of correlation between 2 predictors X sample size, and estimate through
#bootstrap the uncertainty associated with the (partitioned) variance explained by each predictor

#the correlation X sample size settings of the 9 scenarios are the following: 

#correlation among predictors (phi): .1, .5, .9
#sample size (n): 50, 500, 5000


#------------------Simulate response variable and predictors

#set seed for reproducibility
set.seed(647)

#Simulate response variable from fixed parameters
#given phi and n

#predictors are simulated from a multivariate standard normal
#using centered and standardized random variables makes easy to set the population correlation between predictors
#although the sample correlation changes a lot across samples for low n
example_cor_cov <- mvrnorm(50, mu = c(0, 0), Sigma = matrix(c(1, .9, .9, 1),
                                                            nrow = 2, byrow = T))

cor(example_cor_cov)

plot(example_cor_cov[, 1], example_cor_cov[, 2], xlab = "X1", ylab = "X2")

#the function below simulates response variables given fixed parameters and n
SimulResp <- function(phi, n, params, sigma) {
  require(MASS)
  #simulate correlated covariates - let's assume they are Z-normal so that var is 1 (and mean 0)
  #basically correlation is the mean of the product of the 2
  Xs <- mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, phi, phi, 1), nrow = 2))
  #center and scale predicors, and set alpha to 0 - mean response is close to 0
  Xs <- scale(Xs)
  Y <- as.numeric(Xs %*% c(params[c("beta1", "beta2")]) + rnorm(n, mean = 0, sd = sigma))
  dtf <- data.frame(Y = Y, X1 = Xs[, 1], X2 = Xs[, 2])
  return(dtf)
}

#check
cor(SimulResp(phi = .5, n = 1000, params = c("beta1" = 2, "beta2" = 3), sigma = 2))
summary(SimulResp(phi = .5, n = 1000, params = c("beta1" = 2, "beta2" = 3), sigma = 2))

#create an example dataset - phi set to .5
A_dtf <- SimulResp(phi = .5, n = 1000, params = c("beta1" = 2, "beta2" = 3), sigma = 2)

Mod_dtf <- glm(as.formula(A_dtf), data = A_dtf, family = "gaussian")

#see estimated partial effects
Mod_coef <- coef(Mod_dtf)

Mod_modmat <- model.matrix(Mod_dtf)

P_eff_x1 <- as.numeric(cbind(Mod_modmat[, c(1, 2)], mean(Mod_modmat[, 3])) %*% Mod_coef)
P_eff_x2 <- as.numeric(cbind(Mod_modmat[, 1], mean(Mod_modmat[, 2]), Mod_modmat[, 3]) %*% Mod_coef)

A_dtf$Fit_x1 <- P_eff_x1
A_dtf$Fit_x2 <- P_eff_x2

plot(Y ~ X1, data = A_dtf)
lines(Fit_x1 ~ X1, data = A_dtf, col = "red")

plot(Y ~ X2, data = A_dtf)
lines(Fit_x2 ~ X2, data = A_dtf, col = "red")

#------------------Vojta's functions to compute rsq (R-squared) and partitioned variance

# explained deviance of a model
rsq <- function(model) 1-model$deviance/model$null.deviance

#https://stackoverflow.com/questions/7666807/anova-test-fails-on-lme-fits-created-with-pasted-formula
dev.expl2 <- function(model, group, dtf) {
  force(dtf)
  alt_mod <- do.call(update, args = list(object = model,
                                         formula. = paste("~.-", paste(group, collapse = "-"))))
  res <- rsq(model) - rsq(alt_mod)
  return(res)
}


#v.2 with a 'safer' handling of update?
dev.expl2(model = Mod_dtf, group = "X2", dtf = A_dtf)

dev.expl2(model = glm(formula(Mod_dtf), data = A_dtf, family = "gaussian"), 
         group = c("X2"), dtf = A_dtf)

#note that the following doesn't work
dev.expl2(model = glm(as.formula(A_dtf), data = A_dtf, family = "gaussian"), 
          group = c("X2"), dtf = A_dtf) #0

#------------------------Fit model and extract variance explained

#the function below fits the model and extract variance partitions
GetExplDev <- function(dtf) {
  #force(dtf)
  mod_obj <- glm(Y ~ X1 + X2, family = gaussian(link = "identity"), data = dtf)
  #this below works
  #return(rsq(mod_obj) - rsq(update(mod_obj, . ~ . - X1)))
  dev_expl_X1 <- dev.expl2(model = mod_obj, group = "X1", dtf = dtf)
  dev_expl_X2 <- dev.expl2(model = mod_obj, group = "X2", dtf = dtf)
  dev_shared <- (rsq(model = mod_obj) - dev_expl_X1 - dev_expl_X2)
  res <- c("DevExpl_X1" = dev_expl_X1,
           "DevExpl_X2" = dev_expl_X2,
           "DevShared" = dev_shared)
  return(res)
  }


#check that it works as expected
GetExplDev(dtf = A_dtf)

#check that the partitions sum up to 'total' R-sq
sum(GetExplDev(dtf = A_dtf))

#yep
rsq(model = Mod_dtf) #fine

#debugonce(GetExplDev)
#(function(dtf) force(dtf))(dtf = 1)


#-----------------------------Run computations

#run simulations for the different combinations of correlation X sample size
#note that seeds are set before each simulation chunk to make them 'independent':
#1) we don't rely on a single seed number;
#2) we prevents small changes in a chunk potentially having an impcat on next chunk

#low correlation

N <- c(50, 500, 5000)

set.seed(5869)

Low_dfs <- lapply(N, SimulResp, phi = .1, params = c("beta1" = 2, "beta2" = 3), sigma = 2)

names(Low_dfs) <- paste0("N_", as.character(N))

lapply(Low_dfs, cor)

Low_cor <- do.call(rbind, lapply(Low_dfs, function(dtf) {
  #dtf <- SimulResp(phi = .1, n = n, params = c("beta1" = 2, "beta2" = 3), sigma = 2)
  res <- as.data.frame(do.call(rbind, replicate(3000, expr = {
    dtf_boot <- dtf[sample(nrow(dtf), replace = T), ]
    var_expl <- GetExplDev(dtf = dtf_boot)
    return(var_expl)
  }, simplify = F)))
  res$Smpl_sz <- as.character(nrow(dtf))
  return(res)
  }))

Low_cor_res <- do.call(rbind, lapply(N, function(n) {
  n <- as.character(n)
  sub_dtf <- as.matrix(Low_cor[Low_cor$Smpl_sz == n, c(1, 2, 3)])
  expl_var_med <- apply(sub_dtf, 2, median)
  expl_var_lqrt <- apply(sub_dtf, 2, quantile, probs = 0.25)
  expl_var_uqrt <- apply(sub_dtf, 2, quantile, probs = 0.75)
  res <- data.frame(Mdn = expl_var_med, 
                    Lq = expl_var_lqrt, Uq = expl_var_uqrt, Pred = c("X1", "X2", "Shared"), Smpl_sz = n)
  }))


Low_cor_res$Smpl_sz <- factor(Low_cor_res$Smpl_sz, levels = as.character(N))

ggplot(Low_cor_res, aes(x = Smpl_sz, y = Mdn)) +
  geom_errorbar(aes(ymin = Lq, ymax = Uq, col = Pred)) +
  geom_point(col = "black")



#med correlation

set.seed(5690)

Med_dfs <- lapply(N, SimulResp, phi = .5, params = c("beta1" = 2, "beta2" = 3), sigma = 2)

names(Med_dfs) <- paste0("N_", as.character(N))

lapply(Med_dfs, cor)

Med_cor <- do.call(rbind, lapply(Med_dfs, function(dtf) {
  #dtf <- SimulResp(phi = .5, n = n, params = c("beta1" = 2, "beta2" = 3), sigma = 2)
  res <- as.data.frame(do.call(rbind, replicate(3000, expr = {
    dtf_boot <- dtf[sample(nrow(dtf), replace = T), ]
    var_expl <- GetExplDev(dtf = dtf_boot)
    return(var_expl)
  }, simplify = F)))
  res$Smpl_sz <- as.character(nrow(dtf))
  return(res)
}))


Med_cor_res <- do.call(rbind, lapply(N, function(n) {
  n <- as.character(n)
  sub_dtf <- as.matrix(Med_cor[Med_cor$Smpl_sz == n, c(1, 2, 3)])
  expl_var_med <- apply(sub_dtf, 2, median)
  expl_var_lqrt <- apply(sub_dtf, 2, quantile, probs = 0.25)
  expl_var_uqrt <- apply(sub_dtf, 2, quantile, probs = 0.75)
  res <- data.frame(Mdn = expl_var_med, 
                    Lq = expl_var_lqrt, Uq = expl_var_uqrt, Pred = c("X1", "X2", "Shared"), Smpl_sz = n)
}))


Med_cor_res$Smpl_sz <- factor(Med_cor_res$Smpl_sz, levels = as.character(N))

ggplot(Med_cor_res, aes(x = Smpl_sz, y = Mdn)) +
  geom_errorbar(aes(ymin = Lq, ymax = Uq, col = Pred)) +
  geom_point(col = "black")



#high correlation

set.seed(5345)

Hgh_dfs <- lapply(N, SimulResp, phi = .9, params = c("beta1" = 2, "beta2" = 3), sigma = 2)

names(Hgh_dfs) <- paste0("N_", as.character(N))

lapply(Hgh_dfs, cor)

Hgh_cor <- do.call(rbind, lapply(Hgh_dfs, function(dtf) {
  #dtf <- SimulResp(phi = .9, n = n, params = c("beta1" = 2, "beta2" = 3), sigma = 2)
  res <- as.data.frame(do.call(rbind, replicate(3000, expr = {
    dtf_boot <- dtf[sample(nrow(dtf), replace = T), ]
    var_expl <- GetExplDev(dtf = dtf_boot)
    return(var_expl)
  }, simplify = F)))
  res$Smpl_sz <- as.character(nrow(dtf))
  return(res)
}))


Hgh_cor_res <- do.call(rbind, lapply(N, function(n) {
  n <- as.character(n)
  sub_dtf <- as.matrix(Hgh_cor[Hgh_cor$Smpl_sz == n, c(1, 2, 3)])
  expl_var_med <- apply(sub_dtf, 2, median)
  expl_var_lqrt <- apply(sub_dtf, 2, quantile, probs = 0.25)
  expl_var_uqrt <- apply(sub_dtf, 2, quantile, probs = 0.75)
  res <- data.frame(Mdn = expl_var_med, 
                    Lq = expl_var_lqrt, Uq = expl_var_uqrt, Pred = c("X1", "X2", "Shared"), Smpl_sz = n)
}))


Hgh_cor_res$Smpl_sz <- factor(Hgh_cor_res$Smpl_sz, levels = as.character(N))

ggplot(Hgh_cor_res, aes(x = Smpl_sz, y = Mdn)) +
  geom_errorbar(aes(ymin = Lq, ymax = Uq, col = Pred)) +
  geom_point(col = "black")

#unique plot

Boot_var_plot <- rbind(data.frame(Low_cor_res, Cor_lev = "low"),
                       data.frame(Med_cor_res, Cor_lev = "medium"),
                       data.frame(Hgh_cor_res, Cor_lev = "high"))

Boot_var_plot$Cor_lev <- factor(Boot_var_plot$Cor_lev, levels = c("low", "medium", "high"))

Boot_var_plot$Smpl_sz <- factor(Boot_var_plot$Smpl_sz, levels = c("50", "500", "5000"))

Boot_var_result <- ggplot(Boot_var_plot, aes(x = Smpl_sz, y = Mdn)) +
  geom_errorbar(aes(ymin = Lq, ymax = Uq, col = Pred), width = 0.2) +
  geom_point(aes(col = Pred)) +
  scale_color_manual(values = c("X1" = "blue", "X2" = "red", "Shared" = "lightgrey"), name = "Predictor") +
  facet_wrap(~ Cor_lev) +
  ylab("Median (1st, 3rd quartile) explained variance") + xlab("Sample size") +
  theme_pubclean() +
  theme(axis.text = element_text(size = 12), legend.text = element_text(size = 14),
        axis.title = element_text(size = 14), legend.title = element_text(size = 14),
        strip.text = element_text(size = 12), strip.background = element_blank(),
        legend.position = "right")

ggsave(plot = Boot_var_result, device = "jpeg", filename = "BootVar_prel_res.jpeg",
       width = 22, height = 18, units = "cm", dpi = 220, bg = "white")
