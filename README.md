# BootVarParty 🥳

This repo includes simulations (in R) to compute bootstrapped intervals of partitioned variance!

The idea is that, as any other quantity (e.g., parameters) estimated from the data, 'partitioned' variance, which is the variance explained by individual or groups of predictors in a linear modelling setting, should be presented along with intervals measuring the uncertainty around this quantity. Such intervals give us a measure of precision with which we estimate the true, yet unknown, portion of variance explained by predictors. 

Using bootstrap, we can get an estimate of the sampling distribution of 'partitioned' variance. The advantage of using bootstrap is that is a non-parametric technique, so we don't have to come out with a (parametric) model for the estimator of 'partitioned' variance.

Simulations presented here demonstrate that the estimated amount of variance explained by predictors (or shared by them) largely changes across samples in case of low degrees of freedom. So, in case of low sample size or overly complex (high number of parameters) models, our guess on the amount of variance explained by different predictors could be pretty unstable.

<br>

Preliminary results showing how the uncertainty around partitioned variance decreses with increasing sample size:

<br>

![Preliminary results](https://github.com/ManueleBazzichetto/BootVarParty/blob/main/BootVar_prel_res.jpeg)
