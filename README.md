HDChangePoint
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
--------

The R package `HDChangePoint` implements two estimation methods using the parametric nonlinear mixed effects model (NLME) and the multi-stage nonparametric approaches, which are described in the following manuscript:

**U.Lee, R.J.Carroll, K.Marder, Y.wang and T.P.Garcia. (2019+) "Estimating Disease Onset from Change Points of Markers Measured with Error"**

Those two methods estimate individual longitudinal trajectories in an "S"-shape and their respective inflection points.

Installation
------------

To install `HDChangePoint` from GitHub,

``` r
devtools::install_github("unkyunglee/HDChangePoint")
```

Usage
-----

``` r

library(HDChangePoint)

# Specify parameters for data generation
n=80;
model="logist";
p=2;
bb0=0.5;
bb=0.1;
x.sd=0.3;
v1=5;
v2=7;
dist="normal";
eps.sd=0.05;
u.sd=0.05;

set.seed(22)
# Generate a dataset under the logistic model
outdat<-mydata(n=n, model=model, p=p, bb0=bb0, bb=bb, x.sd=x.sd, 
               v1=v1, v2=v2, dist=dist, eps.sd=eps.sd, u.sd=u.sd)


# Specify other parameters for the multi-stage nonparametric procedure
num.interp=45;
newl=45;
k1=20;
k2=20;
time.length=20;
mean.diff=1;
tolerance=0.009
iter=0;

# Multi-stage nonparametric estimation
results<-sim.nonpara(n=n, model="logist", dist="normal", k1=k1, k2=k2, num.interp=num.interp, 
                     newl=newl, eps.sd=eps.sd, mean.diff=1, tolerance=tolerance,
                     itermax=50, iter=iter, time.length=time.length, dat=outdat)
```

Example 1.
----------

We give an example to reproduce the part of the simulation results in Table 1 in our manuscript. We consider the `nonpara_logist_visit1` data available from R package `HDChangePoint`. It contains simluation results from the multi-stage nonparametric estimation when the true data are generated from the logistic model with the total number of visit numbers *m* = 5, 6, or 7.

``` r
library(HDChangePoint)
data(nonpara_logist_visit1_nsim500)
```

We specify parameters and obtain the results as the part of the Table 1 and plots (c) and (d) in Figure 2. Figures will be automatrically saved as .eps files in your working directory.

``` r
# Specify parameters in the summary function
nsim=500;
n=80;
num.interp=45;
newl=45;
time.length=20;
bb.true=0.1;
b0.true=0.5;
true.sigma.u=0.05;
true.sigma.eps=0.05;

# Produce the results 
simus.out1<-simus.results(nsim=nsim, n=n, res=nonpara_logist_visit1_nsim500,
                          model="logist", num.interp=num.interp, newl=newl,
                          time.length=time.length,
                          bb.true=bb.true, b0.true=b0.true, 
                          true.sigma.u=true.sigma.u, true.sigma.eps=true.sigma.eps,
                          file="logist_visit1")




simus.out1$res.table
#>       Bias  ESD  ASE   CP
#> beta0 0.02 0.14 0.16 0.97
#> beta  0.00 0.00 0.00 0.98
simus.out1$relative.avr.est.Z
#>                      Rel.Bias Rel.ESD Rel.Boot.SD Boot.CI.95.lower
#> relative.mean.inf.pt     0.02    0.03        0.02             4.56
#>                      Boot.CI.95.upper
#> relative.mean.inf.pt             4.89
```

Example 2.
----------

We give an example to show reproducibility for data analysis study in our manuscript. We consider the `PSEUDO_PREDICT_HD` data available from R package `HDChangePoint`. This data is a simulated dataset using the data generation function `pseudo.HD()`, see `data-raw\PSEUDO_PREDICT_HD.R`.

``` r
library(HDChangePoint)
data(PSEUDO_PREDICT_HD)
```

The data consist of 80 subjects' information with 8 variables. Each subject made different number of visits *m* = 5, 6 or 7 at clinics. In the `PSEUDO_PREDICT_HD` data, we consider two subject specific covariates: CAG repeats and gender, which may be associated with inflection points.

``` r
head(PSEUDO_PREDICT_HD)
#>   SUBJID event TOTAL_MOTOR_SCORE TRUE_INFL_POINT       AGE      CAG gender
#> 1      1     1        0.09839752        2.750496  8.298965 38.92753      0
#> 2      1     2        0.19024203        2.750496 13.147519 38.92753      0
#> 3      1     3        0.54985245        2.750496 16.249701 38.92753      0
#> 4      1     4        0.55897855        2.750496 16.879371 38.92753      0
#> 5      1     5        0.65459256        2.750496 17.004619 38.92753      0
#> 6      1     6        0.88050864        2.750496 20.820894 38.92753      0
#>     logAGE
#> 1 2.116131
#> 2 2.576233
#> 3 2.788074
#> 4 2.826092
#> 5 2.833485
#> 6 3.035957
?PSEUDO_PREDICT_HD # this gives you more information on the dataset
```

We fit our methods to the `PSEUDO_PREDICT_HD` data. First, we specify the parameters and run the function `hd.study()`.

``` r

# Specify the parameters
simu.dat<-PSEUDO_PREDICT_HD
n=80;
m=42;
num.interp=42;
newl=42;
mean.diff=1;
tolerance=0.01;
itermax=20;
iter=0;

# produce two tables to reproduce similar results as in Table 4 and Figure 1 of our manuscript
simu.analysis.results<-hd.study(simu.data=simu.data, m=m, num.interp=num.interp, 
                                n=n, newl=newl, mean.diff=mean.diff, 
                                tolerance=tolerance, itermax=itermax, iter=iter,     
                                boot.ci=TRUE)
```

We obtain the following results for the multi-stage nonparametric estimates and the parametric NLME estimates, similar to Table 4 in our manuscript. Figure will be automatically saved in your working directory, which is similar to Figure 1 in our manuscript. You can also find the result figure in the following folder: man \\ examples\_figures

``` r
# produce results
simu.analysis.results$nonpara_summary_table4 # multi-stage nonparametric estimates
#>             Estimate Std. Error t value Pr(>|t|)
#> beta0          3.021      0.027 113.590     0.00
#> beta_CAG       0.258      0.023  11.248     0.00
#> beta_gender    0.012      0.050   0.241     0.81
#> sigma_u        0.199         NA      NA       NA
#> sigma_eps      0.173         NA      NA       NA
simu.analysis.results$para_summary_table4    # parametric NLME estimates
#>             Value Std.Error  DF t-value p-value
#> theta1      5.919     0.159 382  37.136   0.000
#> theta2      0.977     0.014 382  68.144   0.000
#> beta0       2.996     0.007 382 454.626   0.000
#> beta_CAG    0.245     0.004 382  64.804   0.000
#> beta_gender 0.001     0.008 382   0.137   0.891
#> sigma_u     0.000        NA  NA      NA      NA
#> sigma_eps   0.057        NA  NA      NA      NA
```
