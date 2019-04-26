Untitled
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
--------

The R package `HDChangePoint` implements two estimation methods using the parametric nonlinear mixed effects model (NLME) and the multi-stage nonparametric approaches, which are described in the following manuscript:

**U.Lee, R.J.Carroll, K.Marder, Y.wang and T.P.Garcia. (2018+) "Estimating Disease Onset from Change Points of Markers Measured with Error"**

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


# Specify parameters for the multi-stage nonparametric procedure
num.interp=45;
newl=45;
k1=20;
k2=20;
mean.diff=1;
tolerance=0.009
itermax=20;
iter=0;

# Multi-stage nonparametric estimation
results<-sim.nonpara(n=n, model="logist", num.interp=num.interp, newl=newl, dist="normal", 
                     k1=k1, k2=k2, eps.sd=eps.sd, 
                     mean.diff=mean.diff, tolerance=tolerance, itermax=itermax, iter=iter, 
                     dat=outdat)
```

Example 1.
----------

We give an example to reproduce the part of the simulation results in Table 1 in our manuscript. We consider the `nonpara_logist_visit1` data available from R package `HDChangePoint`. It contains simluation results from the multi-stage nonparametric estimation when the true data are generated from the logistic model with the total number of visit numbers *m* = 5, 6, or 7.

``` r
library(HDChangePoint)
data(nonpara_logist_visit1)
```

We specify parameters and obtain the results as the part of the Table 1 and plots (c) and (d) in Figure 2. Figures will be automatrically saved as .eps files in your working directory.

``` r
# Specify parameters in the summary function
nsim=100;
n=80;
num.interp=45;
newl=45;
bb.true=0.1;
b0.true=0.5;
true.sigma.u=0.05;
true.sigma.eps=0.05;

# Produce the results 
simus.out1<-simus.results(nsim=nsim, n=n, res=nonpara_logist_visit1, model="logist", 
                          num.interp=num.interp, newl=newl, 
                          bb.true=bb.true, b0.true=b0.true, 
                          true.sigma.u=true.sigma.u, true.sigma.eps=true.sigma.eps,
                          file="logist_visit1")
simus.out1$res.table
#>        Bias   ESD   ASE   CP
#> beta0 0.006 0.115 0.153 0.99
#> beta  0.000 0.003 0.004 0.99
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
#>   SUBJID event TOTAL_MOTOR_SCORE TRUE_INFL_POINT      AGE      CAG gender
#> 1      1     1         0.1496728        3.095996 14.00652 43.33975      0
#> 2      1     2         0.1183109        3.095996 16.14418 43.33975      0
#> 3      1     3         0.1441487        3.095996 16.84745 43.33975      0
#> 4      1     4         0.6248139        3.095996 24.58382 43.33975      0
#> 5      1     5         0.7144454        3.095996 26.90517 43.33975      0
#> 6      1     6         0.9317825        3.095996 31.83840 43.33975      0
#>     logAGE
#> 1 2.639523
#> 2 2.781560
#> 3 2.824199
#> 4 3.202088
#> 5 3.292318
#> 6 3.460673
?PSEUDO_PREDICT_HD # this gives you more information on the dataset
#> starting httpd help server ... done
```

We fit our methods to the `PSEUDO_PREDICT_HD` data. First, we specify the parameters and run the function `hd.study()`.

``` r

# Specify the parameters
simu.dat<-PSEUDO_PREDICT_HD
n=80;
m=45;
num.interp=45;
newl=45;
mean.diff=1;
tolerance=0.01;
itermax=20;
iter=0;

# produce two tables to reproduce similar results as in Table 4 and Figure 1 of our manuscript
simu.analysis.results<-hd.study(simu.data=simu.data, m=m, num.interp=num.interp, 
                                n=n, newl=newl, mean.diff=mean.diff, 
                                tolerance=tolerance, itermax=itermax, iter=iter)
```

We obtain the following results for the multi-stage nonparametric estimates and the parametric NLME estimates, similar to Table 4 and Figure 1 in our manuscript. Figure will be automatically saved in your working directory.

``` r
# produce results
simu.analysis.results$nonpara_summary_table4 # multi-stage nonparametric estimates
#>             Estimate Std. Error t value Pr(>|t|)
#> beta0          2.967      0.018 162.094    0.000
#> beta_CAG       0.246      0.016  15.588    0.000
#> beta_gender    0.012      0.035   0.351    0.727
#> sigma_u        0.139         NA      NA       NA
#> sigma_eps      0.133         NA      NA       NA
simu.analysis.results$para_summary_table4    # parametric NLME estimates
#>              Value Std.Error  DF t-value p-value
#> theta1       5.360     0.158 380  33.968   0.000
#> theta2       1.043     0.017 380  60.052   0.000
#> beta0        3.019     0.007 380 409.613   0.000
#> beta_CAG     0.244     0.004 380  58.877   0.000
#> beta_gender -0.019     0.009 380  -2.001   0.046
#> sigma_u      0.000        NA  NA      NA      NA
#> sigma_eps    0.053        NA  NA      NA      NA
```
