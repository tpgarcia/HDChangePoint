    ## This is examples to produce Tables and Figures in the paper
    ## U.Lee, R.J.Carroll, K.Marder, Y.Wang, T.P.Garcia. (2018+). Estimating Disease Onset from Change Points of Markers Measured with Error.

    ####################################################################################
    ##Example 1. Produce the multi-stage nonparametric simulation results in Table 1.  #
    ####################################################################################

    library(HDChangePoint)


    ## simulation results data when the visit numbers are 5, 6 and 7.
    data(nonpara_logist_visit1)

    ## Specify parameters in the summary function
    nsim=100;
    n=80;
    num.interp=45;
    newl=45;
    bb.true=0.1;
    b0.true=0.5;
    true.sigma.u=0.05;
    true.sigma.eps=0.05;


    ## produce the results and plots (c), (d) in Figure 2 will be automatrically saved in your working directory.
    simus.out1<-HDChangePoint::simus.results(nsim, n, res=nonpara_logist_visit1, model="logist", num.interp, newl, bb.true,
                                  b0.true, true.sigma.u, true.sigma.eps, file="logist_visit1")
    simus.out1$res.table



    ## simulation results data when the visit numbers are 8, 9 and 10.
    data(nonpara_logist_visit2)

    ## produce result
    simus.out2<-simus.results(nsim, n, res=nonpara_logist_visit2, model="logist", num.interp, newl, bb.true,
                                  b0.true, true.sigma.u, true.sigma.eps, file="logist_visit2")
    simus.out2$res.table



    ## load simulated data when the visit numbers are 18, 19 and 20.
    data(nonpara_logist_visit3)

    ## produce result
    simus.out3<-simus.results(nsim, n, res=nonpara_logist_visit3, model="logist", num.interp, newl, bb.true,
                                   b0.true, true.sigma.u, true.sigma.eps, file="logist_visit3")
    simus.out3$res.table

    ####################################################################################
    ##Example 2. Produce the multi-stage nonparametric simulation results in Table 2.  #
    ####################################################################################

    library(HDChangePoint)

    ## load simulated data
    data(nonpara_logist_visit1)

    ## Specify parameters in the summary function
    nsim=100;
    n=80;
    num.interp=45;
    newl=45;
    bb.true=0.1;
    b0.true=0.5;
    true.sigma.u=0.05;
    true.sigma.eps=0.05;


    ## load simulated data when the visit numbers are 5, 6 and 7.
    simus.out1<-simus.results(nsim, n, res=results, model="logist", num.interp, newl, bb.true,
                                  b0.true, true.sigma.u, true.sigma.eps, file="logist_visit1")
    simus.out1$avr.est.Z




    ## load simulated data when the visit numbers are 8, 9 and 10.
    data(nonpara_logist_visit2)

    ## produce result
    simus.out2<-simus.results(nsim, n, res=results, model="logist", num.interp, newl, bb.true,
                                 b0.true, true.sigma.u, true.sigma.eps, file="logist_visit2")
    simus.out2$avr.est.Z




    ## generate simulation results data when the visit numbers are 18, 19 and 20.
    data(nonpara_logist_visit3)

    ## produce result
    simus.out3<-simus.results(nsim, n, res=results, model="logist", num.interp, newl, bb.true,
                                   b0.true, true.sigma.u, true.sigma.eps, file="logist_visit3")
    simus.out3$res.table


  #############################################################################################
  ##Example 3. Produce Figure 1. Note that we cannot produce the real data,                  ##
  ##                             but we simulated a data set                                 ##
  ##                             to show that we can replicate the results from our methods. ##
  #############################################################################################

  library(HDChangePoint)

   ## Load a pseudo-HD data set constructed by using the simulated data.
   data("PSEUDO_PREDICT_HD")
   head(PSEUDO_PREDICT_HD)
   simu.dat<-PSEUDO_PREDICT_HD

   ## Specify the parameters to obtained the analysis results from the simulated dataset.
   n=80;
   m=45;
   num.interp=45;
   newl=45;
   mean.diff=1;
   tolerance=0.01;
   itermax=20;
   iter=0;

   ## produce two tables to reproduce Table 4 in the paper and Figure 1 (automatically saved in your working directory).
   ## we used two covariates: CAG repeats and gender, which are associated with inflection points.
   simu.analysis.results<-hd.study(simu.data, m, num.interp, n, newl,
                                  mean.diff, tolerance, itermax, iter)

   ## produce results
   simu.analysis.results$nonpara_summary_table4
   simu.analysis.results$para_summary_table4

   ## check your working directory for Figure 1.
