#' Summary results of the multi-stage nonparametric procedure
#'
#' @param nsim number of simulation runs.
#' @param n  number of sample size.
#' @param res a nsim-length of list from the simulated results.
#' @param model a character string for a nonlinear model: \code{"logist"} or \code{"arctan"}.
#' @param num.interp  number of pseudo-data generated in the estimation procedure.
#' @param newl length of time points (log transformed ages) at which predictors are required for each individual longitudinal trajectory.
#' @param time.length length of time points for each individual longitudinal trajectory graph to be plotted, should be less than \code{newl}.
#' @param bb.true a parameter for the (p-1)-length of true coefficient vector corresponding to subject specific covariates.
#' @param b0.true a parameter for the true intercept of the log-normal model for the inflection point.
#' @param true.sigma.u a true scale parameter for the random error term in the inflection point model.
#' @param true.sigma.eps a true scale parameter for the within-subject error term in the longitudinal model.
#' @param file a character string of file name to generate plots with .eps file extention.
#'
#''
#' @return A list of the summarized simulation results from the multi-stage nonparametric procedure
#'
#'        \item{res.table}{a data frame of the absolute biases, estimated standard deviations, average of the estimated standard errors and 95\% coverage probabilities for the fixed effects (\code{bb.true}, \code{bb0.true}).}
#'        \item{relative.avr.est.Z}{a data frame of the relative average performance of inflection points for all subjects including the relative absolute biases, the relative empirical standard errors, the relative bootstrap standard deviations and 95\% bootstrap confidence intervals.}
#'
#'
#' @export
#'
#' @examples
#'
#' library(HDChangePoint)
#' data(nonpara_logist_visit1_nsim500)
#'
#' ## Specify parameters to obtain the summary results when the true data were generated from the logistic model.
#'
#' nsim=500;
#' n=80;
#' num.interp=45;
#' newl=45;
#' bb.true=0.1;
#' b0.true=0.5;
#' true.sigma.u=0.05;
#' true.sigma.eps=0.05;
#' time.length=20;
#' file="logistic_visit1";
#'
#' simus.out<-simus.results(nsim=nsim, n=n, res=combi.res, model="logist", num.interp=num.interp, newl=newl, time.length=time.length, bb.true=bb.true,
#'                          b0.true=b0.true, true.sigma.u=true.sigma.u, true.sigma.eps=true.sigma.eps, file=file)
#'
#'

simus.results<-function(nsim=500, n=80, res=combi.res, model="logist", num.interp=45, newl=45, time.length=20, bb.true=0.1, b0.true=0.5,
                        true.sigma.u=0.05, true.sigma.eps=0.05,  file="logist"){



  #################################
  # initialization for true values
  #################################
  true.z<-matrix(0, nrow=nsim, ncol=n, dimnames = list(paste("sim",1:nsim,sep=""),paste("id",1:n,sep="")))
  initial.z<-true.z



  #predicted tms in the population level
  glob.logS<-matrix(0, nrow=nsim, ncol=newl, dimnames = list(paste("nsim",1:nsim,sep=""),paste("x",1:newl,sep="")))
  glob.pred.tms<-glob.logS
  glob.se.pred.tms<-glob.logS
  glob.sigma.eps<-matrix(0, nrow=nsim, ncol=1, dimnames=list(paste("nsim",1:nsim,sep=""),c("sigma_eps1")))
  #glob.str.sigma.eps<-glob.sigma.eps
  glob.T<-matrix(0, nrow=nsim, ncol=1, dimnames = list(paste("nsim",1:nsim,sep=""),c("glob.T")))
  glob.boot.mean<-glob.T
  glob.boot.sd<-glob.T
  glob.boot.ind.cp<-array(0, dim=c(n, 2, nsim), dimnames=list(paste("id",1:n,sep=""),paste(c("95%lower", "95%upper"))
                                                              ,paste("sim",1:nsim,sep="")))


  ## initial and estimated inflection points ##
  update.z<-true.z
  Ta<-true.z
  boot.mean<-true.z
  boot.sd<-true.z
  ind.boot.sd<-true.z
  init.boot.sd<-true.z

  ## estimated beta
  hat.beta<-matrix(0, nrow=nsim, ncol=1, dimnames = list(paste("nsim",1:nsim,sep=""),c("est.beta")))
  hat.str.beta<-hat.beta
  hat.beta0<-hat.beta
  hat.str.beta0<-hat.beta
  hat.sigma.u<-hat.beta
  #hat.str.sigma.u<-hat.beta



  ## convergence info
  #it<-matrix(0, nrow=nsim, ncol=1, dimnames = list(paste("nsim",1:nsim,sep=""),c("iteration")))
  #meandiff<-it


  for (sim in 1:nsim){

    ##{{ true

    ## true inflection points
    true.z[sim,]<-res$combined_ture_z[[sim]]
    ## estimated initial inflection points
    initial.z[sim,]<-res$combined_zstar0[[sim]]


    ##{{ global

    ## global estimated inflection points and tms
    glob.logS[sim,]<-res$combined_newlogS[[sim]]
    glob.pred.tms[sim,]<-res$combined_tms_pred[[sim]]
    glob.se.pred.tms[sim,]<-res$combined_tms_se_pred[[sim]]

    glob.boot.ind.cp[,,sim]<-res$combined.ind.est.cp[[sim]]

    ## inflection points of estiamted omega based on all subjects information
    glob.T[sim]<-res$combined.global.T[[sim]]
    ##  estimated sigma sub epsilon in longitudinal model
    glob.sigma.eps[sim, ]<-res$combined.est.sigma.epsi[[sim]]

    ##global}}


    ##{{ estimated individual inflection points

    ## new inflection points
    Ta[sim,]<-res$combined_indvidual.Ta[[sim]] ## inflection point around 0 (since shifted logAge axes)
    update.z[sim,]<-res$combined_zstar[[sim]]  ## updated inflection points to the original scale

    ## bootstrap mean and bootstrap sd

    ind.boot.sd[sim,]<-res$combined.ind.b.sd[[sim]]
    init.boot.sd[sim,]<-res$combined.init.b.sd[[sim]]

    ## estimated individual inflection points}}


    ##{{ estimated beta
    hat.beta0[sim]<-res$combined.est.beta0[[sim]]
    hat.str.beta0[sim]<-res$combined.str.beta0[[sim]]


    hat.beta[sim]<-res$combined.est.beta[[sim]]
    hat.str.beta[sim]<-res$combined.str.beta[[sim]]

    hat.sigma.u[sim]<-res$combined.est.sigma.ei[[sim]]


    ## estimated beta}}

  }


  newl2<-time.length
  true.xx<-array(0, dim=c(nsim,newl2,n),
                 dimnames=list(paste("nsim", 1:nsim, sep=""),
                               paste("x",1:newl2,sep=""), paste("id",1:n,sep="")))
  true.ind.omega<-true.xx        # ind.org.tms
  pred.ind.tms<- true.xx         #ind.est.tms
  err.ind.omega<-true.xx
  true.ind.omega.deriv1<- true.xx #trueD1
  pred.ind.tms.deriv1<- true.xx  #predD1
  true.ind.omega.deriv2<- true.xx #trueD2
  pred.ind.tms.deriv2<- true.xx  #predD2
  err.ind.omega.deriv1<-true.xx
  err.ind.omega.deriv2<-true.xx

  for (id in 1:n){

    for (sim in 1:500){


      true.xx[sim, ,id]<-res$combined.xx.time[[sim]][,id]

      true.ind.omega[sim, ,id]<-res$combined.ind.org.tms[[sim]][,id]
      pred.ind.tms[sim, ,id]<-res$combined.ind.est.tms[[sim]][,id]

      err.ind.omega[sim,,id]<-abs(pred.ind.tms[sim, ,id]-true.ind.omega[sim, ,id])
      true.ind.omega.deriv1[sim, ,id]<-res$combined.trueD1[[sim]][,id]
      true.ind.omega.deriv2[sim, ,id]<-res$combined.trueD2[[sim]][,id]


      pred.ind.tms.deriv1[sim, ,id]<-res$combined.predD1[[sim]][,id]
      pred.ind.tms.deriv2[sim, ,id]<-res$combined.predD2[[sim]][,id]

      err.ind.omega.deriv1[sim,,id]<-abs(true.ind.omega.deriv1[sim, ,id]- pred.ind.tms.deriv1[sim, ,id])
      err.ind.omega.deriv2[sim,,id]<-abs(true.ind.omega.deriv2[sim, ,id]- pred.ind.tms.deriv2[sim, ,id])

    }
  }

  #
  #     for (sim in 101:200){
  #
  #
  #       true.xx[sim, ,id]<-res$combined.xx.time[[sim]][,id]
  #
  #       true.ind.omega[sim, ,id]<-res$combined.ind.org.tms[[sim]][,id]
  #       pred.ind.tms[sim, ,id]<-res$combined.ind.est.tms[[sim]][,id]
  #
  #       err.ind.omega[sim,,id]<-abs(pred.ind.tms[sim, ,id]-true.ind.omega[sim, ,id])
  #       true.ind.omega.deriv1[sim, ,id]<-res$combined.trueD1[[sim]][,id]
  #       true.ind.omega.deriv2[sim, ,id]<-res$combined.trueD2[[sim]][,id]
  #
  #
  #       pred.ind.tms.deriv1[sim, ,id]<-res$combined.predD1[[sim]][,id]
  #       pred.ind.tms.deriv2[sim, ,id]<-res$combined.predD2[[sim]][,id]
  #
  #       err.ind.omega.deriv1[sim,,id]<-abs(true.ind.omega.deriv1[sim, ,id]- pred.ind.tms.deriv1[sim, ,id])
  #       err.ind.omega.deriv2[sim,,id]<-abs(true.ind.omega.deriv2[sim, ,id]- pred.ind.tms.deriv2[sim, ,id])
  #
  #     }

  # for (sim in 201:300){
  #
  #
  #   true.xx[sim, ,id]<-res$combined.xx.time[[sim]][,id]
  #
  #   true.ind.omega[sim, ,id]<-res$combined.ind.org.tms[[sim]][,id]
  #   pred.ind.tms[sim, ,id]<-res$combined.ind.est.tms[[sim]][,id]
  #
  #   err.ind.omega[sim,,id]<-abs(pred.ind.tms[sim, ,id]-true.ind.omega[sim, ,id])
  #   true.ind.omega.deriv1[sim, ,id]<-res$combined.trueD1[[sim]][,id]
  #   true.ind.omega.deriv2[sim, ,id]<-res$combined.trueD2[[sim]][,id]
  #
  #
  #   pred.ind.tms.deriv1[sim, ,id]<-res$combined.predD1[[sim]][,id]
  #   pred.ind.tms.deriv2[sim, ,id]<-res$combined.predD2[[sim]][,id]
  #
  #   err.ind.omega.deriv1[sim,,id]<-abs(true.ind.omega.deriv1[sim, ,id]- pred.ind.tms.deriv1[sim, ,id])
  #   err.ind.omega.deriv2[sim,,id]<-abs(true.ind.omega.deriv2[sim, ,id]- pred.ind.tms.deriv2[sim, ,id])
  #
  # }
  #
  #     for (sim in 301:400){
  #
  #
  #       true.xx[sim, ,id]<-res$combined.xx.time[[sim]][,id]
  #
  #       true.ind.omega[sim, ,id]<-res$combined.ind.org.tms[[sim]][,id]
  #       pred.ind.tms[sim, ,id]<-res$combined.ind.est.tms[[sim]][,id]
  #
  #       err.ind.omega[sim,,id]<-abs(pred.ind.tms[sim, ,id]-true.ind.omega[sim, ,id])
  #       true.ind.omega.deriv1[sim, ,id]<-res$combined.trueD1[[sim]][,id]
  #       true.ind.omega.deriv2[sim, ,id]<-res$combined.trueD2[[sim]][,id]
  #
  #
  #       pred.ind.tms.deriv1[sim, ,id]<-res$combined.predD1[[sim]][,id]
  #       pred.ind.tms.deriv2[sim, ,id]<-res$combined.predD2[[sim]][,id]
  #
  #       err.ind.omega.deriv1[sim,,id]<-abs(true.ind.omega.deriv1[sim, ,id]- pred.ind.tms.deriv1[sim, ,id])
  #       err.ind.omega.deriv2[sim,,id]<-abs(true.ind.omega.deriv2[sim, ,id]- pred.ind.tms.deriv2[sim, ,id])
  #
  #     }

  # for (sim in 401:500){
  #
  #
  #   true.xx[sim, ,id]<-res$combined.xx.time[[sim]][,id]
  #
  #   true.ind.omega[sim, ,id]<-res$combined.ind.org.tms[[sim]][,id]
  #   pred.ind.tms[sim, ,id]<-res$combined.ind.est.tms[[sim]][,id]
  #
  #   err.ind.omega[sim,,id]<-abs(pred.ind.tms[sim, ,id]-true.ind.omega[sim, ,id])
  #   true.ind.omega.deriv1[sim, ,id]<-res$combined.trueD1[[sim]][,id]
  #   true.ind.omega.deriv2[sim, ,id]<-res$combined.trueD2[[sim]][,id]
  #
  #
  #   pred.ind.tms.deriv1[sim, ,id]<-res$combined.predD1[[sim]][,id]
  #   pred.ind.tms.deriv2[sim, ,id]<-res$combined.predD2[[sim]][,id]
  #
  #   err.ind.omega.deriv1[sim,,id]<-abs(true.ind.omega.deriv1[sim, ,id]- pred.ind.tms.deriv1[sim, ,id])
  #   err.ind.omega.deriv2[sim,,id]<-abs(true.ind.omega.deriv2[sim, ,id]- pred.ind.tms.deriv2[sim, ,id])
  #
  # }

  #}

  ## mean absolute errors over simulation runs for each subject

  mean.abs.err.deriv1<-array(0, dim=c(n,newl2),
                             dimnames = list(paste("ID",1:n,sep=""), paste("visit",1:newl2,sep="")
                             ))
  mean.abs.err.deriv2<-mean.abs.err.deriv1
  mean.abs.err<-mean.abs.err.deriv1


  for (id in 1:n){
    mean.abs.err[id,]<-apply(err.ind.omega[,,id], 2, mean)
    mean.abs.err.deriv1[id,]<-apply(err.ind.omega.deriv1[,,id], 2, mean)
    mean.abs.err.deriv2[id,]<-apply(err.ind.omega.deriv2[,,id], 2, mean)

  }

  # sum
  avr.err.omg<-mean(apply(mean.abs.err,1,mean))
  avr.err.deriv1<-mean(apply(mean.abs.err.deriv1,1,mean))
  avr.err.deriv2<-mean(apply(mean.abs.err.deriv2,2,mean))




  ## plot

  mean.true.omega<-array(0, dim=c(n,newl2),
                         dimnames = list(paste("ID",1:n,sep=""), paste("visit",1:newl2,sep="")
                         ))
  mean.pred.omega<-mean.true.omega
  mean.true.deriv1<-mean.true.omega
  mean.pred.deriv1<-mean.true.omega
  mean.true.deriv2<-mean.true.omega
  mean.pred.deriv2<-mean.true.omega
  mean.true.xx<-mean.true.omega

  for (id in 1:n){
    mean.true.xx[id,]<-apply(true.xx[,,id],2,mean)

    mean.true.omega[id,]<-apply( true.ind.omega[,,id], 2, mean)
    mean.pred.omega[id,]<-apply(pred.ind.tms[,,id], 2, mean)

    mean.true.deriv1[id,]<-apply(true.ind.omega.deriv1[,,id], 2, mean)
    mean.pred.deriv1[id,]<-apply(pred.ind.tms.deriv1[,,id], 2, mean)
    mean.true.deriv2[id,]<-apply(true.ind.omega.deriv2[,,id], 2, mean)
    mean.pred.deriv2[id,]<-apply(pred.ind.tms.deriv2[,,id], 2, mean)
  }

  avr.true.xx<-apply(mean.true.xx,2,mean)

  avr.true.omega<-apply(mean.true.omega,2,mean)
  avr.pred.omega<-apply(mean.pred.omega,2,mean)


  avr.true.d1<-apply(mean.true.deriv1,2,mean)
  avr.pred.d1<-apply(mean.pred.deriv1,2,mean)

  avr.true.d2<-apply(mean.true.deriv2,2,mean)
  avr.pred.d2<-apply(mean.pred.deriv2,2,mean)




  ##############################
  #### summary for results  ####
  ##############################

  ##{{ covergence

  #average.diff<-mean(meandiff)
  #average.iter<-mean(it)

  ## data frame ##

  #converg.info<-data.frame(avr.diff=average.diff, num.iter=average.iter)
  #converg.info<-round(converg.info, digits=3)
  #rownames(converg.info)<-"info"

  ## covergence}}



  ##{{ individual  inflection points

  ##true inflection points ##
  true.infl.z<-true.z[1,]
  #colnames(true.infl.z)<-c("true.z", "true.z2", "true.z3", "true.z4", "true.z5")
  #avr.true.infl<-apply(true.infl.z,2, mean)
  avr.true.infl<-mean(true.infl.z)

  ##Estimated initial inflection point ##
  mean.initial.z<-apply(initial.z,2,mean)#cbind( apply(initial.z[1:100, ],2, mean), apply(initial.z[101:200, ],2, mean), apply(initial.z[201:300, ],2, mean)
  #        , apply(initial.z[301:400, ],2, mean), apply(initial.z[401:500, ],2, mean))
  #colnames(mean.initial.z)<-c("init.z1", "init.z2", "init.z3", "init.z4", "init.z5")
  avr.init.infl<-mean(mean.initial.z) #mean(apply(mean.initial.z,2, mean))



  ##############################
  # Individual Inflection points
  ##############################

  ## average individual inflection points over simulation runs
  mean.update.z<-apply(update.z,2, mean) #cbind(apply(update.z[1:100, ],2, mean), apply(update.z[101:200, ],2, mean), apply(update.z[201:300, ],2, mean)
  #     , apply(update.z[301:400, ],2, mean), apply(update.z[401:500, ],2, mean))
  avr.update.z<-mean(mean.update.z) #apply(combine.update.z,2, mean)
  #mean.update.z<-mean(combined.mean.update.z)

  ## mean absolute errors across subjects
  mae.subj<-matrix(0, nrow=nsim, ncol=1, dimnames = list(paste("nsim",1:nsim,sep=""),c("mae.subj")))
  mae.Ta<-matrix(0, nrow=nsim, ncol=1, dimnames = list(paste("nsim",1:nsim,sep=""),c("mae.Ta")))


  for(sim in 1:500){

    mae.subj[sim]<-mean(abs(update.z[sim, ]-true.infl.z)) #true.infl.z[,"true.z1"]))
    mae.Ta[sim]<-mean(abs(Ta[sim, ]-rep(0, n)))
  }


  # for(sim in 101:200){
  #
  #   mae.subj[sim]<-mean(abs(update.z[sim, ]-true.infl.z[,"true.z2"]))
  #   mae.Ta[sim]<-mean(abs(Ta[sim, ]-rep(0, n)))
  # }
  #
  #
  # for(sim in 201:300){
  #
  #   mae.subj[sim]<-mean(abs(update.z[sim, ]-true.infl.z[,"true.z3"]))
  #   mae.Ta[sim]<-mean(abs(Ta[sim, ]-rep(0, n)))
  # }
  #
  #
  # for(sim in 301:400){
  #
  #   mae.subj[sim]<-mean(abs(update.z[sim, ]-true.infl.z[,"true.z4"]))
  #   mae.Ta[sim]<-mean(abs(Ta[sim, ]-rep(0, n)))
  # }
  #
  #
  #
  # for(sim in 401:500){
  #
  #   mae.subj[sim]<-mean(abs(update.z[sim, ]-true.infl.z[,"true.z5"]))
  #   mae.Ta[sim]<-mean(abs(Ta[sim, ]-rep(0, n)))
  # }



  ## individual inflection points}}

  ##{{ esimation for Z over all subjects:

  ## (1) bias for Z: average mae.subj over simulation runs
  bias.updated.z<-mean(mae.subj)
  bias.Ta<-mean(mae.Ta)

  ## (2) empirical SD for Z: standard deviation over simulations runs for each subject and average out them across subjects
  sd.update.z<-apply(update.z, 2, sd) #cbind(apply(update.z[1:100, ],2, sd), apply(update.z[101:200, ],2, sd), apply(update.z[201:300, ],2, sd), apply(update.z[301:400, ],2, sd),
  #apply(update.z[401:500, ],2, sd))
  #colnames(sd.update.z)<-c("z.sd1", "z.sd2", "z.sd3", "z.sd4", "z.sd5")
  sd.update.z.all<-mean(sd.update.z) #apply(sd.update.z,2, mean))



  ## (3) Boostrap sd for Z
  boot.sd.Z<-apply(ind.boot.sd, 2, mean) #cbind(apply(ind.boot.sd[1:100, ], 2, mean), apply(ind.boot.sd[101:200, ], 2, mean), apply(ind.boot.sd[201:300, ], 2, mean),
  #apply(ind.boot.sd[301:400, ], 2, mean), apply(ind.boot.sd[401:500, ], 2, mean))
  #colnames(boot.sd.Z)<-c("z.boot.sd1", "z.boot.sd2", "z.boot.sd3", "z.boot.sd4", "z.boot.sd5")
  avr.boot.sd.Z<-mean(boot.sd.Z) #(apply(boot.sd.Z,2, mean))


  ### Coverage probability for Z
  avr.boot.ind.cp<-array(0, dim=c(n, 2), dimnames=list(paste("id",1:n,sep=""),paste(c("95%lower", "95%upper"))))
  #avr.boot.ind.cp2<-array(0, dim=c(n, 2), dimnames=list(paste("id",1:n,sep=""),paste(c("95%lower", "95%upper"))))
  #avr.boot.ind.cp3<-array(0, dim=c(n, 2), dimnames=list(paste("id",1:n,sep=""),paste(c("95%lower", "95%upper"))))
  #avr.boot.ind.cp4<-array(0, dim=c(n, 2), dimnames=list(paste("id",1:n,sep=""),paste(c("95%lower", "95%upper"))))
  #avr.boot.ind.cp5<-array(0, dim=c(n, 2), dimnames=list(paste("id",1:n,sep=""),paste(c("95%lower", "95%upper"))))

  for (id in 1:n){
    avr.boot.ind.cp[id,]<-apply(glob.boot.ind.cp[id,,],1,mean) #apply(glob.boot.ind.cp[id,,1:100],1,mean)
    #avr.boot.ind.cp2[id,]<-apply(glob.boot.ind.cp[id,,101:200],1,mean)
    #avr.boot.ind.cp3[id,]<-apply(glob.boot.ind.cp[id,,201:300],1,mean)
    #avr.boot.ind.cp4[id,]<-apply(glob.boot.ind.cp[id,,301:400],1,mean)
    #avr.boot.ind.cp5[id,]<-apply(glob.boot.ind.cp[id,,401:500],1,mean)
  }

  #avr.boot.cp<-rbind(apply(avr.boot.ind.cp1,2,mean), apply(avr.boot.ind.cp2,2,mean), apply(avr.boot.ind.cp3,2,mean),
  #                     apply(avr.boot.ind.cp4,2,mean), apply(avr.boot.ind.cp5,2,mean))
  avr.boot.cp<-colMeans(avr.boot.ind.cp)
  #ind.logT.cp<-cbind( true.infl.z,  mean.initial.z, avr.boot.cp)

  ## esimation for Z}}


  ### data frame ###

  avr.est.Z<-data.frame(abs.bias=bias.updated.z, sd.est= sd.update.z.all, Boot.sd.est= avr.boot.sd.Z,  Boot.cp.est=rbind(avr.boot.cp)) #, cp=cp.z)
  rownames( avr.est.Z)<-"mean.inf.pt"
  avr.est.Z<-round(avr.est.Z, 2)

  relative.avr.est.Z<-data.frame(Rel.Bias=bias.updated.z/avr.true.infl, Rel.ESD= sd.update.z.all/avr.true.infl,
                                 Rel.Boot.SD= avr.boot.sd.Z/avr.true.infl,  Boot.CI=rbind(avr.boot.cp)) #, cp=cp.z)
  rownames(relative.avr.est.Z)<-"relative.mean.inf.pt"
  relative.avr.est.Z<-round(relative.avr.est.Z, 2)

  ##{{ estimated beta

  mean.est.beta<-apply(hat.beta,2, mean)
  bias.beta<-abs(mean.est.beta-bb.true)
  sd.est.beta<-apply(hat.beta,2, sd)
  avr.str.est.beta<-apply(hat.str.beta,2, mean)

  ##Coverage probabilities
  qval<-qt(c(.025, .975), df=78)
  UBbeta<-hat.beta+qval[2]*hat.str.beta
  LBbeta<-hat.beta+qval[1]*hat.str.beta


  #beta.true
  ind1b<- do.call("rbind",lapply(UBbeta,function(x) 1*(x-(bb.true)>0)))
  ind2b<- do.call("rbind",lapply(LBbeta,function(x) 1*(x-(bb.true)<0)))
  indb<-ind1b*ind2b

  cpb<-mean(indb)


  #beta0.true
  mean.est.beta0<-apply(hat.beta0,2, mean)
  bias.beta0<-abs(mean.est.beta0-b0.true)
  sd.est.beta0<-apply(hat.beta0,2, sd)
  avr.str.est.beta0<-apply(hat.str.beta0,2, mean)

  ##Coverage probabilities
  UBbeta0<-hat.beta0+qval[2]*hat.str.beta0
  LBbeta0<-hat.beta0+qval[1]*hat.str.beta0

  ind1b0<- do.call("rbind",lapply(UBbeta0,function(x) 1*(x-(b0.true)>0)))
  ind2b0<- do.call("rbind",lapply(LBbeta0,function(x) 1*(x-(b0.true)<0)))
  indb0<-ind1b0*ind2b0

  cpb0<-mean(indb0)


  ## data.frame
  beta.est<-data.frame(Bias=bias.beta, ESD=sd.est.beta, ASE=avr.str.est.beta, CP=cpb )
  rownames(beta.est)<-"beta"

  beta0.est<-data.frame( Bias=bias.beta0, ESD=sd.est.beta0, ASE=avr.str.est.beta0, CP=cpb0 )
  rownames(beta0.est)<-"beta0"


  ## estimated beta}}


  ## combine all results
  res.table<-rbind(beta0.est, beta.est) #, sigma.u.est, sigma.eps.est)
  res.table<-round(res.table, digits=2)





  ########################
  ## global estimation ###
  ########################

  postscript(paste(file,"_nonpara_average_plot.eps",sep=""))
  #pdf(file="nonpara_plots.pdf")

  ## shift location by the average estimated inflection points (to compare to parametric NLME estimates later)
  shift.true.xx<-avr.true.xx+avr.update.z


  ## Create a blank plot
  par(mar=c(6.1,5.1,2.5,1.1))


  plot(shift.true.xx, avr.true.omega, type="l",  ylim=c(0,2), xlim=c(4.2,5.4),
       xlab="x", ylab=expression("Nonparam" (omega(x))),
       cex.lab=1.3, cex.axis=1.0, lwd=2, #main="Averaged Trajectory",
       col="black", lty=1, font=1.2)
  lines(shift.true.xx, avr.pred.omega, col="red", lty=2, lwd=2)
  segments(avr.update.z, 0, avr.update.z, 1.3, col="orange", lty=4, lwd=2)
  legend(4.2, 2, c("True Logistic Model", "Estimates"), col=c("black", "red"),
         lty=c(1,2), bty="n",cex=1.1,lwd=rep(2,2))
  mtext('(c)', outer=F,side=1,line=4.5, cex=1.2)

  dev.off()




  postscript(paste(file,"_nonpara_second_deriv.eps",sep=""))
  #pdf(file="nonpara_plots_second_deriv.pdf")

  par(mar=c(6.1,5.1,2.5,1.1))
  ####################################
  # Plot of Second Derivatives of Omega
  ####################################

  #average of trajectories of all subjects
  plot(shift.true.xx,   avr.true.d2,
       type="l", ylim=c(-6, 7), xlim=c(4.0, 5.4),
       xlab="x", ylab=expression("Nonparam"(partialdiff^{2}~omega/ partialdiff~x^2)),
       cex.lab=1.3, cex.axis=1.0, lwd=2,
       col="black", lty=1, font=1.3)
  legend(4.0, 7, c("True Logistic Model", "Estimates"), col=c("black", "red"),
         lty=c(1,2), bty="n",cex=1.1,lwd=rep(2,2))
  lines(shift.true.xx,  avr.pred.d2, col="red", lty=2, lwd=2)
  mtext('(d)', outer=F,side=1,line=4.5, cex=1.2)



  dev.off()


  return(list(res.table=res.table,  relative.avr.est.Z=relative.avr.est.Z #avr.true.infl=avr.true.infl, avr.est.Z= avr.est.Z,
              #shift.true.xx=shift.true.xx, avr.true.omega=avr.true.omega, avr.pred.omega=avr.pred.omega, mean.update.z=mean.update.z,
              #avr.true.d2=avr.true.d2, avr.pred.d2=avr.pred.d2
              #avr.bias.updated.z= avr.bias.updated.z, bias.update.initial=bias.update.initial, bias.initial=bias.initial
  ))


}


