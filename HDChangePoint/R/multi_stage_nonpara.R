################################################################
# main simulation for the multi-stage nonparametric approach  ##
################################################################

#' Title
#'
#' @param n
#' @param model
#' @param num.interp
#' @param newl
#' @param dist
#' @param k1
#' @param k2
#' @param eps.sd
#' @param mean.diff
#' @param eps
#' @param itermax
#' @param iter
#' @param dat
#'
#' @return main function to generate simulation results for the multi-stage nonparametric approach including the estimated individual longitudinal
#'         trajectories and their standard errors, estimated inflection points, estimatimated longitudinal model parameters and their standard errors.
#' @export
#'
#' @examples  #library(nlme); library(scam);  library(ShapeChange);
#'
#'             library(registry); library(pkgmaker); library(rngtools); library(iterators);
#'             library(foreach); library(parallel);
#'
#'             libaray(doParallel); nworkers<-detectCores()
#'             cl<-makePSOCKcluster(nworkers);
#'             library(doMC); registerDoMC(nworkers)
#'             library(doRNG)
#'
#'             set.seed(22)
#'             outdat<-mydata(n=80, model="logist", p=2, bb0=0.5, bb=0.1, x.sd=0.3, v1=5, v2=7, dist="normal", eps.sd=0.05, u.sd=0.05)
#'
#'             #parallel computing
#'             results<-foreach(i=1:nsim)%dorng%{
#'
#'                                 sim.nonpara(n=80, model="logist", num.interp=45, newl=45, dist="normal, k1=20, k2=20, eps.sd=0.05,
#'                                 mean.diff=1, esp=0.009, itermax=20, iter=0, dat=outdat)
#'                                 }
#'             stopCluster(cl)
#'

sim.nonpara<-function(n=80, model="arctan",  num.interp=25,  newl=25,  dist="normal", k1=10, k2=10,
                      eps.sd=0.05, mean.diff=1, eps=0.009, itermax=20, iter=0, dat=outdat){



  #########################################################
  ## initialization for pseudo data by the interpolation ##
  #########################################################

  ## log transformed ages, longitudinal resposne
  xx.pseudo<-array(0, dim=c(num.interp,n),
                   dimnames=list(paste("x",1:num.interp,sep=""),paste("id",1:n,sep="")
                   ))

  omega.pseudo<-xx.pseudo
  yy.pseudo<-xx.pseudo
  epstar2<-xx.pseudo


  ## gender
  pseudos<-array(0, dim=c(num.interp, n),
                 dimnames=list(paste("x",1:num.interp,sep=""),paste("id",1:n,sep="")
                 ))


  ###############################################
  ## create pseudo data by linear interpolation #
  ###############################################

  for(id in 1:n){

    #define subdata per ID
    #subdat<-subset(dat, subject==id) #original.dat, subjID==id)
    subdat<-dat[dat$subj.id==id, ]
    ## peudo gender
    id.gend<-as.numeric(subdat[1, "sex"])
    pseudos[,id]<-rep(id.gend, num.interp)

    ## peudo visit times and corresponding motor scores
    approx.points<-approx(subdat[,"xx"], subdat[,"omega"], n=num.interp)

    xx.pseudo[,id]<-approx.points$x
    omega.pseudo[,id]<-approx.points$y

    if(dist=="normal"){
      l.norm2<-length(xx.pseudo[, id])
      epstar2[, id]<-rnorm(l.norm2, 0, eps.sd)
    }

    yy.pseudo[, id]<-omega.pseudo[,id]+epstar2[,id]
  }


  pseudox<-xx.pseudo;
  pseudoy<-yy.pseudo;
  pseudow<-omega.pseudo;
  #plot(xx.pseudo[,5], omega.pseudo[,5], type='l')
  #######################
  ### initializatioin ###
  #######################

  ## num.interp length of log transformed age: xxstar
  xxstar<-array(0, dim=c(num.interp, n),
                dimnames=list(paste("x",1:num.interp,sep=""),paste("id",1:n,sep="")
                ))

  ## nonlinear function omega estaimted at xxstar, where all subjects have the same inflection points at 0
  omg.star<-xxstar

  ## num.interp length of total motor score: yystar
  yystar<-xxstar

  ## estimated nonlinear function omega: when inflection points at intial inflection point, zstar0
  #est.omega.interp<-xxstar
  eps.star<-xxstar
  pseduo.eps<-xxstar

  ###############################################################################
  ### individual data to estimate inflection points based on the fitted curve ###
  ###############################################################################

  ## newd is newdata of log transformed ages at which predictions are required with newl length
  ##                                    for each individual longitudinal trajectories

  newdlogS<-array(0, dim=c(newl, n),
                  dimnames=list(paste("x",1:newl, sep=""),paste("id",1:n,sep="")
                  ))

  ind.predict.tms<-newdlogS
  ind.str.predict.tms<-newdlogS


  ##############################################################################################
  ### initialization for inflection points zstar0
  ###                    individual.Ta, which is
  ###                     difference of previous inflection points and updated inflection points
  ##############################################################################################
  zstar0<-rep(100,n);
  indvidual.Ta<-zstar0;

  b.sd<-zstar0
  init.b.sd<-zstar0
  ind.b.sd<-zstar0
  ## To see convergence behaviors of inflection points in the while loop iterations
  ## Those values are not returned from this simulation function


  ######################################################################
  ### Find initial inflection points: zstar0                         ###
  ######################################################################
  ## intial coverage probability
  init.cp.boot<-array(0, dim=c(n, 2),
                      dimnames=list(paste("id",1:n,sep=""),c("95% lower", "95% upper"))
  )


  for (id in  1:n){
    #pseud.omega<-pseudoy[,id]
    pseud.logS<-pseudox[,id]
    pseud.yy<-pseudoy[,id]

    #####################
    #### ShapeChange ####
    #####################
    change.fit<-changept(pseud.yy~ip(pseud.logS, sh=1), fir=TRUE, ci=TRUE) #
    initial.T<-change.fit$chpt
    initial.z.boot<-change.fit$msbt

    init.cp.boot[id,]<-quantile(sort(initial.z.boot), prob=c(0.025, 0.975))



    ## initial inflection point zstar0 by choosing method of inflection points
    zstar0[id]<- initial.T
    init.b.sd[id]<-sd(initial.z.boot)

  }

  ## From the longforamat of data, choose true individual inflection point
  true_z<-unique(dat$true_z)
  ## average of absolute bias between true Z and initial Z
  mean.diff0<-mean(abs(true_z-zstar0))
  ## intial inflection points for while loop
  zstar<-zstar0;



  ##############################
  # initalization         ######
  # algorithm starts      ######
  ##############################
  interp.pt <-array(0, dim=c(num.interp, 4, n),
                    dimnames=list(paste("interp",1:num.interp,sep=""),
                                  c("AGE","TMS", "logDiff", "zstar"), paste("ID",1:n,sep="")))

  pseudo.data<-vector("list", n)


  ##############################
  # initalization         ######
  # derivatives           ######
  ##############################

  if ((num.interp %% 2)==0){
    org.time.length<-num.interp/2
  }else{
    org.time.length<-ceiling(num.interp/2)
  }

  # org.time.length<-num.interp
  trueD1<-array(0, dim=c(org.time.length,n),
                dimnames=list(paste("visit",1:org.time.length,sep=""),paste("ID",1:n,sep="")))
  trueD2<-trueD1
  predD1<-trueD1
  predD2<-trueD1
  xx.time<-trueD1
  ind.org.tms<-trueD1
  ind.est.tms<-trueD1

  ## intial coverage probability
  ind.cp.boot<-array(0, dim=c(n, 2),
                     dimnames=list(paste("id",1:n,sep=""),c("95% lower", "95% upper"))
  )

  ind.est.cp<-ind.cp.boot

  while((mean.diff>eps)&&(iter<itermax)){

    iter<-iter+1
    old.zstar<-zstar;


    ####################################################################
    # aligned inflection points:                                      ##
    # axis of logitudinal data is shifted by their inflection points  ##
    ####################################################################

    for (id in  1:n){

      ## pseudox: pseudo log transformed age for subject i at jth visit with length of num.interp
      ## zstar: initial inflection points
      ## xxstar: shifted log transformed age aligining by their inflection points for each subject
      ## everyone has same inflection points at 0

      xxstar[,id]<-pseudox[,id]-old.zstar[id]

      if (dist=="normal"){

        ## error term will be divide into two parts
        eps.star[,id]<-rnorm(num.interp, 0, eps.sd)

      }

      ## longitudinal model for total motor scores, yystar
      omg.star[,id]<-w(xxstar[,id],0, model)
      yystar[,id]<- omg.star[,id]+eps.star[,id]

      ### save interpolation points ####
      interp.pt[,,id]<-cbind(pseudox[,id], yystar[,id], xxstar[,id], old.zstar[id])


      #######################################################################################
      ##  get subinterp.data: data set produced by interpolation with indicator of
      ##  negative or positive parts for time points
      #######################################################################################


      ## subdata per each subject ##

      subinterp.data<-data.frame(interp.pt[,,id])


      subinterp.data$fac1<-ifelse(subinterp.data$logDiff<0, 1,0)
      subinterp.data$fac2<-1-subinterp.data$fac1
      subinterp.data$fac<-c(rep(1, sum(subinterp.data$fac1)),rep(2, sum(subinterp.data$fac2)))


      pseudo.data[[id]]<-subinterp.data


    }



    ##############################################
    ## data for monoton increasing  assumption  ##
    ##############################################
    pseud.gen<-do.call("rbind.data.frame", pseudo.data)


    #################################
    ## scam estimation for omega  ##
    #################################
    outscam<- scam(TMS~fac+s(logDiff,k=k1,bs="micx",m=2, by=fac1)+s(logDiff,k=k2,bs="micv",m=2, by=fac2),
                   family=gaussian(link="identity"), data=pseud.gen)

    #scam.check(outscam)

    ## scale parameter for the within-subject errors
    est.sigma.epsi<-sqrt(outscam$sig2)

    # variance of scale parameter
    res.var.eps<-2*(est.sigma.epsi)^{4}/(n-2);
    est.str.sigma.epsi<-sqrt(res.var.eps);

    ####################################################
    ## new points at which predictions are estimated  ##
    ####################################################

    newlogS<-seq(min(pseud.gen$logDiff),max(pseud.gen$logDiff), length.out=newl);
    newfac1<-ifelse(newlogS<0, 1,0)
    newfac2<-1-newfac1
    newfac<-c(rep(1, sum(newfac1)), rep(2, sum(newfac2)))

    ## combined log transfromed Age after shifting x axis of longitudinal process
    mynewlogS<-data.frame(logDiff=newlogS, fac=newfac, fac1=newfac1, fac2=newfac2);

    ######################
    ## predicted values ##
    ######################
    predscam<-predict(outscam, mynewlogS, se.fit=TRUE);

    #########################################################
    ## estimated values and their estimate standard errors ##
    #########################################################
    tms.pred<-predscam$fit;
    tms.se.pred<-predscam$se.fit;

    ############################################
    ### Find global inflection points by scam ##
    ### Three different inflection points     ##
    ############################################

    newlogS<-mynewlogS[,"logDiff"]
    s<-length(newlogS)
    tms.pred<-as.vector(tms.pred)


    #####################
    #### ShapeChange ####
    #####################

    glob.change.fit<-changept(tms.pred~ip(newlogS, sh=1),fir=TRUE,   ci = TRUE) #pen=TRUE, gcv=TRUE,
    glob.T<-glob.change.fit$chpt
    global.T<-glob.T;


    ##########################################
    ### Individual inflection points     ####
    ###########################################

    for (id in  1:n){

      xnew<-as.vector(xxstar[,id]);
      xnewfac1<-ifelse(xnew<0, 1,0)
      xnewfac2<-1-xnewfac1
      xnewfac<-c(rep(1, sum(xnewfac1)), rep(2, sum(xnewfac2)))
      ind.logS<-data.frame(subj.id=id,logDiff=xnew, fac=xnewfac, fac1=xnewfac1, fac2=xnewfac2)


      ind.predscam<-predict(outscam, ind.logS, se.fit=TRUE)

      ind.tms.pred<-ind.predscam$fit
      ind.tms.se.pred<-ind.predscam$se.fit


      ### length of newlogS to be predicted
      fac.ind.logS<-ind.logS$logDiff*ind.logS$fac1+ind.logS$logDiff*ind.logS$fac2
      q<-length(fac.ind.logS)

      ## return individual tms and their standard errors
      ind.predict.tms[,id]<-ind.tms.pred
      ind.str.predict.tms[,id]<-ind.tms.se.pred

      ## return the time points when time points are moved by inflection points
      newdlogS[,id]<-fac.ind.logS


      ######################################
      # methods to find inflection points ##
      ######################################


      #####################
      #### ShapeChange ####
      #####################
      ind.tms.pred<-as.vector(ind.tms.pred)
      fac.ind.logS<-as.vector(fac.ind.logS)

      ## estimate individual inflection pt
      ind.change.fit<-changept(ind.tms.pred~ip(fac.ind.logS, sh=1),fir=TRUE, ci = TRUE)
      Ta<-ind.change.fit$chpt
      indvidual.Ta[id]<-Ta

      ## bootstrap distributions of inflection pt
      boot.dist<-ind.change.fit$msbt
      b.sd[id]<-sd(boot.dist)


      ind.cp.boot[id,]<-quantile(sort(boot.dist), prob=c(0.025, 0.975))


      ######################################
      # Estimates numerical derivatives ##
      ######################################

      ## obtain true tms

      true.ind.omega<-w(fac.ind.logS, Ta, model)

      ## choose points where we want to estimate numerical derivatives for the true tms
      xx.time[,id]<-fac.ind.logS[seq(1,num.interp,2)]
      ind.org.tms[,id]<-true.ind.omega[seq(1,num.interp,2)]

      ## choose points where we want to esimate numerical derivatives for the estimated tms
      ind.est.tms[,id]<-as.vector(ind.tms.pred[seq(1,num.interp,2)])



      ## obtain numerival derivatives: first D1 and Second D2
      trueD1[,id]<-finite.differences(xx.time[,id], ind.org.tms[,id])
      predD1[,id]<-finite.differences(xx.time[,id], ind.est.tms[,id])


      trueD2[,id]<-sec.der.mid(xx.time[,id],ind.org.tms[,id])
      predD2[,id]<-sec.der.mid(xx.time[,id],ind.est.tms[,id])



    }



    ##############################################################
    ### Return to the original scale for inflection points   #####
    ##############################################################

    zstar<-old.zstar+indvidual.Ta;

    ###########################
    ### convergence criteiron #
    ###########################
    mean.diff<-mean(abs(old.zstar-zstar));

  }

  ## variability of update.zstar

  for(id in 1:n){
    ind.b.sd[id]<-init.b.sd[id]+b.sd[id]
    ind.est.cp[id,]<-init.cp.boot[id,]+ind.cp.boot[id,]
  }


  #######################################################
  # Print convergence criterion when it is simulated. ###
  #######################################################
  print(c(iter, mean.diff))


  ######################################################################
  #### Fit linear model for log transformed inflection point ###########
  ######################################################################

  ## initialization for covariates

  W.cov<-array(0, dim=c(n, 2),
               dimnames=list(paste("id",1:n, sep=""),c( "cons", "cov")
               ))

  for (j in 1:n){
    subj.cov<-unique(dat[,"cov.W"])
    W.cov[j, ]<-c(1, subj.cov[j])
  }



  ## standardized covariates ##
  str.cov<-(W.cov[,"cov"]-mean( W.cov[,"cov"]))/sd( W.cov[,"cov"])
  str.zstar<-(zstar-mean(zstar))/sd(zstar)
  ## construct a data frame

  fit.data<-data.frame(cbind(W.cov,zstar,zstar0,str.cov, str.zstar))


  ## fitting linear model for log transformed inflection point
  fit<-lm(zstar~cov, data=fit.data);
  #summary(fit)


  ## Constant term  (baseline coefficient)

  est.beta0<-summary(fit)$coefficients["(Intercept)","Estimate"]
  str.beta0<-summary(fit)$coefficients["(Intercept)","Std. Error"]


  ## Coefficients
  est.beta<-summary(fit)$coefficients["cov","Estimate"]
  str.beta<-summary(fit)$coefficients["cov","Std. Error"]   # same as  sqrt(diag(vcov(fit)))

  ## scale parameter in random error term of inflection point  model
  est.sigma.ei<-summary(fit)$sigma #sigma(fit) #summary(fit)$sigma

  # variance of scale parameter
  res.var<-2*(est.sigma.ei)^{4}/(n-2)
  est.str.sigma.ei<-sqrt(res.var)



  return(list(true_z=true_z,   zstar0=zstar0,   zstar=zstar,  indvidual.Ta=indvidual.Ta,
              newlogS= newlogS, tms.pred=tms.pred, tms.se.pred=tms.se.pred, global.T=global.T,  ## global estimated
              est.sigma.epsi=est.sigma.epsi, est.str.sigma.epsi=est.str.sigma.epsi,
              newdlogS=newdlogS,  ind.predict.tms=ind.predict.tms, ind.str.predict.tms=ind.str.predict.tms, ##  individual estimated
              init.b.sd=init.b.sd, ind.b.sd=ind.b.sd, iter=iter,  mean.diff= mean.diff,    ## interation info
              est.beta0=est.beta0, str.beta0=str.beta0, est.beta=est.beta, str.beta=str.beta, est.sigma.ei=est.sigma.ei, est.str.sigma.ei=est.str.sigma.ei,   # estimated beta
              trueD1=trueD1, trueD2=trueD2, predD1=predD1,predD2=predD2,xx.time=xx.time,ind.org.tms=ind.org.tms,ind.est.tms=ind.est.tms,ind.est.cp=ind.est.cp
  ))

}


################################################################
# summary results for the multi-stage nonparametric approach  ##
################################################################
#' Title
#'
#' @param nsim
#' @param n
#' @param res
#' @param model
#' @param num.interp
#' @param newl
#' @param bb.true
#' @param b0.true
#' @param true.sigma.u
#' @param true.sigma.eps
#' @param num.interval
#' @param file
#'
#' @return summarize simulation results from the multi-stage nonparametric procedure including the performance of longituidinal model parameters
#'         , inflection points (absolute biases, estiamted standard deviations, average of the estimated standard errors and 95% coverage probabilities)
#'         and plots for the average trends of estimated longitudinal trajectories and their sencond derivatives.
#' @export
#'
#' @examples  simus.out<-simus.results(nsim=100, n=80, res=results, model="logist", num.interp=45, newl=45, bb.true=0.1,
#'                                   b0.true=0.5, true.sigma.u=0.05, true.sigma.eps=0.05, num.interval=18, file="logist")
#'
#'
simus.results<-function(nsim=100, n=80, res=results, model="logist", num.interp=45, newl=45, bb.true=0.1, b0.true=0.5,
                        true.sigma.u=0.05, true.sigma.eps=0.05, num.interval=18, file="logist"){


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
  glob.str.sigma.eps<-glob.sigma.eps
  glob.T<-matrix(0, nrow=nsim, ncol=1, dimnames = list(paste("nsim",1:nsim,sep=""),c("glob.T")))
  glob.boot.mean<-glob.T
  glob.boot.sd<-glob.T
  glob.boot.ind.cp<-array(0, dim=c(n, 2, nsim), dimnames=list(paste("id",1:n,sep=""),paste(c("95%lower", "95%upper"))
                                                              ,paste("sim",1:nsim,sep="")))


  ## initial and estimated inflection points ##
  update.z<-true.z
  Ta<-true.z
  #boot.mean<-true.z
  boot.sd<-true.z
  ind.boot.sd<-true.z
  init.boot.sd<-true.z
  ## estimated beta
  hat.beta<-matrix(0, nrow=nsim, ncol=1, dimnames = list(paste("nsim",1:nsim,sep=""),c("est.beta")))
  hat.str.beta<-hat.beta
  hat.beta0<-hat.beta
  hat.str.beta0<-hat.beta
  hat.sigma.u<-hat.beta
  hat.str.sigma.u<-hat.beta



  ## convergence info
  it<-matrix(0, nrow=nsim, ncol=1, dimnames = list(paste("nsim",1:nsim,sep=""),c("iteration")))
  meandiff<-it


  ## individual estimated tms
  ind.newdlogS<-array(0, dim=c(nsim,newl,n),
                      dimnames=list(paste("nsim", 1:nsim, sep=""),
                                    paste("x",1:newl,sep=""), paste("id",1:n,sep="")))


  ind.pred.tms<-ind.newdlogS
  ind.se.pred.tms<-ind.newdlogS


  for (sim in 1:nsim){

    ##{{ true

    ## true inflection points
    true.z[sim,]<-res[[1]]$true_z

    ## estimated initial inflection points
    initial.z[sim,]<-res[[sim]]$zstar0


    ##{{ global

    ## global estimated inflection points and tms
    glob.logS[sim,]<-res[[1]]$newlogS
    glob.pred.tms[sim,]<-res[[sim]]$tms.pred
    glob.se.pred.tms[sim,]<-res[[sim]]$tms.se.pred

    glob.boot.ind.cp[,,sim]<-res[[sim]]$ind.est.cp

    ## inflection points of estiamted omega based on all subjects information
    glob.T[sim]<-res[[sim]]$global.T
    ##  estimated sigma sub epsilon in longitudinal model
    glob.sigma.eps[sim, ]<-res[[sim]]$est.sigma.epsi
    glob.str.sigma.eps[sim, ]<-res[[sim]]$est.str.sigma.epsi
    ##global}}

    ##{{ estimated individual inflection points

    ## new inflection points
    Ta[sim,]<-res[[sim]]$indvidual.Ta  ## inflection point around 0 (since shifted logAge axes)
    update.z[sim,]<-res[[sim]]$zstar   ## updated inflection points to the original scale

    ## bootstrap mean and bootstrap sd

    ind.boot.sd[sim,]<-res[[sim]]$ind.b.sd
    init.boot.sd[sim,]<-res[[sim]]$init.b.sd

    ## estimated individual inflection points}}


    ##{{ estimated beta
    hat.beta0[sim]<-res[[sim]]$est.beta0
    hat.str.beta0[sim]<-res[[sim]]$str.beta0


    hat.beta[sim]<-res[[sim]]$est.beta
    hat.str.beta[sim]<-res[[sim]]$str.beta

    hat.sigma.u[sim]<-res[[sim]]$est.sigma.ei
    hat.str.sigma.u[sim]<-res[[sim]]$est.str.sigma.ei

    ## estimated beta}}

  }


  ## estimated individual inflection points based on global estimated logitudinal response
  ## summarize by subject

  for (id in 1:n){

    for (sim in 1:nsim){


      ind.newdlogS[sim, ,id]<-res[[sim]]$newdlogS[,id]
      ind.pred.tms[sim, ,id]<-res[[sim]]$ind.predict.tms[,id]
      ind.se.pred.tms[sim, ,id]<-res[[sim]]$ind.str.predict.tms[,id]

    }

  }



  newl2<-ceiling(newl/2)
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

    for (sim in 1:nsim){


      true.xx[sim, ,id]<-res[[sim]]$xx.time[,id]

      true.ind.omega[sim, ,id]<-res[[sim]]$ind.org.tms[,id]
      pred.ind.tms[sim, ,id]<-res[[sim]]$ind.est.tms[,id]

      err.ind.omega[sim,,id]<-abs(  pred.ind.tms[sim, ,id]-true.ind.omega[sim, ,id])
      true.ind.omega.deriv1[sim, ,id]<-res[[sim]]$trueD1[,id]
      true.ind.omega.deriv2[sim, ,id]<-res[[sim]]$trueD2[,id]


      pred.ind.tms.deriv1[sim, ,id]<-res[[sim]]$predD1[,id]
      pred.ind.tms.deriv2[sim, ,id]<-res[[sim]]$predD2[,id]

      err.ind.omega.deriv1[sim,,id]<-abs(true.ind.omega.deriv1[sim, ,id]- pred.ind.tms.deriv1[sim, ,id])
      err.ind.omega.deriv2[sim,,id]<-abs(true.ind.omega.deriv2[sim, ,id]- pred.ind.tms.deriv2[sim, ,id])

    }

  }

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

  avr.true.d1<-apply(mean.true.deriv1,2,mean)
  avr.pred.d1<-apply(mean.pred.deriv1,2,mean)

  avr.true.d2<-apply(mean.true.deriv2,2,mean)
  avr.pred.d2<-apply(mean.pred.deriv2,2,mean)




  ##############################
  #### summary for results  ####
  ##############################

  ##{{ covergence

  average.diff<-mean(meandiff)
  average.iter<-mean(it)

  ## data frame ##

  converg.info<-data.frame(avr.diff=average.diff, num.iter=average.iter)
  converg.info<-round(converg.info, digits=3)
  rownames(converg.info)<-"info"

  ## covergence}}



  ##{{ individual  inflection points

  ##true inflection points ##
  true.infl.z<-true.z[1,]
  avr.true.infl<-mean(true.infl.z)

  ##Estimated initial inflection point ##
  mean.initial.z<-apply(initial.z,2, mean)

  ##############################
  # Individual Inflection points
  ##############################

  ## average individual inflection points over simulation runs
  mean.update.z<-mean(apply(update.z,2, mean))

  ## mean absolute errors across subjects
  mae.subj<-matrix(0, nrow=nsim, ncol=1, dimnames = list(paste("nsim",1:nsim,sep=""),c("mae.subj")))
  mae.Ta<-matrix(0, nrow=nsim, ncol=1, dimnames = list(paste("nsim",1:nsim,sep=""),c("mae.Ta")))

  for(sim in 1:nsim){

    mae.subj[sim]<-mean(abs(update.z[sim, ]-true.infl.z))
    mae.Ta[sim]<-mean(abs(Ta[sim, ]-rep(0, n)))
  }

  ## individual inflection points}}

  ##{{ esimation for Z over all subjects:

  ## (1) bias for Z: average mae.subj over simulation runs
  bias.updated.z<-mean(mae.subj)
  bias.Ta<-mean(mae.Ta)

  ## (2) empirical SD for Z: standard deviation over simulations runs for each subject and average out them across subjects
  sd.update.z<-apply(update.z,2, sd)
  sd.update.z.all<-mean(sd.update.z)

  ## (3) Boostrap sd for Z
  boot.sd.Z<-apply(ind.boot.sd, 2, mean)
  avr.boot.sd.Z<-mean(boot.sd.Z)


  ### Coverage probability for Z
  avr.boot.ind.cp<-array(0, dim=c(n, 2), dimnames=list(paste("id",1:n,sep=""),paste(c("95%lower", "95%upper"))))

  for (id in 1:n){
    avr.boot.ind.cp[id,]<-apply(glob.boot.ind.cp[id,,1:nsim],1,mean)
  }

  avr.boot.cp<-colMeans(avr.boot.ind.cp)
  ind.logT.cp<-cbind(  true.infl.z,  mean.initial.z, avr.boot.ind.cp)

  ## esimation for Z}}


  ### data frame ###
  avr.inflex<-data.frame(mean.true.T= avr.true.infl, mean.est.T=mean.update.z, abs.bias.T=bias.updated.z, sd.est.T= sd.update.z.all, str.est.T= avr.boot.sd.Z) #, cp=cp.z)
  rownames(avr.inflex)<-"mean.inf.pt"



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
  res.table<-round(res.table, digits=3)





  ########################
  ## global estimation ###
  ########################

  postscript(paste(file,"_nonpara_plots.eps",sep=""))
  #pdf(file="nonpara_plots.pdf")

  shift.true.xx<-avr.true.xx+mean.update.z
  max.xx<-round(max(shift.true.xx)+0.2,1)
  min.xx<-round(min(shift.true.xx)-0.2,1)
  xx.interval<-(max.xx-min.xx)/num.interval
  xlabs<-round(seq(min.xx, max.xx, by=xx.interval),2)


  max.yy<-max(avr.true.omega)+0.2
  min.yy<-min(avr.true.omega)-0.2
  yy.interval<-(max.yy-min.yy)/num.interval
  ylabs<-round(seq(min.yy, max.yy, by=yy.interval),2)


  ## Create a blank plot
  par(mfcol=c(1,2),  mar=c(6.1,4.5,2.1,0))


  plot(shift.true.xx, avr.true.omega, type="l", ylim=c(min.yy,max.yy), xlim=c(min.xx, max.xx),
       xlab="x", ylab=expression("Nonparam" (omega(x))),
       cex.lab=1.3, cex.axis=1.0, lwd=2, #main="Averaged Trajectory",
       col="black", lty=1, font=1.2, frame.plot=F, axes=F)
  lines(shift.true.xx, avr.pred.omega, col="red", lty=2, lwd=2)
  segments(mean.update.z, min.yy, mean.update.z, ( median(ylabs)+0.3), col="orange", lty=4, lwd=2)
  legend(min.xx, max.yy, c("True Logistic Model", "Estimates"), col=c("black", "red"),
         lty=c(1,2), bty="n",cex=1.1,lwd=rep(2,2))
  axis(1, at=c(xlabs), labels=c(xlabs))
  axis(2, at=c(ylabs), labels=c(ylabs))
  mtext('(c)', outer=F,side=1,line=4.5, cex=1.2)


  ####################################
  # Plot of Second Derivatives of Omega
  ####################################

  max.yy2<-max( avr.true.d2)+0.2
  min.yy2<-min( avr.true.d2)-0.2
  yy2.interval<- (max.yy2-min.yy2)/num.interval
  y2labs<-round(seq(min.yy2, max.yy2, by=yy2.interval), 2)


  #average of trajectories of all subjects
  plot(shift.true.xx,   avr.true.d2,
         type="l", ylim=c(min.yy2, max.yy2), xlim=c(min.xx, max.xx),
         xlab="x", ylab=expression("Nonparam"(partialdiff^{2}~omega/ partialdiff~x^2)),
         cex.lab=1.3, cex.axis=1.0, lwd=2,
         col="black", lty=1, font=1.3, frame.plot = F, axes = F)
  lines(shift.true.xx,  avr.pred.d2, col="red", lty=2, lwd=2)
  mtext('(d)', outer=F,side=1,line=4.5, cex=1.2)
  axis(1, at=c(xlabs), labels=c(xlabs))
  axis(2, at=c(y2labs), labels=c(y2labs))


  dev.off()

  grid.plot<-cbind(avr.true.xx, avr.true.omega, avr.pred.omega,
                   avr.true.d1, avr.pred.d1,avr.true.d2, avr.pred.d2)


  return(list(res.table=res.table, avr.inflex=avr.inflex,  avr.err.omg= avr.err.omg,
              avr.err.deriv1=avr.err.deriv1, avr.err.deriv2=avr.err.deriv2,
              avr.boot.cp= avr.boot.cp, grid.plot=grid.plot))


}

