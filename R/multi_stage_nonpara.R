#' Multi-Stage Nonparametric Procedure
#'
#' @param n number of sample size.
#' @param model a character string for a nonlinear model: \code{"logist"} or \code{"arctan"}.
#' @param num.interp number of pseudo-data points, which are unifromly generated within the observed data points.
#' @param newl a length of time points (log transformed ages) at which predictors are required for each individual longitudinal trajectory.
#' @param dist a character string for the distribution of within-subject error term in the longitudinal model. Default is \code{"normal"}.
#' @param k1 a parameter for the dimension of the basis functions in order to represent the monotone increasing and convex smooth function.
#' @param k2 a parameter for the dimension of the basis functions in order to represent the monotone increasing and concave smooth function.
#' @param eps.sd a true scale parameter of the within-subject error term in the longitudinal model.
#' @param mean.diff an initial tolerence for the convergence criterion. Default is 1.
#' @param tolerance a parameter for the tolerence for the convergence criterion. Default is 0.009.
#' @param itermax a parameter for the maximum number of iterations in the multi-stage nonparametric algorithm. Default is 20.
#' @param iter initial number of iteration in the multi-stage nonparametric algorithm. The default is 0.
#' @param dat a data frame of the generated data set.
#'
#' @return A list of the simulation results includes
#'        \item{true_z}{a n-length of the true individual inflection points vector.}
#'        \item{zstar0}{a n-length of the estimated initial individual inflection points vector.}
#'        \item{zstar}{a n-length of the estimated individual inflection points vector.}
#'        \item{indvidual.Ta}{a n-length of estimated individual inflection points vector after all trajectories were shifted by their inflection points.}
#'        \item{newlogS}{a newl-length of the new time points vector at which predictors are required for all subjects. }
#'        \item{tms.pred}{a newl-length of the predcited longitudinal trajectory corresponding \code{newlogS} from all subjects' longitudinal data. }
#'        \item{tms.se.pred}{a newl-length of the estimated standard errors for \code{tms.pred}.}
#'        \item{global.T}{estimated inflection point on the shifted time points for all subjects.}
#'        \item{est.sigma.epsi}{an estimated scale parameter of the within-subject error in the longitudinal model.}
#'        \item{newdlogS}{a newl x n matrix with new time points, where each column is the newl-length of individual time points vector at which predcitions are required for each subject.}
#'        \item{ind.predict.tms}{a newl x n matrix with predicted longitudinal trajectory, where each column is the newl-length of the predicted longitudinal trajectory at each column of \code{newdlogS} for each subject.}
#'        \item{ind.str.predict.tms}{a newl x n  matrix with estimated standard errors of \code{ind.predict.tms}, where each column is the newl-length of the estimated standard errors of each column of \code{ind.predict.tms} for each subject.}
#'        \item{ind.b.sd}{a n-length of the boostrap standard deviations vector of \code{zstar}, where each component is the bootstrap standard deviaion of the estimated inflection point for each subject. }
#'        \item{ind.est.cp}{A n x 2 95\% bootstrap confidence intervals matrix for \code{zstar}, i.e., each row has the lower bound and the upper bound of the 95\% bootstrap confidence intervals for each individual's estimated inflection point. }
#'        \item{est.beta0}{estimated intercept parameter of the log-normal model for the inflection point.}
#'        \item{str.beta0}{estimated standard error of \code{est.beta0}.}
#'        \item{est.beta}{a (p-1)-length of estimated coefficient vector of subject specific covariates in the log-normal model for the inflection point.}
#'        \item{str.beta}{a (p-1)-length of vector with the estimated standard errors of \code{est.beta}. }
#'        \item{est.sigma.ei}{an estimated scale parameter of the random error in the log-normal model for the inflection point. }
#'        \item{xx.time}{A (num.interp/2) x n time points matrix, where each column is a subset of the each column of \code{newdlogS}, i.e., the individual time points at which estimated numerical derivatives are required for each subject.}
#'        \item{ind.org.tms}{A (num.interp/2) x n true longitudinal trajectories matrix, where  each column is the true individual longitudinal trajectory at each column of \code{xx.time} for each subject.}
#'        \item{ind.est.tms}{A (num.interp/2) x n estimated longitudinal trajectories matrix, where  each column is the subset of \code{ind.predict.tms}, i.e., estimated individual longitudinal trajectory at each column of \code{xx.time} for each subject. }
#'        \item{trueD1}{A (num.interp/2) x n numerically estimated first derivatives matrix of \code{ind.org.tms}, where each column is the numerically estimated first derivative of each column of \code{ind.org.tms} at each column of \code{xx.time} for each subject. }
#'        \item{trueD2}{A (num.interp/2) x n numerically estimated second derivatives matrix of \code{ind.org.tms}, where each column is the numerically estimated second derivative of each column of \code{ind.org.tms} at each column of \code{xx.time} for each subject. }
#'        \item{predD1}{A (num.interp/2) x n numerically estimated first derivatives matrix of \code{ind.est.tms}, where each column is the numerically estimated first derivative of each column of \code{ind.est.tms} at each column of \code{xx.time} for each subject.}
#'        \item{predD2}{A (num.interp/2) x n numerically estimated second derivatives matrix of \code{ind.est.tms}, where each column is the numerically estimated second derivative of each column of \code{ind.est.tms} at each column of \code{xx.time} for each subject.}
#'        \item{iter}{number of interations until convergence. }
#'        \item{mean.diff}{estimated convergence criterion, which is less than a \code{tolerance} if the algorithm convergences.}
#'
#' @export
#'
#' @examples
#'
#' library(HDChangePoint)
#'
#' ## Specify parameters to generate true data
#' n=80;
#' model="logist";
#' p=2;
#' bb0=0.5;
#' bb=0.1;
#' x.sd=0.3;
#' v1=5;
#' v2=7;
#' dist="normal";
#' eps.sd=0.05;
#' u.sd=0.05;
#'
#' ## Specify parameters for the multi-stage nonparametric procedure
#' num.interp=45;
#' newl=45;
#' k1=20;
#' k2=20;
#' mean.diff=1;
#' tolerance=0.009
#' itermax=20;
#' iter=0;
#'
#'
#' set.seed(22)
#'
#' ## Data generation under the logistic model
#' outdat<-mydata(n=n, model=model, p=p, bb0=bb0, bb=bb, x.sd=x.sd, v1=v1, v2=v2, dist=dist, eps.sd=eps.sd, u.sd=u.sd)
#'
#'
#' ## Multi-stage nonparametric estimation
#' results<-sim.nonpara(n=n, model=model, num.interp=num.interp, newl=newl, dist=dist,
#'                      k1=k1, k2=k2, eps.sd=eps.sd, mean.diff=mean.diff, tolerance=tolerance,
#'                      itermax=itermax, iter=iter, dat=outdat)
#'
#'

sim.nonpara<-function(n=80, model="arctan",  num.interp=25,  newl=25,  dist="normal", k1=10, k2=10,
                      eps.sd=0.05, mean.diff=1, tolerance=0.009, itermax=20, iter=0, dat=outdat){



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

  ## boostrap standard deviations for inflection points
  b.sd<-zstar0
  init.b.sd<-zstar0
  ind.b.sd<-zstar0


  ######################################################################
  ### Find initial inflection points: zstar0                         ###
  ######################################################################
  ## intial coverage probability
  init.cp.boot<-array(0, dim=c(n, 2),
                      dimnames=list(paste("id",1:n,sep=""),c("95% lower", "95% upper"))
  )


  for (id in  1:n){

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


  while((mean.diff>tolerance)&&(iter<itermax)){

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


    ## scale parameter for the within-subject errors
    est.sigma.epsi<-sqrt(outscam$sig2)

    # variance of scale parameter
    #res.var.eps<-2*(est.sigma.epsi)^{4}/(n-2);
    #est.str.sigma.epsi<-sqrt(res.var.eps);

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

  ## Constant term  (baseline coefficient)

  est.beta0<-summary(fit)$coefficients["(Intercept)","Estimate"]
  str.beta0<-summary(fit)$coefficients["(Intercept)","Std. Error"]


  ## Coefficients
  est.beta<-summary(fit)$coefficients["cov","Estimate"]
  str.beta<-summary(fit)$coefficients["cov","Std. Error"]   # same as  sqrt(diag(vcov(fit)))

  ## scale parameter in random error term of inflection point  model
  est.sigma.ei<-summary(fit)$sigma

  # variance of scale parameter
  #res.var<-2*(est.sigma.ei)^{4}/(n-2)
  #est.str.sigma.ei<-sqrt(res.var)



  return(list(true_z=true_z,   zstar0=zstar0,   zstar=zstar,  indvidual.Ta=indvidual.Ta,
              newlogS= newlogS, tms.pred=tms.pred, tms.se.pred=tms.se.pred, global.T=global.T,  ## global estimated
              est.sigma.epsi=est.sigma.epsi, #est.str.sigma.epsi=est.str.sigma.epsi,
              newdlogS=newdlogS,  ind.predict.tms=ind.predict.tms, ind.str.predict.tms=ind.str.predict.tms, ##  individual estimated
              #init.b.sd=init.b.sd,
              ind.b.sd=ind.b.sd, iter=iter,  mean.diff= mean.diff,    ## interation info
              est.beta0=est.beta0, str.beta0=str.beta0, est.beta=est.beta, str.beta=str.beta, est.sigma.ei=est.sigma.ei, #est.str.sigma.ei=est.str.sigma.ei,   # estimated beta
              trueD1=trueD1, trueD2=trueD2, predD1=predD1,predD2=predD2,xx.time=xx.time,ind.org.tms=ind.org.tms,ind.est.tms=ind.est.tms,ind.est.cp=ind.est.cp
  ))

}

