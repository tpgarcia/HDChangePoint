#' Analysis example using the simulated HD dataset
#'
#' @param simu.dat a simulated data set, see \code{DATA.R} for details.
#' @param m number of time points at which predictors are required for the longitudinal responses in the parametric NLME procedure
#' @param num.interp number of pseudo-data points, which are unifromly generated within the observed data points.
#' @param n number of sample size.
#' @param newl a length of time points (log transformed ages) at which predictors are required for each individual longitudinal trajectory.
#' @param mean.diff an initial tolerence for the convergence criterion. Default is 1.
#' @param tolerance a parameter for the tolerence for the convergence criterion. Default is 0.009.
#' @param itermax a parameter for the maximum number of iterations in the multi-stage nonparametric algorithm. Default is 20.
#' @param iter initial number of iteration in the multi-stage nonparametric algorithm. The default is 0.
#' @param boot.ci logical value: TRUE if the 95% bootstrap confidence intervals across all subjects are calcuated and plotted on the graph. Default is TRUE.
#'
#' @return a list of
#'        \item{nonpara_summary_table4}{a 5 x 4 array of estimates, standard errors, t-values and p-values for the fixed effects of \code{beta0}, CAG repeates and gender in the inflection point model and scale parameters of the random errors in the inflection point and the longitudinal models.}
#'        \item{para_summary_table4}{a 7 x 4 array of estimates, standard errors, t-values and p-values for the fixed effects of \code{theta1}, \code{theta2}, \code{beta0}, CAG repeates and gender and scale parameters of the random errors in the inflection point and the longitudinal models. The parametric NLME assumes the logistic model.}
#'        , which are similar to Table 4 and return Figure 1 in the paper (See the \code{reference}).
#' @references U.Lee, R.J.Carroll, K.Marder, Y.Wang, T.P.Garcia. (2019+). Estimating Disease Onset from Change Points of Markers Measured with Error.
#' @export
#'
#' @examples
#'
#' library(HDChangePoint)
#'
#' ## a pseudo data set constructed by using the simulated data
#' data("PSEUDO_PREDICT_HD")
#' head(PSEUDO_PREDICT_HD)
#'
#' ## Specify the parameters to obtained the analysis results from the simulated dataset.
#' simu.dat<-PSEUDO_PREDICT_HD
#' n=80;
#' m=45;
#' num.interp=45;
#' newl=45;
#' mean.diff=1;
#' tolerance=0.01;
#' itermax=20;
#' iter=0;
#'
#' simu.analysis.results<-hd.study(simu.data=simu.data, m=m, num.interp=num.interp, n=n, newl=newl,
#'                                mean.diff=mean.diff, tolerance=tolerance, itermax=itermax, iter=iter, boot.ci=TRUE)
#'
#'
#'

hd.study<-function(simu.data=simu.data, m=30, num.interp=30, n=80,
                   newl=30,  mean.diff=1, tolerance=0.005, itermax=20, iter=0, boot.ci=TRUE){



  #simu.data<-data("PSEUDO_PREDICT_HD_DATA")
  ## 80 subjects
  obsID<-unique(simu.dat$SUBJID)

  #########################
  ## estimation starts ####
  #########################

  suborg<-vector("list", n)
  scamfit.data<-vector("list", n)
  merge.two.vectors<-vector("list", n)

  ## initialization for the baseline covariates
  base.cov <-matrix(100,nrow=n, ncol=4,
                    dimnames=list(paste("ID",1:n,sep=""),
                                  c("AGE", "CAG","gender","logAGE"))) #"LCau","LPut","RCau","RPut",



  ## initialization to generate pseudo visit numbers in the range of original scaled logAGE
  xxstar<-array(0, dim=c(num.interp, n),
                dimnames=list(paste("x",1:num.interp,sep=""),paste("id",1:n,sep="")
                ))

  pseudox<-xxstar
  pseudoy<-xxstar

  zstar0<-rep(100,n);
  indvidual.Ta<-zstar0;


  #####################################################################################
  ## newd is newdata of logAGE at which predictions are required with length newl    ##
  ##                                    for each individual trajectories of tms      ##
  #####################################################################################

  newdlogS<-array(0, dim=c(newl, n),
                  dimnames=list(paste("x",1:newl, sep=""),paste("id",1:n,sep="")
                  ))

  ind.predict.tms<-newdlogS
  ind.str.predict.tms<-newdlogS


  for(id in 1:length(obsID)){


    #######################
    #   Individual data   #
    #######################

    subind<-which(simu.dat$SUBJID==obsID[id])
    subdata<-simu.dat[subind,]


    ### save original data ##
    suborg[[id]]<-subdata


    ### baseline covariate ###
    base<-subdata[subdata$event==1,]
    base.cov[id,]<-c(base$AGE, base$CAG, base$gender, base$logAGE)



    ###########################################################
    # interpolation: creating more points between given data ##
    ###########################################################
    ## Time points are AGE

    subdata<-data.frame(subdata)

    interp.pts<-approx(subdata$logAGE, subdata$TOTAL_MOTOR_SCORE, method="linear", n=num.interp, rule=1)
    xxinterp<-interp.pts$x
    yyinterp<-interp.pts$y

    pseudox[,id]<-xxinterp
    pseudoy[,id]<-yyinterp

    ## pseudo dataset
    pseudo.data<-data.frame(ages=xxinterp, tms=yyinterp)

    ####################################################
    # Find inital inflection point using SHAPE CHANGE ##
    ####################################################


    initial.T<-changept(yyinterp~ip(xxinterp, sh=1),fir=TRUE)$chpt

    zstar0[id]<- initial.T

  }

  ## Estimated Initial inflection point
  zstar<-zstar0;

  ############################################################################
  ## while loop to estimate average trajectory and inflection point         ##
  ##            as well as individual trajecotries and inflection points    ##
  ############################################################################
  indvidual.Ta<-rep(100,n);
  ind.b.sd<-rep(100,n);
  ind.cp.boot<-array(0, dim=c(n, 2),
                     dimnames=list(paste("id",1:n,sep=""),c("95% lower", "95% upper"))
  )



  while(( mean.diff > tolerance)&&(iter<itermax)){

    iter<-iter+1
    old.zstar<-zstar;

    for(id in 1:length(obsID)){

      ## all inflection points are aligned at 0

      #id<-1;rm(id)
      subsubj<-which(simu.dat$SUBJID==obsID[id])
      subdat<-simu.dat[subsubj,]

      sub.xx.dat<-subdat$logAGE-old.zstar[id]
      #xxstar[,id]<-pseudox[,id]-old.zstar[id]
      sub.yy.dat<-subdat$TOTAL_MOTOR_SCORE

      ### save interpolation points ####
      scam.dat<-cbind( subdat$logAGE,  sub.yy.dat, sub.xx.dat)
      colnames(scam.dat)<-c("logAGE", "TMS","logDiff")

      ## Separate Data frame compare to logT
      subscam.data<-data.frame(scam.dat)



      subscam.data$fac1<-ifelse(subscam.data$logDiff<0, 1,0)
      subscam.data$fac2<-1-subscam.data$fac1
      subscam.data$fac<-c(rep(1, sum(subscam.data$fac1)),rep(2, sum(subscam.data$fac2)))

      scamfit.data[[id]]<- subscam.data

    }


    ########################################
    ## data generated with pseudo visits  ##
    ########################################
    scam.gen<-do.call("rbind.data.frame",  scamfit.data)


    ################################################################
    ## apply scam function to estimate longitudinal trajectory    ##
    ################################################################

    outscam<- scam(TMS~fac1+fac2+s(logDiff,k=20,bs="micx",m=2, by=fac1)+s(logDiff,k=20,bs="micv",m=2, by=fac2),family=gaussian(link="identity"), data= scam.gen)
    #scam.check(outscam)

    ## standard deviation for the within-subject error
    est.sigma.epsi<-sqrt(outscam$sig2)


    ####################################################
    ## new points at which predictions are estimated  ##
    ####################################################

    for (id in  1:n){

      ## pseudox: pseudo log transformed age for subject i at jth visit with length of num.interp
      ## zstar: initial inflection points

      xxstar[,id]<-pseudox[,id]-old.zstar[id]

      combine.two.vectors<-c(scamfit.data[[id]]$logDiff, xxstar[,id])
      merge.two.vectors[[id]]<-as.vector(unique(sort(combine.two.vectors)))

    }

    pseudo.time<-unlist(merge.two.vectors)
    newlogS<-seq(min(pseudo.time),max(pseudo.time), length.out=newl);

    newfac1<-ifelse(newlogS<0, 1,0)
    newfac2<-1-newfac1

    ### combined log transfromed Age after shifting x axis of longitudinal process
    mynewlogS<-data.frame(logDiff=newlogS, fac1=newfac1, fac2=newfac2) #, fac=newfac);


    ######################
    ## predicted values ##
    ######################
    predscam<-predict(outscam, mynewlogS, se.fit=TRUE) ;

    #########################################################
    ## estimated values and their estimate standard errors ##
    #########################################################
    tms.pred<-predscam$fit;
    tms.se.pred<-predscam$se.fit;


    ############################################
    ### Find global inflection points by       #
    ### ShapeChange                            #
    ############################################
    newlogS<-mynewlogS[,"logDiff"]
    s<-length(newlogS)
    tms.pred<-as.vector(tms.pred)



    ## global inflection point for all subject
    #glob.T<-changept(tms.pred~ip(newlogS, sh=1),fir=TRUE)$chpt
    #global.T<-glob.T;




    if(boot.ci==TRUE){

      glob.change.fit<-changept(tms.pred~ip(newlogS, sh=1),fir=TRUE,  ci = TRUE)
      glob.T<- glob.change.fit$chpt
      boot.dist<-glob.change.fit$msbt

      ## boostrap standarad deviation and 95% confidence interval
      glob.b.sd<-sd(boot.dist)
      glob.cp.boot<-quantile(sort(boot.dist), prob=c(0.025, 0.975))

    }else{

      glob.change.fit<-changept(tms.pred~ip(newlogS, sh=1),fir=TRUE)$chpt
      glob.T<-glob.change.fit$chpt


    }

    global.T<-glob.T;


    #########################################################################
    ## inflection points for all subjects                                  ##
    ## bootstrap 95% confidence interval for the inflection points         ##
    #########################################################################

    if(boot.ci==TRUE){

      glob.lower.cov.boot<-as.vector(glob.cp.boot["2.5%"])
      glob.upper.cov.boot<-as.vector(glob.cp.boot["97.5%"])

    }else{

      glob.lower.cov.boot<-NULL
      glob.upper.cov.boot<-NULL

    }

    #########################################
    ### Individual inflection points     ####
    #########################################

    for (id in  1:n){

      ### combined log transfromed Age after shifting x axis of longitudinal process
      #xnew<-as.vector(xxstar[,id]);
      xnew<-as.vector(merge.two.vectors[[id]]);
      xnewfac1<-ifelse(xnew<0, 1,0)
      xnewfac2<-1-xnewfac1
      xnewfac<-c(rep(1, sum(xnewfac1)), rep(2, sum(xnewfac2)))
      ind.logS<-data.frame(logDiff=xnew, fac1=xnewfac1, fac2=xnewfac2);# fac=xnewfac,


      ind.predscam<-predict(outscam, ind.logS, se.fit=TRUE)
      ind.tms.pred<-ind.predscam$fit
      ind.tms.se.pred<-ind.predscam$se.fit


      ### length of newlogS to be predicted
      fac.ind.logS<-ind.logS$logDiff*ind.logS$fac1+ind.logS$logDiff*ind.logS$fac2
      q<-length(fac.ind.logS)

      #ind.predict.tms[,id]<-ind.tms.pred
      #ind.str.predict.tms[,id]<-ind.tms.se.pred

      #newdlogS[,id]<-fac.ind.logS


      ######################################
      # methods to find inflection points ##
      ######################################

      ###################
      ### ShapeChange  ##
      ###################
      ind.tms.pred<-as.vector(ind.tms.pred)
      fac.ind.logS<-as.vector(fac.ind.logS)

      #Ta<-changept(ind.tms.pred~ip(fac.ind.logS, sh=1),fir=TRUE)$chpt


      if(boot.ci==TRUE){

        ind.Ta<-changept(ind.tms.pred~ip(fac.ind.logS, sh=1),fir=TRUE,  ci = TRUE)
        indvidual.Ta[id]<- ind.Ta$chpt
        ind.boot.dist<-ind.Ta$msbt
<<<<<<< HEAD

        ## boostrap standarad deviation and 95% confidence interval
        ind.b.sd[id]<-sd(ind.boot.dist)
        ind.cp.boot[id,]<-quantile(sort(ind.boot.dist), prob=c(0.025, 0.975))
=======

        ## boostrap standarad deviation and 95% confidence interval
        ind.b.sd[id]<-sd(ind.boot.dist)
        ind.cp.boot[id,]<-quantile(sort(ind.boot.dist), prob=c(0.025, 0.975))

      }else{

        ind.Ta<-changept(tms.pred~ip(newlogS, sh=1),fir=TRUE)$chpt
        indvidual.Ta[id]<-ind.Ta$chpt
>>>>>>> 1e019c4ccae424fd5a5f056475f3487fd881ab67

      }else{

<<<<<<< HEAD
        ind.Ta<-changept(tms.pred~ip(newlogS, sh=1),fir=TRUE)$chpt
        indvidual.Ta[id]<-ind.Ta$chpt


=======
>>>>>>> 1e019c4ccae424fd5a5f056475f3487fd881ab67
      }


      #indvidual.Ta[id]<-Ta

    }



    ##################################################
    ### Gethering results for inflectio points   #####
    ##################################################

    zstar<-old.zstar+indvidual.Ta;

    ## convergence criterion
    mean.diff<-mean(abs(old.zstar-zstar));

  }

  ###########################################################
  ## To pring number of iteration and the convergence info ##
  ## print(c(iter, mean.diff))                             ##
  ###########################################################

  #################################
  #### fit inflection point model #
  #################################

  fit.data<-data.frame(cbind(zstar, base.cov))

  ## baseline covariate are standardized ##
  fit.data$std.age<-(fit.data$AGE-mean(fit.data$AGE))/sd(fit.data$AGE)
  fit.data$std.logage<-(fit.data$logAGE-mean(fit.data$logAGE))/sd(fit.data$logAGE)
  fit.data$std.CAG<-(fit.data$CAG-mean(fit.data$CAG))/sd(fit.data$CAG)

  ## fit linear model to estimate fixed effects ##
  fit0<-lm(zstar~std.CAG+factor(gender), data=fit.data)

  ## summarize estimation
  summary.fit0<-round(summary(fit0)$coefficients, 3)

  ## estimate random error term in the inflection point
  est.sigma.ei0<-summary(fit0)$sigma  # same as "sqrt(sum(resid(fit)^2)/(n-2))"
  res.var0<-2*(est.sigma.ei0)^{4}/(n-2)
  est.str.sigma.ei0<-sqrt(res.var0) ## str of sigma


  ###########################
  # Produce table 4         #
  ###########################

  nonpara_summary_table4<-array(0, dim=c((dim(summary.fit0)[1]+2), dim(summary.fit0)[2]),
                                dimnames=list(c("beta0", "beta_CAG", "beta_gender","sigma_u", "sigma_eps"), c(colnames(summary.fit0))
                                ))

  nonpara_summary_table4[1,]<-summary.fit0[1,]
  nonpara_summary_table4[2,]<-summary.fit0[2,]
  nonpara_summary_table4[3,]<-summary.fit0[3,]
  nonpara_summary_table4[4,]<-c(est.sigma.ei0, NA,NA,NA)
  nonpara_summary_table4[5,]<-c(est.sigma.epsi, NA,NA,NA)

  ## round up to 3 decimal numbers
  nonpara_summary_table4<-round(nonpara_summary_table4, 3)


  ########################################################################
  # Recall suborg: list of dataset for the paprametric NLME approach  ####
  ########################################################################


  #########################
  # STANDARDIZED CAG   ####
  #########################
  stdardized.age<-(base.cov[,"logAGE"]-mean(base.cov[,"logAGE"]))/sd(base.cov[,"logAGE"])
  stdardized.CAG<-(base.cov[,"CAG"]-mean(base.cov[,"CAG"]))/sd(base.cov[,"CAG"])
  #stdardized.SDMT<-(base.cov[,"SDMT"]-mean(base.cov[,"SDMT"]))/sd(base.cov[,"SDMT"])
  #stdardized.SC<-(base.cov[,"SC"]-mean(base.cov[,"SC"]))/sd(base.cov[,"SC"])
  #stdardized.SW<-(base.cov[,"SW"]-mean(base.cov[,"SW"]))/sd(base.cov[,"SW"])
  #stdardized.SI<-(base.cov[,"SI"]-mean(base.cov[,"SI"]))/sd(base.cov[,"SI"])

  for(i in 1:length(obsID)){
    #i<-1; rm(i)
    ll<-dim(suborg[[i]])[1]
    suborg[[i]]$stdz_bage<-rep(as.numeric(stdardized.age[i]),ll)
    suborg[[i]]$stdz_CAG<-rep(as.numeric(stdardized.CAG[i]),ll)

  }




  ## create data.frame from list type of data
  gendat.para<-do.call("rbind.data.frame", suborg)

  ## Create grouped Data
  newgrdata.para<-groupedData(TOTAL_MOTOR_SCORE~logAGE|SUBJID/gender, data=gendat.para, order.groups=FALSE)



  ###################################################################################################################################################
  ## Fit nonlinear mixed model using nlme() in NLME package to estimate fixed effects and random effects
  ## in longitudinal model, which assumped nonlinear logistic function
  ## tms: longitudinal response vector for a subject i and jth visit.
  ## logistf: nonlinear model that data follow. This function is defined in logistf().
  ## theta1: fixed effect. parametrized in logistf
  ## theta2: log transformed inflection points, which is modeled by linear relationship with
  ##         subject specific covariate W.cov, containig fixed effect beta and random effects u
  ##         We do not assume there is constant term.
  ## logS: main covariate in the logitudinal model
  ## model: tms ~logistf(theta1, theta2,  logS)
  ## fixed: two sided linear model in the form of f1~x1, where f1: names of parameters, x1: covariates in the linear relationship with f1
  ##       eg. theat1 ~1
  ## random: two sided formula in the form of r1~x1, where r1: names of parameters, x1 specifies the random effects model for the parameter r1
  ## groups: ~g1 or ~g1/g2../gQ, specify the partitions of the data over which the random effects vary
  ## start: list of initial estimates for the fixed effects and random effects, eg. theta1=5, W.COV=2.2
  ## method: "REML" or "ML", modle is fit by maximizing the restricted log liklihood or log-likelihood.
  ## verbose: TRUE means information on the evoluation of the interative algorithm is printed. Default is FALSE.
  ###################################################################################################################################################



  tms.nlme1<-nlme(TOTAL_MOTOR_SCORE~logistft(theta1, theta2, theta3, logAGE), fixed=list(theta1~1,  theta2~1+stdz_CAG+gender, theta3~1) , random=pdDiag(theta1+theta2+theta3~1),
                  data= newgrdata.para, groups=~SUBJID,  start=list(fixed=c(6, 3.5,0,0,1)), method="ML", verbose=FALSE, na.action=na.omit) #+gender+log(hdage_nobase)+SDMT+STROOP_INTERFERENCE


  ###################
  ### Fixed Effects #
  ###################
  summary.fixed1<-summary(tms.nlme1)$tTable


  #####################
  ##  scale parameter #
  #####################

  ## scale parameter for the within-individual in the longitudinal model
  sig.eps.para<-tms.nlme1$sigma

  ## scale parameter for the inflection point
  sig.u.para<-as.numeric(VarCorr(tms.nlme1)["theta2.(Intercept)", "StdDev"])



  ###############################################################
  ## Extract Random effects with augmented data from groupedData#
  ###############################################################
  randef1<-random.effects(tms.nlme1)

  ##################################################
  ## estimate inflection points (random effects)  ##
  ## estimate bias for inflection points          ##
  ##################################################
  beta0.nlme1<-summary(tms.nlme1)$tTable[,"Value"]["theta2.(Intercept)"]
  beta1.nlme1<-summary(tms.nlme1)$tTable[,"Value"]["theta2.stdz_CAG"]
  beta3.nlme1<-summary(tms.nlme1)$tTable[,"Value"]["theta2.gender"]

  ## estimate inflection points
  est.logT.nlme1<-beta0.nlme1+beta1.nlme1*as.vector(stdardized.CAG)+beta3.nlme1*base.cov[,"gender"]+randef1[,"theta2.(Intercept)"]
  ## estimate average of inflection points
  avr.logT.nlme1<-mean(est.logT.nlme1)



  ## prediction ft
  num.ID<-unique(newgrdata.para[,"SUBJID"])

  #m=30;
  predage<-seq(2.5, 4.5, length=m)
  para.pred <-matrix(100,nrow=length(predage), ncol=n,
                     dimnames=list(paste("visit",1:length(predage),sep=""),paste("ID",1:n,sep="")))

  #logDiff.para <-vector("list", n)
  for(id in 1:n){
    #rm(id); id<-1
    subdata.para<-newgrdata.para[newgrdata.para[,"SUBJID"]==num.ID[id], ]

    newage1<-data.frame(logAGE=predage,stdz_CAG=unique(subdata.para[,"stdz_CAG"]), gender=unique(subdata.para[,"gender"]), SUBJID=num.ID[id]) #,stdz_SDMT=unique(subdata[,"stdz_SDMT"]),stdz_SC=unique(subdata[,"stdz_SC"]),
    #stdz_SW=unique(subdata[,"stdz_SW"]),stdz_SI=unique(subdata[,"stdz_SI"]), gender=unique(subdata[,"gender"]), SUBJID=num.ID[id])
    #logDiff.para[[id]]<-newage1[, 1]-est.logT.nlme1[id]
    pred.sub1<-predict(tms.nlme1, newage1, level=0:1)
    para.pred[,id]<-pred.sub1[,"predict.fixed"]

  }

  avr.para.pred<-rowMeans(para.pred)


  ###########################
  # Produce table 3         #
  ###########################

  para_summary_table4<-array(0, dim=c((dim(summary.fixed1)[1]+2), dim(summary.fixed1)[2]),
                             dimnames=list(c("theta1", "theta2", "beta0", "beta_CAG", "beta_gender","sigma_u", "sigma_eps"), c(colnames(summary.fixed1))
                             ))

  para_summary_table4[1,]<-summary.fixed1[1,]
  para_summary_table4[2,]<-summary.fixed1[5,]
  para_summary_table4[3,]<-summary.fixed1[2,]
  para_summary_table4[4,]<-summary.fixed1[3,]
  para_summary_table4[5,]<-summary.fixed1[4,]
  para_summary_table4[6,]<-c(sig.u.para, NA,NA,NA, NA)
  para_summary_table4[7,]<-c(sig.eps.para, NA,NA,NA, NA)

  ## round up to 3 decimal numbers
  para_summary_table4<-round(para_summary_table4, 3)







  #############################################
  # To produce Figure 1 in manuscript        ##
  #############################################


  group<-1; aa<-list();

  ## obtain data ##

  for (id in 1:80){

    aa[[id]]<-scamfit.data[[id]]
    aa[[id]]$person<-id
    #aa[[id]]$para.logDiff<-logDiff.para[[id]]

  }

  aa.gen<-do.call("rbind.data.frame",aa)

  ###########################################################
  ### To plot average trajectory with confidence interval ###
  ###########################################################

  population.pred<-data.frame(AGE=NA, TMS=tms.pred,  logDiff=newlogS, para.TMS=avr.para.pred, fac1=NA, fac2=NA, fac=NA, person=NA, se=tms.se.pred)


  #data from aa.gen
  p<-ggplot(data=aa.gen, aes(x =logDiff, y = TMS, group = person))

  p<-p+geom_line(color="lightgray", alpha=1)+
    ## spaghetti plot
    stat_smooth(data=aa.gen, method="loess",aes(group=1,colour="blue"), fill="lightblue", linetype="solid", size=1,alpha=0.6)+ #, , se=FALSE, linetype="dashed")
    ## parametric approach
    geom_line(data=population.pred, aes(x =logDiff, y=para.TMS, colour="orange"), size=1, alpha=0.6,  linetype="dotdash")+
    ## scam smooth plot with 95% confidence interval
    geom_ribbon(data=population.pred, aes(ymin =TMS-qnorm(0.975)*se, ymax = TMS+qnorm(0.975)*se), fill = "lightpink", alpha=0.6)+
    geom_line(data=population.pred, aes(x =logDiff, y =TMS, colour="red"),  linetype="dashed", size=1, alpha=0.6)






  ## simple spaghetti plot
  size <- 12
  p +geom_vline(aes(xintercept=global.T),  linetype="dashed")+
    geom_segment(aes(x = glob.lower.cov.boot, y = 0, xend =glob.lower.cov.boot, yend = 0.5))+
    annotate("rect", xmin = glob.lower.cov.boot, xmax = glob.upper.cov.boot, ymin = 0, ymax = 1,
             alpha = .5, colour="lightgrey")+
    scale_colour_manual(name  ="Method",
                        values=c( "blue","orange", "red"), labels =c("Lowess","Parametric NLME","Multi-Stage Nonparametric Method"))+
    scale_linetype_manual(values=c("solid", "dotdash", "dashed" ))+
    guides(color=guide_legend(override.aes=list(linetype=c("solid","dotdash","dashed"), fill=NA)))+
    xlab("Shifted logAges") +
    ylab("Total Motor Score") +
    ggtitle("Motor Sign Decline")+
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          #panel.background = element_rect(fill = "white"),panel.grid.minor=element_blank(),
          #title="Motor-Sign Decline",
          axis.text.x=element_text(size=size),
          axis.text.y=element_text(size=size),
          axis.title.x=element_text(size=size),
          axis.title.y=element_text(size=size,angle=90),
          plot.title=element_text(size=15,hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(colour = 'black', angle = 0, size = 12), #hjust = 3, vjust = 3, face = 'bold')
          legend.position = c(.46, .95),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6))

  ggsave("HDtrjectory.pdf",width=8,height=8)
  #graphics.off()

  glob.boot.ci<-cbind(glob.lower.cov.boot, glob.upper.cov.boot)

  return(list(nonpara_summary_table4=nonpara_summary_table4,
              para_summary_table4=para_summary_table4 #glob.boot.ci=glob.boot.ci

  ))

}




<<<<<<< HEAD
=======



>>>>>>> 1e019c4ccae424fd5a5f056475f3487fd881ab67
