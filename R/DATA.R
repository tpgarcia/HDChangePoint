#' Generate a simulated dataset to reproduce results in the analysis study.
#'
#' A pseudo dataset containing subject's id, number of clinical visits, total motor scores,
#' time of the clinical diagnosis, ages at clinical visits, CAG repeats, Gender, and log transformed ages at clinical visits.
#'
#'
#' @format A dataframe including 80 subjects' longitudinal data with 8 variables:
#' \describe{
#'   \item{n}{number of sample size}
#'   \item{model}{a character string to choose a nonlinear model \code{w()}, \code{"logist"} or \code{"arctan"}.}
#'     \item{p}{dimension of parameters corresponding to subject specific covariates including a constant term. The default is 2.}
#'   \item{bb0}{an intercept term for the inflection point model.}
#'     \item{bb}{(p-1)-length of parameters of subject specific covariates in the inflection point model.}
#'   \item{x.sd}{standard deviation of log transformed ages at clinical visits to generate time points.}
#'     \item{v1}{the least number of visits for all subjects.}
#'   \item{v2}{the largest number of visits for all subjects.}
#'    \item{dist}{distribution for the error terms for the longitudinal model.}
#'   \item{eps.sd}{standard deviation for the error term for the longitudinal model.}
#'   \item{u.sd}{standard deviation for the error term for the inflection point model.}
#' }
#'
#' @source see data-raw/PSEUDO_PREDICT_HD.R
#'
#'
#'
#'
#'
"PSEUDO_PREDICT_HD"

#' Generated simulation results from the multi-stage nonparametric approach.
#'
#' nonpara_logist_visit1 is the simulation results data from the multi-stage nonparametric method. We use this data to reproduce simulation results in our manuscript.
#' The true data were generated from the logistic model and the total number of visit numbers \eqn{m=5, 6,} or 7.
#' The generated data are based on 100 simulation runs of the function \code{sim.nonpara()} in multi_stage_nonpara.R.
#'
#'
#' @format A 100 lists of results and each list includes :
#' \describe{
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
#' }
#'
#' @source see returns of the function \code{sim.nonpara()} in R/multi_stage_nonpara.R
#'
#'
#'
#'
#'
"nonpara_logist_visit1"


#' Generated simulation results from the multi-stage nonparametric approach.
#'
#' nonpara_logist_visit2 is the simulation results data from the multi-stage nonparametric method. We use this data to reproduce simulation results in our manuscript.
#' The true data were generated from the logistic model and the total number of visit numbers \eqn{m=8, 9,} or 10.
#' The generated data are based on 100 simulation runs of the function \code{sim.nonpara()} in multi_stage_nonpara.R.
#'
#'
#' @format A 100 lists of results and each list includes :
#' \describe{
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
#' }
#'
#' @source see returns of the function \code{sim.nonpara()} in R/multi_stage_nonpara.R
#'
#'
#'
#'
#'
"nonpara_logist_visit2"


#' Generated simulation results from the multi-stage nonparametric approach.
#'
#' nonpara_logist_visit3 is the simulation results data from the multi-stage nonparametric method. We use this data to reproduce simulation results in our manuscript.
#' The true data were generated from the logistic model and the total number of visit numbers \eqn{m=8, 9,} or 10.
#' The generated data are based on 100 simulation runs of the function \code{sim.nonpara()} in multi_stage_nonpara.R.
#'
#'
#' @format A 100 lists of results and each list includes :
#' \describe{
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
#' }
#'
#' @source see returns of the function \code{sim.nonpara()} in R/multi_stage_nonpara.R
#'
#'
#'
#'
#'
"nonpara_logist_visit3"
