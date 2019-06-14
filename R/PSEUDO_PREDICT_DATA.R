#' Generate a simulated dataset to reproduce results in the analysis study.
#'
#' A pseudo dataset containing subject's id, number of clinical visits, total motor scores,
#' time of the clinical diagnosis, ages at clinical visits, CAG repeats, Gender, and log transformed ages at clinical visits.
#'
#'
#' @format A dataframe including 80 subjects' longitudinal data with 8 variables:
#' \describe{
#'   \item{SUBJID}{identification number of a subject.}
#'   \item{event}{number of events (visits) for a subject.}
#'     \item{TOTAL_MOTOR_SCORE}{total motor score measured at each visit for a subject.}
#'   \item{TRUE_INFL_POINT}{simulated true inflection point for each subject. In practice, this is time of the clinically diagnosed disease onset.}
#'     \item{AGE}{age at each visit for a subject.}
#'   \item{CAG}{length of CAG repeats for a subject.}
#'     \item{gender}{male for 0; female for 1.}
#'   \item{logAGE}{log-transformed ages corresponding \code{AGE} at each visit for a subject.}
#' }
#'
#' @source see data-raw/PSEUDO_PREDICT_HD.R
#'
#'
#'
#'
#'
"PSEUDO_PREDICT_HD"
