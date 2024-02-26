#' pmvalsampsize
#' - Sample Size for External Validation of a Prediction Model
#'
#' @description \code{pmvalsampsize} computes the minimum sample size required for the external validation of an existing multivariable prediction model using the criteria proposed by Archer \emph{et al}. and Riley \emph{et al}.  \code{pmvalsampsize} can currently be used to calculate the minimum sample size for the external validation of models with binary outcomes. \cr \cr \emph{Continuous and survival (time-to-event) outcome model calculations are a work in progress}. \cr \cr
#'
#' @return A list including a matrix of calculated sample size requirements for each criteria defined under \emph{'Details'}, and a series of vectors of parameters used in the calculations as well as the the final recommended minimum sample size and number of events required for external validation.
#' 
#' @details A series of criteria define the sample size needed to ensure precise estimation of key measures of prediction model performance, allowing conclusions to be drawn about whether the model is potentially accurate and useful in a given population of interest.\cr \cr
#'
# For \bold{continuous outcome models}, Archer \emph{et al}. specify four criteria to calculate the sample size (N) needed for: \cr
# i) precise estimation of the R-squared performance, \cr
# ii) precise estimation of the calibration-in-the-large (CITL), \cr
# iii) precise estimation of the calibration slope (c-slope), and \cr
# iv) precise estimation of the residual variance for CITL & c-slope. \cr \cr
#
# The sample size calculation requires the user to pre-specify (e.g. based on previous evidence from the development study or information known from the validation sample if available) the following;
# - the anticipated R-squared of the model
# - the target precision (CI width) for the R-squared
# - the anticipated CITL performance
# - the target precision (CI width) for the CITL
# - the anticipated c-slope
# - the target precision (CI width) for the c-slope
# - the anticipated variance of observed values in the validation sample \cr \cr
#'
#' For \bold{binary outcomes}, there are three criteria to calculate the sample size (N) needed for: \cr
#' i) precise estimation of the Observed/Expected (O/E) statistic, \cr
#' ii) precise estimation of the calibration slope (c-slope), and \cr
#' iii) precise estimation of the c-statistic.
#'
#' The sample size calculation requires the user to pre-specify the following; \cr
#' - the outcome event proportion \cr
#' - the anticipated O/E performance \cr
#' - the target precision (CI width) for the O/E \cr
#' - the anticipated c-slope \cr
#' - the target precision (CI width) for the c-slope \cr
#' - the anticipated c-statisitic performance \cr
#' - the target precision (CI width) for the c-statisitic \cr
#' - the distribution of estimated probabilities from the model, ideally specified on the log-odds scale - AKA the Linear Predictor (LP) \cr \cr
#'
#' With thanks to Richard Riley for helpful feedback
#'
#'
#' @author Joie Ensor (University of Birmingham, j.ensor@bham.ac.uk)
#'
#'
#' @param type specifies the type of analysis for which sample size is being calculated
#'      \itemize{
# #'          \item \code{"c"} specifies sample size calculation for a prediction model with a continuous outcome
#'          \item \code{"b"} specifies sample size calculation for a prediction model with a binary outcome
#'      }
#'
#' @param cslope specifies the anticipated c-slope performance in the validation sample.  Default conservatively assumes perfect c-slope=1.  The value could alternatively be based on a previous validation study for example.  For \bold{binary outcomes} the c-slope calculation requires the user to specify a distribution for the assumed LP in the validation sample (or alternatively the distribution of predicted probabilities in the validation sample). See \code{lp*()} options below.
#'
#' @param csciwidth specifies the target CI width (acceptable precision) for the c-slope performance. Default assumes CI width=0.2.
#'
# #' \bold{Binary outcome options}
#'
#' @param prevalence specifies the overall outcome proportion (for a prognostic model) or overall prevalence (for a diagnostic model) expected within the model validation sample.  This is a \bold{required input}.  This should be derived based on previous studies in the same population or directly from the validation sample if to hand.
#'
#' @param simobs specifies the number of observations to use when simulating the LP distribution for c-slope calculation in criteria 2.  Default observations=1,000,000.  Higher \code{simobs()} values will reduce random variation further.
#'
#' @param cstatistic specifies the anticipated c-statistic performance in the validation sample.  This is a \bold{required input}. May be based on the optimism-adjusted c-statistic reported in the development study for the existing prediction model. Ideally, this should be an optimism-adjusted c-statistic.  NB: This input is also used when using the \code{lpcstat()} option.
#'
#' @param cstatciwidth specifies the target CI width (acceptable precision) for the c-statistic performance.  Default assumes CI width=0.1.
#'
#' @param oe specifies the anticipated O/E performance in the validation sample.  Default conservatively assumes perfect O/E=1.
#'
#' @param oeciwidth specifies the target CI width (acceptable precision) for the E/O performance.  Default assumes CI width=0.2.  The choice of CI width is context specific, and depends on the event probability in the population.  See Riley et al. for further details.
#'
#' @param oeseincrement sets the increment by which to iterate when identifying the SE(ln(OE)) to meet the target CI width specified for OE.  The default iteration increment=0.0001.  In the majority of cases this will be suitably small to ensure a precise SE is identified.  The user should check the output table to ensure that the target CI width has been attained and adjust the increment if necessary.
#'
#' @param lpnormal defines parameters to simulate the LP distribution for criteria 2 from a normal distribution.  The user must specify the mean and standard deviation (in this order) of the LP distribution.
#'
# #' @param lpskewednormal defines parameters to simulate the LP distribution for criteria 2 from a skewed normal distribution.  The user must specify the mean, variance, skew and kurtosis parameters (in this order) of the LP distribution.  NB: \code{lpskewednormal()} option can take a little longer than other distributional assumptions.
#'
#' @param lpbeta defines parameters to simulate the distribution of predicted probabilities for criteria 2 from a beta distribution.  The user must specify the alpha and beta parameters (in this order) of the probability distribution. The LP distribution is then generated internally using this probability distribution.
#'
#' @param lpcstat defines parameters to simulate the LP distribution for criteria 2 assuming that the distribution of events and non-events are normal with a common variance.  The user specifies a single input value - the expected mean for the non-events distribution. This could be informed by clinical guidance.  However, this input is taken as a starting value and an iterative process is used to identify the most appropriate values for the event and non-event distributions so as to closely match the anticipated prevalence in the validation sample.  NB: this approach makes strong assumptions of normality and equal variances in each outcome group, which may be unrealistic in most situations.
#'
#' @param tolerance for use with \code{lpcstat} option.  Sets the tolerance for agreement between the simulated and expected event proportion during the iterative procedure for calculating the mean for the non-events distribution.
#'
#' @param increment for use with \code{lpcstat} option.  Sets increment by which to iterate the value of the mean for the non-events distribution.  Trial and error may be necessary as it is dependent on how close the initial input for the non-event mean in \code{lpcstat} is to the required value.  If the process takes a particularly long time then the user could try an alternative increment value, or an alternative non-event mean value in \code{lpcstat}.  The \code{trace} option may be useful in such circumstances.
#'
#' @param trace for use with \code{lpcstat} option.  Specifies that a trace of the values obtained in each iteration when identifying the non-event mean is output.  Useful when finding the appropriate values for \code{lpcstat} & \code{increment()} is proving difficult!
#'
#' @param seed specifies the initial value of the random-number seed used by the random-number functions when simulating data to approximate the LP distribution for criteria 2.
#'
#' @param graph specifies that a histogram of the simulated LP distribution for criteria 2 is produced.  The graph also details summary statistics for the simulated distribution.  Useful option for checking the simulated LP distribution against the source of input parameters.  Also useful for reporting at publication.
#'
#' @param sensitivity specifies the anticipated sensitivity performance in the validation sample at the chosen risk threshold (specified using \code{threshold}). If sensitivity and specificity are not provided then \code{pmvalsampsize} uses the simulated LP distribution from criteria 2 and the user-specified risk threshold to estimate the anticipated sensitivity and specificity to be used in calculation of net benefit. NB: net benefit criteria is not calculated if either i) \code{sensitivity}, \code{specificity} and \code{threshold} or ii) \code{threshold} option are not specified.
#'
#' @param specificity specifies the anticipated specificity performance in the validation sample at the chosen risk threshold (specified using \code{threshold}). If sensitivity and specificity are not provided then \code{pmvalsampsize} uses the simulated LP distribution from criteria 2 and the user-specified risk threshold to estimate the anticipated sensitivity and specificity to be used in calculation of net benefit. NB: net benefit criteria is not calculated if either i) \code{sensitivity}, \code{specificity} and \code{threshold} or ii) \code{threshold} option are not specified.
#'
#' @param threshold specifies the risk threshold to be used for calculation of net benefit performance of the model in the validation sample. If sensitivity and specificity are not provided then \code{threshold} must be given in order for \code{pmvalsampsize} to assess sample size requirements for net benefit. NB: net benefit criteria is not calculated if either i) \code{sensitivity}, \code{specificity} and \code{threshold} or ii) \code{threshold} option are not specified.
#'
#' @param nbciwidth specifies the target CI width (acceptable precision) for the standardised net benefit performance.  Default assumes CI width=0.2.  The choice of CI width is context specific.  See Riley et al. for further details.
#'
#' @param nbseincrement sets the increment by which to iterate when identifying the SE(standardised net benefit) to meet the target CI width specified for standardised net benefit.  The default iteration increment=0.0001.  In the majority of cases this will be suitably small to ensure a precise SE is identified.  The user should check the output table to ensure that the target CI width has been attained and adjust the increment if necessary.
#'
# #' \bold{Continuous outcome options}
#'
# #' @param rsquared specifies the anticipated value of R-squared in the validation sample.  This is a required input.  May be based on the optimism-adjusted R-squared reported in the development study for the existing prediction model. Ideally, this should be an \emph{optimism-adjusted} R-squared value, not the apparent R-squared value, as the latter is optimistic (biased).
#'
# #' @param rsqciwidth specifies the target CI width (acceptable precision) for the R-squared performance.  Default assumes CI width=0.1.
#'
# #' @param citl specifies the anticipated value of CITL in the validation sample.  Default conservatively assumes perfect CITL=0.
#'
# #' @param citlciwidth specifies the target CI width (acceptable precision) for the CITL performance.  The choice of CI width is context specific, and depends on the scale of the outcome.  See Archer et al. for further details.
#'
# #' @param varobs the anticipated variance of observed values in the validation sample.  This is a required input.  This could again be based on the previous development study, or on clinical knowledge.
#'
# #' @param mmoe multiplicative margin of error (MMOE) acceptable for calculation of the variance for c-slope and CITL in the calibration model.  The default is a MMOE of 10\%.  See Archer et al. reference below for further information.
#'
#' @importFrom stats qnorm rbeta rbinom rnorm
#' @import graphics
#' @import utils
#' @importFrom pROC roc
#'
#' @examples
#' # Examples based on those included in the papers referenced below by
#' # Riley et al. & Archer et al. published in Statistics in Medicine.
#' 
#' # Note that the examples below use a very low simulation sample for criteria
#' # 2 for brevity. Sizes matching the default or higher are recommended to 
#' # minimise uncertainty in the calculated sample size requirements.
#'
#' # Binary outcomes (Logistic prediction models)
#' # Use pmvalsampsize to calculate the minimum sample size required to
#' # externally validate an existing multivariable prediction model for a
#' # binary outcome (e.g. mechanical heart valve failure).  Based on previous
#' # evidence, the outcome prevalence is anticipated to be 0.018 (1.8%) and the
#' # reported c-statistic was 0.8.  The LP distribution was published and
#' # appeared normally distributed with mean(SD) of -5(2.5).  We target default
#' # CI widths for all but O/E CI width=1 (see Riley et al. for details). We can
#' # use the graph option to check the simulated distribution is appropriate.
#'
#' pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8,
#' lpnorm = c(-5,2.5), oeciwidth = 1, simobs = 1000)
#'
#'
#' # Alternatively, lets assume that the authors provided a distribution of
#' # predicted probabilities (e.g. as part of a calibration plot).  We can use
#' # this to specify parameters for a beta distribution to simulate the LP
#' # distribution as below:
#' pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8,
#' lpbeta = c(0.5,0.5), oeciwidth = 1, simobs = 1000)
#'
#' # Finally, we can use the anticipated c-statistic to simulate the event and
#' # non-event distributions assuming normality and common variances.  We input
#' # a starting value for the mean of the non-events distribution as below:
#'
#' pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8,
#' lpcstat = -4.7, oeciwidth = 1, seed = 1234, simobs = 10000)
#'
# #' ## Continuous outcomes (Linear prediction models)
# #' # Use pmvalsampsize to calculate the minimum sample size required for validation of a multivariable prediction model for a continuous outcome (here, fat-free mass in children).  We know the existing prediction model has been previously validated in a small sample with an R-squared of 0.9.  Riley et al. compute the anticipated variance of outcome values in the validation sample based on reported upper and lower quartile values deriving var(Y)=0.089 (see Riley et al. for details).  We conservatively assume perfect calibration performance at validation.  If we target a R-sq CI width=0.1, CITL CI width=0.08 and c-slope CI width=0.2, then we can run the following command:
#'
# #' pmvalsampsize(type = "c", rsquared = 0.9, varobs = 0.089, citlciwidth = 0.08)
#'
#' @references  Archer L, Snell K, Ensor J, Hudda MT, Collins GS, Riley RD.  Minimum sample size for external validation of a clinical prediction model with a continuous outcome. \emph{Statistics in Medicine. 2020.}
#'
#' @references Riley RD, Debray TPA, Collins G, Archer L, Ensor J, van Smeden M, Snell KIE.  Minimum sample size for external validation of a clinical prediction model with a binary outcome. \emph{Statistics in Medicine. 2021.}
#'
#'
#' @export
pmvalsampsize <- function(type,
                       cslope = 1,
                       csciwidth = 0.2,
                       oe = 1,
                       oeciwidth = 0.2,
                       cstatistic = NA,
                       cstatciwidth = 0.1,
                       simobs = 1000000,
                       lpnormal = NA,
                       #lpskewednormal = NA,
                       lpbeta = NA,
                       lpcstat = NA,
                       tolerance = 0.0005,
                       increment = 0.1,
                       oeseincrement = 0.0001,
                       graph = FALSE,
                       trace = FALSE,
                       prevalence = NA,
                       seed = 123456,
                       sensitivity = NA,
                       specificity = NA,
                       threshold = NA,
                       nbciwidth = 0.2,
                       nbseincrement = 0.0001) {
                      # rsquared = 0,
                      # rsqciwidth = 0.1,
                      # citl = 0,
                      # citlciwidth = 0,
                     #  varobs = 0,
                     #  mmoe=1.1

  # error checking
  pmvalsampsize_errorcheck(type=type,prevalence=prevalence,cstatistic=cstatistic,oe=oe,oeciwidth=oeciwidth,cslope=cslope,csciwidth=csciwidth,cstatciwidth=cstatciwidth,simobs=simobs,lpnormal=lpnormal,lpbeta=lpbeta,lpcstat=lpcstat,tolerance=tolerance,increment=increment,oeseincrement=oeseincrement,seed=seed,graph=graph,trace=trace,sensitivity=sensitivity,specificity=specificity,threshold=threshold,nbciwidth=nbciwidth,nbseincrement=nbseincrement)

  # choose function based on analysis type
  # if (type == "c") out <- pmvalsampsize_cont()

  if (type == "b") out <- pmvalsampsize_bin(prevalence=prevalence,cstatistic=cstatistic,oe=oe,oeciwidth=oeciwidth,cslope=cslope,csciwidth=csciwidth,cstatciwidth=cstatciwidth,simobs=simobs,lpnormal=lpnormal,lpbeta=lpbeta,lpcstat=lpcstat,tolerance=tolerance,increment=increment,oeseincrement=oeseincrement,seed=seed,graph=graph,trace=trace,sensitivity=sensitivity,specificity=specificity,threshold=threshold,nbciwidth=nbciwidth,nbseincrement=nbseincrement)


  est <- out
  class(est) <- "pmvalsampsize"
  est

}


#' @export
print.pmvalsampsize <- function(x, ...) {

  print(x$results_table)

  # if (x$type == "continuous") {
  #   cat(" \n Minimum sample size required for new model development based on user inputs =", x$sample_size, " \n \n ")
  #   cat("* 95% CI for intercept = (",x$int_lci,", ", x$int_uci, "), for sample size n = ",x$sample_size,sep = "")
  # }
  if (x$type == "binary") {
    cat(" \n Minimum sample size required for model validation based on user inputs = ", x$sample_size, ", \n ",sep = "")
    cat("with ", ceiling(x$events), " events (assuming an outcome prevalence = ", x$prevalence, ") \n \n ",sep = "")

  }
}

#' @export
summary.pmvalsampsize <- function(object, ...) {

  print(object$results_table)
#
#   if (object$type == "continuous") {
#     cat(" \n Minimum sample size required for new model development based on user inputs =", object$sample_size, " \n \n ")
#     cat("* 95% CI for intercept = (",object$int_lci,", ", object$int_uci, "), for sample size n = ",object$sample_size,sep = "")
#   }
  if (object$type == "binary") {
    cat(" \n Minimum sample size required for new model development based on user inputs = ", object$sample_size, ", \n ",sep = "")
    cat("with", ceiling(object$events), "events (assuming an outcome prevalence =", object$prevalence, ") and an EPP =", object$EPP, " \n \n ")

  }
}
