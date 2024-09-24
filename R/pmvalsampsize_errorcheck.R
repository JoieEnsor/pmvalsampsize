pmvalsampsize_errorcheck <- function(type,prevalence,cstatistic,oe,oeciwidth,cslope,csciwidth,cstatciwidth,simobs,lpnormal,lpskewednormal,lpbeta,lpcstat,tolerance,increment,oeseincrement,seed,graph,trace,sensitivity,specificity,threshold,nbciwidth,nbseincrement) {

  # syntax checks
  if (type != "c" & type != "b" & type != "s") stop('type must be "c" for continuous, "b" for binary or "s" for survival')
  if (is.numeric(simobs) == F) stop('simobs must be an integer')
  if (simobs != round(simobs)) stop('simobs must be an integer')


  # parameter restrictions
  #if (shrinkage <= 0 | shrinkage >= 1) stop('shrinkage must be between 0 and 1')

  # parameters for continuous
  # if (type == "c") {
  #   # syntax
  #   if (is.na(rsquared)) stop('rsquared must be specified for continuous outcome models')
  #   # parameters needed
  #   if (is.na(sd)) stop('sd must be specified for continuous sample size')
  #   if (is.na(intercept)) stop('intercept must be specified for continuous sample size')
  #   # parameter conditions
  #   if (is.numeric(rsquared) == F) stop('rsquared must be numeric')
  #   if (is.numeric(intercept) == F) stop('intercept must be numeric')
  #   if (is.numeric(sd) == F) stop('sd must be numeric')
  #   # parameters not needed
  #   if (is.na(prevalence) == F) stop('prevalence not required for continuous sample size')
  #   if (is.na(rate) == F) stop('rate not required for continuous sample size')
  #   if (is.na(timepoint) == F) stop('timepoint not required for continuous sample size')
  #   if (is.na(meanfup) == F) stop('meanfup not required for continuous sample size')
  # }


  # parameters for binary
  if (type == "b") {
    # parameters needed
    if (is.na(prevalence)) stop('prevalence must be specified for binary sample size')
    if (is.na(cstatistic)) stop('cstatistic must be specified for binary outcome models')

    if (is.na(lpnormal[1])==F) {

      if (is.na(lpbeta[1])==F | is.na(lpcstat[1])==F) {
        stop('Only one LP distribution option can be specified')

      }

    }
    else if (is.na(lpbeta[1])==F) {

      if (is.na(lpcstat[1])==F) {
        stop('Only one LP distribution option can be specified')

      }

    }
    else if (is.na(lpnormal[1]) & is.na(lpbeta[1]) & is.na(lpcstat[1])) {

      stop('An LP distribution must be specified')

    }

    # parameter conditions
    if (is.numeric(prevalence) == F) stop('prevalence must be numeric')

    if (cstatistic < 0 | cstatistic > 1) stop('cstatistic must be between 0 and 1')
    if (is.numeric(cstatistic) == F) stop('cstatistic must be numeric')

    if (is.numeric(cslope) == F) stop('cslope must be numeric')

    if (prevalence <= 0 | prevalence >= 1) stop('prevalence must be between 0 and 1')

    # parameters not needed

  }

}

