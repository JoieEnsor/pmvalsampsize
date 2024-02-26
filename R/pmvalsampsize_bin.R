# library(pROC)

pmvalsampsize_bin <- function(prevalence,cstatistic,oe,oeciwidth,cslope,csciwidth,cstatciwidth,simobs,lpnormal,lpbeta,lpcstat,tolerance,increment,oeseincrement,seed,graph,trace,sensitivity,specificity,threshold,nbciwidth,nbseincrement) {

#############################

  set.seed(seed)

#############################
  # criteria 1 - o/e

  # prevalence <-0.018
  # oeciwidth <- 1
  # oeseincrement <- 0.0001
  # oe <- 1

  width_oe <- 0
  se_oe <- 0
  while (width_oe < oeciwidth) {
     se_oe = se_oe + oeseincrement
     ub_oe = exp(log(oe) + (1.96*se_oe))
     lb_oe = exp(log(oe) - (1.96*se_oe))
     width_oe = ub_oe - lb_oe
  }

  n1 <- ceiling((1-prevalence)/(prevalence*se_oe^2))

  E1 <- n1*prevalence

#############################
# criteria 2 - c-slope

  # cslope <- 1
  # csciwidth <- 0.2
  # lpnormal <- c(-5,2.5)
  # simobs <-1000000
  # is.vector(lpnormal)
  # # graph <- "TRUE"
  # mean <- 5.8
  # sd <-  5
  # skew <- -0.5

  if (!is.na(lpnormal[1])) {
    lpdist <- "normal"
    message("Normal LP distribution with parameters - mean=",lpnormal[1],", sd=",lpnormal[2],"\n")
    LP <- rnorm(n=simobs, mean=lpnormal[1], sd=lpnormal[2])
  }
  # else if (!is.na(lpskewednormal[1])) {
  #   lpdist <- "skewed normal"
  #   message("Skewed normal LP distribution with parameters - mean=",lpskewednormal[1],", sd=",lpskewednormal[2],", skew=",lpskewednormal[3],"\n")
  #
  #   mean <- lpskewednormal[1]
  #   sd <-  lpskewednormal[2]
  #   skew <- lpskewednormal[3]
  #
  #   if (skew<0) {
  #     sign <- -1
  #   }
  #   else if (skew>0) {
  #     sign <- 1
  #   }
  #
  #   delta <- ((abs(skew)^(2/3)*(pi/2))/(abs(skew)^(2/3))+((4-pi)/2)^(2/3))^0.5
  #   delta <- sign*delta
  #   # delta <- (((pi/2)*((abs(skew)^(2/3))))/((abs(skew)^(2/3))+((4-pi)/2)^(2/3)))^0.5
  #   shape <- delta/(1-delta^2)^0.5
  #   scale <- sd/((1-2*delta^2)/pi)^0.5
  #   location <- mean-scale*delta*(2/pi)^0.5
  #
  #   LP <- rsn(n=simobs, xi = location, omega = scale, alpha = shape)
  # }
  else if (!is.na(lpbeta[1])) {
    lpdist <- "beta"
    message("Beta P distribution with parameters - alpha=",lpbeta[1],", beta=",lpbeta[2],"\n")
    P <- rbeta(n=simobs, shape1=lpbeta[1], shape2=lpbeta[2])
    LP <- log(P/(1-P))
  }
  else if (!is.na(lpcstat[1])) {
    lpdist <- "cstat"
    m2 <- lpcstat[1]
    var <- 2*(qnorm(cstatistic)^2)

    # m2<-0.8
    # cstatistic<-0.8
    # var<- 2*(qnorm(cstatistic)^2)
    # simobs<-1000
    # prevalence <- 0.2
    # tolerance <- 0.0005
    # increment <- 0.1

    outcome <- rbinom(n=simobs, size=1, prob=prevalence)

    LP <- rnorm(n=simobs, mean=m2, sd=sqrt(var))
    LP2 <- rnorm(n=simobs, mean=m2+var, sd=sqrt(var))

    LP[outcome==1] <- LP2[outcome==1]

    # df <- data.frame(outcome, LP2, LP)
    # df$LP[df$outcome==1] <- df$LP2[df$outcome==1]

    P <- exp(LP)/(1+exp(LP))
    outcome_test <- rbinom(n=simobs, size=1, prob=P)
    diff <- abs(prevalence - mean(outcome_test))

    if (trace==TRUE) {
      if (diff>tolerance) {
          message("Proportion of observed outcome events does not match input prevalence","\n","Beginning iterative approach ...","\n")

        message("-------------- TRACE ON --------------","\n")

        n <- 1
        diff_new <- diff      # ask em about a name including a local
        m2 <- m2+increment

        outcome <- rbinom(n=simobs, size=1, prob=prevalence)
        LP <- rnorm(n=simobs, mean=m2, sd=sqrt(var))
        LP2 <- rnorm(n=simobs, mean=m2+var, sd=sqrt(var))

        LP[outcome==1] <- LP2[outcome==1]

        P <- exp(LP)/(1+exp(LP))
        outcome_test <- rbinom(n=simobs, size=1, prob=P)

        #n <- n+1

        diff_new <- abs(prevalence - mean(outcome_test))

        if (diff_new < diff) {
          while (diff_new > tolerance) {
            m2 <- m2+increment

            outcome <- rbinom(n=simobs, size=1, prob=prevalence)
            LP <- rnorm(n=simobs, mean=m2, sd=sqrt(var))
            LP2 <- rnorm(n=simobs, mean=m2+var, sd=sqrt(var))

            LP[outcome==1] <- LP2[outcome==1]

            P <- exp(LP)/(1+exp(LP))
            outcome_test <- rbinom(n=simobs, size=1, prob=P)

            diff_new <- abs(prevalence - mean(outcome_test))

            message("Proportion of outcome events under simulation =",mean(outcome_test),"\n","Target prevalence =",prevalence,"\n","Mean in non-event group = ",m2,"\n")
          }
        }
        else {
          while (diff_new > tolerance) {
            m2 <- m2-increment

            outcome <- rbinom(n=simobs, size=1, prob=prevalence)
            LP <- rnorm(n=simobs, mean=m2, sd=sqrt(var))
            LP2 <- rnorm(n=simobs, mean=m2+var, sd=sqrt(var))

            LP[outcome==1] <- LP2[outcome==1]

            P <- exp(LP)/(1+exp(LP))
            outcome_test <- rbinom(n=simobs, size=1, prob=P)

            diff_new <- abs(prevalence - mean(outcome_test))

            message("Proportion of outcome events under simulation =",mean(outcome_test),"\n","Target prevalence =",prevalence,"\n","Mean in non-event group = ",m2,"\n")
          }
        }

        message("-------------- TRACE OFF --------------","\n")
        message("Proportion of observed outcome events is within tolerance","\n","Proportion of outcome events under simulation =",mean(outcome_test),"\n","Target prevalence =",prevalence,"\n","Mean in non-event group = ",m2,"\n","\n")
      }
      else {
        message("Proportion of observed outcome events is within tolerance","\n","Proportion of outcome events under simulation =",mean(outcome_test),"\n","Target prevalence =",prevalence,"\n","Mean in non-event group = ",m2,"\n","\n")
      }

    }
    else {
      if (diff > tolerance) {
          message("Proportion of observed outcome events does not match input prevalence","\n","Beginning iterative approach ...","\n")

          n <- 1
          diff_new <- diff
          m2 <- m2+increment

          outcome <- rbinom(n=simobs, size=1, prob=prevalence)
          LP <- rnorm(n=simobs, mean=m2, sd=sqrt(var))
          LP2 <- rnorm(n=simobs, mean=m2+var, sd=sqrt(var))

          LP[outcome==1] <- LP2[outcome==1]

          P <- exp(LP)/(1+exp(LP))
          outcome_test <- rbinom(n=simobs, size=1, prob=P)

          #n <- n+1

          diff_new <- abs(prevalence - mean(outcome_test))

          if (diff_new < diff) {
            while (diff_new > tolerance) {
              m2 <- m2+increment

              outcome <- rbinom(n=simobs, size=1, prob=prevalence)
              LP <- rnorm(n=simobs, mean=m2, sd=sqrt(var))
              LP2 <- rnorm(n=simobs, mean=m2+var, sd=sqrt(var))

              LP[outcome==1] <- LP2[outcome==1]

              P <- exp(LP)/(1+exp(LP))
              outcome_test <- rbinom(n=simobs, size=1, prob=P)

              diff_new <- abs(prevalence - mean(outcome_test))


            }
          }
          else {
            while (diff_new > tolerance) {
              m2 <- m2-increment

              outcome <- rbinom(n=simobs, size=1, prob=prevalence)
              LP <- rnorm(n=simobs, mean=m2, sd=sqrt(var))
              LP2 <- rnorm(n=simobs, mean=m2+var, sd=sqrt(var))

              LP[outcome==1] <- LP2[outcome==1]

              P <- exp(LP)/(1+exp(LP))
              outcome_test <- rbinom(n=simobs, size=1, prob=P)

              diff_new <- abs(prevalence - mean(outcome_test))


            }
          }


          message("Proportion of observed outcome events is within tolerance","\n","Proportion of outcome events under simulation =",mean(outcome_test),"\n","Target prevalence =",prevalence,"\n","Mean in non-event group = ",m2,"\n","\n")
        }
        else {
          message("Proportion of observed outcome events is within tolerance","\n","Proportion of outcome events under simulation =",mean(outcome_test),"\n","Target prevalence =",prevalence,"\n","Mean in non-event group = ",m2,"\n","\n")
        }

      }

    # check C-statistic is correct
    c <- suppressMessages(roc(outcome_test~P,ci=TRUE))
    simulated_data_cstat <- c$auc

  }

  if (graph==TRUE) {
    hist(LP, freq=FALSE)
  }

  # input assumed parameters of calibration model (in future vr these could be options)
   beta0 <- 0
   beta1 <- 1

  # caclulate elements of I matrix
   Borenstein_00 <- exp(beta0 + (beta1*LP))/((1+ exp(beta0 + (beta1*LP)))^2)
	 Borenstein_01 <- LP*exp(beta0 + (beta1*LP))/((1+ exp(beta0 + (beta1*LP)))^2)
  Borenstein_11 <- LP*LP*exp(beta0 + (beta1*LP))/((1+ exp(beta0 + (beta1*LP)))^2)

	 I_00 <- mean(Borenstein_00)

	 I_01 <- mean(Borenstein_01)

	 I_11 <- mean(Borenstein_11)

	# calculate SE from input target CI width
	 se_cslope <- csciwidth/(2*1.96) # rounding

  # calculate sample size
  n2 <- ceiling(I_00/(se_cslope*se_cslope*((I_00*I_11)-(I_01*I_01))))

	E2 <- n2*prevalence

#############################
# criteria 3 - c-statistic

# cstatistic <-0.8
# prevalence <-0.018
# simobs <-1000000
# cstatciwidth <- 0.1

size <- 1:1000000

	suppressWarnings(se_cstatsq <- cstatistic*(1-cstatistic)*(1+(((size/2)-1)*((1-cstatistic)/(2-cstatistic))) +((((size/2)-1)*cstatistic)/(1+cstatistic)))/(size*size*prevalence*(1-prevalence)))

se_cstat <- se_cstatsq^0.5

CIwidth <- 2*1.96*se_cstat

cstat_df <- data.frame(size = size,se_cstatsq=se_cstatsq,se_cstat=se_cstat,CIwidth=CIwidth)

cstat_df2 <- cstat_df[cstat_df$CIwidth<=cstatciwidth,]

n3 <- cstat_df2$size[1]

suppressWarnings(se_cstat <- sqrt(cstatistic*(1-cstatistic)*(1+(((n3/2)-1)*((1-cstatistic)/(2-cstatistic))) +((((n3/2)-1)*cstatistic)/(1+cstatistic)))/(n3*n3*prevalence*(1-prevalence))))

#############################
# criteria 4 - net benefit

# prevalence <-0.018
# sensitivity <-0.5
# specificity <- 0.5
# threshold <- 0.08
# nbciwidth <- 0.2
# nbseincrement <- 0.0001


if (!is.na(sensitivity) & !is.na(specificity)) {

  nb <- (sensitivity*prevalence) - ((1-specificity)*(1-prevalence)*(threshold/(1-threshold)))
  standardised_nb <- nb/prevalence

  w <- ((1-prevalence)/prevalence)*(threshold/(1-threshold))

  # calc se for target ci width
  width_nb <- 0
  se_nb <- 0
  while (width_nb<nbciwidth) {
    se_nb <- se_nb+nbseincrement
    ub_nb <- (standardised_nb) + (1.96*se_nb)
    lb_nb <- (standardised_nb) - (1.96*se_nb)
    width_nb <- ub_nb-lb_nb
  }

  se_nb <- round(se_nb,digits = 3)

  # calculate sample size
  n4 <- ceiling((1/(se_nb^2))*(((sensitivity*(1-sensitivity))/prevalence)+(w*w*specificity*(1-specificity)/(1-prevalence))+(w*w*(1-specificity)*(1-specificity)/(prevalence*(1-prevalence)))))

  no_nb <- FALSE
}
else if (!is.na(threshold)) {
  nb_p <- exp(LP)/(1+exp(LP))
  nb_outcome <- rbinom(n=simobs, size=1, prob=nb_p)

  nb_df <- data.frame(nb_outcome,nb_p)
  nb_df$classification <- 1
  nb_df$classification[nb_df$nb_p<threshold] <- 0

	sensitivity = round(sum(nb_df$classification[nb_df$nb_outcome==1]) / sum(nb_df$nb_outcome==1), digits = 2)

	specificity = round((sum(nb_df$nb_outcome==0) - sum(nb_df$classification[nb_df$nb_outcome==0])) / sum(nb_df$nb_outcome==0), digits = 2)

	nb <- (sensitivity*prevalence) - ((1-specificity)*(1-prevalence)*(threshold/(1-threshold)))
  standardised_nb <- nb/prevalence

  w <- ((1-prevalence)/prevalence)*(threshold/(1-threshold))

  # calc se for target ci width
  width_nb <- 0
  se_nb <- 0
  while (width_nb<nbciwidth) {
    se_nb <- se_nb+nbseincrement
    ub_nb <- (standardised_nb) + (1.96*se_nb)
    lb_nb <- (standardised_nb) - (1.96*se_nb)
    width_nb <- ub_nb-lb_nb
  }

  se_nb <- round(se_nb,digits = 3)

  # calculate sample size
  n4 <- ceiling((1/(se_nb^2))*(((sensitivity*(1-sensitivity))/prevalence)+(w*w*specificity*(1-specificity)/(1-prevalence))+(w*w*(1-specificity)*(1-specificity)/(prevalence*(1-prevalence)))))

  no_nb <- FALSE
}
else {
  no_nb <- TRUE
}

################################### summary

if (no_nb==TRUE) {
  # minimum n
  nfinal <- max(n1,n2,n3)
  efinal <- nfinal*prevalence

  # create output table
  res <- matrix(NA,4,4)
  colnames(res) <- c("Samp_size","Perf","SE","CI width")
  rownames(res) <- c("Criteria 1 - O/E", "Criteria 2 - C-slope", "Criteria 3 - C statistic", "Final SS")
  res[,1] <- c(n1,n2,n3,1)
  res[,2] <- round(c(oe,cslope,cstatistic,1),digits = 3)
  res[,3] <- round(c(se_oe,se_cslope,se_cstat,1),digits = 3)
  res[,4] <- round(c(width_oe,csciwidth,cstatciwidth,1),digits = 3)

  res_sort <- order(res[,1],decreasing=FALSE)
  res[4,1] <- nfinal
  res[4,2] <- res[res_sort[4],2]
  res[4,3] <- res[res_sort[4],3]
  res[4,4] <- res[res_sort[4],4]

  # returns
  if (!is.na(lpcstat[1])) {
    out <- list(results_table = res,
                sample_size = nfinal,
                events = efinal,
                prevalence = prevalence,
                type = "binary",
                cstatistic = cstatistic,
                oe  = oe,
                se_oe  = se_oe,
                lb_oe  =lb_oe ,
                ub_oe  = ub_oe,
                width_oe  = width_oe,
                cslope  = cslope,
                se_cslope  = se_cslope,
                csciwidth = csciwidth ,
                se_cstat  = se_cstat,
                cstatciwidth = cstatciwidth ,
                simulated_data_cstat = simulated_data_cstat)
  }
  else {
    out <- list(results_table = res,
                sample_size = nfinal,
                events = efinal,
                prevalence = prevalence,
                type = "binary",
                cstatistic = cstatistic,
                oe  = oe,
                se_oe  = se_oe,
                lb_oe  =lb_oe ,
                ub_oe  = ub_oe,
                width_oe  = width_oe,
                cslope  = cslope,
                se_cslope  = se_cslope,
                csciwidth = csciwidth ,
                se_cstat  = se_cstat,
                cstatciwidth = cstatciwidth)
  }
}
else {
 # minimum n
  nfinal <- max(n1,n2,n3,n4)
  efinal <- nfinal*prevalence

  # create output table
  res <- matrix(NA,5,4)
  colnames(res) <- c("Samp_size","Perf","SE","CI width")
  rownames(res) <- c("Criteria 1 - O/E", "Criteria 2 - C-slope", "Criteria 3 - C statistic", "Criteria 4 - St Net Benefit", "Final SS")
  res[,1] <- c(n1,n2,n3,n4,1)
  res[,2] <- round(c(oe,cslope,cstatistic,standardised_nb,1),digits = 3)
  res[,3] <- round(c(se_oe,se_cslope,se_cstat,se_nb,1),digits = 3)
  res[,4] <- round(c(width_oe,csciwidth,cstatciwidth,nbciwidth,1),digits = 3)

  res_sort <- order(res[,1],decreasing=FALSE)
  res[5,1] <- nfinal
  res[5,2] <- res[res_sort[5],2]
  res[5,3] <- res[res_sort[5],3]
  res[5,4] <- res[res_sort[5],4]

  # returns
  if (!is.na(lpcstat[1])) {
  out <- list(results_table = res,
              sample_size = nfinal,
              events = efinal,
              prevalence = prevalence,
              type = "binary",
              cstatistic = cstatistic,
              oe  = oe,
              se_oe  = se_oe,
              lb_oe  =lb_oe ,
              ub_oe  = ub_oe,
              width_oe  = width_oe,
              cslope  = cslope,
              se_cslope  = se_cslope,
              csciwidth = csciwidth ,
              se_cstat  = se_cstat,
              cstatciwidth = cstatciwidth ,
              simulated_data_cstat = simulated_data_cstat,
              standardised_nb = standardised_nb,
              net_benefit = nb,
              se_st_nb = se_nb,
              nbciwidth = nbciwidth,
              sensitivity = sensitivity,
              specificity = specificity,
              threshold = threshold)
  }
  else {
    out <- list(results_table = res,
                sample_size = nfinal,
                events = efinal,
                prevalence = prevalence,
                type = "binary",
                cstatistic = cstatistic,
                oe  = oe,
                se_oe  = se_oe,
                lb_oe  =lb_oe ,
                ub_oe  = ub_oe,
                width_oe  = width_oe,
                cslope  = cslope,
                se_cslope  = se_cslope,
                csciwidth = csciwidth ,
                se_cstat  = se_cstat,
                cstatciwidth = cstatciwidth,
                standardised_nb = standardised_nb,
                net_benefit = nb,
                se_st_nb = se_nb,
                nbciwidth = nbciwidth,
                sensitivity = sensitivity,
                specificity = specificity,
                threshold = threshold)
  }
}

}
