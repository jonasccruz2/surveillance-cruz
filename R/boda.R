boda <- function(sts, control=list(range=NULL, X=NULL, trend=FALSE, season=FALSE,
                                   prior=c('iid','rw1','rw2'), alpha=0.05,
                                   alphaExclusion=1e-4, mc.munu=100,
                                   mc.y=10, verbose=FALSE,
                                   samplingMethod=c('joint','marginals'),
                                   quantileMethod=c("MC","MM"),
                                   offset=FALSE)) {
  
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("The boda function requires the INLA package to be installed.")
  }
  
  if (ncol(sts) > 1) {
    stop("boda currently only handles univariate sts objects.")
  }
  
  if (is.null(control[["quantileMethod", exact=TRUE]])) {
    control$quantileMethod <- "MC"
  } else {
    control$quantileMethod <- match.arg(control$quantileMethod,
                                        c("MC","MM"))
  }
  
  observed_values <- as.vector(observed(sts))
  time_index <- 1:length(observed_values)

  # offset = TRUE: the population of the sts object (e.g. live births) is
  # used as the exposure E in INLA, i.e. mean = population * exp(eta).
  # offset = FALSE (default): original behaviour, E = mean of the counts.
  if (is.null(control[["offset", exact=TRUE]])) {
    control$offset <- FALSE
  }

  if (!is.logical(control$offset) || length(control$offset) != 1 || is.na(control$offset)) {
    stop("offset has to be TRUE or FALSE.")
  }

  population_values <- as.vector(sts@populationFrac)

  if (control$offset) {
    if (length(population_values) != length(observed_values) ||
        any(is.na(population_values)) || any(population_values < 0)) {
      stop("offset = TRUE needs a non-negative population for every time point of the sts object.")
    }
    if (length(unique(population_values)) == 1) {
      warning("offset = TRUE but the population is constant; was 'population' given to sts()?")
    }
  }

  if (is.null(control[["range", exact=TRUE]])) {
    warning("No range given. Range is defined as time from second period until end of time series.")
    control$range <- (sts@freq + 1):length(observed_values)
  }
  
  if (!all(control$range %in% time_index)) {
    stop("Evaluation period 'range' has to be vector of time series indices.")
  }
  
  control$range <- sort(control$range)
  
  if (is.null(control[["verbose", exact=TRUE]])) {
    control$verbose <- FALSE
  }
  
  if (is.null(control[["trend", exact=TRUE]])) {
    control$trend <- FALSE
  }
  
  if (is.null(control[["season", exact=TRUE]])) {
    control$season <- FALSE
  }
  
  if (!is.logical(control$trend) || !is.logical(control$season)) {
    stop("trend and season are logical parameters.")
  }
  
  prior <- match.arg(control$prior, c("iid","rw1","rw2"))
  
  if (is.vector(control$X)) {
    control$X <- as.matrix(control$X, ncol=1)
  }
  
  samplingMethod <- match.arg(control$samplingMethod,
                              c("joint","marginals"))
  
  if (is.null(control[["alpha", exact=TRUE]])) {
    control$alpha <- 0.05
  }
  
  if (control$alpha <= 0 || control$alpha >= 1) {
    stop("The significance level 'alpha' has to be between 0 and 1.")
  }
  
  if (is.null(control[["alphaExclusion", exact=TRUE]])) {
    control$alphaExclusion <- 1e-4
  }
  
  if (control$alphaExclusion <= 0 || control$alphaExclusion >= 1) {
    stop("The exclusion level 'alphaExclusion' has to be between 0 and 1.")
  }
  
  if (is.null(control[["mc.munu", exact=TRUE]])) {
    control$mc.munu <- 100
  }
  
  if (is.null(control[["mc.y", exact=TRUE]])) {
    control$mc.y <- 10
  }
  
  if (!control$mc.munu > 0 ||
      control$mc.munu != round(control$mc.munu,0) ||
      !control$mc.y > 0 ||
      control$mc.y != round(control$mc.y,0)) {
    stop("Number of Monte Carlo trials has to be an integer larger than zero")
  }
  
  model_formula <- paste(
    "observed_values ~ f(time_index, model='",
    prior,
    "', cyclic=FALSE)",
    sep=""
  )
  
  data_frame <- data.frame(
    observed_values=observed_values,
    time_index=time_index
  )

  if (control$offset) {
    data_frame$exposure <- population_values
  }

  if (control$trend) {
    model_formula <- paste(
      model_formula,
      "+ f(time_trend, model='linear')",
      sep=""
    )
    
    data_frame <- data.frame(
      data_frame,
      time_trend=time_index
    )
  }
  
  if (control$season) {
    model_formula <- paste(
      model_formula,
      "+ f(time_season, model='seasonal', season.length=",
      sts@freq,
      ")",
      sep=""
    )
    
    data_frame <- data.frame(
      data_frame,
      time_season=time_index
    )
  }
  
  covariable_formula <- NULL
  
  if (!is.null(control$X)) {
    
    if (nrow(control$X) != length(observed_values)) {
      stop("Argument for covariates 'X' has to have the same length like the time series")
    }
    
    for (i in 1:ncol(control$X)) {
      covariable_formula <- paste(
        covariable_formula,
        "+",
        colnames(control$X)[i]
      )
    }
    
    model_formula <- paste(
      model_formula,
      covariable_formula,
      sep=""
    )
    
    data_frame <- data.frame(
      data_frame,
      control$X
    )
  }
  
  model_formula <- as.formula(model_formula)
  
  use_progress_bar <- length(control$range) > 1
  
  if (use_progress_bar) {
    progress_bar <- txtProgressBar(
      min=min(control$range),
      max=max(control$range),
      initial=0,
      style=if (interactive()) 3 else 1
    )
  }
  
  upperbound_values <- rep(NA, length(observed_values))
  lowerbound_values <- rep(NA, length(observed_values))
  predicted_values <- rep(NA, length(observed_values))
  
  observed_clean <- observed_values
  
  for (i in control$range) {
    
    exclusion_data <- data_frame[1:i, ]
    exclusion_data$observed_values <- observed_clean[1:i]
    exclusion_data$observed_values[i] <- NA
    
    exclusion_fit <- bodaFit(
      dat=exclusion_data,
      samplingMethod=samplingMethod,
      modelformula=model_formula,
      alpha=control$alphaExclusion,
      mc.munu=control$mc.munu,
      mc.y=control$mc.y,
      quantileMethod=control$quantileMethod,
      verbose=control$verbose
    )
    
    if (!is.na(observed_values[i]) &&
        !is.na(exclusion_fit$upperbound) &&
        observed_values[i] > exclusion_fit$upperbound) {
      
      observed_clean[i] <- NA
    }
    
    final_data <- data_frame[1:i, ]
    final_data$observed_values <- observed_clean[1:i]
    final_data$observed_values[i] <- NA
    
    fit_result <- bodaFit(
      dat=final_data,
      samplingMethod=samplingMethod,
      modelformula=model_formula,
      alpha=control$alpha,
      mc.munu=control$mc.munu,
      mc.y=control$mc.y,
      quantileMethod=control$quantileMethod,
      verbose=control$verbose
    )
    
    upperbound_values[i] <- fit_result$upperbound
    lowerbound_values[i] <- fit_result$lowerbound
    predicted_values[i] <- fit_result$predicted
    
    if (use_progress_bar) {
      setTxtProgressBar(progress_bar, i)
    }
  }
  
  if (use_progress_bar) {
    close(progress_bar)
  }
  
  sts@alarm[,1] <- observed_values > upperbound_values
  sts@upperbound[,1] <- upperbound_values
  
  control$predicted <- predicted_values
  control$lowerbound <- lowerbound_values
  control$name <- paste("boda(prior=", prior, ")", sep="")
  
  sts@control <- control
  
  return(sts[control$range, ])
}

bodaFit <- function(dat, modelformula, alpha, mc.munu, mc.y,
                    samplingMethod, quantileMethod, verbose = FALSE) {
  
  T1 <- nrow(dat)
  
  environment(modelformula) <- environment()
  
  if (!is.null(dat$exposure)) {
    # offset: population of each time point as exposure E
    exposure_value <- dat$exposure

    # no population at the predicted time point: nothing expected
    if (exposure_value[T1] <= 0) {
      return(list(
        upperbound=0,
        lowerbound=0,
        predicted=0
      ))
    }

    # time points without population carry no information for the fit
    no_population <- exposure_value <= 0
    dat$observed_values[no_population] <- NA

    # Rescale E so that its mean equals the mean count. With a large
    # population (e.g. ~1e5 live births per period) the linear predictor sits
    # far from 0 (around log(rate) ~ -8) and INLA becomes numerically
    # unstable: some fits fail and the posterior samples explode (absurdly
    # high or low thresholds). Multiplying E by a constant does not change
    # the model -- the constant is absorbed by the intercept -- and keeps the
    # linear predictor near 0, as with the original scaling by the mean count.
    with_data <- !is.na(dat$observed_values) & !no_population
    if (!any(with_data)) with_data <- !no_population
    target <- mean(dat$observed_values[with_data], na.rm=TRUE)
    if (!is.finite(target) || target <= 0) target <- 1
    exposure_value <- exposure_value * target / mean(exposure_value[with_data])
    exposure_value[no_population] <- target
  } else {
    # original behaviour: constant scaling by the mean of the counts
    exposure_value <- mean(dat$observed_values, na.rm=TRUE)
  }

  # exposure of the predicted time point (used for the predictive draws)
  exposure_T1 <- if (length(exposure_value) > 1) exposure_value[T1] else exposure_value

  model <- INLA::inla(
    modelformula,
    data=dat,
    family='nbinomial',
    E=exposure_value,
    verbose=verbose,
    control.predictor=list(compute=TRUE, link=1),
    control.compute=c(
      list(cpo=FALSE, config=TRUE),
      if (packageVersion("INLA") >= "21.07.10")
        list(return.marginals.predictor=TRUE)
    ),
    control.inla=list(
      int.strategy="grid",
      dz=1,
      diff.logdens=10
    )
  )
  
  if (is.null(model)) {
    warning("NULL result from INLA at t = ", T1)
    return(list(
      upperbound=NA_real_,
      lowerbound=NA_real_,
      predicted=NA_real_
    ))
  }
  
  if (samplingMethod == "marginals") {
    
    marginal_distribution <- model$marginals.fitted.values[[T1]]
    
    mean_values <- try(
      INLA::inla.rmarginal(
        n=mc.munu,
        marginal_distribution
      ),
      silent=TRUE
    )
    
    if (inherits(mean_values, "try-error")) {
      warning("degenerate marginal posterior at t = ", T1)
      
      return(list(
        upperbound=NA_real_,
        lowerbound=NA_real_,
        predicted=NA_real_
      ))
    }
    
    hyperparameter_distribution <- model$internal.marginals.hyperpar[[1]]
    
    theta_values <- exp(
      INLA::inla.rmarginal(
        n=mc.munu,
        hyperparameter_distribution
      )
    )
  }
  
  if (samplingMethod == "joint") {
    
    posterior_samples <- INLA::inla.posterior.sample(
      n=mc.munu,
      result=model,
      intern=TRUE,
      skew.corr=FALSE
    )
    
    theta_values <- exp(
      t(
        sapply(
          posterior_samples,
          function(x) x$hyperpar[[1]]
        )
      )
    )
    
    mean_values <- exp(
      t(
        sapply(
          posterior_samples,
          function(x) x$latent[[T1]]
        )
      )
    )
  }
  
  valid_values <- mean_values >= 0 & theta_values > 0
  
  if (any(!valid_values)) {
    
    warning("degenerate posterior sampling at t = ", T1)
    
    return(list(
      upperbound=NA_real_,
      lowerbound=NA_real_,
      predicted=NA_real_
    ))
  }
  
  if (quantileMethod == "MC") {
    
    predicted_distribution <- unlist(
      mapply(
        rnbinom,
        size=theta_values,
        mu=exposure_T1 * mean_values,
        MoreArgs=list(n=mc.y),
        SIMPLIFY=FALSE
      )
    )
    
    upperbound_value <- quantile(
      predicted_distribution,
      probs=(1 - alpha),
      type=3,
      na.rm=TRUE
    )
    
    lowerbound_value <- quantile(
      predicted_distribution,
      probs=alpha,
      type=3,
      na.rm=TRUE
    )
    
    predicted_value <- mean(
      predicted_distribution,
      na.rm=TRUE
    )
  }
  
  if (quantileMethod == "MM") {
    
    minimum_bracket <- qnbinom(
      p=(1 - alpha),
      mu=exposure_T1 * min(mean_values),
      size=max(theta_values)
    )
    
    maximum_bracket <- qnbinom(
      p=(1 - alpha),
      mu=exposure_T1 * max(mean_values),
      size=min(theta_values)
    )
    
    upperbound_value <- qmix(
      p=(1 - alpha),
      mu=exposure_T1 * mean_values,
      size=theta_values,
      bracket=c(minimum_bracket, maximum_bracket)
    )
    
    lowerbound_value <- qmix(
      p=alpha,
      mu=exposure_T1 * mean_values,
      size=theta_values,
      bracket=c(0, maximum_bracket)
    )
    
    predicted_value <- mean(
      exposure_T1 * mean_values,
      na.rm=TRUE
    )
  }
  
  return(list(
    upperbound=as.numeric(upperbound_value),
    lowerbound=as.numeric(lowerbound_value),
    predicted=as.numeric(predicted_value)
  ))
}