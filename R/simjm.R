#' Simulate data for a univariate or multivariate joint model
#'
#' Returns a data frame containing data simulated from a joint model for
#' longitudinal and time-to-event data.
#'
#' @export
#'
#' @param n Number of subjects
#' @param M Number of longitudinal markers
#' @param random_trajectory The desired type of trajectory in the random effects
#'   part of the longitudinal model
#' @param assoc The desired type of association structure for the joint model.
#'   Only one association structure type can be specified. The options are:
#'   \code{"etavalue"}, \code{"etaslope"}, \code{"shared_b(1)"},
#'   \code{"shared_coef(1)"}, \code{"shared_b(2)"}, \code{"shared_coef(2)"}.
#' @param basehaz The desired baseline hazard for the event submodel. Can be
#'   \code{"weibull"}.
#' @param betaLong_intercept True intercept in the longitudinal submodel. Can
#'   be a scalar or a vector of length M.
#' @param betaLong_binary True coefficient for the binary covariate in the
#'   longitudinal submodel. Can be a scalar or a vector of length M.
#' @param betaLong_continuous True coefficient for the continuous covariate in the
#'   longitudinal submodel. Can be a scalar or a vector of length M.
#' @param betaLong_slope True coefficient for the fixed effect slope in the
#'   longitudinal submodel. Can be a scalar or a vector of length M.
#' @param betaLong_aux True parameter value for the auxiliary parameter in the
#'   longitudinal submodel (sigma for Gaussian models, number of trials for
#'   binomial models, size for negative binomial models, shape for Gamma models,
#'   lambda for inverse Gaussian models). Can be a scalar or a vector of length M.
#' @param betaEvent_intercept True intercept term (log hazard scale) in the
#'   event submodel. Only required for Weibull models.
#' @param betaEvent_binary True coefficient (log hazard ratio) for the binary
#'   covariate in the event submodel.
#' @param betaEvent_continuous True coefficient (log hazard ratio) for the
#'   continuous covariate in the event submodel.
#' @param betaEvent_assoc True association parameter (log hazard ratio) in the
#'   event submodel. Can be a scalar or a vector of length M.
#' @param betaEvent_aux True parameter value(s) for the auxiliary parameter(s) in
#'   the event submodel (shape parameter for Weibull models). Can be a scalar
#'   or a vector, depending on how many auxiliary parameters are required for
#'   specifying the baseline hazard.
#' @param b_sd Vector of standard deviations for the random effects.
#' @param b_rho Correlation between the random effects.
#' @param prob_Z1 Probability for the binary covariate included in each of the
#'   submodels.
#' @param mean_Z2 Mean of the (normally distributed) continuous covariate included
#'   in each of the submodels.
#' @param sd_Z2 Standard deviation of the (normally distributed) continuous
#'   covariate included in each of the submodels.
#' @param max_yobs The maximum allowed number of longitudinal measurements. The
#'   actual number of observed measurements will depend on the individuals event time.
#' @param max_fuptime The maximum follow up time in whatever the desired time
#'   units are. This time will also be used as the censoring time (i.e. for subjects
#'   who have a simulated survival time that is after \code{max_fuptime})
#' @param family A family for the the longitudinal submodel, or for a multivariate
#'   joint model this can be a list of families. See \code{\link[stats]{glm}} and
#'   \code{\link[stats]{family}}.
#' @param seed An optional \code{\link[=set.seed]{seed}}.
#' @param interval The interval over which to search for the
#'   \code{\link[stats]{uniroot}} corresponding to each simulated event time.
#'
#' @return A data frame in long format. The data frame also has a number of
#'   attributes that record the values specified for the arguments in the
#'   \code{simjm} call (for example the "true" parameter values, the number
#'   of individuals, the number of longitudinal markers, and so on).
#'
#' @examples
#' #####
#' # By default simjm will simulate data for three longitudinal biomarkers.
#' # However, if we only want one or two longitudinal biomarkers then we need to
#' # adjust the "true" parameters accordingly so that we still get realistic
#' # survival times. The easiest option is to change the "true" association
#' # parameter value. Some suggested values that lead to realistic survival
#' # times are given below.
#'
#' # For one longitudinal marker:
#' simdat1 <- simjm(M = 1, betaEvent_assoc = 0.03)
#'
#' # For two longitudinal markers:
#' simdat2 <- simjm(M = 2, betaEvent_assoc = 0.015)
#'
#' # For three longitudinal markers, we can just use the defaults:
#' simdat3 <- simjm()
#'
#' # Simulate three markers for just 100 individuals and then return
#' # the true parameter values (which are stored as an attribute):
#' simdat4 <- simjm(n = 100)
#' attr(simdat4, "params")
#'
#' # Simulate three longitudinal markers, using "etaslope"
#' # association structure
#' simdat5 <- simjm(assoc = "etaslope", betaEvent_assoc = 0.3)
#'
#' # For one longitudinal marker, with a bernoulli outcome.
#' # (Note that 'betaLong_aux' specifies the number of trials
#' # for the binomial outcome).
#' simdat6 <- simjm(M = 1,
#'                  betaLong_intercept = 1.3,
#'                  betaLong_binary = -0.6,
#'                  betaLong_continuous = -0.03,
#'                  betaLong_slope = 0.05,
#'                  b_sd = c(1, 0.05),
#'                  betaEvent_intercept = -9,
#'                  betaEvent_assoc = 0.05,
#'                  family = binomial(),
#'                  betaLong_aux = 1)
#'
simjm <- function(n = 200, M = 1,
                  random_trajectory = "linear",
                  assoc = "etavalue",
                  basehaz = c("weibull"),
                  betaLong_intercept = 0,
                  betaLong_binary = 1,
                  betaLong_continuous = 1,
                  betaLong_slope = 1,
                  betaLong_aux = 1,
                  betaEvent_intercept = -4,
                  betaEvent_binary = 1,
                  betaEvent_continuous = 0,
                  betaEvent_assoc = 0.2,
                  betaEvent_aux = 1.2,
                  b_sd = c(2,1), b_rho = -0.2,
                  prob_Z1 = 0.5,
                  mean_Z2 = 0, sd_Z2 = 1,
                  max_yobs = 10,
                  max_fuptime = 5,
                  family = gaussian,
                  seed = sample.int(.Machine$integer.max, 1),
                  interval = c(0, 200))
{

  #----- Preliminaries

  set.seed(seed)

  basehaz <- match.arg(basehaz)

  ok_assocs <- c("etavalue", "etaslope", "etaauc", "muvalue",
                 "shared_b(1)", "shared_coef(1)",
                 "shared_b(2)", "shared_coef(2)")
  assoc <- maybe_broadcast(assoc, M)
  if (!all(assoc %in% ok_assocs))
    stop("'assoc' must be one of: ", paste(ok_assocs, collapse = ", "))

  ok_trajs <- c("linear", "none", "poly")
  random_trajectory <- maybe_broadcast(random_trajectory, M)
  if (!all(random_trajectory %in% ok_trajs))
    stop("'assoc' must be one of: ", paste(ok_assocs, collapse = ", "))

  for (m in 1:M) {
    shared_slope <- (assoc[m] %in% c("shared_b(2)", "shared_coef(2)"))
    if (shared_slope && (!random_trajectory[m] == "linear"))
      stop("Can only use 'shared_b(2)' or 'shared_coef(2)' when ",
           "'random_trajectory' is linear.")
  }

  if (!is(family, "list"))
    family <- list(family)
  family <- maybe_broadcast(family, M)
  family <- lapply(family, validate_family)
  ok_families <- c("gaussian", "binomial")
  lapply(family, function(x) if (!x$family %in% ok_families)
    stop("'family' must be one of: ", paste(ok_families, collapse = ", ")))

  #----- Parameters

  # Broadcast user-specified parameters
  betaLong_intercept  <- maybe_broadcast(betaLong_intercept,  M)
  betaLong_binary     <- maybe_broadcast(betaLong_binary,     M)
  betaLong_continuous <- maybe_broadcast(betaLong_continuous, M)
  betaLong_slope      <- maybe_broadcast(betaLong_slope,      M)
  betaLong_aux        <- maybe_broadcast(betaLong_aux,        M)
  betaEvent_assoc     <- maybe_broadcast(betaEvent_assoc,     M)

  # Generate subject-specific random effects
  b_dim <- sapply(random_trajectory, function(x)
    switch(x,
           none   = 1L, # random intercepts model
           linear = 2L, # random slopes model
           poly   = 3L) # random poly (degree = 2) model
  )
  b_dim_total <- sum(b_dim) # total num of random effects
  if (length(unique(b_dim)) == 1L) { # same num of reffs in each submodel
    if (length(b_sd) == b_dim[1]) {
      b_sd <- rep(b_sd, times = M)
    }
  }
  if (!length(b_sd) == b_dim_total)
    stop("b_sd appears to be the wrong length.")
  corr_mat <- matrix(rep(b_rho, b_dim_total ^ 2), ncol = b_dim_total)
  diag(corr_mat) <- 1
  dd <- MASS::mvrnorm(n = n, mu = rep(0, b_dim_total), Sigma = corr_mat)
  b <- sapply(1:length(b_sd), function(x) b_sd[x] * dd[,x])
  colnames(b) <- paste0("b", 1:b_dim_total)

  # Construct data frame of parameters
  pars <- data.frame(
    id = 1:n,
    betaEvent_intercept  = rep(betaEvent_intercept,  n),
    betaEvent_binary     = rep(betaEvent_binary,     n),
    betaEvent_continuous = rep(betaEvent_continuous, n),
    betaEvent_aux        = rep(betaEvent_aux,        n)
  )
  for (m in 1:M) {
    # fixed effect parameters
    pars[[paste0("betaEvent_assoc", m)]] <- rep(betaEvent_assoc[m], n)
    pars[[paste0("betaLong_intercept", m)]] <- rep(betaLong_intercept[m], n)
    pars[[paste0("betaLong_binary", m)]] <- rep(betaLong_binary[m], n)
    pars[[paste0("betaLong_continuous", m)]] <- rep(betaLong_continuous[m], n)
    if (random_trajectory[m] %in% c("none", "linear")) { # fixed effect slope
      pars[[paste0("betaLong_slope", m)]] <- rep(betaLong_slope[m], n)
    } else if (random_trajectory[m] == "poly") { # fixed effect quadratic
      pars[[paste0("betaLong_poly1", m)]] <- rep(betaLong_poly1[m], n)
      pars[[paste0("betaLong_poly2", m)]] <- rep(betaLong_poly2[m], n)
    }
    # add on subject-specific intercept
    shift <- if (m == 1) 0 else sum(b_dim[1:(m - 1)])
    b_idx <- shift + 1
    pars[[paste0("betaLong_intercept",  m)]] <-
      pars[[paste0("betaLong_intercept",  m)]] + b[, b_idx]
    # add on subject-specific linear slope
    if (random_trajectory[m] == "linear") {
      b_idx <- shift + 2
      pars[[paste0("betaLong_slope",  m)]] <-
        pars[[paste0("betaLong_slope",  m)]] + b[, b_idx]
    }
    # add on subject-specific quadratic terms
    if (random_trajectory[m] == "poly") {
      b_idx <- shift + 2 # index of first quadratic term
      pars[[paste0("betaLong_poly1",  m)]] <-
        pars[[paste0("betaLong_poly1",  m)]] + b[, b_idx]
      b_idx <- shift + 3 # index of second quadratic term
      pars[[paste0("betaLong_poly2",  m)]] <-
        pars[[paste0("betaLong_poly2",  m)]] + b[, b_idx]
    }
  }

  #----- Data

  # Generate baseline covariates - binary
  Z1 <- rbinom(n, 1, prob_Z1)      # covariate value for each subject

  # Generate baseline covariates - continuous
  Z2 <- rnorm(n, mean_Z2, sd_Z2)   # covariate value for each subject

  # Construct data frame of baseline covariates
  covs <- data.frame(id = 1:n, Z1, Z2)

  # Generate survival times
  ss <- simsurv::simsurv(hazfn = jm_hazfn, x = covs, pars = pars,
                         maxt = max_fuptime, basehaz = basehaz, M = M,
                         random_trajectory = random_trajectory,
                         assoc = assoc, family = family, interval = interval)

  # Randomly generate observation times between 0 and max_fuptime
  tij <- runif(n * max_yobs, 0, max_fuptime)

  # Construct data frame of longitudinal data
  dat.id <- data.frame(pars, covs, ss)              # single row per subject
  dat <- dat.id[rep(row.names(dat.id), max_yobs), ] # multiple row per subject
  dat <- data.frame(dat, tij)                       # merge on measurement times and errors
  dat <- dat[order(dat$id, dat$tij), ]              # sort on ID and time

  # Calculate longitudinal outcomes
  for (m in 1:M) {
    dat[[paste0("Xij_", m)]] <-
      dat[[paste0("betaLong_intercept", m)]] +
      dat[[paste0("betaLong_binary", m)]] * dat[["Z1"]] +
      dat[[paste0("betaLong_continuous", m)]] * dat[["Z2"]]
    if (random_trajectory[m] %in% c("none", "linear")) {
      # if no random slope, then still use fixed linear slope
      dat[[paste0("Xij_", m)]] <- dat[[paste0("Xij_", m)]] +
        dat[[paste0("betaLong_slope", m)]] * dat[["tij"]]
    } else if (random_trajectory[m] == "poly") {
      dat[[paste0("Xij_", m)]] <- dat[[paste0("Xij_", m)]] +
        dat[[paste0("betaLong_poly1",  m)]] * dat[["tij"]] +
        dat[[paste0("betaLong_poly2",  m)]] * dat[["tij"]] * dat[["tij"]]
    }
    fam <- family[[m]]$family
    invlink <- family[[m]]$linkinv
    mu <- invlink(dat[[paste0("Xij_", m)]])
    if (fam == "gaussian") {
      sigma <- betaLong_aux[m]
      dat[[paste0("Yij_", m)]] <- rnorm(length(mu), mu, sigma)
    } else if (fam == "binomial") {
      trials <- betaLong_aux[m]
      dat[[paste0("Yij_", m)]] <- rbinom(length(mu), trials, mu)
    } else if (fam == "poisson") {
      dat[[paste0("Yij_", m)]] <- rpois(length(mu), mu)
    } else if (fam == "neg_binomial_2") {
      size <- betaLong_aux[m]
      dat[[paste0("Yij_", m)]] <- rnbinom(length(mu), size, mu)
    } else if (fam == "Gamma") {
      shape <- betaLong_aux[m]
      dat[[paste0("Yij_", m)]] <- rgamma(length(mu), shape, rate = shape / mu)
    } else if (fam == "inverse.gaussian") {
      lambda <- betaLong_aux[m]
      dat[[paste0("Yij_", m)]] <- .rinvGauss(length(mu), mu, lambda)
    }
  }

  # Final dataset
  ret <- dat[dat$tij <= dat$eventtime, ]  # only keep rows before event time
  sel <- grep("^id$|^Z|^tij|^Yij|eventtime|status", colnames(ret))
  ret <- ret[, sel, drop = FALSE]
  rownames(ret) <- NULL

  # Store 'true' parameter values
  long_params <- nlist(
    betaLong_intercept,
    betaLong_binary,
    betaLong_continuous,
    betaLong_aux
  )
  if (any(random_trajectory %in% c("none", "linear")))
    long_params$betaLong_slope <- betaLong_slope
  if (any(random_trajectory == "poly")) {
    long_params$betaLong_poly1 <- betaLong_poly1
    long_params$betaLong_poly2 <- betaLong_poly2
  }
  event_params <- nlist(
    betaEvent_intercept,
    betaEvent_binary,
    betaEvent_continuous,
    betaEvent_assoc,
    betaEvent_aux
  )
  re_params <- nlist(
    b_sd,
    b_corr = corr_mat
  )

  # Return object
  structure(ret, params = c(long_params, event_params, re_params), n = n, M = M,
            max_yobs = max_yobs, max_fuptime = max_fuptime, assoc = assoc,
            family = family, random_trajectory = random_trajectory, seed = seed)
}


#--------------------- internal

# Return the hazard at time t, for a shared parameter joint model
#
# @param t The time variable
# @param x A named vector of covariate values
# @param pars A named vector of parameter values
# @param basehaz The type of baseline hazard
# @param M The number of longitudinal submodels
# @param assoc A vector of character strings, indicating the association
#   structure to use for each longitudinal submodel
# @param A list of family objects, providing the family (and link function)
#   for each longitudinal submodel
# @return A scalar
jm_hazfn <- function(t, x, pars, basehaz = "weibull", M = 1,
                     random_trajectory = "linear",
                     assoc = "etavalue", family = list(gaussian())) {

  if (t == 0)
    return(0) # boundary condition

  # Baseline hazard
  if (basehaz == "weibull") {
    h0 <- pars[["betaEvent_aux"]] * (t ^ (pars[["betaEvent_aux"]] - 1))
  }

  # Time-invariant part of event submodel eta
  etaevent <-
    pars[["betaEvent_intercept"]] +
    pars[["betaEvent_binary"]] * x[["Z1"]] +
    pars[["betaEvent_continuous"]] * x[["Z2"]]

  # Association structure
  for (m in 1:M) {

    # eta value
    if (assoc[m] == "etavalue") {
      etabaseline_m <- etavalue_m <-
        pars[[paste0("betaLong_intercept", m)]] +
        pars[[paste0("betaLong_binary", m)]] * x[["Z1"]] +
        pars[[paste0("betaLong_continuous", m)]] * x[["Z2"]]
      if (random_trajectory[m] %in% c("none", "linear")) {
        # if no random slope, then still use fixed linear slope
        etavalue_m <- etavalue_m +
          pars[[paste0("betaLong_slope", m)]] * t
      } else if (random_trajectory[m] == "poly") {
        etavalue_m <- etavalue_m +
          pars[[paste0("betaLong_poly1",  m)]] * t +
          pars[[paste0("betaLong_poly2",  m)]] * (t * t)
      }
      etaevent <- etaevent + pars[[paste0("betaEvent_assoc", m)]] * etavalue_m
    }

    # eta slope
    if (assoc[m] == "etaslope") {
      if (random_trajectory[m] %in% c("none", "linear")) {
        # if no random slope, then still use fixed linear slope
        etaslope_m <-
          pars[[paste0("betaLong_slope", m)]]
      } else if (random_trajectory[m] == "poly") {
        etaslope_m <-
          pars[[paste0("betaLong_poly1",  m)]] +
          pars[[paste0("betaLong_poly2",  m)]] * (2 * t)
      }
      etaevent <- etaevent + pars[[paste0("betaEvent_assoc", m)]] * etaslope_m
    }

    # eta auc
    if (assoc[m] == "etaauc") {
      if (random_trajectory[m] %in% c("none", "linear")) {
        etaauc_m <-
          etabaseline_m * t +
          0.5 * (etavalue_m - etabaseline_m) * t
      } else if (random_trajectory[m] == "poly") {
        stop("Calculation of hazard function for 'etaauc' with polynomial ",
             "trajectory not yet implemented.", call. = FALSE)
      }
      etaevent <- etaevent + pars[[paste0("betaEvent_assoc", m)]] * etaauc_m
    }

    # mu value
    if (assoc[m] == "etaauc") {
      invlink <- family[[m]]$invlink
      muvalue_m <- invlink(etavalue_m)
      etaevent <- etaevent + pars[[paste0("betaEvent_assoc", m)]] * muvalue_m
    }

  }
  # Calculate and return hazard
  h <- h0 * exp(etaevent)
  return(h)
}
