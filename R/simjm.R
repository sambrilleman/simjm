#' Simulate data for a univariate or multivariate joint model
#'
#' Returns a list of data frames containing data simulated from a joint model for
#' longitudinal and time-to-event data.
#'
#' @export
#'
#' @param n Number of individuals
#' @param M Number of longitudinal markers
#' @param random_trajectory The desired type of trajectory in the random effects
#'   part of the longitudinal model.
#' @param assoc A character string, or a character vector of length M,
#'   specifying the desired type of association structure for
#'   linking each longitudinal outcome to the hazard of the event.
#'   Only one association structure type can be specified for each longitudinal
#'   outcome. The options are: \code{"etavalue"}, \code{"etaslope"},
#'   \code{"shared_b(1)"}, \code{"shared_coef(1)"}, \code{"shared_b(2)"},
#'   \code{"shared_coef(2)"}.
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
#' @param prob_Z1 Probability for the binary covariate included in each of the
#'   submodels.
#' @param mean_Z2 Mean of the (normally distributed) continuous covariate included
#'   in each of the submodels.
#' @param sd_Z2 Standard deviation of the (normally distributed) continuous
#'   covariate included in each of the submodels.
#' @param b_sd A list, with each element of the list containing the vector of
#'   standard deviations for the random effects relating to one grouping factor.
#' @param b_rho Correlation between the random effects. This is assumed to be
#'   the correlation between all random effects within a grouping factor, and
#'   the same correlation is used for each grouping factor.
#' @param max_yobs The maximum allowed number of longitudinal measurements. The
#'   actual number of observed measurements will depend on the individuals event time.
#' @param max_fuptime The maximum follow up time in whatever the desired time
#'   units are. This time will also be used as the censoring time (i.e. for
#'   individuals who have a simulated survival time that is after \code{max_fuptime}).
#' @param family A family for the the longitudinal submodel, or for a multivariate
#'   joint model this can be a list of families. See \code{\link[stats]{glm}} and
#'   \code{\link[stats]{family}}.
#' @param clust_control A named list providing the arguments that are
#'   necessary for simulating data where the first longitudinal submodel has
#'   lower level clustering within individuals. The named list should contain
#'   the following elements:
#'   \describe{
#'   \item{L}{Integer specifying the maximum number of lower level units
#'   clustered within an individual. The actual number of lower level units
#'   clustered within each individual is taken to be a random uniform (truncated
#'   to an integer) on the range \code{[1,L+1]}.}
#'   \item{assoc}{Character string specifying the method for combining information
#'   across the lower level units clustered within an individual when forming the
#'   association structure. Can be \code{"sum"} for specifying which indicates
#'   the association structure should be based on a summation across the lower
#'   level units clustered within an individual. Can be \code{"mean"} which
#'   indicates that the association structure should be based on the mean
#'   (i.e. average) taken across the lower level units clustered within an
#'   individual.}
#'   \item{u_sd}{Numeric vector providing the standard deviations of the random
#'   effect at the cluster level.}
#'   \item{random_trajectory}{The desired type of trajectory in the random effects
#'   part of the longitudinal model at the cluster level. Can be \code{"none"} to
#'   only include a random intercept at the cluster level, or otherwise
#'   \code{"linear"} or \code{"poly"}.}
#'   }
#' @param seed An optional \code{\link[=set.seed]{seed}}.
#' @param interval The interval over which to search for the
#'   \code{\link[stats]{uniroot}} corresponding to each simulated event time.
#'
#' @details The \code{simjm} function returns data simulated under a joint
#'   longitudinal and time-to-event model. The joint model can be univariate
#'   (i.e. one longitudinal outcome) or multivariate (i.e. more than one
#'   longitudinal outcome). If more than one longitudinal outcome is specified
#'   then the longitudinal outcomes are assumed to be correlated via a joint
#'   multivariate normal distribution for the individual-level random effects.
#'   Each longitudinal outcome may have a different family and link function.
#'   The time-to-event model is assumed to be a parametric proportional hazards
#'   model with the baseline hazard specified via the \code{basehaz} argument.
#'   The (log) hazard of the event is assumed to be related to each longitudinal
#'   outcome via the association structure specified in the \code{assoc} argument.
#'   A different association may be specified for linking each longitudinal
#'   outcome to the hazard of the event.
#'
#'   Note that by default \code{simjm} will simulate data for one longitudinal
#'   outcome (i.e. a univariate joint model). If we want to simulate data for
#'   a multivariate joint model, or we wish to use one of the non-default
#'   association structures (for example "etaslope"), then we may need to adjust
#'   the "true" parameters accordingly so that we still get realistic survival
#'   times (for example, by changing the "true" association parameter value).
#'
#' @return A list of data frames, one for each longitudinal biomarker (in long
#'   format) and one for the event time data. The returned object also has a
#'   number of attributes that record the values specified for the arguments in
#'   the \code{simjm} call (for example the "true" parameter values, the number
#'   of individuals, the number of longitudinal markers, and so on).
#'
#' @examples
#' # For one longitudinal marker, we can just use the defaults:
#' simdat1 <- simjm()
#'
#' # For two longitudinal markers:
#' simdat2 <- simjm(M = 2, betaEvent_assoc = 0.1)
#'
#' # Simulate three markers for just 100 individuals and then return
#' # the true parameter values (which are stored as an attribute):
#' simdat3 <- simjm(M = 3, betaEvent_assoc = 0.1, n = 100)
#' attr(simdat3, "params")
#'
#' # Simulate one longitudinal marker, using "etaslope"
#' # association structure
#' simdat4 <- simjm(assoc = "etaslope", betaEvent_assoc = 0.8)
#'
#' # For one longitudinal marker, with a bernoulli outcome.
#' # (Note that 'betaLong_aux' specifies the number of trials
#' # for the binomial outcome).
#' simdat5 <- simjm(M = 1,
#'                  betaLong_intercept = 1.3,
#'                  betaLong_binary = -0.6,
#'                  betaLong_continuous = -0.03,
#'                  betaLong_slope = 0.05,
#'                  b_sd = c(1, 0.05),
#'                  betaEvent_intercept = -9,
#'                  betaEvent_assoc = 0.3,
#'                  family = binomial(),
#'                  betaLong_aux = 1)
#'
#' # For one longitudinal marker, with lower level clustering
#' # within individuals. The model includes only a random intercept
#' # at the individual-level, and then a random intercept and
#' # random slope at the lower cluster level. The association
#' # structure is based on a summation of the expected values for
#' # each of the lower level units.
#' simdat6 <- simjm(M = 1,
#'                  random_trajectory = "none",
#'                  betaLong_intercept = -1,
#'                  b_sd = 1,
#'                  clust_control = list(
#'                    L = 6,
#'                    assoc = "sum",
#'                    random_trajectory = "linear",
#'                    u_sd = c(1,1)
#'                  ))
#'
simjm <- function(n = 200, M = 1,
                  random_trajectory = "linear",
                  assoc = "etavalue",
                  basehaz = c("weibull"),
                  betaLong_intercept = 0,
                  betaLong_binary = 1,
                  betaLong_continuous = 1,
                  betaLong_slope = 1,
                  betaLong_poly1 = 1,
                  betaLong_poly2 = 0.2,
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
                  clust_control = list(),
                  seed = sample.int(.Machine$integer.max, 1),
                  interval = c(1E-8, 200)) {

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

  has_clust <- !missing(clust_control)
  if (has_clust) { # has lower-level clustering within the individual
    if (!is(clust_control, "list"))
      stop("'clust_control' should be a named list.")
    clust_nms <- names(clust_control)
    req_clust_args <- c("L", "assoc", "random_trajectory", "u_sd")
    if (!identical(sort(clust_nms), sort(req_clust_args)))
      stop("'clust_control' should include the following named arguments: ",
           paste(req_clust_args, collapse = ", "))
    if (clust_control$L < 1)
      stop("In clust_control, 'L' should be a positive integer.")
    if (!is.numeric(clust_control$u_sd) || any(clust_control$u_sd < 0))
      stop("In clust_control, 'u_sd' must be a positive scalar.")
    ok_clust_assocs <- c("sum", "mean")
    if (!clust_control$assoc %in% ok_clust_assocs)
      stop("In clust_control, 'assoc' must be one of: ",
           paste(ok_clust_assocs, collapse = ", "))
    marker1_traj_types <- c(random_trajectory[1L], clust_control$random_trajectory)
    if (!any(marker1_traj_types == "none")) {
      stop("If lower-level clustering within the individual is specified, ",
           "then the random effects must be limited to a random intercept ",
           "(ie. 'random_trajectory = none') at either the individual-level ",
           "or the clustering-level within individuals.")
    }
  }

  #----- Parameters

  # Broadcast user-specified parameters
  betaLong_intercept  <- maybe_broadcast(betaLong_intercept,  M)
  betaLong_binary     <- maybe_broadcast(betaLong_binary,     M)
  betaLong_continuous <- maybe_broadcast(betaLong_continuous, M)
  betaLong_slope      <- maybe_broadcast(betaLong_slope,      M)
  betaLong_poly1      <- maybe_broadcast(betaLong_poly1,      M)
  betaLong_poly2      <- maybe_broadcast(betaLong_poly2,      M)
  betaLong_aux        <- maybe_broadcast(betaLong_aux,        M)
  betaEvent_assoc     <- maybe_broadcast(betaEvent_assoc,     M)

  # Draw individual-level random effects
  b_dim <- sapply(random_trajectory, function(x)
    switch(x,
           none   = 1L, # random intercepts model
           linear = 2L, # random slopes model
           poly   = 3L) # random poly (degree = 2) model
  )
  b_dim_total <- sum(b_dim) # total num of individual-level random effects
  if (length(unique(b_dim)) == 1L) { # same num of reffs in each submodel
    if (length(b_sd) == b_dim[1]) {
      b_sd <- rep(b_sd, times = M)
    }
  }
  if (!length(b_sd) == b_dim_total)
    stop("'b_sd' appears to be the wrong length.")
  b_corr_mat <- matrix(rep(b_rho, b_dim_total ^ 2), ncol = b_dim_total)
  diag(b_corr_mat) <- 1
  b_dd <- MASS::mvrnorm(n = n, mu = rep(0, b_dim_total), Sigma = b_corr_mat)
  b <- sapply(1:length(b_sd), function(x) b_sd[x] * b_dd[,x])
  colnames(b) <- paste0("b", 1:b_dim_total)

  # Draw random effects if lower level clustering within individuals
  # NB lower level clustering only applies to the first longitudinal submodel
  if (has_clust) {
    L <- clust_control$L
    u_sd <- clust_control$u_sd
    u_dim <- switch(clust_control$random_trajectory,
                    none   = 1L, # random intercepts model
                    linear = 2L, # random slopes model
                    poly   = 3L) # random poly (degree = 2) model
    if (!length(u_sd) == u_dim)
      stop("In clust_control, 'u_sd' appears to be the wrong length. ",
           "Should be length ", u_dim, ".")
    Li <- as.integer(runif(n, 1, L+1)) # num units within each individual
    u_corr_mat <- matrix(rep(b_rho, u_dim ^ 2), ncol = u_dim)
    diag(u_corr_mat) <- 1
    u_dd <- MASS::mvrnorm(n = sum(Li), mu = rep(0, u_dim), Sigma = u_corr_mat)
    u <- sapply(1:length(clust_control$u_sd), function(x) u_sd[x] * u_dd[,x])
    colnames(u) <- paste0("u", 1:u_dim)
  }

  # Construct data frame of parameters
  betas <- list()
  for (m in 1:M) {
    nm <- paste0("Long", m)
    betas[[nm]] <- data.frame(id = 1:n)
    # fixed effect parameters
    betas[[nm]][[paste0("betaLong_intercept", m)]] <- rep(betaLong_intercept[m], n)
    betas[[nm]][[paste0("betaLong_binary", m)]] <- rep(betaLong_binary[m], n)
    betas[[nm]][[paste0("betaLong_continuous", m)]] <- rep(betaLong_continuous[m], n)
    if (random_trajectory[m] %in% c("none", "linear")) { # fixed effect slope
      betas[[nm]][[paste0("betaLong_slope", m)]] <- rep(betaLong_slope[m], n)
    } else if (random_trajectory[m] == "poly") { # fixed effect quadratic
      betas[[nm]][[paste0("betaLong_poly1", m)]] <- rep(betaLong_poly1[m], n)
      betas[[nm]][[paste0("betaLong_poly2", m)]] <- rep(betaLong_poly2[m], n)
    }
    # add on subject-specific intercept
    shift <- if (m == 1) 0 else sum(b_dim[1:(m - 1)])
    b_idx <- shift + 1
    betas[[nm]][[paste0("betaLong_intercept",  m)]] <-
      betas[[nm]][[paste0("betaLong_intercept",  m)]] + b[, b_idx]
    # add on subject-specific linear slope
    if (random_trajectory[m] == "linear") {
      b_idx <- shift + 2
      betas[[nm]][[paste0("betaLong_slope",  m)]] <-
        betas[[nm]][[paste0("betaLong_slope",  m)]] + b[, b_idx]
    }
    # add on subject-specific quadratic terms
    if (random_trajectory[m] == "poly") {
      b_idx <- shift + 2 # index of first quadratic term
      betas[[nm]][[paste0("betaLong_poly1",  m)]] <-
        betas[[nm]][[paste0("betaLong_poly1",  m)]] + b[, b_idx]
      b_idx <- shift + 3 # index of second quadratic term
      betas[[nm]][[paste0("betaLong_poly2",  m)]] <-
        betas[[nm]][[paste0("betaLong_poly2",  m)]] + b[, b_idx]
    }
  }
  betas[["Event"]] <- data.frame(
    id = 1:n,
    betaEvent_intercept  = rep(betaEvent_intercept,  n),
    betaEvent_binary     = rep(betaEvent_binary,     n),
    betaEvent_continuous = rep(betaEvent_continuous, n),
    betaEvent_aux        = rep(betaEvent_aux,        n)
  )
  for (m in 1:M)
    betas[["Event"]][[paste0("betaEvent_assoc", m)]] <- rep(betaEvent_assoc[m], n)
  if (has_clust) { # expand rows and add on cluster-specific random effects
    betas[["Long1"]] <-
      betas[["Long1"]][rep(row.names(betas[["Long1"]]), Li), , drop = FALSE]
    # cluster id within each individual
    betas[["Long1"]]$clust <- sequence(Li)
    # unique cluster id
    betas[["Long1"]]$clust_id <-
      paste(betas[["Long1"]][["id"]], betas[["Long1"]]$clust, sep = "_")
    # add on cluster-specific intercept
    betas[["Long1"]][["betaLong_intercept1"]] <-
      betas[["Long1"]][["betaLong_intercept1"]] + u[, 1]
    # add on cluster-specific linear slope
    if (clust_control$random_trajectory == "linear") {
      betas[["Long1"]][["betaLong_slope1"]] <-
        betas[["Long1"]][["betaLong_slope1"]] + u[, 2]
    }
    # add on cluster-specific quadratic terms
    if (clust_control$random_trajectory == "poly") {
      betas[["Long1"]][["betaLong_poly11"]] <-
        betas[["Long1"]][["betaLong_poly11"]] + u[, 2]
      betas[["Long1"]][["betaLong_poly21"]] <-
        betas[["Long1"]][["betaLong_poly21"]] + u[, 3]
    }
  }

  # Determine ultimate shape of the trajectory (used to determine which
  # parameters to build y_eta from). If there is no lower level clustering
  # within individuals, then this will just be equal to random_trajectory.
  # Otherwise, for the first longitudinal submodel it will be set equal to the
  # value of either random_trajectory[1L] or clust_control$random_trajectory
  # if one of those isn't set equal to "none".
  trajectory <- random_trajectory
  if (has_clust) {
    sel <- which(!marker1_traj_types == "none")
    trajectory[1L] <- marker1_traj_types[sel]
  }

  #----- Data

  # Generate baseline covariates - binary
  Z1 <- rbinom(n, 1, prob_Z1)      # covariate value for each subject

  # Generate baseline covariates - continuous
  Z2 <- rnorm(n, mean_Z2, sd_Z2)   # covariate value for each subject

  # Construct data frame of baseline covariates
  covs <- data.frame(id = 1:n, Z1, Z2)

  # Generate survival times
  ss <- simsurv::simsurv(hazard = jm_hazard, x = covs, betas = betas,
                         idvar = "id", ids = covs$id, trajectory = trajectory,
                         maxt = max_fuptime, basehaz = basehaz, M = M,
                         assoc = assoc, family = family, interval = interval)

  # Construct data frame of event data
  dat <- list(
    Event = data.frame(betas[["Event"]], covs, ss) # single row per subject
  )

  # Construct data frame of longitudinal data
  for (m in 1:M) {
    nm <- paste0("Long", m)
    dat[[nm]] <- merge(betas[[nm]], dat[["Event"]])
    dat[[nm]] <- merge(dat[[nm]], covs)
    dat[[nm]] <- dat[[nm]][rep(row.names(dat[[nm]]), max_yobs), ] # multiple row per subject
    dat[[nm]]$tij <- runif(nrow(dat[[nm]]), 0, max_fuptime)       # create observation times
    if (has_clust && m == 1) { # sort on ID, cluster ID and time
      dat[[nm]] <- dat[[nm]][order(dat[[nm]]$id, dat[[nm]]$clust_id, dat[[nm]]$tij), ]
    } else { # sort on ID and time
      dat[[nm]] <- dat[[nm]][order(dat[[nm]]$id, dat[[nm]]$tij), ]
    }
    dat[[nm]][[paste0("Xij_", m)]] <-
      dat[[nm]][[paste0("betaLong_intercept", m)]] +
      dat[[nm]][[paste0("betaLong_binary", m)]] * dat[[nm]][["Z1"]] +
      dat[[nm]][[paste0("betaLong_continuous", m)]] * dat[[nm]][["Z2"]]
    if (random_trajectory[m] %in% c("none", "linear")) {
      # if no random slope, then still use fixed linear slope
      dat[[nm]][[paste0("Xij_", m)]] <- dat[[nm]][[paste0("Xij_", m)]] +
        dat[[nm]][[paste0("betaLong_slope", m)]] * dat[[nm]][["tij"]]
    } else if (random_trajectory[m] == "poly") {
      dat[[nm]][[paste0("Xij_", m)]] <- dat[[nm]][[paste0("Xij_", m)]] +
        dat[[nm]][[paste0("betaLong_poly1",  m)]] * dat[[nm]][["tij"]] +
        dat[[nm]][[paste0("betaLong_poly2",  m)]] * dat[[nm]][["tij"]] * dat[[nm]][["tij"]]
    }
    fam <- family[[m]]$family
    invlink <- family[[m]]$linkinv
    mu <- invlink(dat[[nm]][[paste0("Xij_", m)]])
    if (fam == "gaussian") {
      sigma <- betaLong_aux[m]
      dat[[nm]][[paste0("Yij_", m)]] <- rnorm(length(mu), mu, sigma)
    } else if (fam == "binomial") {
      trials <- betaLong_aux[m]
      dat[[nm]][[paste0("Yij_", m)]] <- rbinom(length(mu), trials, mu)
    } else if (fam == "poisson") {
      dat[[nm]][[paste0("Yij_", m)]] <- rpois(length(mu), mu)
    } else if (fam == "neg_binomial_2") {
      size <- betaLong_aux[m]
      dat[[nm]][[paste0("Yij_", m)]] <- rnbinom(length(mu), size, mu)
    } else if (fam == "Gamma") {
      shape <- betaLong_aux[m]
      dat[[nm]][[paste0("Yij_", m)]] <- rgamma(length(mu), shape, rate = shape / mu)
    } else if (fam == "inverse.gaussian") {
      lambda <- betaLong_aux[m]
      dat[[nm]][[paste0("Yij_", m)]] <- .rinvGauss(length(mu), mu, lambda)
    }
  }

  # Prepare final datasets
  ret <- lapply(dat, function(x) {
    if ("tij" %in% colnames(x))
      x <- x[x$tij <= x$eventtime, ] # only keep rows before event time
    sel <- grep("^id$|^clust|^Z|^tij|^Yij|eventtime|status", colnames(x))
    x <- x[, sel, drop = FALSE]
    rownames(x) <- NULL
    return(x)
  })
  commonids <- ret[[1L]]$id
  for (i in 1:length(ret))
    commonids <- intersect(commonids, ret[[i]]$id)
  ret <- lapply(ret, function(x) x[x$id %in% commonids, , drop = FALSE])

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
    b_corr = b_corr_mat
  )

  # Return object
  structure(ret, params = c(long_params, event_params, re_params),
            n = length(unique(ret$Event$id)), M = M,
            max_yobs = max_yobs, max_fuptime = max_fuptime, assoc = assoc,
            family = family, random_trajectory = random_trajectory,
            clust_control = clust_control, seed = seed)
}


#--------------------- internal

# Return the hazard at time t, for a shared parameter joint model
#
# @param t The time variable
# @param x A named vector of covariate values
# @param betas A named vector of parameter values
# @param basehaz The type of baseline hazard
# @param M The number of longitudinal submodels
# @param assoc A vector of character strings, indicating the association
#   structure to use for each longitudinal submodel
# @param A list of family objects, providing the family (and link function)
#   for each longitudinal submodel
# @return A scalar
jm_hazard <- function(t, x, betas, basehaz = "weibull", M = 1,
                      trajectory = "linear", assoc = "etavalue",
                      family = list(gaussian()), grp_assoc = NULL) {

  if (t == 0)
    return(0) # boundary condition

  # Baseline hazard
  if (basehaz == "weibull") {
    h0 <- betas[["Event"]][["betaEvent_aux"]] *
      (t ^ (betas[["Event"]][["betaEvent_aux"]] - 1))
  }

  # Time-invariant part of event submodel eta
  etaevent <-
    betas[["Event"]][["betaEvent_intercept"]] +
    betas[["Event"]][["betaEvent_binary"]] * x[["Z1"]] +
    betas[["Event"]][["betaEvent_continuous"]] * x[["Z2"]]

  # Association structure
  for (m in 1:M) {

    nm <- paste0("Long", m)

    etabaseline_m <- etavalue_m <-
      betas[[nm]][[paste0("betaLong_intercept", m)]] +
      betas[[nm]][[paste0("betaLong_binary", m)]] * x[["Z1"]] +
      betas[[nm]][[paste0("betaLong_continuous", m)]] * x[["Z2"]]
    if (trajectory[m] %in% c("none", "linear")) {
      # if no random slope, then still use fixed linear slope
      etavalue_m <- etavalue_m +
        betas[[nm]][[paste0("betaLong_slope", m)]] * t
    } else if (trajectory[m] == "poly") {
      etavalue_m <- etavalue_m +
        betas[[nm]][[paste0("betaLong_poly1",  m)]] * t +
        betas[[nm]][[paste0("betaLong_poly2",  m)]] * (t * t)
    }
    if (!is.null(grp_assoc)) {
      if (grp_assoc == "sum") {
        res_etavalue_m <- sum(etavalue_m)
      } else if (grp_assoc == "mean") {
        res_etavalue_m <- mean(etavalue_m)
      }
    } else res_etavalue_m <- etavalue_m

    # eta value
    if (assoc[m] == "etavalue") {
      etaevent <- etaevent +
        betas[["Event"]][[paste0("betaEvent_assoc", m)]] * res_etavalue_m
    }

    # eta slope
    if (assoc[m] == "etaslope") {
      if (trajectory[m] %in% c("none", "linear")) {
        # if no random slope, then still use fixed linear slope
        etaslope_m <-
          betas[[nm]][[paste0("betaLong_slope", m)]]
      } else if (trajectory[m] == "poly") {
        etaslope_m <-
          betas[[nm]][[paste0("betaLong_poly1",  m)]] +
          betas[[nm]][[paste0("betaLong_poly2",  m)]] * (2 * t)
      }
      if (!is.null(grp_assoc)) {
        if (grp_assoc == "sum") {
          res_etaslope_m <- sum(etaslope_m)
        } else if (grp_assoc == "mean") {
          res_etaslope_m <- mean(etaslope_m)
        }
      } else res_etaslope_m <- etaslope_m
      etaevent <- etaevent +
        betas[["Event"]][[paste0("betaEvent_assoc", m)]] * res_etaslope_m
    }

    # eta auc
    if (assoc[m] == "etaauc") {
      if (trajectory[m] %in% c("none", "linear")) {
        etaauc_m <-
          etabaseline_m * t +
          0.5 * (etavalue_m - etabaseline_m) * t
      } else if (trajectory[m] == "poly") {
        stop("Calculation of hazard function for 'etaauc' with polynomial ",
             "trajectory not yet implemented.", call. = FALSE)
      }
      if (!is.null(grp_assoc))
        stop("'etaauc' association structure cannot currently be used ",
             "when there is lower level clustering within individuals.")
      etaevent <- etaevent +
        betas[["Event"]][[paste0("betaEvent_assoc", m)]] * etaauc_m
    }

    # mu value
    if (assoc[m] == "muvalue") {
      invlink <- family[[m]]$invlink
      muvalue_m <- invlink(etavalue_m)
      if (!is.null(grp_assoc)) {
        if (grp_assoc == "sum") {
          res_muvalue_m <- sum(muvalue_m)
        } else if (grp_assoc == "mean") {
          res_muvalue_m <- mean(muvalue_m)
        }
      } else res_muvalue_m <- muvalue_m
      etaevent <- etaevent +
        betas[["Event"]][[paste0("betaEvent_assoc", m)]] * muvalue_m
    }

  }
  # Calculate and return hazard
  h <- h0 * exp(etaevent)
  return(h)
}
