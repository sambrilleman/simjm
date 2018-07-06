# Copyright (C) 2017,2018 Sam Brilleman

#' Simulate data for a univariate or multivariate joint model
#'
#' Returns a list of data frames containing data simulated from a joint model for
#' longitudinal and time-to-event data.
#'
#' @export
#' @importFrom stats gaussian
#' @importFrom methods is
#'
#' @param n Number of individuals
#' @param M Number of longitudinal markers
#' @param fixed_trajectory The desired type of trajectory in the fixed effects
#'   part of the longitudinal model. Can be \code{"none"} (intercept only),
#'   \code{"linear"} (the default, intercept and linear slope),
#'   \code{"quadratic"} (intercept, linear slope, and quadratic term),
#'   or \code{"cubic"} (intercept, linear slope, quadratic term, and
#'   cubic term). Can be a single character string, or a character
#'   vector of length M (if a different trajectory type is to be used
#'   for each longitudinal submodel). Note that in addition to these time
#'   effects, two baseline covariates (one binary and one continuous) are
#'   always included in the longitudinal submodel as well.
#' @param random_trajectory The desired type of trajectory in the random effects
#'   part of the longitudinal model. Can be \code{"none"} (random intercept only),
#'   \code{"linear"} (the default, random intercept and linear slope),
#'   \code{"quadratic"} (random intercept, linear slope, and quadratic term),
#'   or \code{"cubic"} (random intercept, linear slope, quadratic term, and
#'   cubic term). Can be a single character string, or a character
#'   vector of length M (if a different trajectory type is to be used
#'   for each longitudinal submodel). Note that the corresponding
#'   \code{fixed_trajectory} argument must be at least as complex as the
#'   random effect structure; for example, you cannot specify
#'   \code{fixed_trajectory = "linear"} and \code{random_trajectory = "cubic"},
#'   because there would be no corresponding fixed effects parameters for
#'   the quadratic and cubic random terms in the model.
#' @param assoc A character string, or a character vector of length M,
#'   specifying the desired type of association structure for
#'   linking each longitudinal outcome to the hazard of the event.
#'   Only one association structure type can be specified for each longitudinal
#'   outcome. The options are: \code{"null"}, \code{"etavalue"},
#'   \code{"etaslope"}, \code{"etaauc"}, \code{"muvalue"},
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
#' @param betaLong_linear True coefficient for the fixed effect linear term in
#'   the longitudinal submodel when \code{fixed_trajectory = "linear"},
#'   \code{fixed_trajectory = "quadratic"} or \code{fixed_trajectory = "cubic"}.
#'   Can be a scalar or a vector of length M.
#' @param betaLong_quadratic True coefficient for the fixed effect quadratic term
#'   in the longitudinal submodel when \code{fixed_trajectory = "quadratic"} or
#'   \code{fixed_trajectory = "cubic"}. Can be a scalar or a vector of length M.
#' @param betaLong_cubic True coefficient for the fixed effect cubic term in
#'   the longitudinal submodel when \code{fixed_trajectory = "cubic"}.
#'   Can be a scalar or a vector of length M.
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
#'   standard deviations for the individual-level random effects.
#' @param b_rho Correlation between the individual-level random effects. This is
#'   only relevant when there is a total of >1 individual-level random effects in
#'   the joint model. This can be specified as a scalar correlation term, which
#'   assumes the same (true) correlation
#'   between each of the individual-level random effects, or it can be a correlation
#'   matrix for the correlation between the individual-level random effects. When
#'   simulating data for a multivariate joint model, the structure of the
#'   correlation matrix is such that the first \eqn{K_1} columns/rows correspond
#'   to the individual-level random effects for the first longitudinal submodel,
#'   and the next \eqn{K_2} columns/rows correspond to the individual-level random
#'   effects for the second longitudinal submodel, and so on.
#' @param max_yobs The maximum allowed number of longitudinal measurements for
#'   each biomarker. The actual number of observed measurements will depend on the
#'   individuals event time. Every individual is forced to have at least a baseline
#'   measurement for each biomarker (i.e. a biomarker measurement at time 0). The
#'   remaining biomarker measurement times will be uniformly distributed between 0
#'   and \code{max_fuptime} if \code{balanced = FALSE}, or evenly spaced between
#'   0 and \code{max_fuptime} if \code{balanced = TRUE}.
#' @param max_fuptime The maximum follow up time in whatever the desired time
#'   units are. This time will also be used as the censoring time (i.e. for
#'   individuals who have a simulated survival time that is after \code{max_fuptime}).
#' @param balanced A logical, specifying whether the timings of the longitudinal
#'   measurements should be balanced across individuals. If \code{FALSE} (the
#'   default), then each individual will have a baseline longitudinal measurement
#'   and the remaining measurement times will be chosen randomly from a uniform
#'   distribution on the range \code{[0,max_fuptime]}. If \code{TRUE}, then each
#'   individual will have a baseline longitudinal measurement and the remaining
#'   measurement times will be at common times for each individual, chosen to be
#'   evenly spaced between \code{[0,max_fuptime]}.
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
#'   effects at the cluster level.}
#'   \item{u_rho}{A scalar or a correlation matrix providing the correlation
#'   between the random effects for the lower level clustering factor. This is
#'   only relevant if there are >1 random effects for the lower level clustering
#'   factor. This is specified in a similar way as the \code{b_rho} argument
#'   described above.}
#'   \item{random_trajectory}{The desired type of trajectory in the random effects
#'   part of the longitudinal model at the cluster level. Can be \code{"none"} to
#'   only include a random intercept at the cluster level, or otherwise
#'   \code{"linear"} to include a random intercept and linear slope term.}
#'   }
#' @param return_eta A logical, if \code{return_eta = TRUE}, then the simulated
#'   data will also include the value of longitudinal submodel's linear predictor
#'   evaluated at the various measurement times. These will be stored in the
#'   variables named \code{Xij_1}, \code{Xij_2}, etc. (Note that if
#'   \code{family = "gaussian"} then these are equivalent to the error-free
#'   values of the biomarker at the measurement times).
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
#'   number of attributes, including a class attribute \code{simjm} and
#'   a number of other attributes that record the arguments that were specified in
#'   the \code{simjm} call; for example the "true" parameter values, the number
#'   of individuals, the number of longitudinal markers, and so on.
#'
#' @examples
#' # Note that throughout the examples we specify 'n = 30' to
#' # ensure a small sample size so that the examples run quickly
#'
#' # Simulate one longitudinal marker (we can just use the defaults):
#' simdat1 <- simjm(n = 30)
#'
#' # Simulate two longitudinal markers and then check the
#' # true parameter values (stored as an attribute):
#' simdat2 <- simjm(M = 2, n = 30, betaEvent_assoc = 0.1)
#' attr(simdat2, "params")
#'
#' # Simulate one longitudinal marker, using "etaslope"
#' # association structure:
#' simdat3 <- simjm(n = 30, assoc = "etaslope", betaEvent_assoc = 0.8)
#'
#' # For one longitudinal marker, with a bernoulli outcome.
#' # (Note that 'betaLong_aux' specifies the number of trials
#' # for the binomial outcome).
#' simdat4 <- simjm(M = 1, n = 30,
#'                  betaLong_intercept = 1.3,
#'                  betaLong_binary = -0.6,
#'                  betaLong_continuous = -0.03,
#'                  betaLong_linear = 0.05,
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
#' simdat5 <- simjm(M = 1, n = 30,
#'                  fixed_trajectory = "linear",
#'                  random_trajectory = "none",
#'                  betaLong_intercept = -1,
#'                  b_sd = 1,
#'                  clust_control = list(
#'                    L = 6,
#'                    assoc = "sum",
#'                    random_trajectory = "linear",
#'                    u_sd = c(1,1),
#'                    u_rho = 0.2
#'                  ))
#'
#' # Simulate two longitudinal markers, where the first longitudinal
#' # marker has a linear trajectory, and the second marker has
#' # quadratic trajectory. We will specify a correlation matrix for
#' # the individual-level random effects (via the 'b_rho' argument).
#' corrmat <- matrix(c(
#'   1.0, 0.2, 0.5, 0.0, 0.0,
#'   0.2, 1.0, 0.0, 0.2, 0.0,
#'   0.5, 0.0, 1.0, 0.3, 0.0,
#'   0.0, 0.2, 0.3, 1.0, -.1,
#'   0.0, 0.0, 0.0, -.1, 1.0), ncol = 5)
#' simdat6 <- simjm(M = 2, n = 30,
#'                  fixed_trajectory = c("linear", "quadratic"),
#'                  random_trajectory = c("linear", "quadratic"),
#'                  b_sd = c(2, 1, 2, 1, 0.5),
#'                  b_rho = corrmat,
#'                  betaEvent_assoc = c(0.1, 0.2))
#'
simjm <- function(n = 200, M = 1,
                  fixed_trajectory = "cubic",
                  random_trajectory = "linear",
                  assoc = "etavalue",
                  basehaz = c("weibull"),
                  betaLong_intercept = 10,
                  betaLong_binary = -1,
                  betaLong_continuous = 1,
                  betaLong_linear = -0.25,
                  betaLong_quadratic = 0.03,
                  betaLong_cubic = -0.0015,
                  betaLong_aux = 0.5,
                  betaEvent_intercept = -7.5,
                  betaEvent_binary = -0.5,
                  betaEvent_continuous = 0.5,
                  betaEvent_assoc = 0.5,
                  betaEvent_aux = 1.2,
                  b_sd = c(1.5, 0.07), b_rho = -0.2,
                  prob_Z1 = 0.5,
                  mean_Z2 = 0, sd_Z2 = 1,
                  max_yobs = 10,
                  max_fuptime = 20,
                  balanced = FALSE,
                  family = gaussian,
                  clust_control = list(),
                  return_eta = FALSE,
                  seed = sample.int(.Machine$integer.max, 1),
                  interval = c(1E-8, 200)) {

  #----- Preliminaries

  set.seed(seed)

  basehaz <- match.arg(basehaz)

  if (max_yobs < 1)
    stop("'max_yobs' must be at least 1.")

  # Check input assoc is valid
  ok_assocs <- c("etavalue", "etaslope", "etaauc", "muvalue",
                 "null", "shared_b(1)", "shared_coef(1)",
                 "shared_b(2)", "shared_coef(2)")
  assoc <- maybe_broadcast(assoc, M)
  if (!all(assoc %in% ok_assocs))
    stop("'assoc' must be one of: ", paste(ok_assocs, collapse = ", "))

  # Check input to trajectory type is valid
  ok_trajs  <- c("none", "linear", "quadratic", "cubic")
  fixed_trajectory  <- maybe_broadcast(fixed_trajectory,  M)
  random_trajectory <- maybe_broadcast(random_trajectory, M)
  if (!all(fixed_trajectory %in% ok_trajs))
    stop("'fixed_trajectory' must be one of: ", paste(ok_trajs, collapse = ", "))
  if (!all(random_trajectory %in% ok_trajs))
    stop("'random_trajectory' must be one of: ", paste(ok_trajs, collapse = ", "))
  if (!length(fixed_trajectory) == M)
    stop("'fixed_trajectory' is the wrong length.")
  if (!length(random_trajectory) == M)
    stop("'random_trajectory' is the wrong length.")

  # Check family is valid
  if (!is(family, "list"))
    family <- list(family)
  family <- maybe_broadcast(family, M)
  family <- lapply(family, validate_family)
  ok_families <- c("gaussian", "binomial")
  lapply(family, function(x) if (!x$family %in% ok_families)
    stop("'family' must be one of: ", paste(ok_families, collapse = ", ")))

  # Error check to ensure that the random effects structure isn't
  # more complex than the corresponding fixed effects structure
  fixed_traj_index  <- match(fixed_trajectory,  ok_trajs)
  random_traj_index <- match(random_trajectory, ok_trajs)
  sel <- which(random_traj_index > fixed_traj_index)
  if (length(sel))
    stop("The 'random_trajectory' cannot be more complex than the ",
         "corresponding 'fixed_trajectory'. This problem was encountered ",
         "for longitudinal submodel(s): ", paste(sel, collapse = ", "))

  # Error checks for assoc type
  for (m in 1:M) {
    shared_slope <- (assoc[m] %in% c("shared_b(2)", "shared_coef(2)"))
    if (shared_slope && (!random_trajectory[m] == "linear"))
      stop("Can only use 'shared_b(2)' or 'shared_coef(2)' when ",
           "'random_trajectory' is linear.")
  }

  # Check specified structure for any lower-level clustering
  has_clust <- !missing(clust_control)
  if (has_clust) { # has lower-level clustering within the individual
    if (!is(clust_control, "list"))
      stop("'clust_control' should be a named list.")
    clust_nms <- names(clust_control)
    ok_clust_args <- c("L", "assoc", "random_trajectory", "u_sd", "u_rho")
    if (!all(clust_nms %in% ok_clust_args))
      stop("'clust_control' should only include the following named arguments: ",
           paste(ok_clust_args, collapse = ", "))
    req_clust_args <- c("L", "assoc", "random_trajectory", "u_sd")
    if (!all(req_clust_args %in% clust_nms))
      stop("'clust_control' must include the following named arguments: ",
           paste(req_clust_args, collapse = ", "))
    if (clust_control$L < 1)
      stop("In clust_control, 'L' should be a positive integer.")
    if (!is.numeric(clust_control$u_sd) || any(clust_control$u_sd < 0))
      stop("In clust_control, 'u_sd' must be a positive scalar.")
    ok_clust_assocs <- c("sum", "mean", "max", "min")
    if (!clust_control$assoc %in% ok_clust_assocs)
      stop("In clust_control, 'assoc' must be one of: ",
           paste(ok_clust_assocs, collapse = ", "))
    ok_clust_trajs  <- c("none", "linear")
    if (!(clust_control$random_trajectory %in% ok_clust_trajs))
      stop("'clust_control$random_trajectory' must be one of: ",
           paste(ok_clust_trajs, collapse = ", "))
    marker1_traj_types <- c(random_trajectory[1L], clust_control$random_trajectory)
    if (!any(marker1_traj_types == "none")) {
      stop("If lower-level clustering within the individual is specified, ",
           "then the random effects must be limited to a random intercept ",
           "(ie. 'random_trajectory = none') at either the individual-level ",
           "or the clustering-level within individuals.")
    }
    grp_assoc <- clust_control$assoc
  } else {
    grp_assoc <- NULL
  }

  #----- Parameters

  # Broadcast user-specified parameters
  betaLong_intercept  <- maybe_broadcast(betaLong_intercept,  M)
  betaLong_binary     <- maybe_broadcast(betaLong_binary,     M)
  betaLong_continuous <- maybe_broadcast(betaLong_continuous, M)
  betaLong_linear     <- maybe_broadcast(betaLong_linear,     M)
  betaLong_quadratic  <- maybe_broadcast(betaLong_quadratic,  M)
  betaLong_cubic      <- maybe_broadcast(betaLong_cubic,      M)
  betaLong_aux        <- maybe_broadcast(betaLong_aux,        M)
  betaEvent_assoc     <- maybe_broadcast(betaEvent_assoc,     M)

  # Draw individual-level REs
  b_dim <- sapply(random_trajectory, function(x)
    switch(x,
           none      = 1L, # random intercepts model
           linear    = 2L, # random slopes model
           quadratic = 3L, # random quadratic terms
           cubic     = 4L) # random cubic terms
  )
  b_dim_total <- sum(b_dim) # total num of individual-level REs

  # Validate b_sd
  if (length(unique(b_dim)) == 1L) { # same num of REs in each submodel
    if (length(b_sd) == b_dim[1]) {
      b_sd <- rep(b_sd, times = M)
    }
  }
  if (!length(b_sd) == b_dim_total)
    stop("'b_sd' appears to be the wrong length.")

  # Validate b_rho
  if (b_dim_total == 1) { # only one RE, no corr matrix needed
    b_corr_mat = matrix(1,1,1)
  } else { # >1 RE, requires corr matrix
    if (is.scalar(b_rho) && b_dim_total > 1L) {
      # user supplied a constant correlation term for REs
      b_corr_mat <- matrix(rep(b_rho, b_dim_total ^ 2), ncol = b_dim_total)
      diag(b_corr_mat) <- 1
    } else {
      # user supplied a correlation matrix for REs
      b_corr_mat <- validate_corr_matrix(b_rho)
    }
  }

  # Draw standardised REs and scale them by b_sd
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
                    linear = 2L) # random slopes model
    if (!length(u_sd) == u_dim)
      stop("In clust_control, 'u_sd' appears to be the wrong length. ",
           "Should be length ", u_dim, ".")
    Li <- as.integer(stats::runif(n, 1, L+1)) # num units within each individual
    if (u_dim == 1) { # only one RE, no corr matrix needed
      u_corr_mat = matrix(1,1,1)
    } else { # >1 RE, requires corr matrix
      u_rho <- clust_control$u_rho
      if (is.null(u_rho))
        stop("'u_rho' must be provided in the clust_control list when there is ",
             "more than one random effect for the lower level clustering factor.")
      if (is.scalar(u_rho) && u_dim > 1L) {
        # user supplied a constant correlation term for REs
        u_corr_mat <- matrix(rep(u_rho, u_dim ^ 2), ncol = u_dim)
        diag(u_corr_mat) <- 1
      } else {
        # user supplied a correlation matrix for REs
        u_corr_mat <- validate_corr_matrix(u_rho)
      }
    }
    u_dd <- MASS::mvrnorm(n = sum(Li), mu = rep(0, u_dim), Sigma = u_corr_mat)
    u <- sapply(1:length(clust_control$u_sd), function(x) u_sd[x] * u_dd[,x])
    colnames(u) <- paste0("u", 1:u_dim)
  }

  # Construct data frame of parameters
  betas <- list()

  # Longitudinal submodel parameters
  for (m in 1:M) {
    nm <- paste0("Long", m)
    betas[[nm]] <- data.frame(id = 1:n)

    # fixed effect parameters
    betas[[nm]][[paste0("betaLong_intercept", m)]] <- rep(betaLong_intercept[m], n)
    betas[[nm]][[paste0("betaLong_binary", m)]] <- rep(betaLong_binary[m], n)
    betas[[nm]][[paste0("betaLong_continuous", m)]] <- rep(betaLong_continuous[m], n)
    if (fixed_trajectory[m] %in% c("linear", "quadratic", "cubic")) {
      # fixed effect linear
      betas[[nm]][[paste0("betaLong_linear", m)]] <- rep(betaLong_linear[m], n)
    }
    if (fixed_trajectory[m] %in% c("quadratic", "cubic")) {
      # fixed effect quadratic
      betas[[nm]][[paste0("betaLong_quadratic", m)]] <- rep(betaLong_quadratic[m], n)
    }
    if (fixed_trajectory[m] %in% c("cubic")) {
      # fixed effect cubic
      betas[[nm]][[paste0("betaLong_cubic", m)]] <- rep(betaLong_cubic[m], n)
    }

    # add on subject-specific intercept
    shift <- if (m == 1) 0 else sum(b_dim[1:(m - 1)])
    b_idx <- shift + 1
    betas[[nm]][[paste0("betaLong_intercept",  m)]] <-
      betas[[nm]][[paste0("betaLong_intercept",  m)]] + b[, b_idx]

    # add on subject-specific linear term
    if (random_trajectory[m] %in% c("linear", "quadratic", "cubic")) {
      b_idx <- shift + 2
      betas[[nm]][[paste0("betaLong_linear",  m)]] <-
        betas[[nm]][[paste0("betaLong_linear",  m)]] + b[, b_idx]
    }

    # add on subject-specific quadratic term
    if (random_trajectory[m] %in% c("quadratic", "cubic")) {
      b_idx <- shift + 3
      betas[[nm]][[paste0("betaLong_quadratic",  m)]] <-
        betas[[nm]][[paste0("betaLong_quadratic",  m)]] + b[, b_idx]
    }

    # add on subject-specific cubic term
    if (random_trajectory[m] %in% c("cubic")) {
      b_idx <- shift + 4
      betas[[nm]][[paste0("betaLong_cubic",  m)]] <-
        betas[[nm]][[paste0("betaLong_cubic",  m)]] + b[, b_idx]
    }
  }

  # Additional clustering factors (only allowed for longitudinal submodel 1)
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
    if (clust_control$random_trajectory %in% "linear") {
      betas[["Long1"]][["betaLong_linear1"]] <-
        betas[["Long1"]][["betaLong_linear1"]] + u[, 2]
    }
  }

  # Event submodel parameters
  betas[["Event"]] <- data.frame(
    id = 1:n,
    betaEvent_intercept  = rep(betaEvent_intercept,  n),
    betaEvent_binary     = rep(betaEvent_binary,     n),
    betaEvent_continuous = rep(betaEvent_continuous, n),
    betaEvent_aux        = rep(betaEvent_aux,        n)
  )
  for (m in 1:M)
    betas[["Event"]][[paste0("betaEvent_assoc", m)]] <- rep(betaEvent_assoc[m], n)

  #----- Data

  # Generate baseline covariates - binary
  Z1 <- stats::rbinom(n, 1, prob_Z1)      # covariate value for each subject

  # Generate baseline covariates - continuous
  Z2 <- stats::rnorm(n, mean_Z2, sd_Z2)   # covariate value for each subject

  # Construct data frame of baseline covariates
  covs <- data.frame(id = 1:n, Z1, Z2)

  # Generate survival times
  ss <- simsurv::simsurv(# the following arguments apply to 'simsurv'
                         hazard = jm_hazard, x = covs, betas = betas,
                         idvar = "id", ids = covs$id,
                         maxt = max_fuptime, interval = interval,
                         # the following arguments apply to 'jm_hazard'
                         basehaz = basehaz, M = M, trajectory = fixed_trajectory,
                         assoc = assoc, family = family, grp_assoc = grp_assoc)

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
    # create observation times
    if (balanced) { # longitudinal observation times balanced across individuals
      tij_seq  <- max_fuptime * 0:(max_yobs - 1) / max_yobs # baseline and post-baseline
      dat[[nm]]$tij <- rep(tij_seq, each = nrow(betas[[nm]]))
    } else { # longitudinal observation times unbalanced across individuals
      tij_seq1 <- rep(0, nrow(betas[[nm]])) # baseline
      tij_seq2 <- stats::runif(nrow(dat[[nm]]) - nrow(betas[[nm]]), 0, max_fuptime) # post-baseline
      dat[[nm]]$tij <- c(tij_seq1, tij_seq2)
    }
    if (has_clust && m == 1) { # sort on ID, cluster ID and time
      dat[[nm]] <- dat[[nm]][order(dat[[nm]]$id, dat[[nm]]$clust_id, dat[[nm]]$tij), ]
    } else { # sort on ID and time
      dat[[nm]] <- dat[[nm]][order(dat[[nm]]$id, dat[[nm]]$tij), ]
    }
    dat[[nm]][[paste0("Xij_", m)]] <-
      dat[[nm]][[paste0("betaLong_intercept", m)]] +
      dat[[nm]][[paste0("betaLong_binary", m)]] * dat[[nm]][["Z1"]] +
      dat[[nm]][[paste0("betaLong_continuous", m)]] * dat[[nm]][["Z2"]]
    if (fixed_trajectory[m] %in% c("linear", "quadratic", "cubic")) {
      dat[[nm]][[paste0("Xij_", m)]] <- dat[[nm]][[paste0("Xij_", m)]] +
        dat[[nm]][[paste0("betaLong_linear", m)]] *
        dat[[nm]][["tij"]]
    }
    if (fixed_trajectory[m] %in% c("quadratic", "cubic")) {
      dat[[nm]][[paste0("Xij_", m)]] <- dat[[nm]][[paste0("Xij_", m)]] +
        dat[[nm]][[paste0("betaLong_quadratic", m)]] *
        dat[[nm]][["tij"]] *
        dat[[nm]][["tij"]]
    }
    if (fixed_trajectory[m] %in% c("cubic")) {
      dat[[nm]][[paste0("Xij_", m)]] <- dat[[nm]][[paste0("Xij_", m)]] +
        dat[[nm]][[paste0("betaLong_cubic", m)]] *
        dat[[nm]][["tij"]] *
        dat[[nm]][["tij"]] *
        dat[[nm]][["tij"]]
    }
    fam <- family[[m]]$family
    invlink <- family[[m]]$linkinv
    mu <- invlink(dat[[nm]][[paste0("Xij_", m)]])
    if (fam == "gaussian") {
      sigma <- betaLong_aux[m]
      dat[[nm]][[paste0("Yij_", m)]] <- stats::rnorm(length(mu), mu, sigma)
    } else if (fam == "binomial") {
      trials <- betaLong_aux[m]
      dat[[nm]][[paste0("Yij_", m)]] <- stats::rbinom(length(mu), trials, mu)
    } else if (fam == "poisson") {
      dat[[nm]][[paste0("Yij_", m)]] <- stats::rpois(length(mu), mu)
    } else if (fam == "neg_binomial_2") {
      size <- betaLong_aux[m]
      dat[[nm]][[paste0("Yij_", m)]] <- stats::rnbinom(length(mu), size, mu)
    } else if (fam == "Gamma") {
      shape <- betaLong_aux[m]
      dat[[nm]][[paste0("Yij_", m)]] <- stats::rgamma(length(mu), shape, rate = shape / mu)
    } else if (fam == "inverse.gaussian") {
      lambda <- betaLong_aux[m]
      dat[[nm]][[paste0("Yij_", m)]] <- .rinvGauss(length(mu), mu, lambda)
    }
  }

  # Prepare final datasets
  ret <- lapply(dat, function(x) {
    if ("tij" %in% colnames(x))
      x <- x[x$tij <= x$eventtime, ] # only keep rows before event time
    if (return_eta) {
      sel <- grep("^id$|^clust|^Z|^tij|^Yij|^Xij|eventtime|status", colnames(x))
    } else {
      sel <- grep("^id$|^clust|^Z|^tij|^Yij|eventtime|status", colnames(x))
    }
    x <- x[, sel, drop = FALSE]
    rownames(x) <- NULL
    return(x)
  })

  # Only keep individuals who have at least one measurement for each biomarker
  # (NB Following the issues around delayed entry, this is now all individuals,
  #     since every individual is forced to have a baseline measurement. This
  #     may change in the future though, if there is a need to test whether
  #     rstanarm is correctly accommodating delayed entry, i.e. the situation
  #     in which we exclude individuals who do not have a baseline measurement.)
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
  if (any(fixed_trajectory %in% c("linear", "quadratic", "cubic")))
    long_params$betaLong_linear <- betaLong_linear
  if (any(fixed_trajectory %in% c("quadratic", "cubic")))
    long_params$betaLong_quadratic <- betaLong_quadratic
  if (any(fixed_trajectory %in% c("cubic")))
    long_params$betaLong_cubic <- betaLong_cubic
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
  cov_params <- nlist(
    prob_Z1, mean_Z2, sd_Z2
  )

  # Return object
  structure(ret,
            params = c(long_params, event_params, re_params, cov_params),
            n = length(unique(ret$Event$id)),
            M = M,
            max_yobs = max_yobs,
            max_fuptime = max_fuptime,
            balanced = balanced,
            assoc = assoc,
            family = family,
            fixed_trajectory = fixed_trajectory,
            random_trajectory = random_trajectory,
            return_eta = return_eta,
            clust_control = clust_control,
            seed = seed,
            class = c("simjm", class(ret)))
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
    if (trajectory[m] %in% c("linear", "quadratic", "cubic")) {
      etavalue_m <- etavalue_m +
        betas[[nm]][[paste0("betaLong_linear", m)]] * t
    }
    if (trajectory[m] %in% c("quadratic", "cubic")) {
      etavalue_m <- etavalue_m +
        betas[[nm]][[paste0("betaLong_quadratic", m)]] * (t * t)
    }
    if (trajectory[m] %in% c("cubic")) {
      etavalue_m <- etavalue_m +
        betas[[nm]][[paste0("betaLong_cubic", m)]] * (t * t * t)
    }

    if (!is.null(grp_assoc)) {
      if (grp_assoc == "sum") {
        res_etavalue_m <- sum(etavalue_m)
      } else if (grp_assoc == "mean") {
        res_etavalue_m <- mean(etavalue_m)
      }
    } else res_etavalue_m <- etavalue_m

    if (assoc[m] == "etavalue") {
      # eta value
      res_etavalue_m <- collapse_across_clusters(etavalue_m, grp_assoc)
      etaevent <- etaevent +
        betas[["Event"]][[paste0("betaEvent_assoc", m)]] * res_etavalue_m
    } else if (assoc[m] == "etaslope") {
      # eta slope
      if (trajectory[m] == "none") {
        etaslope_m <- 0
      } else if (trajectory[m] == "linear") {
        etaslope_m <-
          betas[[nm]][[paste0("betaLong_linear", m)]]
      } else if (trajectory[m] == "quadratic") {
        etaslope_m <-
          betas[[nm]][[paste0("betaLong_linear",  m)]] +
          betas[[nm]][[paste0("betaLong_quadratic",  m)]] * 2 * t
      } else if (trajectory[m] == "cubic") {
        etaslope_m <-
          betas[[nm]][[paste0("betaLong_linear",  m)]] +
          betas[[nm]][[paste0("betaLong_quadratic",  m)]] * 2 * t +
          betas[[nm]][[paste0("betaLong_cubic",  m)]] * 3 * I(t ^ 2)
      }
      res_etaslope_m <- collapse_across_clusters(etaslope_m, grp_assoc)
      etaevent <- etaevent +
        betas[["Event"]][[paste0("betaEvent_assoc", m)]] * res_etaslope_m
    } else if (assoc[m] == "etaauc") {
      # eta auc
      if (trajectory[m] == "none") {
        etaauc_m <-
          etabaseline_m * t
      } else if (trajectory[m] == "linear") {
        etaauc_m <-
          etabaseline_m * t +
          I(1/2) * betas[[nm]][[paste0("betaLong_linear", m)]] * I(t ^ 2)
      } else if (trajectory[m] == "quadratic") {
        etaauc_m <-
          etabaseline_m * t +
          I(1/2) * betas[[nm]][[paste0("betaLong_linear", m)]] * I(t ^ 2) +
          I(1/3) * betas[[nm]][[paste0("betaLong_linear", m)]] * I(t ^ 3)
      } else if (trajectory[m] == "cubic") {
        etaauc_m <-
          etabaseline_m * t +
          I(1/2) * betas[[nm]][[paste0("betaLong_linear", m)]] * I(t ^ 2) +
          I(1/3) * betas[[nm]][[paste0("betaLong_linear", m)]] * I(t ^ 3) +
          I(1/4) * betas[[nm]][[paste0("betaLong_linear", m)]] * I(t ^ 4)
      }
      res_etaauc_m <- collapse_across_clusters(etaauc_m, grp_assoc)
      etaevent <- etaevent +
        betas[["Event"]][[paste0("betaEvent_assoc", m)]] * res_etaauc_m
    } else if (assoc[m] == "muvalue") {
      # mu value
      invlink <- family[[m]]$invlink
      muvalue_m <- invlink(etavalue_m)
      res_muvalue_m <- collapse_across_clusters(muvalue_m, grp_assoc)
      etaevent <- etaevent +
        betas[["Event"]][[paste0("betaEvent_assoc", m)]] * res_muvalue_m
    }

  }

  # Calculate hazard
  h <- h0 * exp(etaevent)

  # Validate and return hazard
  if (!length(h) == 1) {
    stop("Bug found: returned hazard should be a scalar.")
  }
  return(h)
}

# Apply summary function in association structure when there is lower-level
# clustering within patients
#
# @param x A scalar or numeric vector with the 'etavalue', 'etaslope', 'etaauc'
#   value (or whatever quantity is being used in the association structure) for
#   one patient.
# @param collapse_fun A function to match and then apply to 'x'. If NULL, then
#   'x' is just returned, in which case 'x' should just be a scalar (i.e. there
#   should only be a patient-level value in the association structure and no
#   lower-level clusters within a patient).
# @return A scalar, used as the association term in the event submodel.
collapse_across_clusters <- function(x, collapse_fun = NULL) {
  if (is.null(collapse_fun)) {
    return(x)
  } else {
    fn <- match.fun(collapse_fun)
    return(fn(x))
  }
}
