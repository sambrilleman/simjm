# Copyright (C) 2017,2018 Sam Brilleman

#' Plot method for simjm simulated data frames
#'
#' Plot the underlying longitudinal trajectory and hazard function
#' for the data generating model underpinning a \code{\link{simjm}}
#' call. The plots are deterministic functions that use the true
#' parameter values from the data generating model. These true
#' parameter values are stored as an attribute in the object returned
#' by a call to \code{\link{simjm}} (which is a list of data frames
#' with class \code{c("simjm", "list")}).
#'
#' @export
#'
#' @param x An object of class \code{"simjm"}, which is a list of
#'   data frames returned by a call to \code{\link{simjm}}.
#' @param type The type of plot to produce. Type \code{"trajectory"}
#'   produces a plot of the mean trajectory for each longitudinal
#'   outcome, or if \code{m} is specified then the mean trajectory
#'   for the \emph{m}th longitudinal outcome. Type \code{"basehaz"}
#'   produces a plot of the hazard function evaluated at the mean
#'   value of the baseline covariates and a zero contribution from
#'   the association structure. Type \code{"hazard"} produces a
#'   plot of the hazard function evaluated at the mean value of the
#'   baseline covariates and at the current mean value (i.e. the
#'   value at time \emph{t}) for each of the longitudinal outcome(s).
#'   The default \code{"all"} produces a plot matrix with all of the
#'   above plot types.
#' @param m A positive integer. If \code{NULL} then plots are shown
#'   for all the longitudinal outcomes (only relevant if the data
#'   was simulated for a multivariate joint model). Otherwise, one
#'   can specify a specific longitudinal outcome to produce the plots
#'   for. Note that is ignored if \code{type = "basehaz"} or
#'   \code{type = "hazard"}.
#'
#' @return A ggplot object, which can be further customised using
#'   the \strong{ggplot2} package.
#'
plot.simjm <- function(x,
                       type = c("all",
                                "trajectory",
                                "basehaz",
                                "hazard"),
                       m = 1) {

  type <- match.arg(type)

  M                <- attr(x, "M")
  max_fuptime      <- attr(x, "max_fuptime")
  params           <- attr(x, "params")
  assoc            <- attr(x, "assoc")
  fixed_trajectory <- attr(x, "fixed_trajectory")

  eta_y_fun <- function(t, params, m) {
    eta_y <-
      params$betaLong_intercept[m] +
      params$betaLong_binary[m]     * params$prob_Z1 +
      params$betaLong_continuous[m] * params$mean_Z2
    if (fixed_trajectory[m] %in% c("linear", "quadratic", "cubic"))
      eta_y <- eta_y + params$betaLong_linear[m] * (t^1)
    if (fixed_trajectory[m] %in% c("quadratic", "cubic"))
      eta_y <- eta_y + params$betaLong_quadratic[m] * (t^2)
    if (fixed_trajectory[m] %in% c("cubic"))
      eta_y <- eta_y + params$betaLong_cubic[m] * (t^3)
    eta_y
  }

  basehaz_fun <- function(t, params) {
    basehaz <- params$betaEvent_aux * (t ^ (params$betaEvent_aux - 1))
    basehaz <- basehaz * exp(
      params$betaEvent_intercept +
      params$betaEvent_binary * params$prob_Z1 +
      params$betaEvent_continuous * params$mean_Z2)
    basehaz
  }

  haz_fun <- function(t, params, M) {
    haz <- params$betaEvent_aux * (t ^ (params$betaEvent_aux - 1))
    haz <- haz * exp(
      params$betaEvent_intercept +
      params$betaEvent_binary * params$prob_Z1 +
      params$betaEvent_continuous * params$mean_Z2)
    for (m in 1:M) {
      if (assoc[m] == "etavalue") {
        eta_y <-
          params$betaLong_intercept[m] +
          params$betaLong_binary[m]     * params$prob_Z1 +
          params$betaLong_continuous[m] * params$mean_Z2
        if (fixed_trajectory[m] %in% c("linear", "quadratic", "cubic"))
          eta_y <- eta_y + params$betaLong_linear[m] * (t^1)
        if (fixed_trajectory[m] %in% c("quadratic", "cubic"))
          eta_y <- eta_y + params$betaLong_quadratic[m] * (t^2)
        if (fixed_trajectory[m] %in% c("cubic"))
          eta_y <- eta_y + params$betaLong_cubic[m] * (t^3)
        haz <- haz * exp(params$betaEvent_assoc[m] * eta_y)
      }
    }
    haz
  }

  gg <-
    ggplot2::ggplot(data.frame(t = c(0, max_fuptime)), ggplot2::aes(t)) +
    ggplot2::theme_set(ggplot2::theme_bw()) +
    ggplot2::xlab("Time")

  gg_eta_y <- lapply(1:M, function(m) {
    gg +
    ggplot2::stat_function(fun = eta_y_fun, args = list(params = params, m = m)) +
    ggplot2::ylab("Mean of longitudinal outcome")
  })

  gg_basehaz <-
    gg +
    ggplot2::stat_function(fun = basehaz_fun, args = list(params = params)) +
    ggplot2::ylab("Baseline hazard rate")

  gg_haz <-
    gg +
    ggplot2::stat_function(fun = haz_fun, args = list(params = params, M = M)) +
    ggplot2::ylab("Hazard rate (as function \n of the longitudinal outcome)")

  if (type == "all") {
    ggpubr::ggarrange(plotlist = c(
      gg_eta_y,
      list(gg_basehaz),
      list(gg_haz)))
  } else if (type == "trajectory") {
    gg_eta_y[[m]]
  } else if (type == "basehaz") {
    gg_basehaz
  } else if (type == "hazard") {
    gg_haz
  }
}



