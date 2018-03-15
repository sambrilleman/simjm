# Copyright (C) 2017,2018 Sam Brilleman

#' Plot observed longitudinal trajectories using facets
#'
#' This is a convenience function that allows the user to easily plot the
#' observed longitudinal biomarker trajectory for each individual. Each individual
#' is shown in a separate facet in the plot. This function can be useful
#' to assess whether the "true" parameter values passed to \code{\link{simjm}}
#' are generating realistic longitudinal trajectories for the simulated
#' biomarkers.
#'
#' @export
#'
#' @param yvar Name of the variable corresponding to the longitudinal biomarker.
#' @param data The data frame with simulated data for the longitudinal
#'   submodel. The variable indicating the individual should be called "id",
#'   as is the case in the data frames returned by \code{\link{simjm}}.
#' @param data_filter Optionally, a condition wrapped in quotes indicating the
#'   subset of individuals to include in the plot (for example we could use
#'   \code{data_filter = "id <= 20"} to include plots the first 20 individuals).
#'   If \code{NULL} then a separate facet is included for every individual in
#'   \code{data}.
#' @param timevar Name of the variable corresponding to time.
#' @param idvar Name of the variable identifying individuals.
#'
#' @return A \code{ggplot} object that can be further customised using the
#'   \pkg{ggplot2} package.
#'
#' @examples
#'   dat <- simjm(n = 20, random_trajectory = "poly", b_sd = c(2,2,2), seed = 123)
#'   plot_traj("Yij_1", data = dat$Long1)
#'   plot_traj("Yij_1", data = dat$Long1, data_filter = "id <= 10")
#'
plot_traj <- function(yvar = "Yij_1", data, data_filter = NULL,
                      timevar = "tij", idvar = "id") {
  if (!requireNamespace("ggplot2"))
    stop("This function requires 'ggplot2' to be installed.")
  if (!requireNamespace("dplyr"))
    stop("This function requires 'dplyr' to be installed.")
  dat <- if (!is.null(data_filter)) dplyr::filter_(data, data_filter) else data
  map <- ggplot2::aes_string(y = substitute(yvar), x = substitute(timevar))
  ggplot2::ggplot(dat, map) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(substitute(idvar), scales = "free")
}
