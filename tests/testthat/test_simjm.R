
library(lme4)

sim_run <- function(nsims = 200, npats = 200) {

  dat1 <- simjm(n = npats, assoc = "null", random_trajectory = "poly", b_sd = c(2,1,1))
  pars <- unlist(attr(dat1, "params"))

  sim_run <- function() {
    simdat <- simjm(n = npats, assoc = "null", random_trajectory = "poly", b_sd = c(2,1,1))
    simdat$Long1$tij2 <- simdat$Long1$tij * simdat$Long1$tij
    m1 <- lme4::lmer(Yij_1 ~ Z1 + Z2 + tij + tij2 + (tij + tij2 | id), data = simdat$Long1)
    fixef(m1)
  }

  sims <- replicate(nsims, sim_run())

  translate <- function(nm) {
    switch(nm,
           "betaLong_intercept" = "^\\(Intercept\\)$",
           "betaLong_binary" = "^Z1$",
           "betaLong_continuous" = "^Z2$",
           "betaLong_poly1" = "^tij$",
           "betaLong_poly2" = "^tij2$")
  }

  for (nm in names(pars)) {
    est_nm <- translate(nm)
    if (length(est_nm)) {
      rownames(sims) <- gsub(est_nm, nm, rownames(sims))
    }
  }

  value   <- rowMeans(sims)
  bias    <- rowMeans(sims) - pars[rownames(sims)]
  relbias <- 100 * (rowMeans(sims) - pars[rownames(sims)]) / pars[rownames(sims)]

  data.frame(value, bias, relbias)
}

res <- sim_run(200, 200)


