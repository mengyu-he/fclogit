# Example file for FCLogit / FClogit prototype
#
# This file gives small examples for fitting conditional logistic regression
# with the R wrapper and Rcpp score functions.

library(Rcpp)
library(RcppArmadillo)
library(nleqslv)
library(MASS)

Rcpp::sourceCpp("fclogit_hg_sp.cpp")
source("fclogit_wrapper.R")

#-----------------------------------------------------------------------------#
#------------------------------ data simulation -------------------------------
#-----------------------------------------------------------------------------#

# simulate.fclogit.data: generate stratified binary data from the conditional
# logistic model with a fixed number of cases in each stratum.

simulate.fclogit.data = function(n.strata = 50, n = 6, m1 = 2,
                                 beta = c(log(2), -log(1.5)),
                                 seed = 1) {
  
  set.seed(seed)
  n.var <- length(beta)
  
  if (length(n) == 1) n <- rep(n, n.strata)
  if (length(m1) == 1) m1 <- rep(m1, n.strata)
  stopifnot(length(n) == n.strata, length(m1) == n.strata)
  stopifnot(all(m1 > 0), all(m1 < n))
  
  n.total <- sum(n)
  y <- rep(0, n.total)
  stratum <- rep(1:n.strata, times = n)
  x <- matrix(rnorm(n.total * n.var), ncol = n.var)
  colnames(x) <- paste0("x", 1:n.var)
  
  first <- c(1, cumsum(n)[-n.strata] + 1)
  last <- cumsum(n)
  
  for (s in 1:n.strata) {
    ind <- first[s]:last[s]
    eta <- as.vector(x[ind,, drop = FALSE] %*% beta)
    prob <- exp(eta - max(eta))
    case.ind <- sample(ind, size = m1[s], prob = prob)
    y[case.ind] <- 1
  }
  
  data.frame(y = y, stratum = stratum, x)
}


#-----------------------------------------------------------------------------#
#-------- Example 1: ordinary fclogit using automatic HG/SP selection --------
#-----------------------------------------------------------------------------#

n.vec <- c(rep(8, 10), rep(700, 10))
m1.vec <- c(rep(2, 10), rep(350, 10))

ex1 <- simulate.fclogit.data(n.strata = length(n.vec), n = n.vec, m1 = m1.vec,
                             beta = c(log(1.5), -log(1.3)), seed = 1001)

fit.auto <- fclogit.fit(y = ex1$y,
                        x = as.matrix(ex1[, c("x1", "x2")]),
                        stratum = ex1$stratum,
                        method = "auto",
                        firth = TRUE)

fit.auto
summary(fit.auto)

# The automatic rule is made per stratum.
table(fit.auto$stratum.method)

#-----------------------------------------------------------------------------#
#---------------------- Example 2: force HG and force SP ----------------------
#-----------------------------------------------------------------------------#

# For small strata, HG is feasible and can be used as an exact benchmark.
ex2 <- simulate.fclogit.data(n.strata = 80, n = 6, m1 = 2,
                             beta = c(log(2), -log(1.5)), seed = 2002)

fit.hg <- fclogit.fit(y = ex2$y,
                      x = as.matrix(ex2[, c("x1", "x2")]),
                      stratum = ex2$stratum,
                      method = "hg",
                      firth = FALSE)

# SP can also be forced for checking and comparison.
fit.sp <- fclogit.fit(y = ex2$y,
                      x = as.matrix(ex2[, c("x1", "x2")]),
                      stratum = ex2$stratum,
                      method = "sp",
                      firth = FALSE)

cbind(HG = coef(fit.hg), SP = coef(fit.sp))

#-----------------------------------------------------------------------------#
#--------------------- Example 3: Firth-corrected fclogit ---------------------
#-----------------------------------------------------------------------------#

# A sparse example where Firth correction may be useful.
ex3 <- simulate.fclogit.data(n.strata = 60, n = 4, m1 = 1,
                             beta = c(log(4)), seed = 2026)

fit.firth <- fclogit.fit(y = ex3$y,
                         x = as.matrix(ex3[, "x1", drop = FALSE]),
                         stratum = ex3$stratum,
                         method = "auto",
                         firth = TRUE)

summary(fit.firth)
exp(coef(fit.firth))

#-----------------------------------------------------------------------------#
#------------------ Example 4: demonstrate hybrid switching ------------------
#-----------------------------------------------------------------------------#

# The default c0 is intentionally permissive for many small examples.
# To demonstrate mixed HG/SP selection, use a smaller c0.
ex4 <- simulate.fclogit.data(n.strata = 20, n = 80, m1 = 40,
                             beta = c(log(1.5), -log(1.3)), seed = 3030)

fit.c0.large <- fclogit.fit(y = ex4$y,
                            x = as.matrix(ex4[, c("x1", "x2")]),
                            stratum = ex4$stratum,
                            method = "auto",
                            firth = FALSE,
                            c0 = 312500)

fit.c0.small <- fclogit.fit(y = ex4$y,
                            x = as.matrix(ex4[, c("x1", "x2")]),
                            stratum = ex4$stratum,
                            method = "auto",
                            firth = FALSE,
                            c0 = 1000)

table(default.c0 = fit.c0.large$stratum.method)
table(smaller.c0 = fit.c0.small$stratum.method)

#-----------------------------------------------------------------------------#
#------------------ Example 5: case/control flipping for HG ------------------
#-----------------------------------------------------------------------------#

# When a stratum has more cases than controls and HG is used, the implementation
# evaluates the HG recursion on the smaller side and converts the contribution
# back to the original beta scale.
ex5 <- simulate.fclogit.data(n.strata = 50, n = 8, m1 = 6,
                             beta = c(log(1.8), -log(1.4)), seed = 4040)

fit.flip <- fclogit.fit(y = ex5$y,
                        x = as.matrix(ex5[, c("x1", "x2")]),
                        stratum = ex5$stratum,
                        method = "hg",
                        firth = FALSE)

summary(fit.flip)
table(fit.flip$hg.flip)

#-----------------------------------------------------------------------------#
#-------- Example 6: extract estimates, standard errors, and intervals -------
#-----------------------------------------------------------------------------#

beta.hat <- coef(fit.auto)
se.hat <- sqrt(diag(vcov(fit.auto)))
ci.beta <- cbind(beta.hat - 1.96 * se.hat,
                 beta.hat + 1.96 * se.hat)
ci.or <- exp(ci.beta)

out <- data.frame(Estimate = beta.hat,
                  SE = se.hat,
                  OR = exp(beta.hat),
                  OR.Lower = ci.or[, 1],
                  OR.Upper = ci.or[, 2])

round(out, 4)


#-----------------------------------------------------------------------------#
#----------------- Example 7: check against survival::clogit -----------------
#-----------------------------------------------------------------------------#

# This check is only for ordinary conditional logistic regression without Firth
# correction. It is optional and requires the survival package.
if (requireNamespace("survival", quietly = TRUE)) {
  fit.clogit <- survival::clogit(y ~ x1 + x2 + strata(stratum), data = ex1)
  cbind(fclogit = coef(fit.hg), clogit = coef(fit.clogit))
}


#-----------------------------------------------------------------------------#
#----------------- Example 8: data in grouped representation -----------------
#-----------------------------------------------------------------------------#

dat = data.frame(
  pair = c(1, 1),
  X = c(1, 0),
  taxon_count = c(10, 20),
  ref_count = c(1000, 500)
)

fit.grouped = fclogit.fit.grouped(
  t1 = dat$taxon_count,
  c = dat$taxon_count + dat$ref_count,
  x = dat[, "X", drop = FALSE],
  stratum = dat$pair,
  method = "sp",
  firth = FALSE
)

summary(fit.grouped)

