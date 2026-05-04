#-----------------------------------------------------------------------------#
#----------------- Functions for Individual Data Representation ---------------
#-----------------------------------------------------------------------------#

# make.fclogit: pre-processing function for fclogit data
#             observations are ordered by stratum and then by outcome,
#             with cases before controls within each stratum.

make.fclogit = function(y, x, stratum = NULL) {
  
  if (class(x)[1] != 'matrix') {
    if (class(x) == 'data.frame') {
      x <- as.matrix(x)
    } else if (class(x) == 'numeric') {
      x <- matrix(x, ncol = 1)
    }
  }
  
  if (is.null(stratum)) stratum <- rep(1, length(y))
  ord <- order(stratum, -y)
  x <- x[ord,, drop = FALSE]
  y <- y[ord]
  stratum <- stratum[ord]
  
  unique.strata <- unique(stratum)
  n.strata <- length(unique.strata)
  n.i <- as.vector(table(factor(stratum, levels = unique.strata)))
  last.i <- cumsum(n.i)
  first.i <- c(1, last.i[-n.strata] + 1)
  
  new.strata <- rep(0, length(y))
  m1 <- n.c <- rep(0, n.strata)
  c <- x.unique <- t1 <- NULL
  
  for (i in 1:n.strata) {
    j1 <- first.i[i]
    j2 <- last.i[i]
    new.strata[j1:j2] <- i
    m1[i] <- sum(y[j1:j2])
    
    x.this <- x[j1:j2,, drop = FALSE]
    unique.x <- unique(x.this, margin = 1)
    index <- match(data.frame(t(x.this)), data.frame(t(unique.x)))
    n.c[i] <- nrow(unique.x)
    
    c <- c(c, as.vector(table(index)))
    x.unique <- rbind(x.unique, unique.x)
    tmp <- factor(y[j1:j2], levels = c(0, 1))
    t1 <- c(t1, as.vector(table(index, tmp)[,"1"]))
  }
  
  last.unique <- cumsum(n.c)
  first.unique <- c(1, last.unique[-n.strata] + 1)
  
  fclogit.data <- list(y = y, x = x, strata = new.strata,
                       n.i = n.i, first.i = first.i, last.i = last.i,
                       first.unique = first.unique, last.unique = last.unique,
                       c = c, n.c = n.c, x.unique = x.unique,
                       n.strata = n.strata, m1 = m1, t1 = t1)
  return(fclogit.data)
}

# (Main function)
# fclogit.fit: wrapper for conditional logistic regression with/without Firth correction
#           the underlying score functions are calculated by Rcpp.
#
# Notes:
#   method = "hg": force exact Howard-Gail recursion for all strata
#   method = "sp": force saddlepoint approximation for all strata
#   method = "auto": choose HG/SP separately for each stratum

fclogit.fit = function(y, x, stratum = NULL, method = c("auto", "hg", "sp"),
                       firth = FALSE, beta0 = NULL, alpha0 = NULL,
                       c0 = 312500, h0 = 500,
                       se = TRUE, eps = 10^-6,
                       control = list(ftol = 10^-8, xtol = 10^-8, maxit = 100)) {
  
  method <- match.arg(method)
  fclogit.data <- make.fclogit(y = y, x = x, stratum = stratum)
  
  n.strata <- fclogit.data$n.strata
  n.var <- ncol(fclogit.data$x)
  
  if (is.null(beta0)) beta0 <- rep(0, n.var)
  beta0 <- as.vector(beta0)
  
  if (is.null(colnames(fclogit.data$x))) {
    x.names <- paste0("x", 1:n.var)
  } else {
    x.names <- colnames(fclogit.data$x)
  }
  
  if (method == "auto") {
    stratum.method <- choose.fclogit.method(beta = beta0, fclogit.data = fclogit.data,
                                            firth = firth, c0 = c0, h0 = h0)
  } else {
    stratum.method <- rep(method, n.strata)
  }
  
  hg.flip <- get.hg.flip(fclogit.data = fclogit.data, stratum.method = stratum.method)
  n.sp <- sum(stratum.method == "sp")
  
  if (n.sp > 0) {
    if (is.null(alpha0)) {
      alpha0.all <- get.alpha0.sp(beta = beta0, fclogit.data = fclogit.data)
      alpha0 <- alpha0.all[stratum.method == "sp"]
    }
    par0 <- c(alpha0, beta0)
  } else {
    par0 <- beta0
  }
  
  stratum.method.num <- ifelse(stratum.method == "hg", 1, 2)
  
  score.fun <- function(par) {
    score_hybrid(par = par, x_all = fclogit.data$x, x_unique = fclogit.data$x.unique,
                 ns = fclogit.data$n.i, m1s = fclogit.data$m1,
                 t1s = fclogit.data$t1, cs = fclogit.data$c,
                 first_unique = fclogit.data$first.unique,
                 last_unique = fclogit.data$last.unique,
                 stratum_method = stratum.method.num,
                 hg_flip = as.integer(hg.flip), firth = firth)
  }
  
  fit <- nleqslv(x = par0, fn = score.fun, control = control)
  par.hat <- as.vector(fit$x)
  beta.ind <- n.sp + 1:n.var
  beta.hat <- as.vector(par.hat[beta.ind])
  names(beta.hat) <- x.names
  
  alpha.hat <- rep(NA, n.strata)
  if (n.sp > 0) alpha.hat[stratum.method == "sp"] <- par.hat[1:n.sp]
  
  vcov <- NULL
  std.err <- rep(NA, n.var)
  names(std.err) <- x.names
  
  if (se) {
    jac <- num.jacobian(fn = score.fun, theta = par.hat, eps = eps)
    vcov.all <- tryCatch(solve(-jac), error = function(e) MASS::ginv(-jac))
    vcov <- vcov.all[beta.ind, beta.ind, drop = FALSE]
    rownames(vcov) <- colnames(vcov) <- x.names
    std.err <- sqrt(diag(vcov))
  }
  
  method.label <- ifelse(all(stratum.method == "hg"), "hg",
                         ifelse(all(stratum.method == "sp"), "sp", "hybrid"))
  
  res <- list(coefficients = beta.hat, alpha = alpha.hat,
              std.err = std.err, vcov = vcov,
              method = method.label, requested.method = method,
              stratum.method = stratum.method, hg.flip = hg.flip,
              n.hg = sum(stratum.method == "hg"), n.sp = n.sp,
              firth = firth, nleqslv = fit, fclogit.data = fclogit.data,
              call = match.call())
  class(res) <- "fclogit"
  return(res)
}


# choose.fclogit.method: per-stratum switching rule for auto method

choose.fclogit.method = function(beta, fclogit.data, firth = FALSE, c0 = 312500, h0 = 500) {
  
  n.strata <- fclogit.data$n.strata
  n.var <- ncol(fclogit.data$x)
  q <- ifelse(firth, 3, 2)
  stratum.method <- rep("hg", n.strata)
  
  for (k in 1:n.strata) {
    n <- fclogit.data$n.i[k]
    m1 <- fclogit.data$m1[k]
    r <- min(m1, n - m1)
    C <- (n - r) * r * n.var^q
    
    if (C > c0) {
      stratum.method[k] <- "sp"
      next
    }
    
    i <- fclogit.data$first.i[k]
    j <- fclogit.data$last.i[k]
    x.this <- fclogit.data$x[i:j,, drop = FALSE]
    beta.this <- if (m1 > n / 2) -beta else beta
    eta <- sort(as.vector(x.this %*% beta.this))
    L.beta <- sum(eta[1:r])
    U.beta <- sum(eta[(n - r + 1):n])
    
    if (L.beta < -h0 || U.beta > h0) stratum.method[k] <- "sp"
  }
  
  return(stratum.method)
}


# get.hg.flip: switch case/control labels for HG when controls are fewer than cases

get.hg.flip = function(fclogit.data, stratum.method) {
  hg.flip <- rep(FALSE, fclogit.data$n.strata)
  ind <- which(stratum.method == "hg")
  hg.flip[ind] <- fclogit.data$m1[ind] > fclogit.data$n.i[ind] / 2
  return(hg.flip)
}


# get.alpha0.sp: initial values for stratum-specific intercepts in SP

get.alpha0.sp = function(beta, fclogit.data) {
  
  n.strata <- fclogit.data$n.strata
  alpha0 <- rep(0, n.strata)
  eta <- as.vector(fclogit.data$x.unique %*% beta)
  
  for (k in 1:n.strata) {
    i <- fclogit.data$first.unique[k]
    j <- fclogit.data$last.unique[k]
    p0 <- fclogit.data$m1[k] / fclogit.data$n.i[k]
    alpha0[k] <- qlogis(p0) - weighted.mean(eta[i:j], w = fclogit.data$c[i:j])
  }
  
  return(alpha0)
}


# num.jacobian: numerical Jacobian for score equations

num.jacobian = function(fn, theta, eps = 10^-6, method = c("central", "forward")) {
  
  method <- match.arg(method)
  theta <- as.vector(theta)
  f0 <- fn(theta)
  p <- length(theta)
  q <- length(f0)
  jac <- matrix(0, q, p)
  
  for (j in 1:p) {
    step <- rep(0, p)
    step[j] <- eps * max(1, abs(theta[j]))
    
    if (method == "central") {
      fp <- fn(theta + step)
      fm <- fn(theta - step)
      jac[, j] <- (fp - fm) / (2 * step[j])
    } else {
      fp <- fn(theta + step)
      jac[, j] <- (fp - f0) / step[j]
    }
  }
  
  return(jac)
}


# print and summary methods

print.fclogit = function(x, digits = max(3, getOption("digits") - 3), ...) {
  
  cat("Conditional logistic regression\n")
  cat("Method:", toupper(x$method), "\n")
  cat("Firth correction:", ifelse(x$firth, "yes", "no"), "\n")
  cat("Strata using HG:", x$n.hg, "\n")
  cat("Strata using SP:", x$n.sp, "\n")
  cat("Convergence code:", x$nleqslv$termcd, "\n\n")
  print(round(x$coefficients, digits))
  invisible(x)
}

summary.fclogit = function(object, ...) {
  
  beta <- object$coefficients
  se <- object$std.err
  z <- beta / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  tab <- cbind(Estimate = beta, `Std. Error` = se, `z value` = z,
               `Pr(>|z|)` = p, OR = exp(beta))
  
  res <- list(call = object$call, method = object$method, firth = object$firth,
              n.hg = object$n.hg, n.sp = object$n.sp,
              convergence = object$nleqslv$termcd, coefficients = tab)
  class(res) <- "summary.fclogit"
  return(res)
}

print.summary.fclogit = function(x, digits = max(3, getOption("digits") - 3), ...) {
  
  cat("Conditional logistic regression\n")
  cat("Method:", toupper(x$method), "\n")
  cat("Firth correction:", ifelse(x$firth, "yes", "no"), "\n")
  cat("Strata using HG:", x$n.hg, "\n")
  cat("Strata using SP:", x$n.sp, "\n")
  cat("Convergence code:", x$convergence, "\n\n")
  printCoefmat(x$coefficients, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
  invisible(x)
}

coef.fclogit = function(object, ...) {
  object$coefficients
}

vcov.fclogit = function(object, ...) {
  object$vcov
}

#-----------------------------------------------------------------------------#
#------------------- Functions for Grouped Data Representation ----------------
#-----------------------------------------------------------------------------#

# make.fclogit.grouped: pre-processing function for grouped fclogit data
#     each row corresponds to a unique covariate pattern within a stratum.
#   t1: number of cases/events in each group
#   c: total number of observations in each group
#   x: grouped covariate matrix
#   stratum: stratum indicator

make.fclogit.grouped = function(t1, c, x, stratum = NULL) {
  
  if (class(x)[1] != 'matrix') {
    if (class(x) == 'data.frame') {
      x <- as.matrix(x)
    } else if (class(x) == 'numeric') {
      x <- matrix(x, ncol = 1)
    }
  }
  
  if (is.null(stratum)) stratum <- rep(1, length(t1))
  
  t1 <- as.vector(t1)
  c <- as.vector(c)
  stratum <- as.vector(stratum)
  
  if (length(t1) != length(c)) stop("t1 and c must have the same length.")
  if (nrow(x) != length(t1)) stop("nrow(x) must equal length(t1).")
  if (length(stratum) != length(t1)) stop("stratum must have the same length as t1.")
  if (any(t1 < 0) || any(c < 0)) stop("t1 and c must be non-negative.")
  if (any(t1 > c)) stop("Each t1 must be less than or equal to c.")
  
  ord <- order(stratum)
  x <- x[ord,, drop = FALSE]
  t1 <- t1[ord]
  c <- c[ord]
  stratum <- stratum[ord]
  
  unique.strata <- unique(stratum)
  n.strata <- length(unique.strata)
  n.c <- as.vector(table(factor(stratum, levels = unique.strata)))
  last.unique <- cumsum(n.c)
  first.unique <- c(1, last.unique[-n.strata] + 1)
  
  new.strata <- rep(0, length(t1))
  n.i <- m1 <- rep(0, n.strata)
  
  for (i in 1:n.strata) {
    j1 <- first.unique[i]
    j2 <- last.unique[i]
    new.strata[j1:j2] <- i
    n.i[i] <- sum(c[j1:j2])
    m1[i] <- sum(t1[j1:j2])
  }
  
  fclogit.data <- list(y = NULL, x = NULL, strata = NULL,
                       n.i = n.i, first.i = NULL, last.i = NULL,
                       first.unique = first.unique, last.unique = last.unique,
                       c = c, n.c = n.c, x.unique = x,
                       n.strata = n.strata, m1 = m1, t1 = t1,
                       grouped = TRUE)
  return(fclogit.data)
}


# expand.fclogit.grouped: expand grouped data to individual-level rows for HG.

expand.fclogit.grouped = function(fclogit.data) {
  
  n.strata <- fclogit.data$n.strata
  n.var <- ncol(fclogit.data$x.unique)
  x.all <- matrix(NA, nrow = sum(fclogit.data$n.i), ncol = n.var)
  y.all <- rep(NA, sum(fclogit.data$n.i))
  strata.all <- rep(NA, sum(fclogit.data$n.i))
  
  row.start <- 1
  first.i <- last.i <- rep(0, n.strata)
  
  for (s in 1:n.strata) {
    j1 <- fclogit.data$first.unique[s]
    j2 <- fclogit.data$last.unique[s]
    ind <- j1:j2
    
    x.s <- fclogit.data$x.unique[ind,, drop = FALSE]
    t1.s <- fclogit.data$t1[ind]
    c.s <- fclogit.data$c[ind]
    t0.s <- c.s - t1.s
    
    x.case <- x.s[rep(seq_along(t1.s), t1.s),, drop = FALSE]
    x.ctrl <- x.s[rep(seq_along(t0.s), t0.s),, drop = FALSE]
    x.s.all <- rbind(x.case, x.ctrl)
    y.s.all <- c(rep(1, sum(t1.s)), rep(0, sum(t0.s)))
    
    row.end <- row.start + nrow(x.s.all) - 1
    first.i[s] <- row.start
    last.i[s] <- row.end
    x.all[row.start:row.end,, drop = FALSE] <- x.s.all
    y.all[row.start:row.end] <- y.s.all
    strata.all[row.start:row.end] <- s
    row.start <- row.end + 1
  }
  
  colnames(x.all) <- colnames(fclogit.data$x.unique)
  fclogit.data$x <- x.all
  fclogit.data$y <- y.all
  fclogit.data$strata <- strata.all
  fclogit.data$first.i <- first.i
  fclogit.data$last.i <- last.i
  return(fclogit.data)
}


# get.hg.bound.grouped: compute lower/upper log-product bounds without expanding counts

get.hg.bound.grouped = function(eta, count, r) {
  
  ord <- order(eta)
  eta <- eta[ord]
  count <- count[ord]
  
  left <- r
  L.beta <- 0
  for (j in seq_along(eta)) {
    a <- min(left, count[j])
    L.beta <- L.beta + a * eta[j]
    left <- left - a
    if (left == 0) break
  }
  
  left <- r
  U.beta <- 0
  for (j in rev(seq_along(eta))) {
    a <- min(left, count[j])
    U.beta <- U.beta + a * eta[j]
    left <- left - a
    if (left == 0) break
  }
  
  return(c(L.beta = L.beta, U.beta = U.beta))
}

# choose.fclogit.method.grouped: per-stratum switching rule for grouped data

choose.fclogit.method.grouped = function(beta, fclogit.data, firth = FALSE,
                                         c0 = 312500, h0 = 500) {
  
  n.strata <- fclogit.data$n.strata
  n.var <- ncol(fclogit.data$x.unique)
  q <- ifelse(firth, 3, 2)
  stratum.method <- rep("hg", n.strata)
  
  for (k in 1:n.strata) {
    n <- fclogit.data$n.i[k]
    m1 <- fclogit.data$m1[k]
    r <- min(m1, n - m1)
    C <- (n - r) * r * n.var^q
    
    if (C > c0) {
      stratum.method[k] <- "sp"
      next
    }
    
    i <- fclogit.data$first.unique[k]
    j <- fclogit.data$last.unique[k]
    x.this <- fclogit.data$x.unique[i:j,, drop = FALSE]
    c.this <- fclogit.data$c[i:j]
    beta.this <- if (m1 > n / 2) -beta else beta
    eta <- as.vector(x.this %*% beta.this)
    bound <- get.hg.bound.grouped(eta = eta, count = c.this, r = r)
    
    if (bound["L.beta"] < -h0 || bound["U.beta"] > h0) stratum.method[k] <- "sp"
  }
  
  return(stratum.method)
}


# fclogit.fit.grouped: wrapper for grouped conditional logistic regression data
#
#  t1: number of y=1 observations in each covariate group.
#  (For microbiome analyses, t1 is the target-taxon count)
#  c: the total number of observations in each covariate group.
#  (For microbiome analyses, c = target-taxon count + reference-taxon count )

fclogit.fit.grouped = function(t1, c, x, stratum = NULL,
                               method = c("auto", "hg", "sp"),
                               firth = FALSE, beta0 = NULL, alpha0 = NULL,
                               c0 = 312500, h0 = 500,
                               se = TRUE, eps = 10^-6,
                               max.expand = 5 * 10^6,
                               control = list(ftol = 10^-8, xtol = 10^-8,
                                              maxit = 100)) {
  
  method <- match.arg(method)
  fclogit.data <- make.fclogit.grouped(t1 = t1, c = c, x = x, stratum = stratum)
  
  n.strata <- fclogit.data$n.strata
  n.var <- ncol(fclogit.data$x.unique)
  
  if (is.null(beta0)) beta0 <- rep(0, n.var)
  beta0 <- as.vector(beta0)
  
  if (is.null(colnames(fclogit.data$x.unique))) {
    x.names <- paste0("x", 1:n.var)
  } else {
    x.names <- colnames(fclogit.data$x.unique)
  }
  
  if (method == "auto") {
    stratum.method <- choose.fclogit.method.grouped(beta = beta0,
                                                    fclogit.data = fclogit.data,
                                                    firth = firth, c0 = c0, h0 = h0)
  } else {
    stratum.method <- rep(method, n.strata)
  }
  
  hg.flip <- get.hg.flip(fclogit.data = fclogit.data,
                         stratum.method = stratum.method)
  n.sp <- sum(stratum.method == "sp")
  n.hg <- sum(stratum.method == "hg")
  
  if (n.hg > 0) {
    if (sum(fclogit.data$n.i) > max.expand) {
      stop("HG strata require individual-level expansion in the current implementation. ",
           "Use method = 'sp', reduce c0, or increase max.expand if expansion is intended.")
    }
    fclogit.data <- expand.fclogit.grouped(fclogit.data)
    x.all <- fclogit.data$x
  } else {
    x.all <- matrix(numeric(0), nrow = 0, ncol = n.var)
    colnames(x.all) <- x.names
  }
  
  if (n.sp > 0) {
    if (is.null(alpha0)) {
      alpha0.all <- get.alpha0.sp(beta = beta0, fclogit.data = fclogit.data)
      alpha0 <- alpha0.all[stratum.method == "sp"]
    }
    par0 <- c(alpha0, beta0)
  } else {
    par0 <- beta0
  }
  
  stratum.method.num <- ifelse(stratum.method == "hg", 1, 2)
  
  score.fun <- function(par) {
    score_hybrid(par = par, x_all = x.all, x_unique = fclogit.data$x.unique,
                 ns = fclogit.data$n.i, m1s = fclogit.data$m1,
                 t1s = fclogit.data$t1, cs = fclogit.data$c,
                 first_unique = fclogit.data$first.unique,
                 last_unique = fclogit.data$last.unique,
                 stratum_method = stratum.method.num,
                 hg_flip = as.integer(hg.flip), firth = firth)
  }
  
  fit <- nleqslv(x = par0, fn = score.fun, control = control)
  par.hat <- as.vector(fit$x)
  beta.ind <- n.sp + seq_len(n.var)
  beta.hat <- as.vector(par.hat[beta.ind])
  names(beta.hat) <- x.names
  
  alpha.hat <- rep(NA, n.strata)
  if (n.sp > 0) alpha.hat[stratum.method == "sp"] <- par.hat[seq_len(n.sp)]
  
  vcov <- NULL
  std.err <- rep(NA, n.var)
  names(std.err) <- x.names
  
  if (se) {
    jac <- num.jacobian(fn = score.fun, theta = par.hat, eps = eps)
    vcov.all <- tryCatch(solve(-jac), error = function(e) MASS::ginv(-jac))
    vcov <- vcov.all[beta.ind, beta.ind, drop = FALSE]
    rownames(vcov) <- colnames(vcov) <- x.names
    std.err <- sqrt(diag(vcov))
  }
  
  method.label <- ifelse(all(stratum.method == "hg"), "hg",
                         ifelse(all(stratum.method == "sp"), "sp", "hybrid"))
  
  res <- list(coefficients = beta.hat, alpha = alpha.hat,
              std.err = std.err, vcov = vcov,
              method = method.label, requested.method = method,
              stratum.method = stratum.method, hg.flip = hg.flip,
              n.hg = n.hg, n.sp = n.sp,
              firth = firth, nleqslv = fit, fclogit.data = fclogit.data,
              grouped = TRUE, call = match.call())
  class(res) <- "fclogit"
  return(res)
}
