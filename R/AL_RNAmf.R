#' object to optimize the point by ALM criterion updating at level 1 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param fit an object of class RNAmf.
#' @return A negative predictive posterior variance at Xcand.
#' @noRd
#'

obj.ALM_level_1 <- function(fit, Xcand) {
  newx <- matrix(Xcand, nrow = 1)
  fit1 <- fit$fits[[1]]
  kernel <- fit$kernel
  
  ### calculate the posterior predictive variance ###
  if (kernel == "sqex") {
    predsig2 <- pred.GP(fit1, newx)$sig2
  } else if (kernel == "matern1.5") {
    predsig2 <- pred.matGP(fit1, newx)$sig2
  } else if (kernel == "matern2.5") {
    predsig2 <- pred.matGP(fit1, newx)$sig2
  }
  
  -predsig2 # to maximize the current variance.
}


#' object to optimize the point by ALM criterion updating at level 2 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param fit an object of class RNAmf.
#' @return A negative predictive posterior variance at Xcand.
#' @noRd
#'

obj.ALM_level_2 <- function(fit, Xcand) {
  newx <- matrix(Xcand, nrow = 1)
  -predict(fit, newx)$sig2[[2]] # to maximize the current variance.
}


#' object to optimize the point by ALM criterion updating at level 3 with three levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param fit an object of class RNAmf.
#' @return A predictive posterior variance at Xcand.
#' @noRd
#'

obj.ALM_level_3 <- function(fit, Xcand) {
  newx <- matrix(Xcand, nrow = 1)
  -predict(fit, newx)$sig2[[3]] # to maximize the current variance.
}


#' object to optimize the point by ALD criterion updating V1 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param fit an object of class RNAmf.
#' @return A negative V1 at Xcand.
#' @noRd
#'

obj.ALD_V1_2level <- function(fit, Xcand) { # low
  newx <- matrix(Xcand, nrow = 1)
  
  kernel <- fit$kernel
  constant <- fit$constant
  fit1 <- fit$fits[[1]]
  fit2 <- fit$fits[[2]]
  
  if (kernel == "sqex") {
    if (constant) {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.GP(fit1, newx)
      x.mu <- pred.fit$mu # mean of f1(u)
      sig2 <- pred.fit$sig2 #* 0
      
      ### calculate the closed form ###
      X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
      w1.x2 <- fit2$X[, d + 1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat
      mu2 <- fit2$mu.hat
      
      Ci <- fit2$Ki
      a <- crossprod(Ci, y2 - mu2)
      
      # mean
      predy <- mu2 + crossprod(t(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))) *
                                   1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                                   exp(-(drop(outer(x.mu, w1.x2, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))), a)
      
      # var
      mat <- drop(crossprod(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))),
                            exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))))) * # common components
        1 / sqrt(1 + 4 * sig2 / theta[d + 1]) *
        exp(-(outer(w1.x2, w1.x2, FUN = "+") / 2 - matrix(x.mu, n, n))^2 / (theta[d + 1] / 2 + 2 * sig2)) *
        exp(-(outer(w1.x2, w1.x2, FUN = "-"))^2 / (2 * theta[d + 1]))
      
      VE <- - (predy - mu2)^2 + drop(crossprod(a, crossprod(mat, a)))
    } else {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.GP(fit1, newx)
      x.mu <- pred.fit$mu # mean of f1(u)
      sig2 <- pred.fit$sig2
      
      ### calculate the closed form ###
      X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
      w1.x2 <- fit2$X[, d + 1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat
      
      Ci <- fit2$Ki
      a <- crossprod(Ci, y2)
      
      # mean
      predy <- crossprod(t(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))) *
                             1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                             exp(-(drop(outer(x.mu, w1.x2, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))), a)
      
      # var
      mat <- drop(crossprod(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))),
                            exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))))) * # common components
        1 / sqrt(1 + 4 * sig2 / theta[d + 1]) *
        exp(-(outer(w1.x2, w1.x2, FUN = "+") / 2 - matrix(x.mu, n, n))^2 / (theta[d + 1] / 2 + 2 * sig2)) *
        exp(-(outer(w1.x2, w1.x2, FUN = "-"))^2 / (2 * theta[d + 1]))
      
      VE <- - predy^2 + drop(crossprod(a, crossprod(mat, a)))
    }
  } else if (kernel == "matern1.5") {
    if (constant) {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
      x.mu <- drop(pred.fit$mu) # mean of f1(u)
      sig2 <- pred.fit$sig2 #* 0
      
      ### calculate the closed form ###
      X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
      w1.x2 <- fit2$X[, d + 1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat
      mu2 <- fit2$mu.hat
      
      Ci <- fit2$Ki
      a <- crossprod(Ci, y2 - mu2)
      
      # mean
      mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
      mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]
      
      lambda11 <- c(1, mua)
      lambda12 <- c(0, 1)
      lambda21 <- c(1, -mub)
      
      e1 <- cbind(matrix(1 - sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])
      e2 <- cbind(matrix(1 + sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])
      
      predy <- mu2 + drop(crossprod(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                                      (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                         (crossprod(t(e1), lambda11) * pnorm((mua - w1.x2) / sqrt(sig2)) +
                                            crossprod(t(e1), lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                                         exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                         (crossprod(t(e2), lambda21) * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                                            crossprod(t(e2), lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2)))), a))
      
      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 1.5, theta = theta[d + 1])
      }
      mat <- drop(crossprod(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5), cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))
      
      VE <- - (predy - mu2)^2 + drop(crossprod(a, crossprod(mat, a)))
    } else {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
      x.mu <- drop(pred.fit$mu) # mean of f1(u)
      sig2 <- pred.fit$sig2
      
      ### calculate the closed form ###
      X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
      w1.x2 <- fit2$X[, d + 1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat
      
      Ci <- fit2$Ki
      a <- crossprod(Ci, y2)
      
      # mean
      mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
      mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]
      
      lambda11 <- c(1, mua)
      lambda12 <- c(0, 1)
      lambda21 <- c(1, -mub)
      
      e1 <- cbind(matrix(1 - sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])
      e2 <- cbind(matrix(1 + sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])
      
      predy <- drop(crossprod(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                                (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e1), lambda11) * pnorm((mua - w1.x2) / sqrt(sig2)) +
                                      crossprod(t(e1), lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                                   exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e2), lambda21) * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                                      crossprod(t(e2), lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2)))), a))
      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 1.5, theta = theta[d + 1])
      }
      mat <- drop(crossprod(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5), cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))
      
      VE <- - predy^2 + drop(crossprod(a, crossprod(mat, a)))
    }
  } else if (kernel == "matern2.5") {
    if (constant) {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
      x.mu <- drop(pred.fit$mu) # mean of f1(u)
      sig2 <- pred.fit$sig2 #* 0
      
      ### calculate the closed form ###
      X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
      w1.x2 <- fit2$X[, d + 1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat
      mu2 <- fit2$mu.hat
      
      Ci <- fit2$Ki
      a <- crossprod(Ci, y2 - mu2)
      
      # mean
      mua <- x.mu - sqrt(5) * sig2 / theta[d + 1]
      mub <- x.mu + sqrt(5) * sig2 / theta[d + 1]
      
      lambda11 <- c(1, mua, mua^2 + sig2)
      lambda12 <- cbind(0, 1, matrix(mua + w1.x2))
      lambda21 <- c(1, -mub, mub^2 + sig2)
      lambda22 <- cbind(0, 1, matrix(-mub - w1.x2))
      
      e1 <- cbind(
        matrix(1 - sqrt(5) * w1.x2 / theta[d + 1] + 5 * w1.x2^2 / (3 * theta[d + 1]^2)),
        matrix(sqrt(5) / theta[d + 1] - 10 * w1.x2 / (3 * theta[d + 1]^2)),
        5 / (3 * theta[d + 1]^2)
      )
      e2 <- cbind(
        matrix(1 + sqrt(5) * w1.x2 / theta[d + 1] + 5 * w1.x2^2 / (3 * theta[d + 1]^2)),
        matrix(sqrt(5) / theta[d + 1] + 10 * w1.x2 / (3 * theta[d + 1]^2)),
        5 / (3 * theta[d + 1]^2)
      )
      
      predy <- mu2 + drop(crossprod(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5)) *
                                      (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                         (crossprod(t(e1), lambda11) * pnorm((mua - w1.x2) / sqrt(sig2)) +
                                            rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                                         exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                         (crossprod(t(e2), lambda21) * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                                            rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2)))), a))
      
      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 2.5, theta = theta[d + 1])
      }
      mat <- drop(crossprod(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5), cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))
      
      VE <- - (predy - mu2)^2 + drop(crossprod(a, crossprod(mat, a)))
    } else {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
      x.mu <- drop(pred.fit$mu) # mean of f1(u)
      sig2 <- pred.fit$sig2
      
      ### calculate the closed form ###
      X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
      w1.x2 <- fit2$X[, d + 1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat
      
      Ci <- fit2$Ki
      a <- crossprod(Ci, y2)
      
      # mean
      mua <- x.mu - sqrt(5) * sig2 / theta[d + 1]
      mub <- x.mu + sqrt(5) * sig2 / theta[d + 1]
      
      lambda11 <- c(1, mua, mua^2 + sig2)
      lambda12 <- cbind(0, 1, matrix(mua + w1.x2))
      lambda21 <- c(1, -mub, mub^2 + sig2)
      lambda22 <- cbind(0, 1, matrix(-mub - w1.x2))
      
      e1 <- cbind(
        matrix(1 - sqrt(5) * w1.x2 / theta[d + 1] + 5 * w1.x2^2 / (3 * theta[d + 1]^2)),
        matrix(sqrt(5) / theta[d + 1] - 10 * w1.x2 / (3 * theta[d + 1]^2)),
        5 / (3 * theta[d + 1]^2)
      )
      e2 <- cbind(
        matrix(1 + sqrt(5) * w1.x2 / theta[d + 1] + 5 * w1.x2^2 / (3 * theta[d + 1]^2)),
        matrix(sqrt(5) / theta[d + 1] + 10 * w1.x2 / (3 * theta[d + 1]^2)),
        5 / (3 * theta[d + 1]^2)
      )
      
      predy <- drop(crossprod(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5)) *
                                (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e1), lambda11) * pnorm((mua - w1.x2) / sqrt(sig2)) +
                                      rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                                   exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e2), lambda21) * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                                      rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2)))), a))
      
      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 2.5, theta = theta[d + 1])
      }
      mat <- drop(crossprod(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5), cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))
      
      VE <- - predy^2 + drop(crossprod(a, crossprod(mat, a)))
    }
  }
  return(-VE) # to maximize the V1.
}


#' object to optimize the point by ALD criterion updating V2 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param fit an object of class RNAmf.
#' @return A negative V2 at Xcand.
#' @noRd
#'

obj.ALD_V2_2level <- function(fit, Xcand) { # high
  newx <- matrix(Xcand, nrow = 1)
  
  kernel <- fit$kernel
  constant <- fit$constant
  fit1 <- fit$fits[[1]]
  fit2 <- fit$fits[[2]]
  
  if (kernel == "sqex") {
    if (constant) {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.GP(fit1, newx)
      x.mu <- pred.fit$mu # mean of f1(u)
      sig2 <- pred.fit$sig2 #* 0
      
      ### calculate the closed form ###
      X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
      w1.x2 <- fit2$X[, d + 1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat
      mu2 <- fit2$mu.hat
      
      Ci <- fit2$Ki
      a <- crossprod(Ci, y2 - mu2)
      
      # mean
      predy <- mu2 + crossprod(t(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))) *
                                   1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                                   exp(-(drop(outer(x.mu, w1.x2, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))), a)
      
      # var
      mat <- drop(crossprod(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))),
                            exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))))) * # common components
        1 / sqrt(1 + 4 * sig2 / theta[d + 1]) *
        exp(-(outer(w1.x2, w1.x2, FUN = "+") / 2 - matrix(x.mu, n, n))^2 / (theta[d + 1] / 2 + 2 * sig2)) *
        exp(-(outer(w1.x2, w1.x2, FUN = "-"))^2 / (2 * theta[d + 1]))
      
      EV <- tau2hat - tau2hat * sum(Ci * mat)
    } else {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.GP(fit1, newx)
      x.mu <- pred.fit$mu # mean of f1(u)
      sig2 <- pred.fit$sig2
      
      ### calculate the closed form ###
      X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
      w1.x2 <- fit2$X[, d + 1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat
      
      Ci <- fit2$Ki
      a <- crossprod(Ci, y2)
      
      # mean
      predy <- crossprod(t(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))) *
                             1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                             exp(-(drop(outer(x.mu, w1.x2, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))), a)
      
      # var
      mat <- drop(crossprod(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))),
                            exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))))) * # common components
        1 / sqrt(1 + 4 * sig2 / theta[d + 1]) *
        exp(-(outer(w1.x2, w1.x2, FUN = "+") / 2 - matrix(x.mu, n, n))^2 / (theta[d + 1] / 2 + 2 * sig2)) *
        exp(-(outer(w1.x2, w1.x2, FUN = "-"))^2 / (2 * theta[d + 1]))
      
      EV <- tau2hat - tau2hat * sum(Ci * mat)
    }
  } else if (kernel == "matern1.5") {
    if (constant) {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
      x.mu <- drop(pred.fit$mu) # mean of f1(u)
      sig2 <- pred.fit$sig2 #* 0
      
      ### calculate the closed form ###
      X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
      w1.x2 <- fit2$X[, d + 1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat
      mu2 <- fit2$mu.hat
      
      Ci <- fit2$Ki
      a <- crossprod(Ci, y2 - mu2)
      
      # mean
      mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
      mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]
      
      lambda11 <- c(1, mua)
      lambda12 <- c(0, 1)
      lambda21 <- c(1, -mub)
      
      e1 <- cbind(matrix(1 - sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])
      e2 <- cbind(matrix(1 + sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])
      
      predy <- mu2 + drop(crossprod(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                                      (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                         (crossprod(t(e1), lambda11) * pnorm((mua - w1.x2) / sqrt(sig2)) +
                                            crossprod(t(e1), lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                                         exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                         (crossprod(t(e2), lambda21) * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                                            crossprod(t(e2), lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2)))), a))
      
      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 1.5, theta = theta[d + 1])
      }
      mat <- drop(crossprod(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5), cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))
      
      EV <- tau2hat - tau2hat * sum(Ci * mat)
    } else {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
      x.mu <- drop(pred.fit$mu) # mean of f1(u)
      sig2 <- pred.fit$sig2
      
      ### calculate the closed form ###
      X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
      w1.x2 <- fit2$X[, d + 1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat
      
      Ci <- fit2$Ki
      a <- crossprod(Ci, y2)
      
      # mean
      mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
      mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]
      
      lambda11 <- c(1, mua)
      lambda12 <- c(0, 1)
      lambda21 <- c(1, -mub)
      
      e1 <- cbind(matrix(1 - sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])
      e2 <- cbind(matrix(1 + sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])
      
      predy <- drop(crossprod(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                                (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e1), lambda11) * pnorm((mua - w1.x2) / sqrt(sig2)) +
                                      crossprod(t(e1), lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                                   exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e2), lambda21) * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                                      crossprod(t(e2), lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2)))), a))
      
      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 1.5, theta = theta[d + 1])
      }
      mat <- drop(crossprod(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5), cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))
      
      EV <- tau2hat - tau2hat * sum(Ci * mat)
    }
  } else if (kernel == "matern2.5") {
    if (constant) {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
      x.mu <- drop(pred.fit$mu) # mean of f1(u)
      sig2 <- pred.fit$sig2 #* 0
      
      ### calculate the closed form ###
      X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
      w1.x2 <- fit2$X[, d + 1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat
      mu2 <- fit2$mu.hat
      
      Ci <- fit2$Ki
      a <- crossprod(Ci, y2 - mu2)
      
      # mean
      mua <- x.mu - sqrt(5) * sig2 / theta[d + 1]
      mub <- x.mu + sqrt(5) * sig2 / theta[d + 1]
      
      lambda11 <- c(1, mua, mua^2 + sig2)
      lambda12 <- cbind(0, 1, matrix(mua + w1.x2))
      lambda21 <- c(1, -mub, mub^2 + sig2)
      lambda22 <- cbind(0, 1, matrix(-mub - w1.x2))
      
      e1 <- cbind(
        matrix(1 - sqrt(5) * w1.x2 / theta[d + 1] + 5 * w1.x2^2 / (3 * theta[d + 1]^2)),
        matrix(sqrt(5) / theta[d + 1] - 10 * w1.x2 / (3 * theta[d + 1]^2)),
        5 / (3 * theta[d + 1]^2)
      )
      e2 <- cbind(
        matrix(1 + sqrt(5) * w1.x2 / theta[d + 1] + 5 * w1.x2^2 / (3 * theta[d + 1]^2)),
        matrix(sqrt(5) / theta[d + 1] + 10 * w1.x2 / (3 * theta[d + 1]^2)),
        5 / (3 * theta[d + 1]^2)
      )
      
      predy <- mu2 + drop(crossprod(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5)) *
                                      (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                         (crossprod(t(e1), lambda11) * pnorm((mua - w1.x2) / sqrt(sig2)) +
                                            rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                                         exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                         (crossprod(t(e2), lambda21) * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                                            rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2)))), a))
      
      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 2.5, theta = theta[d + 1])
      }
      mat <- drop(crossprod(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5), cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))
      
      EV <- tau2hat - tau2hat * sum(Ci * mat)
    } else {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
      x.mu <- drop(pred.fit$mu) # mean of f1(u)
      sig2 <- pred.fit$sig2
      
      ### calculate the closed form ###
      X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
      w1.x2 <- fit2$X[, d + 1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat
      
      Ci <- fit2$Ki
      a <- crossprod(Ci, y2)
      
      # mean
      mua <- x.mu - sqrt(5) * sig2 / theta[d + 1]
      mub <- x.mu + sqrt(5) * sig2 / theta[d + 1]
      
      lambda11 <- c(1, mua, mua^2 + sig2)
      lambda12 <- cbind(0, 1, matrix(mua + w1.x2))
      lambda21 <- c(1, -mub, mub^2 + sig2)
      lambda22 <- cbind(0, 1, matrix(-mub - w1.x2))
      
      e1 <- cbind(
        matrix(1 - sqrt(5) * w1.x2 / theta[d + 1] + 5 * w1.x2^2 / (3 * theta[d + 1]^2)),
        matrix(sqrt(5) / theta[d + 1] - 10 * w1.x2 / (3 * theta[d + 1]^2)),
        5 / (3 * theta[d + 1]^2)
      )
      e2 <- cbind(
        matrix(1 + sqrt(5) * w1.x2 / theta[d + 1] + 5 * w1.x2^2 / (3 * theta[d + 1]^2)),
        matrix(sqrt(5) / theta[d + 1] + 10 * w1.x2 / (3 * theta[d + 1]^2)),
        5 / (3 * theta[d + 1]^2)
      )
      
      predy <- drop(crossprod(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5)) *
                                (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e1), lambda11) * pnorm((mua - w1.x2) / sqrt(sig2)) +
                                      rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                                   exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e2), lambda21) * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                                      rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2)))), a))
      
      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 2.5, theta = theta[d + 1])
      }
      mat <- drop(crossprod(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5), cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))
      
      EV <- tau2hat - tau2hat * sum(Ci * mat)
    }
  }
  return(-EV) # to maximize the V2.
}


#' object to optimize the point by ALD criterion updating V1 with three levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param fit an object of class RNAmf.
#' @param mc.sample a number of mc samples generated for this approach. Default is 100.
#' @param parallel logical indicating whether to run parallel or not. Default is FALSE.
#' @param ncore the number of core for parallel. Default is 1.
#' @return A negative V1 at Xcand.
#' @importFrom stats var
#' @noRd
#'

obj.ALD_V1_3level <- function(fit, Xcand, mc.sample, parallel = FALSE, ncore = 1) { # low
  
  kernel <- fit$kernel
  constant <- fit$constant
  fit1 <- fit$fits[[1]]
  fit2 <- fit$fits[[2]]
  fit3 <- fit$fits[[3]]
  
  d <- ncol(fit1$X)
  newx <- matrix(Xcand, nrow = 1)
  
  if (kernel == "sqex") {
    y1.sample <- rnorm(mc.sample, mean = pred.GP(fit1, newx)$mu, sd = sqrt(pred.GP(fit1, newx)$sig2))
  } else if (kernel == "matern1.5") {
    y1.sample <- rnorm(mc.sample, mean = pred.matGP(fit1, newx)$mu, sd = sqrt(pred.matGP(fit1, newx)$sig2))
  } else if (kernel == "matern2.5") {
    y1.sample <- rnorm(mc.sample, mean = pred.matGP(fit1, newx)$mu, sd = sqrt(pred.matGP(fit1, newx)$sig2))
  }
  
  ### V1 MC approximation ###
  if (parallel) {
    VEE.out <- foreach(i = 1:mc.sample, .combine = c) %dorng% {
      if (kernel == "sqex") {
        pred2 <- pred.GP(fit2, cbind(newx, y1.sample[i]))
      } else if (kernel == "matern1.5") {
        pred2 <- pred.matGP(fit2, cbind(newx, y1.sample[i]))
      } else if (kernel == "matern2.5") {
        pred2 <- pred.matGP(fit2, cbind(newx, y1.sample[i]))
      }
      x.mu <- pred2$mu
      sig2 <- pred2$sig2
      
      if (kernel == "sqex") {
        ### calculate the closed form ###
        X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
        w2.x3 <- fit3$X[, d + 1]
        y3 <- fit3$y
        n <- length(y3)
        theta <- fit3$theta
        tau2hat <- fit3$tau2hat
        Ci <- fit3$Ki
        
        if (constant) {
          mu3 <- fit3$mu.hat
        } else {
          mu3 <- 0
        }
        a <- crossprod(Ci, y3 - mu3)
        # mean
        VEE <- crossprod(t(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)])))) *
                             1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                             exp(-(drop(outer(x.mu, w2.x3, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))), a)
      } else if (kernel == "matern1.5") {
        ### calculate the closed form ###
        X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
        w2.x3 <- fit3$X[, d + 1]
        y3 <- fit3$y
        n <- length(y3)
        theta <- fit3$theta
        tau2hat <- fit3$tau2hat
        Ci <- fit3$Ki
        
        if (constant) {
          mu3 <- fit3$mu.hat
        } else {
          mu3 <- 0
        }
        a <- crossprod(Ci, y3 - mu3)
        # mean
        mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
        mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]
        
        lambda11 <- c(1, mua)
        lambda12 <- c(0, 1)
        lambda21 <- c(1, -mub)
        
        e1 <- cbind(matrix(1 - sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])
        e2 <- cbind(matrix(1 + sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])
        
        VEE <- drop(crossprod(t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                                (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e1), lambda11) * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                      crossprod(t(e1), lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                                   exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e2), lambda21) * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                      crossprod(t(e2), lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2)))), a))
      } else if (kernel == "matern2.5") {
        ### calculate the closed form ###
        X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
        w2.x3 <- fit3$X[, d + 1]
        y3 <- fit3$y
        n <- length(y3)
        theta <- fit3$theta
        tau2hat <- fit3$tau2hat
        Ci <- fit3$Ki
        
        if (constant) {
          mu3 <- fit3$mu.hat
        } else {
          mu3 <- 0
        }
        a <- crossprod(Ci, y3 - mu3)
        # mean
        mua <- x.mu - sqrt(5) * sig2 / theta[d + 1]
        mub <- x.mu + sqrt(5) * sig2 / theta[d + 1]
        
        lambda11 <- c(1, mua, mua^2 + sig2)
        lambda12 <- cbind(0, 1, matrix(mua + w2.x3))
        lambda21 <- c(1, -mub, mub^2 + sig2)
        lambda22 <- cbind(0, 1, matrix(-mub - w2.x3))
        
        e1 <- cbind(
          matrix(1 - sqrt(5) * w2.x3 / theta[d + 1] + 5 * w2.x3^2 / (3 * theta[d + 1]^2)),
          matrix(sqrt(5) / theta[d + 1] - 10 * w2.x3 / (3 * theta[d + 1]^2)),
          5 / (3 * theta[d + 1]^2)
        )
        e2 <- cbind(
          matrix(1 + sqrt(5) * w2.x3 / theta[d + 1] + 5 * w2.x3^2 / (3 * theta[d + 1]^2)),
          matrix(sqrt(5) / theta[d + 1] + 10 * w2.x3 / (3 * theta[d + 1]^2)),
          5 / (3 * theta[d + 1]^2)
        )
        
        VEE <- drop(crossprod(t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 2.5)) *
                                (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e1), lambda11) * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                      rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                                   exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e2), lambda21) * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                      rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2)))), a))
      }
      return(VEE) # to maximize the V1.
    }
    attr(VEE.out, "rng") <- NULL
    attr(VEE.out, "doRNG_version") <- NULL
  } else {
    VEE.out <- rep(0, mc.sample)
    for (i in 1:mc.sample) {
      if (kernel == "sqex") {
        pred2 <- pred.GP(fit2, cbind(newx, y1.sample[i]))
      } else if (kernel == "matern1.5") {
        pred2 <- pred.matGP(fit2, cbind(newx, y1.sample[i]))
      } else if (kernel == "matern2.5") {
        pred2 <- pred.matGP(fit2, cbind(newx, y1.sample[i]))
      }
      x.mu <- pred2$mu
      sig2 <- pred2$sig2
      
      if (kernel == "sqex") {
        ### calculate the closed form ###
        X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
        w2.x3 <- fit3$X[, d + 1]
        y3 <- fit3$y
        n <- length(y3)
        theta <- fit3$theta
        tau2hat <- fit3$tau2hat
        Ci <- fit3$Ki
        
        if (constant) {
          mu3 <- fit3$mu.hat
        } else {
          mu3 <- 0
        }
        a <- crossprod(Ci, y3 - mu3)
        # mean
        VEE <- crossprod(t(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)])))) *
                             1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                             exp(-(drop(outer(x.mu, w2.x3, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))), a)
      } else if (kernel == "matern1.5") {
        ### calculate the closed form ###
        X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
        w2.x3 <- fit3$X[, d + 1]
        y3 <- fit3$y
        n <- length(y3)
        theta <- fit3$theta
        tau2hat <- fit3$tau2hat
        Ci <- fit3$Ki
        
        if (constant) {
          mu3 <- fit3$mu.hat
        } else {
          mu3 <- 0
        }
        a <- crossprod(Ci, y3 - mu3)
        # mean
        mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
        mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]
        
        lambda11 <- c(1, mua)
        lambda12 <- c(0, 1)
        lambda21 <- c(1, -mub)
        
        e1 <- cbind(matrix(1 - sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])
        e2 <- cbind(matrix(1 + sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])
        
        VEE <- drop(crossprod(t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                                (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e1), lambda11) * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                      crossprod(t(e1), lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                                   exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e2), lambda21) * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                      crossprod(t(e2), lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2)))), a))
      } else if (kernel == "matern2.5") {
        ### calculate the closed form ###
        X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
        w2.x3 <- fit3$X[, d + 1]
        y3 <- fit3$y
        n <- length(y3)
        theta <- fit3$theta
        tau2hat <- fit3$tau2hat
        Ci <- fit3$Ki
        
        if (constant) {
          mu3 <- fit3$mu.hat
        } else {
          mu3 <- 0
        }
        a <- crossprod(Ci, y3 - mu3)
        # mean
        mua <- x.mu - sqrt(5) * sig2 / theta[d + 1]
        mub <- x.mu + sqrt(5) * sig2 / theta[d + 1]
        
        lambda11 <- c(1, mua, mua^2 + sig2)
        lambda12 <- cbind(0, 1, matrix(mua + w2.x3))
        lambda21 <- c(1, -mub, mub^2 + sig2)
        lambda22 <- cbind(0, 1, matrix(-mub - w2.x3))
        
        e1 <- cbind(
          matrix(1 - sqrt(5) * w2.x3 / theta[d + 1] + 5 * w2.x3^2 / (3 * theta[d + 1]^2)),
          matrix(sqrt(5) / theta[d + 1] - 10 * w2.x3 / (3 * theta[d + 1]^2)),
          5 / (3 * theta[d + 1]^2)
        )
        e2 <- cbind(
          matrix(1 + sqrt(5) * w2.x3 / theta[d + 1] + 5 * w2.x3^2 / (3 * theta[d + 1]^2)),
          matrix(sqrt(5) / theta[d + 1] + 10 * w2.x3 / (3 * theta[d + 1]^2)),
          5 / (3 * theta[d + 1]^2)
        )
        
        VEE <- drop(crossprod(t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 2.5)) *
                                (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e1), lambda11) * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                      rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                                   exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                   (crossprod(t(e2), lambda21) * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                      rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2)))), a))
      }
      VEE.out[i] <- VEE # to maximize the V1.
    }
  }
  return(-var(VEE.out)) # to maximize the V1.
}


#' object to optimize the point by ALD criterion updating V2 with three levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param fit an object of class RNAmf.
#' @param mc.sample a number of mc samples generated for this approach. Default is 100.
#' @param parallel logical indicating whether to run parallel or not. Default is FALSE.
#' @param ncore the number of core for parallel. Default is 1.
#' @return A negative V2 at Xcand.
#' @noRd
#'

obj.ALD_V2_3level <- function(fit, Xcand, mc.sample, parallel = FALSE, ncore = 1) { # med
  
  V <- predict(fit, Xcand)$sig2[[3]]
  EVE <- V +
    obj.ALD_V1_3level(fit, Xcand, mc.sample, parallel = FALSE, ncore = 1) +# -V1
    obj.ALD_V3_3level(fit, Xcand) # -V3
  
  return(-EVE) # to maximize the V2.
}


#' object to optimize the point by ALD criterion updating V3 with three levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param fit an object of class RNAmf.
#' @return A negative V2 at Xcand.
#' @noRd
#'

obj.ALD_V3_3level <- function(fit, Xcand) { # high
  newx <- matrix(Xcand, nrow = 1)
  
  kernel <- fit$kernel
  constant <- fit$constant
  fit1 <- fit$fits[[1]]
  fit2 <- fit$fits[[2]]
  fit3 <- fit$fits[[3]]
  
  if (kernel == "sqex") {
    pred.RNAmf_two_level <- predict(fit, newx)
    x.mu <- pred.RNAmf_two_level$mu[[2]]
    sig2 <- pred.RNAmf_two_level$sig2[[2]]
    
    d <- ncol(fit1$X)
    newx <- matrix(newx, ncol = d)
    
    ### calculate the closed form ###
    X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
    w2.x3 <- fit3$X[, d + 1]
    y3 <- fit3$y
    n <- length(y3)
    theta <- fit3$theta
    tau2hat <- fit3$tau2hat
    Ci <- fit3$Ki
    
    if (constant) {
      mu3 <- fit3$mu.hat
    } else {
      mu3 <- 0
    }
    a <- crossprod(Ci, y3 - mu3)
    
    # mean
    predy <- mu3 + crossprod(t(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)])))) *
                                 1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                                 exp(-(drop(outer(x.mu, w2.x3, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))), a)
    # var
    mat <- drop(crossprod(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)])))),
                          exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)])))))) * # common components
      1 / sqrt(1 + 4 * sig2 / theta[d + 1]) *
      exp(-(outer(w2.x3, w2.x3, FUN = "+") / 2 - matrix(x.mu, n, n))^2 / (theta[d + 1] / 2 + 2 * sig2)) *
      exp(-(outer(w2.x3, w2.x3, FUN = "-"))^2 / (2 * theta[d + 1]))
    
    EEV <- tau2hat - tau2hat * sum(Ci * mat)
  } else if (kernel == "matern1.5") {
    d <- ncol(fit1$X)
    newx <- matrix(newx, ncol = d)
    pred.RNAmf_two_level <- predict(fit, newx)
    x.mu <- pred.RNAmf_two_level$mu[[2]]
    sig2 <- pred.RNAmf_two_level$sig2[[2]]
    
    ### calculate the closed form ###
    X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
    w2.x3 <- fit3$X[, d + 1]
    y3 <- fit3$y
    n <- length(y3)
    theta <- fit3$theta
    tau2hat <- fit3$tau2hat
    Ci <- fit3$Ki
    
    if (constant) {
      mu3 <- fit3$mu.hat
    } else {
      mu3 <- 0
    }
    a <- crossprod(Ci, y3 - mu3)
    
    # mean
    mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
    mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]
    
    lambda11 <- c(1, mua)
    lambda12 <- c(0, 1)
    lambda21 <- c(1, -mub)
    
    e1 <- cbind(matrix(1 - sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])
    e2 <- cbind(matrix(1 + sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])
    
    predy <- mu3 + drop(crossprod(t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                                    (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                       (crossprod(t(e1), lambda11) * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                          crossprod(t(e1), lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                                       exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                       (crossprod(t(e2), lambda21) * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                          crossprod(t(e2), lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2)))), a))
    # var
    zeta <- function(x, y) {
      zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 1.5, theta = theta[d + 1])
    }
    mat <- drop(crossprod(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 1.5), cor.sep(t(newx), X3, theta[-(d + 1)], nu = 1.5)))* # constant depends on kernel
      outer(w2.x3, w2.x3, FUN = Vectorize(zeta))
    
    EEV <- tau2hat - tau2hat * sum(Ci * mat)
  } else if (kernel == "matern2.5") {
    d <- ncol(fit1$X)
    newx <- matrix(newx, ncol = d)
    pred.RNAmf_two_level <- predict(fit, newx)
    x.mu <- pred.RNAmf_two_level$mu[[2]]
    sig2 <- pred.RNAmf_two_level$sig2[[2]]
    
    ### calculate the closed form ###
    X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
    w2.x3 <- fit3$X[, d + 1]
    y3 <- fit3$y
    n <- length(y3)
    theta <- fit3$theta
    tau2hat <- fit3$tau2hat
    Ci <- fit3$Ki
    
    if (constant) {
      mu3 <- fit3$mu.hat
    } else {
      mu3 <- 0
    }
    a <- crossprod(Ci, y3 - mu3)
    
    # mean
    mua <- x.mu - sqrt(5) * sig2 / theta[d + 1]
    mub <- x.mu + sqrt(5) * sig2 / theta[d + 1]
    
    lambda11 <- c(1, mua, mua^2 + sig2)
    lambda12 <- cbind(0, 1, matrix(mua + w2.x3))
    lambda21 <- c(1, -mub, mub^2 + sig2)
    lambda22 <- cbind(0, 1, matrix(-mub - w2.x3))
    
    e1 <- cbind(
      matrix(1 - sqrt(5) * w2.x3 / theta[d + 1] + 5 * w2.x3^2 / (3 * theta[d + 1]^2)),
      matrix(sqrt(5) / theta[d + 1] - 10 * w2.x3 / (3 * theta[d + 1]^2)),
      5 / (3 * theta[d + 1]^2)
    )
    e2 <- cbind(
      matrix(1 + sqrt(5) * w2.x3 / theta[d + 1] + 5 * w2.x3^2 / (3 * theta[d + 1]^2)),
      matrix(sqrt(5) / theta[d + 1] + 10 * w2.x3 / (3 * theta[d + 1]^2)),
      5 / (3 * theta[d + 1]^2)
    )
    
    predy <- mu3 + drop(crossprod(t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 2.5)) *
                                    (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                       (crossprod(t(e1), lambda11) * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                          rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                                       exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                       (crossprod(t(e2), lambda21) * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                          rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2)))), a))
    
    # var
    zeta <- function(x, y) {
      zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 2.5, theta = theta[d + 1])
    }
    mat <- drop(crossprod(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 2.5), cor.sep(t(newx), X3, theta[-(d + 1)], nu = 2.5))) * # constant depends on kernel
      outer(w2.x3, w2.x3, FUN = Vectorize(zeta))
    
    EEV <- tau2hat - tau2hat * sum(Ci * mat)
  }
  return(-EEV) # to maximize the V3.
}


#' object to optimize the next point by ALC criterion updating at level 1 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param Xref vector or matrix of reference data.
#' @param fit an object of class RNAmf.
#' @param mc.sample a number of mc samples generated for this approach. Default is 100.
#' @param parallel logical indicating whether to run parallel or not. Default is FALSE.
#' @param ncore the number of core for parallel. Default is 1.
#' @return A mean of the deduced variance at Xref.
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm
#' @noRd
#'

obj.ALC_level_1 <- function(Xcand, Xref, fit, MC = FALSE, mc.sample, parallel = FALSE, ncore = 1) {
  fit1 <- f1 <- fit$fits[[1]]
  fit2 <- f2 <- fit$fits[[2]]
  constant <- fit$constant
  kernel <- fit$kernel
  g <- fit1$g
  
  Xcand <- matrix(Xcand, nrow = 1)
  if(MC){
    if (kernel == "sqex") {
      y1.sample <- rnorm(mc.sample, mean = pred.GP(f1, Xcand)$mu, sd = sqrt(pred.GP(f1, Xcand)$sig2))
    } else if (kernel == "matern1.5") {
      y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
    } else if (kernel == "matern2.5") {
      y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
    }
  }else{
    if (kernel == "sqex") {
      y1.sample <- pred.GP(f1, Xcand)$mu
    } else if (kernel == "matern1.5") {
      y1.sample <- pred.matGP(f1, Xcand)$mu
    } else if (kernel == "matern2.5") {
      y1.sample <- pred.matGP(f1, Xcand)$mu
    }
  }
  Xcand1 <- Xcand
  
  ### Choose level 1 ###
  ### update Ki1
  if (kernel == "sqex") {
    cov.newx1 <- covar.sep(X1 = Xcand1, d = f1$theta, g = g)
    cov.Xnewx1 <- covar.sep(X1 = f1$X, X2 = Xcand1, d = f1$theta, g = 0)
  } else if (kernel == "matern1.5") {
    cov.newx1 <- cor.sep(X = Xcand1, theta = f1$theta, nu = 1.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand1, theta = f1$theta, nu = 1.5)
  } else if (kernel == "matern2.5") {
    cov.newx1 <- cor.sep(X = Xcand1, theta = f1$theta, nu = 2.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand1, theta = f1$theta, nu = 2.5)
  }
  v.next1 <- drop(cov.newx1 - crossprod(t(crossprod(cov.Xnewx1, f1$Ki)), cov.Xnewx1))
  g.next1 <- -1 / drop(v.next1) * crossprod(f1$Ki, cov.Xnewx1)
  
  
  fit1$Ki <- rbind(
    cbind(f1$Ki + crossprod(t(g.next1), t(g.next1) * v.next1), g.next1),
    cbind(t(g.next1), 1 / drop(v.next1))
  )
  
  fit1$X <- rbind(f1$X, Xcand1)
  
  if(MC){
    ALC.out <- rep(0, mc.sample)
    if (parallel) {
      ALC.out <- foreach(i = 1:mc.sample, .combine = c) %dorng% {
        fit.tmp <- fit
        fit1$y <- c(f1$y, y1.sample[i])
        fit1$tau2hat <- drop(crossprod(t(crossprod(fit1$y - fit1$mu.hat, fit1$Ki)), (fit1$y - fit1$mu.hat)) / length(fit1$y))
        
        fit.tmp$fits[[1]] <- fit1
        
        return(mean(predict(fit.tmp, Xref)$sig2[[2]])) # to minimize the deduced variance. To maximize, -mean
      }
      attr(ALC.out, "rng") <- NULL
      attr(ALC.out, "doRNG_version") <- NULL
    } else {
      for (i in 1:mc.sample) {
        fit.tmp <- fit
        fit1$y <- c(f1$y, y1.sample[i])
        fit1$tau2hat <- drop(crossprod(t(crossprod(fit1$y - fit1$mu.hat, fit1$Ki)), (fit1$y - fit1$mu.hat)) / length(fit1$y))
        
        fit.tmp$fits[[1]] <- fit1
        
        ALC.out[i] <- mean(predict(fit.tmp, Xref)$sig2[[2]]) # to minimize the deduced variance. To maximize, -mean
      }
    }
  }else{
    fit.tmp <- fit
    fit1$y <- c(f1$y, y1.sample)
    fit1$tau2hat <- drop(crossprod(t(crossprod(fit1$y - fit1$mu.hat, fit1$Ki)), (fit1$y - fit1$mu.hat)) / length(fit1$y))
    
    fit.tmp$fits[[1]] <- fit1
    
    ALC.out <- mean(predict(fit.tmp, Xref)$sig2[[2]])
  }
  return(mean(ALC.out))
}


#' object to optimize the next point by ALC criterion updating at level 2 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param Xref vector or matrix of reference data.
#' @param fit an object of class RNAmf.
#' @param mc.sample a number of mc samples generated for this approach. Default is 100.
#' @param parallel logical indicating whether to run parallel or not. Default is FALSE.
#' @param ncore the number of core for parallel. Default is 1.
#' @return A mean of the deduced variance at Xref.
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm
#' @noRd
#'

obj.ALC_level_2 <- function(Xcand, Xref, fit, MC = FALSE, mc.sample, parallel = FALSE, ncore = 1) {
  fit1 <- f1 <- fit$fits[[1]]
  fit2 <- f2 <- fit$fits[[2]]
  constant <- fit$constant
  kernel <- fit$kernel
  g <- fit1$g
  
  Xcand <- matrix(Xcand, nrow = 1)
  if(MC){
    if (kernel == "sqex") {
      y1.sample <- rnorm(mc.sample, mean = pred.GP(f1, Xcand)$mu, sd = sqrt(pred.GP(f1, Xcand)$sig2))
    } else if (kernel == "matern1.5") {
      y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
    } else if (kernel == "matern2.5") {
      y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
    }
  }else{
    if (kernel == "sqex") {
      y1.sample <- pred.GP(f1, Xcand)$mu
    } else if (kernel == "matern1.5") {
      y1.sample <- pred.matGP(f1, Xcand)$mu
    } else if (kernel == "matern2.5") {
      y1.sample <- pred.matGP(f1, Xcand)$mu
    }
  }
  Xcand1 <- Xcand
  
  ### Choose level 1 ###
  ### update Ki1
  if (kernel == "sqex") {
    cov.newx1 <- covar.sep(X1 = Xcand1, d = f1$theta, g = g)
    cov.Xnewx1 <- covar.sep(X1 = f1$X, X2 = Xcand1, d = f1$theta, g = 0)
  } else if (kernel == "matern1.5") {
    cov.newx1 <- cor.sep(X = Xcand1, theta = f1$theta, nu = 1.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand1, theta = f1$theta, nu = 1.5)
  } else if (kernel == "matern2.5") {
    cov.newx1 <- cor.sep(X = Xcand1, theta = f1$theta, nu = 2.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand1, theta = f1$theta, nu = 2.5)
  }
  v.next1 <- drop(cov.newx1 - crossprod(t(crossprod(cov.Xnewx1, f1$Ki)), cov.Xnewx1))
  g.next1 <- -1 / drop(v.next1) * crossprod(f1$Ki, cov.Xnewx1)
  
  
  fit1$Ki <- rbind(
    cbind(f1$Ki + crossprod(t(g.next1), t(g.next1) * v.next1), g.next1),
    cbind(t(g.next1), 1 / drop(v.next1))
  )
  
  fit1$X <- rbind(f1$X, Xcand1)
  
  if(MC){
    if (parallel) {
      ALC.out <- foreach(i = 1:mc.sample, .combine = c) %dorng% {
        fit.tmp <- fit
        fit1$y <- c(f1$y, y1.sample[i])
        fit1$tau2hat <- drop(crossprod(t(crossprod(fit1$y - fit1$mu.hat, fit1$Ki)), (fit1$y - fit1$mu.hat)) / length(fit1$y))
        
        fit.tmp$fits[[1]] <- fit1
        
        ### Choose level 2 ###
        if (kernel == "sqex") {
          pred2 <- pred.GP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        } else if (kernel == "matern1.5") {
          pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        } else if (kernel == "matern2.5") {
          pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        }
        
        ### update Ki2
        newx2 <- cbind(Xcand, y1.sample[i])
        
        if (kernel == "sqex") {
          cov.newx2 <- covar.sep(X1 = newx2, d = f2$theta, g = g)
          cov.Xnewx2 <- covar.sep(X1 = f2$X, X2 = newx2, d = f2$theta, g = 0)
        } else if (kernel == "matern1.5") {
          cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 1.5)
          cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 1.5)
        } else if (kernel == "matern2.5") {
          cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 2.5)
          cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 2.5)
        }
        v.next2 <- drop(cov.newx2 - crossprod(t(crossprod(cov.Xnewx2, f2$Ki)), cov.Xnewx2))
        g.next2 <- -1 / drop(v.next2) * crossprod(f2$Ki, cov.Xnewx2)
        
        fit2$Ki <- rbind(
          cbind(f2$Ki + crossprod(t(g.next2), t(g.next2) * v.next2), g.next2),
          cbind(t(g.next2), 1 / drop(v.next2))
        )
        
        fit2$X <- rbind(f2$X, newx2)
        fit2$y <- c(f2$y, y2.sample)
        fit2$tau2hat <- drop(crossprod(t(crossprod(fit2$y - fit2$mu.hat, fit2$Ki)), (fit2$y - fit2$mu.hat)) / length(fit2$y))
        
        fit.tmp$fits[[2]] <- fit2
        
        return(mean(predict(fit.tmp, Xref)$sig2[[2]])) # to minimize the deduced variance. To maximize, -mean
      }
      attr(ALC.out, "rng") <- NULL
      attr(ALC.out, "doRNG_version") <- NULL
    } else {
      ALC.out <- rep(0, mc.sample)
      for (i in 1:mc.sample) {
        fit.tmp <- fit
        fit1$y <- c(f1$y, y1.sample[i])
        fit1$tau2hat <- drop(crossprod(t(crossprod(fit1$y - fit1$mu.hat, fit1$Ki)), (fit1$y - fit1$mu.hat)) / length(fit1$y))
        
        fit.tmp$fits[[1]] <- fit1
        
        ### Choose level 2 ###
        if (kernel == "sqex") {
          pred2 <- pred.GP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        } else if (kernel == "matern1.5") {
          pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        } else if (kernel == "matern2.5") {
          pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        }
        
        ### update Ki2
        newx2 <- cbind(Xcand, y1.sample[i])
        
        if (kernel == "sqex") {
          cov.newx2 <- covar.sep(X1 = newx2, d = f2$theta, g = g)
          cov.Xnewx2 <- covar.sep(X1 = f2$X, X2 = newx2, d = f2$theta, g = 0)
        } else if (kernel == "matern1.5") {
          cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 1.5)
          cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 1.5)
        } else if (kernel == "matern2.5") {
          cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 2.5)
          cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 2.5)
        }
        v.next2 <- drop(cov.newx2 - crossprod(t(crossprod(cov.Xnewx2, f2$Ki)), cov.Xnewx2))
        g.next2 <- -1 / drop(v.next2) * crossprod(f2$Ki, cov.Xnewx2)
        
        fit2$Ki <- rbind(
          cbind(f2$Ki + crossprod(t(g.next2), t(g.next2) * v.next2), g.next2),
          cbind(t(g.next2), 1 / drop(v.next2))
        )
        
        fit2$X <- rbind(f2$X, newx2)
        fit2$y <- c(f2$y, y2.sample)
        fit2$tau2hat <- drop(crossprod(t(crossprod(fit2$y - fit2$mu.hat, fit2$Ki)), (fit2$y - fit2$mu.hat)) / length(fit2$y))
        
        fit.tmp$fits[[2]] <- fit2
        
        ALC.out[i] <- mean(predict(fit.tmp, Xref)$sig2[[2]]) # to minimize the deduced variance. To maximize, -mean
      }
    }
  }else{
    fit.tmp <- fit
    fit1$y <- c(f1$y, y1.sample)
    fit1$tau2hat <- drop(crossprod(t(crossprod(fit1$y - fit1$mu.hat, fit1$Ki)), (fit1$y - fit1$mu.hat)) / length(fit1$y))
    
    fit.tmp$fits[[1]] <- fit1
    
    ### Choose level 2 ###
    if (kernel == "sqex") {
      pred2 <- pred.GP(fit2, cbind(Xcand, y1.sample))
    } else if (kernel == "matern1.5") {
      pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample))
    } else if (kernel == "matern2.5") {
      pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample))
    }
    y2.sample <- pred2$mu
    
    ### update Ki2
    newx2 <- cbind(Xcand, y1.sample)
    
    if (kernel == "sqex") {
      cov.newx2 <- covar.sep(X1 = newx2, d = f2$theta, g = g)
      cov.Xnewx2 <- covar.sep(X1 = f2$X, X2 = newx2, d = f2$theta, g = 0)
    } else if (kernel == "matern1.5") {
      cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 1.5)
      cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 1.5)
    } else if (kernel == "matern2.5") {
      cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 2.5)
      cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 2.5)
    }
    v.next2 <- drop(cov.newx2 - crossprod(t(crossprod(cov.Xnewx2, f2$Ki)), cov.Xnewx2))
    g.next2 <- -1 / drop(v.next2) * crossprod(f2$Ki, cov.Xnewx2)
    
    fit2$Ki <- rbind(
      cbind(f2$Ki + crossprod(t(g.next2), t(g.next2) * v.next2), g.next2),
      cbind(t(g.next2), 1 / drop(v.next2))
    )
    
    fit2$X <- rbind(f2$X, newx2)
    fit2$y <- c(f2$y, y2.sample)
    fit2$tau2hat <- drop(crossprod(t(crossprod(fit2$y - fit2$mu.hat, fit2$Ki)), (fit2$y - fit2$mu.hat)) / length(fit2$y))
    
    fit.tmp$fits[[2]] <- fit2
    
    ALC.out <- mean(predict(fit.tmp, Xref)$sig2[[2]])
  }
  return(mean(ALC.out))
}


#' object to optimize the next point by ALC criterion updating at level 1 with three levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param Xref vector or matrix of reference data.
#' @param fit an object of class RNAmf.
#' @param mc.sample a number of mc samples generated for this approach. Default is 100.
#' @param parallel logical indicating whether to run parallel or not. Default is FALSE.
#' @param ncore the number of core for parallel. Default is 1.
#' @return A mean of the deduced variance at Xref.
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm
#' @noRd
#'

obj.ALC_level_3_1 <- function(Xcand, Xref, fit, MC=FALSE, mc.sample, parallel = FALSE, ncore = 1) {
  fit1 <- f1 <- fit$fits[[1]]
  fit2 <- f2 <- fit$fits[[2]]
  fit3 <- f3 <- fit$fits[[3]]
  constant <- fit$constant
  kernel <- fit$kernel
  g <- fit1$g
  
  Xcand <- matrix(Xcand, nrow = 1)
  if(MC){
    if (kernel == "sqex") {
      y1.sample <- rnorm(mc.sample, mean = pred.GP(f1, Xcand)$mu, sd = sqrt(pred.GP(f1, Xcand)$sig2))
    } else if (kernel == "matern1.5") {
      y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
    } else if (kernel == "matern2.5") {
      y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
    }
  }else{
    if (kernel == "sqex") {
      y1.sample <- pred.GP(f1, Xcand)$mu
    } else if (kernel == "matern1.5") {
      y1.sample <- pred.matGP(f1, Xcand)$mu
    } else if (kernel == "matern2.5") {
      y1.sample <- pred.matGP(f1, Xcand)$mu
    }
  }
  Xcand1 <- Xcand
  
  ### Choose level 1 ###
  ### update Ki1
  if (kernel == "sqex") {
    cov.newx1 <- covar.sep(X1 = Xcand1, d = f1$theta, g = g)
    cov.Xnewx1 <- covar.sep(X1 = f1$X, X2 = Xcand1, d = f1$theta, g = 0)
  } else if (kernel == "matern1.5") {
    cov.newx1 <- cor.sep(X = Xcand1, theta = f1$theta, nu = 1.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand1, theta = f1$theta, nu = 1.5)
  } else if (kernel == "matern2.5") {
    cov.newx1 <- cor.sep(X = Xcand1, theta = f1$theta, nu = 2.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand1, theta = f1$theta, nu = 2.5)
  }
  v.next1 <- drop(cov.newx1 - crossprod(t(crossprod(cov.Xnewx1, f1$Ki)), cov.Xnewx1))
  g.next1 <- -1 / drop(v.next1) * crossprod(f1$Ki, cov.Xnewx1)
  
  
  fit1$Ki <- rbind(
    cbind(f1$Ki + crossprod(t(g.next1), t(g.next1) * v.next1), g.next1),
    cbind(t(g.next1), 1 / drop(v.next1))
  )
  
  fit1$X <- rbind(f1$X, Xcand1)
  
  if(MC){
    if (parallel) {
      ALC.out <- foreach(i = 1:mc.sample, .combine = c) %dorng% {
        fit.tmp <- fit
        fit1$y <- c(f1$y, y1.sample[i])
        fit1$tau2hat <- drop(crossprod(t(crossprod(fit1$y - fit1$mu.hat, fit1$Ki)), (fit1$y - fit1$mu.hat)) / length(fit1$y))
        
        fit.tmp$fits[[1]] <- fit1
        
        return(mean(predict(fit.tmp, Xref)$sig2[[3]])) # to minimize the deduced variance. To maximize, -mean
      }
      attr(ALC.out, "rng") <- NULL
      attr(ALC.out, "doRNG_version") <- NULL
    } else {
      ALC.out <- rep(0, mc.sample)
      for (i in 1:mc.sample) {
        fit.tmp <- fit
        fit1$y <- c(f1$y, y1.sample[i])
        fit1$tau2hat <- drop(crossprod(t(crossprod(fit1$y - fit1$mu.hat, fit1$Ki)), (fit1$y - fit1$mu.hat)) / length(fit1$y))
        
        fit.tmp$fits[[1]] <- fit1
        
        ALC.out[i] <- mean(predict(fit.tmp, Xref)$sig2[[3]]) # to minimize the deduced variance. To maximize, -mean
      }
    }
  }else{
    fit.tmp <- fit
    fit1$y <- c(f1$y, y1.sample)
    fit1$tau2hat <- drop(crossprod(t(crossprod(fit1$y - fit1$mu.hat, fit1$Ki)), (fit1$y - fit1$mu.hat)) / length(fit1$y))
    
    fit.tmp$fits[[1]] <- fit1
    
    ALC.out <- mean(predict(fit.tmp, Xref)$sig2[[3]])
  }
  return(mean(ALC.out))
}


#' object to optimize the next point by ALC criterion updating at level 2 with three levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param Xref vector or matrix of reference data.
#' @param fit an object of class RNAmf.
#' @param mc.sample a number of mc samples generated for this approach. Default is 100.
#' @param parallel logical indicating whether to run parallel or not. Default is FALSE.
#' @param ncore the number of core for parallel. Default is 1.
#' @return A mean of the deduced variance at Xref.
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm
#' @noRd
#'

obj.ALC_level_3_2 <- function(Xcand, Xref, fit, MC=FALSE, mc.sample, parallel = FALSE, ncore = 1) {
  fit1 <- f1 <- fit$fits[[1]]
  fit2 <- f2 <- fit$fits[[2]]
  fit3 <- f3 <- fit$fits[[3]]
  constant <- fit$constant
  kernel <- fit$kernel
  g <- fit1$g
  
  Xcand <- matrix(Xcand, nrow = 1)
  if(MC){
    if (kernel == "sqex") {
      y1.sample <- rnorm(mc.sample, mean = pred.GP(f1, Xcand)$mu, sd = sqrt(pred.GP(f1, Xcand)$sig2))
    } else if (kernel == "matern1.5") {
      y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
    } else if (kernel == "matern2.5") {
      y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
    }
  }else{
    if (kernel == "sqex") {
      y1.sample <- pred.GP(f1, Xcand)$mu
    } else if (kernel == "matern1.5") {
      y1.sample <- pred.matGP(f1, Xcand)$mu
    } else if (kernel == "matern2.5") {
      y1.sample <- pred.matGP(f1, Xcand)$mu
    }
  }
  Xcand1 <- Xcand
  
  ### Choose level 1 ###
  ### update Ki1
  if (kernel == "sqex") {
    cov.newx1 <- covar.sep(X1 = Xcand1, d = f1$theta, g = g)
    cov.Xnewx1 <- covar.sep(X1 = f1$X, X2 = Xcand1, d = f1$theta, g = 0)
  } else if (kernel == "matern1.5") {
    cov.newx1 <- cor.sep(X = Xcand1, theta = f1$theta, nu = 1.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand1, theta = f1$theta, nu = 1.5)
  } else if (kernel == "matern2.5") {
    cov.newx1 <- cor.sep(X = Xcand1, theta = f1$theta, nu = 2.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand1, theta = f1$theta, nu = 2.5)
  }
  v.next1 <- drop(cov.newx1 - crossprod(t(crossprod(cov.Xnewx1, f1$Ki)), cov.Xnewx1))
  g.next1 <- -1 / drop(v.next1) * crossprod(f1$Ki, cov.Xnewx1)
  
  
  fit1$Ki <- rbind(
    cbind(f1$Ki + crossprod(t(g.next1), t(g.next1) * v.next1), g.next1),
    cbind(t(g.next1), 1 / drop(v.next1))
  )
  
  fit1$X <- rbind(f1$X, Xcand1)
  
  if(MC){
    if (parallel) {
      ALC.out <- foreach(i = 1:mc.sample, .combine = c) %dorng% {
        fit.tmp <- fit
        fit1$y <- c(f1$y, y1.sample[i])
        fit1$tau2hat <- drop(crossprod(t(crossprod(fit1$y - fit1$mu.hat, fit1$Ki)), (fit1$y - fit1$mu.hat)) / length(fit1$y))
        
        fit.tmp$fits[[1]] <- fit1
        
        ### Choose level 2 ###
        if (kernel == "sqex") {
          pred2 <- pred.GP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        } else if (kernel == "matern1.5") {
          pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        } else if (kernel == "matern2.5") {
          pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        }
        
        ### update Ki2
        newx2 <- cbind(Xcand, y1.sample[i])
        
        if (kernel == "sqex") {
          cov.newx2 <- covar.sep(X1 = newx2, d = f2$theta, g = g)
          cov.Xnewx2 <- covar.sep(X1 = f2$X, X2 = newx2, d = f2$theta, g = 0)
        } else if (kernel == "matern1.5") {
          cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 1.5)
          cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 1.5)
        } else if (kernel == "matern2.5") {
          cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 2.5)
          cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 2.5)
        }
        v.next2 <- drop(cov.newx2 - crossprod(t(crossprod(cov.Xnewx2, f2$Ki)), cov.Xnewx2))
        g.next2 <- -1 / drop(v.next2) * crossprod(f2$Ki, cov.Xnewx2)
        
        fit2$Ki <- rbind(
          cbind(f2$Ki + crossprod(t(g.next2), t(g.next2) * v.next2), g.next2),
          cbind(t(g.next2), 1 / drop(v.next2))
        )
        
        fit2$X <- rbind(f2$X, newx2)
        fit2$y <- c(f2$y, y2.sample)
        fit2$tau2hat <- drop(crossprod(t(crossprod(fit2$y - fit2$mu.hat, fit2$Ki)), (fit2$y - fit2$mu.hat)) / length(fit2$y))
        
        fit.tmp$fits[[2]] <- fit2
        
        return(mean(predict(fit.tmp, Xref)$sig2[[3]])) # to minimize the deduced variance. To maximize, -mean
      }
      attr(ALC.out, "rng") <- NULL
      attr(ALC.out, "doRNG_version") <- NULL
    } else {
      ALC.out <- rep(0, mc.sample)
      for (i in 1:mc.sample) {
        fit.tmp <- fit
        fit1$y <- c(f1$y, y1.sample[i])
        fit1$tau2hat <- drop(crossprod(t(crossprod(fit1$y - fit1$mu.hat, fit1$Ki)), (fit1$y - fit1$mu.hat)) / length(fit1$y))
        
        fit.tmp$fits[[1]] <- fit1
        
        ### Choose level 2 ###
        if (kernel == "sqex") {
          pred2 <- pred.GP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        } else if (kernel == "matern1.5") {
          pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        } else if (kernel == "matern2.5") {
          pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        }
        
        ### update Ki2
        newx2 <- cbind(Xcand, y1.sample[i])
        
        if (kernel == "sqex") {
          cov.newx2 <- covar.sep(X1 = newx2, d = f2$theta, g = g)
          cov.Xnewx2 <- covar.sep(X1 = f2$X, X2 = newx2, d = f2$theta, g = 0)
        } else if (kernel == "matern1.5") {
          cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 1.5)
          cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 1.5)
        } else if (kernel == "matern2.5") {
          cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 2.5)
          cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 2.5)
        }
        v.next2 <- drop(cov.newx2 - crossprod(t(crossprod(cov.Xnewx2, f2$Ki)), cov.Xnewx2))
        g.next2 <- -1 / drop(v.next2) * crossprod(f2$Ki, cov.Xnewx2)
        
        fit2$Ki <- rbind(
          cbind(f2$Ki + crossprod(t(g.next2), t(g.next2) * v.next2), g.next2),
          cbind(t(g.next2), 1 / drop(v.next2))
        )
        
        fit2$X <- rbind(f2$X, newx2)
        fit2$y <- c(f2$y, y2.sample)
        fit2$tau2hat <- drop(crossprod(t(crossprod(fit2$y - fit2$mu.hat, fit2$Ki)), (fit2$y - fit2$mu.hat)) / length(fit2$y))
        
        fit.tmp$fits[[2]] <- fit2
        
        ALC.out[i] <- mean(predict(fit.tmp, Xref)$sig2[[3]]) # to minimize the deduced variance. To maximize, -mean
      }
    }
  }else{
    fit.tmp <- fit
    fit1$y <- c(f1$y, y1.sample)
    fit1$tau2hat <- drop(crossprod(t(crossprod(fit1$y - fit1$mu.hat, fit1$Ki)), (fit1$y - fit1$mu.hat)) / length(fit1$y))
    
    fit.tmp$fits[[1]] <- fit1
    
    ### Choose level 2 ###
    if (kernel == "sqex") {
      pred2 <- pred.GP(fit2, cbind(Xcand, y1.sample))
    } else if (kernel == "matern1.5") {
      pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample))
    } else if (kernel == "matern2.5") {
      pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample))
    }
    y2.sample <- pred2$mu
    
    ### update Ki2
    newx2 <- cbind(Xcand, y1.sample)
    
    if (kernel == "sqex") {
      cov.newx2 <- covar.sep(X1 = newx2, d = f2$theta, g = g)
      cov.Xnewx2 <- covar.sep(X1 = f2$X, X2 = newx2, d = f2$theta, g = 0)
    } else if (kernel == "matern1.5") {
      cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 1.5)
      cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 1.5)
    } else if (kernel == "matern2.5") {
      cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 2.5)
      cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 2.5)
    }
    v.next2 <- drop(cov.newx2 - crossprod(t(crossprod(cov.Xnewx2, f2$Ki)), cov.Xnewx2))
    g.next2 <- -1 / drop(v.next2) * crossprod(f2$Ki, cov.Xnewx2)
    
    fit2$Ki <- rbind(
      cbind(f2$Ki + crossprod(t(g.next2), t(g.next2) * v.next2), g.next2),
      cbind(t(g.next2), 1 / drop(v.next2))
    )
    
    fit2$X <- rbind(f2$X, newx2)
    fit2$y <- c(f2$y, y2.sample)
    fit2$tau2hat <- drop(crossprod(t(crossprod(fit2$y - fit2$mu.hat, fit2$Ki)), (fit2$y - fit2$mu.hat)) / length(fit2$y))
    
    fit.tmp$fits[[2]] <- fit2
    
    ALC.out <- mean(predict(fit.tmp, Xref)$sig2[[3]])
  }
  return(mean(ALC.out))
}


#' object to optimize the next point by ALC criterion updating at level 3 with three levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param Xref vector or matrix of reference data.
#' @param fit an object of class RNAmf.
#' @param mc.sample a number of mc samples generated for this approach. Default is 100.
#' @param parallel logical indicating whether to run parallel or not. Default is FALSE.
#' @param ncore the number of core for parallel. Default is 1.
#' @return A mean of the deduced variance at Xref.
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm
#' @noRd
#'

obj.ALC_level_3_3 <- function(Xcand, Xref, fit, MC=FALSE, mc.sample, parallel = FALSE, ncore = 1) {
  fit1 <- f1 <- fit$fits[[1]]
  fit2 <- f2 <- fit$fits[[2]]
  fit3 <- f3 <- fit$fits[[3]]
  constant <- fit$constant
  kernel <- fit$kernel
  g <- fit1$g
  
  Xcand <- matrix(Xcand, nrow = 1)
  if(MC){
    if (kernel == "sqex") {
      y1.sample <- rnorm(mc.sample, mean = pred.GP(f1, Xcand)$mu, sd = sqrt(pred.GP(f1, Xcand)$sig2))
    } else if (kernel == "matern1.5") {
      y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
    } else if (kernel == "matern2.5") {
      y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
    }
  }else{
    if (kernel == "sqex") {
      y1.sample <- pred.GP(f1, Xcand)$mu
    } else if (kernel == "matern1.5") {
      y1.sample <- pred.matGP(f1, Xcand)$mu
    } else if (kernel == "matern2.5") {
      y1.sample <- pred.matGP(f1, Xcand)$mu
    }
  }
  Xcand1 <- Xcand
  
  ### Choose level 1 ###
  ### update Ki1
  if (kernel == "sqex") {
    cov.newx1 <- covar.sep(X1 = Xcand1, d = f1$theta, g = g)
    cov.Xnewx1 <- covar.sep(X1 = f1$X, X2 = Xcand1, d = f1$theta, g = 0)
  } else if (kernel == "matern1.5") {
    cov.newx1 <- cor.sep(X = Xcand1, theta = f1$theta, nu = 1.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand1, theta = f1$theta, nu = 1.5)
  } else if (kernel == "matern2.5") {
    cov.newx1 <- cor.sep(X = Xcand1, theta = f1$theta, nu = 2.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand1, theta = f1$theta, nu = 2.5)
  }
  v.next1 <- drop(cov.newx1 - crossprod(t(crossprod(cov.Xnewx1, f1$Ki)), cov.Xnewx1))
  g.next1 <- -1 / drop(v.next1) * crossprod(f1$Ki, cov.Xnewx1)
  
  
  fit1$Ki <- rbind(
    cbind(f1$Ki + crossprod(t(g.next1), t(g.next1) * v.next1), g.next1),
    cbind(t(g.next1), 1 / drop(v.next1))
  )
  
  fit1$X <- rbind(f1$X, Xcand1)
  
  if(MC){
    if (parallel) {
      ALC.out <- foreach(i = 1:mc.sample, .combine = c) %dorng% {
        fit.tmp <- fit
        fit1$y <- c(f1$y, y1.sample[i])
        fit1$tau2hat <- drop(crossprod(t(crossprod(fit1$y - fit1$mu.hat, fit1$Ki)), (fit1$y - fit1$mu.hat)) / length(fit1$y))
        
        fit.tmp$fits[[1]] <- fit1
        
        ### Choose level 2 ###
        if (kernel == "sqex") {
          pred2 <- pred.GP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        } else if (kernel == "matern1.5") {
          pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        } else if (kernel == "matern2.5") {
          pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        }
        
        ### update Ki2
        newx2 <- cbind(Xcand, y1.sample[i])
        
        if (kernel == "sqex") {
          cov.newx2 <- covar.sep(X1 = newx2, d = f2$theta, g = g)
          cov.Xnewx2 <- covar.sep(X1 = f2$X, X2 = newx2, d = f2$theta, g = 0)
        } else if (kernel == "matern1.5") {
          cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 1.5)
          cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 1.5)
        } else if (kernel == "matern2.5") {
          cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 2.5)
          cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 2.5)
        }
        v.next2 <- drop(cov.newx2 - crossprod(t(crossprod(cov.Xnewx2, f2$Ki)), cov.Xnewx2))
        g.next2 <- -1 / drop(v.next2) * crossprod(f2$Ki, cov.Xnewx2)
        
        fit2$Ki <- rbind(
          cbind(f2$Ki + crossprod(t(g.next2), t(g.next2) * v.next2), g.next2),
          cbind(t(g.next2), 1 / drop(v.next2))
        )
        
        fit2$X <- rbind(f2$X, newx2)
        fit2$y <- c(f2$y, y2.sample)
        fit2$tau2hat <- drop(crossprod(t(crossprod(fit2$y - fit2$mu.hat, fit2$Ki)), (fit2$y - fit2$mu.hat)) / length(fit2$y))
        
        fit.tmp$fits[[2]] <- fit2
        
        ### Choose level 3 ###
        if (kernel == "sqex") {
          pred3 <- pred.GP(fit2, cbind(Xcand, y2.sample))
          y3.sample <- rnorm(1, pred3$mu, sqrt(pred3$sig2))
        } else if (kernel == "matern1.5") {
          pred3 <- pred.matGP(fit2, cbind(Xcand, y2.sample))
          y3.sample <- rnorm(1, pred3$mu, sqrt(pred3$sig2))
        } else if (kernel == "matern2.5") {
          pred3 <- pred.matGP(fit2, cbind(Xcand, y2.sample))
          y3.sample <- rnorm(1, pred3$mu, sqrt(pred3$sig2))
        }
        
        ### update Ki2
        newx3 <- cbind(Xcand, y2.sample)
        
        if (kernel == "sqex") {
          cov.newx3 <- covar.sep(X1 = newx3, d = f3$theta, g = g)
          cov.Xnewx3 <- covar.sep(X1 = f3$X, X2 = newx3, d = f3$theta, g = 0)
        } else if (kernel == "matern1.5") {
          cov.newx3 <- cor.sep(X = newx3, theta = f3$theta, nu = 1.5)
          cov.Xnewx3 <- cor.sep(X = f3$X, x = newx3, theta = f3$theta, nu = 1.5)
        } else if (kernel == "matern2.5") {
          cov.newx3 <- cor.sep(X = newx3, theta = f3$theta, nu = 2.5)
          cov.Xnewx3 <- cor.sep(X = f3$X, x = newx3, theta = f3$theta, nu = 2.5)
        }
        v.next3 <- drop(cov.newx3 - crossprod(t(crossprod(cov.Xnewx3, f3$Ki)), cov.Xnewx3))
        g.next3 <- -1 / drop(v.next3) * crossprod(f3$Ki, cov.Xnewx3)
        
        fit3$Ki <- rbind(
          cbind(f3$Ki + crossprod(t(g.next3), t(g.next3) * v.next3), g.next3),
          cbind(t(g.next3), 1 / drop(v.next3))
        )
        
        fit3$X <- rbind(f3$X, newx3)
        fit3$y <- c(f3$y, y3.sample)
        fit3$tau2hat <- drop(crossprod(t(crossprod(fit3$y - fit3$mu.hat, fit3$Ki)), (fit3$y - fit3$mu.hat)) / length(fit3$y))
        
        fit.tmp$fits[[3]] <- fit3
        
        return(mean(predict(fit.tmp, Xref)$sig2[[3]])) # to minimize the deduced variance. To maximize, -mean
      }
      attr(ALC.out, "rng") <- NULL
      attr(ALC.out, "doRNG_version") <- NULL
    } else {
      ALC.out <- rep(0, mc.sample)
      for (i in 1:mc.sample) {
        fit.tmp <- fit
        fit1$y <- c(f1$y, y1.sample[i])
        fit1$tau2hat <- drop(crossprod(t(crossprod(fit1$y - fit1$mu.hat, fit1$Ki)), (fit1$y - fit1$mu.hat)) / length(fit1$y))
        
        fit.tmp$fits[[1]] <- fit1
        
        ### Choose level 2 ###
        if (kernel == "sqex") {
          pred2 <- pred.GP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        } else if (kernel == "matern1.5") {
          pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        } else if (kernel == "matern2.5") {
          pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
          y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
        }
        
        ### update Ki2
        newx2 <- cbind(Xcand, y1.sample[i])
        
        if (kernel == "sqex") {
          cov.newx2 <- covar.sep(X1 = newx2, d = f2$theta, g = g)
          cov.Xnewx2 <- covar.sep(X1 = f2$X, X2 = newx2, d = f2$theta, g = 0)
        } else if (kernel == "matern1.5") {
          cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 1.5)
          cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 1.5)
        } else if (kernel == "matern2.5") {
          cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 2.5)
          cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 2.5)
        }
        v.next2 <- drop(cov.newx2 - crossprod(t(crossprod(cov.Xnewx2, f2$Ki)), cov.Xnewx2))
        g.next2 <- -1 / drop(v.next2) * crossprod(f2$Ki, cov.Xnewx2)
        
        fit2$Ki <- rbind(
          cbind(f2$Ki + crossprod(t(g.next2), t(g.next2) * v.next2), g.next2),
          cbind(t(g.next2), 1 / drop(v.next2))
        )
        
        fit2$X <- rbind(f2$X, newx2)
        fit2$y <- c(f2$y, y2.sample)
        fit2$tau2hat <- drop(crossprod(t(crossprod(fit2$y - fit2$mu.hat, fit2$Ki)), (fit2$y - fit2$mu.hat)) / length(fit2$y))
        
        fit.tmp$fits[[2]] <- fit2
        
        ### Choose level 3 ###
        if (kernel == "sqex") {
          pred3 <- pred.GP(fit2, cbind(Xcand, y2.sample))
          y3.sample <- rnorm(1, pred3$mu, sqrt(pred3$sig2))
        } else if (kernel == "matern1.5") {
          pred3 <- pred.matGP(fit2, cbind(Xcand, y2.sample))
          y3.sample <- rnorm(1, pred3$mu, sqrt(pred3$sig2))
        } else if (kernel == "matern2.5") {
          pred3 <- pred.matGP(fit2, cbind(Xcand, y2.sample))
          y3.sample <- rnorm(1, pred3$mu, sqrt(pred3$sig2))
        }
        
        ### update Ki2
        newx3 <- cbind(Xcand, y2.sample)
        
        if (kernel == "sqex") {
          cov.newx3 <- covar.sep(X1 = newx3, d = f3$theta, g = g)
          cov.Xnewx3 <- covar.sep(X1 = f3$X, X2 = newx3, d = f3$theta, g = 0)
        } else if (kernel == "matern1.5") {
          cov.newx3 <- cor.sep(X = newx3, theta = f3$theta, nu = 1.5)
          cov.Xnewx3 <- cor.sep(X = f3$X, x = newx3, theta = f3$theta, nu = 1.5)
        } else if (kernel == "matern2.5") {
          cov.newx3 <- cor.sep(X = newx3, theta = f3$theta, nu = 2.5)
          cov.Xnewx3 <- cor.sep(X = f3$X, x = newx3, theta = f3$theta, nu = 2.5)
        }
        v.next3 <- drop(cov.newx3 - crossprod(t(crossprod(cov.Xnewx3, f3$Ki)), cov.Xnewx3))
        g.next3 <- -1 / drop(v.next3) * crossprod(f3$Ki, cov.Xnewx3)
        
        fit3$Ki <- rbind(
          cbind(f3$Ki + crossprod(t(g.next3), t(g.next3) * v.next3), g.next3),
          cbind(t(g.next3), 1 / drop(v.next3))
        )
        
        fit3$X <- rbind(f3$X, newx3)
        fit3$y <- c(f3$y, y3.sample)
        fit3$tau2hat <- drop(crossprod(t(crossprod(fit3$y - fit3$mu.hat, fit3$Ki)), (fit3$y - fit3$mu.hat)) / length(fit3$y))
        
        fit.tmp$fits[[3]] <- fit3
        
        ALC.out[i] <- mean(predict(fit.tmp, Xref)$sig2[[3]]) # to minimize the deduced variance. To maximize, -mean
      }
    }
  }else{
    fit.tmp <- fit
    fit1$y <- c(f1$y, y1.sample)
    fit1$tau2hat <- drop(crossprod(t(crossprod(fit1$y - fit1$mu.hat, fit1$Ki)), (fit1$y - fit1$mu.hat)) / length(fit1$y))
    
    fit.tmp$fits[[1]] <- fit1
    
    ### Choose level 2 ###
    if (kernel == "sqex") {
      pred2 <- pred.GP(fit2, cbind(Xcand, y1.sample))
    } else if (kernel == "matern1.5") {
      pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample))
    } else if (kernel == "matern2.5") {
      pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample))
    }
    y2.sample <- pred2$mu
    
    ### update Ki2
    newx2 <- cbind(Xcand, y1.sample)
    
    if (kernel == "sqex") {
      cov.newx2 <- covar.sep(X1 = newx2, d = f2$theta, g = g)
      cov.Xnewx2 <- covar.sep(X1 = f2$X, X2 = newx2, d = f2$theta, g = 0)
    } else if (kernel == "matern1.5") {
      cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 1.5)
      cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 1.5)
    } else if (kernel == "matern2.5") {
      cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 2.5)
      cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 2.5)
    }
    v.next2 <- drop(cov.newx2 - crossprod(t(crossprod(cov.Xnewx2, f2$Ki)), cov.Xnewx2))
    g.next2 <- -1 / drop(v.next2) * crossprod(f2$Ki, cov.Xnewx2)
    
    fit2$Ki <- rbind(
      cbind(f2$Ki + crossprod(t(g.next2), t(g.next2) * v.next2), g.next2),
      cbind(t(g.next2), 1 / drop(v.next2))
    )
    
    fit2$X <- rbind(f2$X, newx2)
    fit2$y <- c(f2$y, y2.sample)
    fit2$tau2hat <- drop(crossprod(t(crossprod(fit2$y - fit2$mu.hat, fit2$Ki)), (fit2$y - fit2$mu.hat)) / length(fit2$y))
    
    fit.tmp$fits[[2]] <- fit2
    
    ### Choose level 3 ###
    if (kernel == "sqex") {
      pred3 <- pred.GP(fit2, cbind(Xcand, y2.sample))
    } else if (kernel == "matern1.5") {
      pred3 <- pred.matGP(fit2, cbind(Xcand, y2.sample))
    } else if (kernel == "matern2.5") {
      pred3 <- pred.matGP(fit2, cbind(Xcand, y2.sample))
    }
    y3.sample <- pred3$mu
    
    ### update Ki2
    newx3 <- cbind(Xcand, y2.sample)
    
    if (kernel == "sqex") {
      cov.newx3 <- covar.sep(X1 = newx3, d = f3$theta, g = g)
      cov.Xnewx3 <- covar.sep(X1 = f3$X, X2 = newx3, d = f3$theta, g = 0)
    } else if (kernel == "matern1.5") {
      cov.newx3 <- cor.sep(X = newx3, theta = f3$theta, nu = 1.5)
      cov.Xnewx3 <- cor.sep(X = f3$X, x = newx3, theta = f3$theta, nu = 1.5)
    } else if (kernel == "matern2.5") {
      cov.newx3 <- cor.sep(X = newx3, theta = f3$theta, nu = 2.5)
      cov.Xnewx3 <- cor.sep(X = f3$X, x = newx3, theta = f3$theta, nu = 2.5)
    }
    v.next3 <- drop(cov.newx3 - crossprod(t(crossprod(cov.Xnewx3, f3$Ki)), cov.Xnewx3))
    g.next3 <- -1 / drop(v.next3) * crossprod(f3$Ki, cov.Xnewx3)
    
    fit3$Ki <- rbind(
      cbind(f3$Ki + crossprod(t(g.next3), t(g.next3) * v.next3), g.next3),
      cbind(t(g.next3), 1 / drop(v.next3))
    )
    
    fit3$X <- rbind(f3$X, newx3)
    fit3$y <- c(f3$y, y3.sample)
    fit3$tau2hat <- drop(crossprod(t(crossprod(fit3$y - fit3$mu.hat, fit3$Ki)), (fit3$y - fit3$mu.hat)) / length(fit3$y))
    
    fit.tmp$fits[[3]] <- fit3
    
    ALC.out <- mean(predict(fit.tmp, Xref)$sig2[[3]])
  }
  return(mean(ALC.out))
}


#' @title Active Learning for Recursive Non-Additive Emulator
#'
#' @description The function acquires the new point and fidelity level by maximizing one of the
#' four active learning criteria: ALM, ALC, ALD, or ALMC.
#'
#' \itemize{
#'   \item \strong{ALM} (Active Learning MacKay): It calculates the ALM criterion \eqn{\frac{\sigma^{*2}_l(\bm{x})}{\sum^l_{j=1}C_j}},
#'   where \eqn{\sigma^{*2}_l(\bm{x})} is the posterior predictive variance
#'   at each fidelity level \eqn{l} and \eqn{C_j} is the simulation cost at level \eqn{j}.
#'
#'   \item \strong{ALD} (Active Learning Decomposition): It calculates the ALD criterion \eqn{\frac{V_l(\bm{x})}{\sum^l_{j=1}C_j}},
#'   where \eqn{V_l(\bm{x})} is the variance contribution of GP emulator
#'   at each fidelity level \eqn{l} and \eqn{C_j} is the simulation cost at level \eqn{j}.
#'
#'   \item \strong{ALC} (Active Learning Cohn): It calculates the ALC criterion
#'   \eqn{\frac{\Delta \sigma_L^{2}(l,\bm{x})}{\sum^l_{j=1}C_j} =
#'   \frac{\int_{\Omega} \sigma_L^{*2}(\bm{\xi})-\tilde{\sigma}_L^{*2}(\bm{\xi};l,\bm{x}){\rm{d}}\bm{\xi}}{\sum^l_{j=1}C_j}},
#'   where \eqn{f_L} is the highest-fidelity simulation code,
#'   \eqn{\sigma_L^{*2}(\bm{\xi})} is the posterior variance of \eqn{f_L(\bm{\xi})},
#'   \eqn{C_j} is the simulation cost at fidelity level \eqn{j},
#'   and \eqn{\tilde{\sigma}_L^{*2}(\bm{\xi};l,\bm{x})} is the posterior variance
#'   based on the augmented design combining the current design and a new input location \eqn{\bm{x}}
#'   at each fidelity level lower than or equal to \eqn{l}.
#'   The integration is approximated by MC integration using uniform reference samples.
#'
#'   \item \strong{ALMC} (Active Learning MacKay-Cohn): A hybrid approach.
#'   It finds the optimal input location \eqn{\bm{x}^*}
#'   by maximizing \eqn{\sigma^{*2}_L(\bm{x})},
#'   the posterior predictive variance at the highest-fidelity level \eqn{L}.
#'   After selecting \eqn{\bm{x}^*},
#'   it finds the optimal fidelity level by maximizing ALC criterion at \eqn{\bm{x}^*},
#'   \eqn{\text{argmax}_{l\in\{1,\ldots,L\}} \frac{\Delta \sigma_L^{2}(l,\bm{x}^*)}{\sum^l_{j=1}C_j}},
#'   where \eqn{C_j} is the simulation cost at level \eqn{j}.
#' }
#'
#' A new point is acquired on \code{Xcand}. If \code{Xcand=NULL} and \code{Xref=NULL}, a new point is acquired on unit hypercube \eqn{[0,1]^d}.
#'
#' For details, see Heo and Sung (2025, <\doi{https://doi.org/10.1080/00401706.2024.2376173}>).
#'
#' @details For \code{"ALC"}, or \code{"ALMC"}, \code{Xref} plays a role of \eqn{\bm{\xi}} to approximate the integration.
#' To impute the posterior variance based on the augmented design \eqn{\tilde{\sigma}_L^{*2}(\bm{\xi};l,\bm{x})},
#' MC approximation is used.
#' Due to the nested assumption, imputing \eqn{y^{[s]}_{n_s+1}} for each \eqn{1\leq s\leq l} by drawing samples
#' from the posterior distribution of \eqn{f_s(\bm{x}^{[s]}_{n_s+1})}
#' based on the current design allows to compute \eqn{\tilde{\sigma}_L^{*2}(\bm{\xi};l,\bm{x})}.
#' Inverse of covariance matrix is computed by the Sherman-Morrison formula.
#'
#' To search for the next acquisition \eqn{\bm{x}^*} by maximizing AL criterion,
#' the gradient-based optimization can be used by \code{optim=TRUE}.
#' Firstly, \eqn{\tilde{\sigma}_L^{*2}(\bm{\xi};l,\bm{x})} is computed at
#' \eqn{5 \times d} number of points.
#' After that, the point minimizing \eqn{\tilde{\sigma}_L^{*2}(\bm{\xi};l,\bm{x})}
#' serves as a starting point of optimization by \code{L-BFGS-B} method.
#' Otherwise, when \code{optim=FALSE}, AL criterion is optimized only on \code{Xcand}.
#'
#' The point is selected by maximizing the ALC criterion:
#' \eqn{\text{argmax}_{l\in\{1,\ldots,L\}; \bm{x} \in \Omega}
#' \frac{\Delta \sigma_L^{2}(l,\bm{x})}{\sum^l_{j=1}C_j}}.
#'
#' @param criterion character string specifying the active learning criterion to use. Must be one of \code{"ALM"}, \code{"ALD"}, \code{"ALC"}, or \code{"ALMC"}. Default is \code{"ALM"}.
#' @param fit object of class \code{RNAmf}.
#' @param Xref vector or matrix of reference locations to approximate the integral of ALC. If \code{Xref=NULL}, \eqn{100 \times d} points from 0 to 1 are generated by Latin hypercube design. Only used when \code{criterion="ALC"} or \code{"ALMC"}. Default is \code{NULL}.
#' @param Xcand vector or matrix of the candidate set for grid-based search. If \code{use_optim=FALSE}, the criterion is evaluated and optimized only on this set. If \code{Xcand=NULL}, \eqn{100 \times d} points from 0 to 1 generated by Latin hypercube design (or \code{Xref} for ALC and ALMC) are used. Default is \code{NULL}.
#' @param MC logical indicating whether to use Monte Carlo approximation to impute the posterior variance (for ALC/ALMC). If \code{FALSE}, posterior means are used. Default is \code{FALSE}.
#' @param mc.sample a number of MC samples generated for the imputation through MC approximation. Default is \code{100}.
#' @param cost vector of the costs for each level of fidelity. If \code{cost=NULL}, total costs at all levels would be 1. \code{cost} is encouraged to have an ascending order of positive values. Default is \code{NULL}.
#' @param use_optim logical indicating whether to optimize the criterion using \code{optim}'s gradient-based \code{L-BFGS-B} method. If \code{TRUE}, \eqn{5 \times d} starting points are generated by Latin hypercube design for optimization. If \code{FALSE}, the point is selected from \code{Xcand}. Default is \code{TRUE}.
#' @param parallel logical indicating whether to use parallel computation. Default is \code{FALSE}.
#' @param ncore integer specifying the number of cores for parallel computation. Used only if \code{parallel=TRUE}. Default is 1.
#' @param trace logical indicating whether to print the computational progress and time. Default is \code{TRUE}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{AL}: The values of the selected criterion at the candidate points (if \code{use_optim=FALSE}) or the optimized value (if \code{use_optim=TRUE}). For \code{ALMC}, it returns the ALC scores for each level at the chosen point.
#'   \item \code{cost}: A copy of the \code{cost} argument.
#'   \item \code{Xcand}: A copy of the \code{Xcand} argument used.
#'   \item \code{chosen}: A list containing the chosen fidelity \code{level} and the new input location \code{Xnext}.
#'   \item \code{time}: The computation time in seconds.
#' }
#'
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm optim
#' @importFrom lhs randomLHS
#' @importFrom foreach foreach
#' @importFrom doRNG %dorng%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#'
#' @usage AL_RNAmf(criterion = c("ALM", "ALC", "ALD", "ALMC"), fit,
#' Xref = NULL, Xcand = NULL, MC = FALSE, mc.sample = 100, cost = NULL,
#' use_optim = TRUE, parallel = FALSE, ncore = 1, trace = TRUE)
#' @export
#' @examples
#' \donttest{
#' ### simulation costs ###
#' cost <- c(1, 3)
#'
#' ### 1-d Perdikaris function in Perdikaris, et al. (2017) ###
#' # low-fidelity function
#' f1 <- function(x) {
#'   sin(8 * pi * x)
#' }
#'
#' # high-fidelity function
#' f2 <- function(x) {
#'   (x - sqrt(2)) * (sin(8 * pi * x))^2
#' }
#'
#' ### training data ###
#' n1 <- 13
#' n2 <- 8
#'
#' ### fix seed to reproduce the result ###
#' set.seed(1)
#'
#' ### generate initial nested design ###
#' X <- NestedX(c(n1, n2), 1)
#' X1 <- X[[1]]
#' X2 <- X[[2]]
#'
#' y1 <- f1(X1)
#' y2 <- f2(X2)
#'
#' ### n=100 uniform test data ###
#' x <- seq(0, 1, length.out = 100)
#'
#' ### fit an RNAmf ###
#' fit.RNAmf <- RNAmf(list(X1, X2), list(y1, y2), kernel = "sqex", constant=TRUE)
#'
#' ### 1. ALM Criterion ###
#' alm.RNAmf <- AL_RNAmf(criterion="ALM",
#'                       Xcand = x, fit=fit.RNAmf, cost = cost,
#'                       use_optim = FALSE, parallel = TRUE, ncore = 2)
#' print(alm.RNAmf$chosen)
#'
#' ### visualize ALM ###
#' oldpar <- par(mfrow = c(1, 2))
#' plot(x, alm.RNAmf$AL$ALM1,
#'      type = "l", lty = 2,
#'      xlab = "x", ylab = "ALM criterion at the low-fidelity level",
#'      ylim = c(min(c(alm.RNAmf$AL$ALM1, alm.RNAmf$AL$ALM2)),
#'               max(c(alm.RNAmf$AL$ALM1, alm.RNAmf$AL$ALM2))))
#' points(alm.RNAmf$chosen$Xnext,
#'        alm.RNAmf$AL$ALM1[which(x == drop(alm.RNAmf$chosen$Xnext))],
#'        pch = 16, cex = 1, col = "red")
#' plot(x, alm.RNAmf$AL$ALM2,
#'      type = "l", lty = 2,
#'      xlab = "x", ylab = "ALM criterion at the high-fidelity level",
#'      ylim = c(min(c(alm.RNAmf$AL$ALM1, alm.RNAmf$AL$ALM2)),
#'               max(c(alm.RNAmf$AL$ALM1, alm.RNAmf$AL$ALM2))))
#' par(oldpar)
#'
#'
#' ### 2. ALD Criterion ###
#' ald.RNAmf <- AL_RNAmf(criterion="ALD",
#'                       Xcand = x, fit=fit.RNAmf, cost = cost,
#'                       use_optim = FALSE, parallel = TRUE, ncore = 2)
#' print(ald.RNAmf$chosen)
#'
#' ### visualize ALD ###
#' oldpar <- par(mfrow = c(1, 2))
#' plot(x, ald.RNAmf$AL$ALD1,
#'      type = "l", lty = 2,
#'      xlab = "x", ylab = "ALD criterion at the low-fidelity level",
#'      ylim = c(min(c(ald.RNAmf$AL$ALD1, ald.RNAmf$AL$ALD2)),
#'               max(c(ald.RNAmf$AL$ALD1, ald.RNAmf$AL$ALD2))))
#' points(ald.RNAmf$chosen$Xnext,
#'        ald.RNAmf$AL$ALD1[which(x == drop(ald.RNAmf$chosen$Xnext))],
#'        pch = 16, cex = 1, col = "red")
#' plot(x, ald.RNAmf$AL$ALD2,
#'      type = "l", lty = 2,
#'      xlab = "x", ylab = "ALD criterion at the high-fidelity level",
#'      ylim = c(min(c(ald.RNAmf$AL$ALD1, ald.RNAmf$AL$ALD2)),
#'               max(c(ald.RNAmf$AL$ALD1, ald.RNAmf$AL$ALD2))))
#' par(oldpar)
#'
#'
#' ### 3. ALC Criterion ###
#' alc.RNAmf <- AL_RNAmf(criterion="ALC",
#'                       Xref = x, Xcand = x, fit=fit.RNAmf, cost = cost,
#'                       use_optim = FALSE, parallel = TRUE, ncore = 2)
#' print(alc.RNAmf$chosen)
#'
#' ### visualize ALC ###
#' oldpar <- par(mfrow = c(1, 2))
#' plot(x, alc.RNAmf$AL$ALC1,
#'      type = "l", lty = 2,
#'      xlab = "x", ylab = "ALC criterion augmented at the low-fidelity level",
#'      ylim = c(min(c(alc.RNAmf$AL$ALC1, alc.RNAmf$AL$ALC2)),
#'               max(c(alc.RNAmf$AL$ALC1, alc.RNAmf$AL$ALC2))))
#' points(alc.RNAmf$chosen$Xnext,
#'        alc.RNAmf$AL$ALC1[which(x == drop(alc.RNAmf$chosen$Xnext))],
#'        pch = 16, cex = 1, col = "red")
#' plot(x, alc.RNAmf$AL$ALC2,
#'      type = "l", lty = 2,
#'      xlab = "x", ylab = "ALC criterion augmented at the high-fidelity level",
#'      ylim = c(min(c(alc.RNAmf$AL$ALC1, alc.RNAmf$AL$ALC2)),
#'               max(c(alc.RNAmf$AL$ALC1, alc.RNAmf$AL$ALC2))))
#' par(oldpar)
#'
#'
#' ### 4. ALMC Criterion ###
#' almc.RNAmf <- AL_RNAmf(criterion="ALMC",
#'                        Xref = x, Xcand = x, fit=fit.RNAmf, cost = cost,
#'                        use_optim = FALSE, parallel = TRUE, ncore = 2)
#' print(almc.RNAmf$chosen)
#' }

AL_RNAmf <- function(criterion = c("ALM", "ALC", "ALD", "ALMC"), fit, Xref = NULL, Xcand = NULL, MC = FALSE, mc.sample = 100, cost = NULL, use_optim = TRUE, parallel = FALSE, ncore = 1, trace = TRUE) {
  t0 <- proc.time()
  ### check the object ###
  if (!inherits(fit, "RNAmf")) stop("The object is not of class \"RNAmf\" \n")
  criterion <- match.arg(criterion)
  L <- fit$level
  d <- ncol(fit$fits[[1]]$X)
  
  ### check the cost ###
  if (is.null(cost)) {
    cost <- c(1, rep(0, L - 1))
  } else {
    if (!is.numeric(cost)) stop("'cost' must be a numeric vector.")
    if (length(cost) != L) stop(paste0("The length of 'cost' must equal the number of fidelity levels (", L, ")."))
    if (any(cost < 0)) stop("'cost' must contain non-negative values.")
    for(l in 1:(L-1)) {
      if (cost[l] >= cost[l+1]) warning(paste("Cost for level", l, "is >= level", l+1, ". Usually, higher fidelity is more expensive."))
    }
  }
  
  ### Xref check ###
  if (!is.null(Xref)) {
    if (is.vector(Xref)) Xref <- matrix(Xref, ncol = 1)
    if (!is.matrix(Xref)) stop("'Xref' must be a matrix or vector.")
    if (ncol(Xref) != d) stop(paste0("'Xref' must have ", d, " columns (matching the model input dimension)."))
  }
  ### Xcand check ###
  if (!is.null(Xcand)) {
    if (is.vector(Xcand)) Xcand <- matrix(Xcand, ncol = 1)
    if (!is.matrix(Xcand)) stop("'Xcand' must be a matrix or vector.")
    if (ncol(Xcand) != d) stop(paste0("'Xcand' must have ", d, " columns (matching the model input dimension)."))
  }
  
  if (!is.logical(MC)) stop("'MC' must be TRUE or FALSE.")
  if (!is.numeric(mc.sample) || length(mc.sample) != 1 || mc.sample < 1) stop("'mc.sample' must be a positive integer.")
  if (!is.logical(use_optim)) stop("'use_optim' must be TRUE or FALSE.")
  if (!is.logical(parallel)) stop("'parallel' must be TRUE or FALSE.")
  if (!is.numeric(ncore) || length(ncore) != 1 || ncore < 1) stop("'ncore' must be a positive integer.")
  if (!is.logical(trace)) stop("'trace' must be TRUE or FALSE.")
  
  if ((criterion %in% c("ALC", "ALMC")) && is.null(Xref)) {
    Xref <- randomLHS(100 * d, d)
  }
  
  ### parallel setup ###
  if (parallel) {
    cl <- makeCluster(ncore)
    registerDoParallel(cl)
    clusterEvalQ(cl, library(RNAmf))
    on.exit(stopCluster(cl), add = TRUE)
  }
  
  ### Xcand setup ###
  if (use_optim) {
    Xcand <- randomLHS(5 * d, d)
  } else {
    if (is.null(Xcand)) {
      Xcand <- randomLHS(100 * d, d)
    } else if (is.null(dim(Xcand))) {
      Xcand <- matrix(Xcand, ncol = 1)
    }
  }
  
  ### Helper: Get Function Name ###
  get_obj_func_name <- function(crit, l, total_L) {
    if (crit == "ALM") return(paste0("obj.ALM_level_", l))
    if (crit == "ALC") {
      if (total_L == 3) return(paste0("obj.ALC_level_3_", l))
      return(paste0("obj.ALC_level_", l))
    }
    if (crit == "ALD") {
      suffix <- if (total_L == 2) "2level" else "3level"
      return(paste0("obj.ALD_V", l, "_", suffix))
    }
    return(NULL)
  }
  
  ### Helper: Current Variance (for ALC/ALMC) ###
  Icurrent <- 0
  if (criterion %in% c("ALC", "ALMC")) {
    pred_sigs <- predict(fit, Xref)$sig2
    Icurrent <- mean(pred_sigs[[L]])
  }
  
  if (criterion == "ALMC") { ### ALMC ###
    ### 1. Find Xnext using ALM at the highest level (L) ###
    obj_fun <- get(get_obj_func_name("ALM", L, L))
    
    # A. Grid Search on Xcand
    if (parallel) {
      vals <- foreach(i = 1:nrow(Xcand), .combine = c) %dorng% {
        -obj_fun(matrix(Xcand[i, ], nrow = 1), fit = fit)
      }
      attr(vals, "rng") <- NULL
      attr(vals, "doRNG_version") <- NULL
    } else {
      vals <- numeric(nrow(Xcand))
      for(i in 1:nrow(Xcand)) vals[i] <- -obj_fun(matrix(Xcand[i, ], nrow = 1), fit = fit)
    }
    
    # B. Optimization on Xcand
    best_idx <- which.max(vals)
    X.start <- matrix(Xcand[best_idx, ], nrow = 1)
    if (use_optim) {
      o <- optim(X.start, obj_fun, method = "L-BFGS-B", lower = 0, upper = 1, fit = fit)
      Xnext <- matrix(o$par, nrow = 1)
    } else {
      Xnext <- X.start
    }
    
    if(trace) cat("Finding the next point (ALMC):", (proc.time() - t0)[3], "seconds\n")
    
    ### 2. Calculate ALC scores for all levels at the chosen point Xnext ###
    ALMC_scores <- numeric(L)
    ALMC_vals <- numeric(L)
    
    for (l in 1:L) {
      alc_func <- get(get_obj_func_name("ALC", l, L))
      val <- alc_func(Xnext, Xref = Xref, fit = fit, MC = MC, mc.sample = mc.sample, parallel = parallel, ncore = ncore)
      ALMC_vals[l] <- Icurrent - val
      ALMC_scores[l] <- (Icurrent - val) / sum(cost[1:l])
    }
    
    if(trace) cat("Calculating deduced variance:", (proc.time() - t0)[3], "seconds\n")
    
    best_level <- which.max(ALMC_scores)
    final_Xnext <- Xnext
    AL_out <- ALMC_vals / cumsum(cost)
    names(AL_out) <- paste0("ALMC", 1:L)
    
  } else { ### ALM, ALD, ALC ###
    AL_out <- vector("list", L)
    level_scores <- numeric(L)
    best_pars <- list()
    
    for (l in 1:L) {
      obj_fun <- get(get_obj_func_name(criterion, l, L))
      
      # 1. Grid Search on Xcand
      if (criterion == "ALC") {
        wrapper <- function(x) obj_fun(x, Xref = Xref, fit = fit, MC = MC, mc.sample = mc.sample)
      } else if (criterion == "ALD" && L == 3 && l < 3) {
        wrapper <- function(x) obj_fun(x, fit = fit, mc.sample = mc.sample)
      } else {
        wrapper <- function(x) obj_fun(x, fit = fit)
      }
      
      if (parallel) {
        vals <- foreach(i = 1:nrow(Xcand), .combine = c, .export = c("pred.GP", "covar.sep", "distance")) %dorng% {
          val <- wrapper(matrix(Xcand[i, ], nrow = 1))
          if (criterion %in% c("ALM", "ALD")) return(-val) else return(val)
        }
        attr(vals, "rng") <- NULL
        attr(vals, "doRNG_version") <- NULL
      } else {
        vals <- numeric(nrow(Xcand))
        for (i in 1:nrow(Xcand)) {
          val <- wrapper(matrix(Xcand[i, ], nrow = 1))
          vals[i] <- if (criterion %in% c("ALM", "ALD")) -val else val
        }
      }
      
      # 2. Optimization on Xcand
      if (criterion %in% c("ALM", "ALD")) best_idx <- which.max(vals) else best_idx <- which.min(vals)
      X.start <- matrix(Xcand[best_idx, ], nrow = 1)
      
      if (use_optim) {
        if (criterion == "ALC") {
          o <- optim(par = X.start, fn = obj_fun, method = "L-BFGS-B", lower = 0, upper = 1, fit = fit, Xref = Xref, MC = MC, mc.sample = mc.sample)
        } else if (criterion == "ALD" && L == 3 && l < 3) { # ALD V1 or V2 when L = 3
          o <- optim(par = X.start, fn = obj_fun, method = "L-BFGS-B", lower = 0, upper = 1, fit = fit, mc.sample = mc.sample)
        } else { # ALM or ALD V3
          o <- optim(par = X.start, fn = obj_fun, method = "L-BFGS-B", lower = 0, upper = 1, fit = fit)
        }
        
        best_val <- if (criterion %in% c("ALM", "ALD")) -o$value else o$value
        best_par <- o$par
      } else {
        best_val <- vals[best_idx]
        best_par <- X.start
      }
      
      # Calculate Score
      if (criterion == "ALC") {
        score <- (Icurrent - best_val) / sum(cost[1:l])
        grid_out <- vals
      } else {
        score <- best_val / sum(cost[1:l])
        grid_out <- vals
      }
      
      level_scores[l] <- score
      best_pars[[l]] <- best_par
      AL_out[[l]] <- grid_out / sum(cost[1:l])
    }
    
    names(AL_out) <- paste0(criterion, 1:L)
    best_level <- which.max(level_scores)
    final_Xnext <- matrix(best_pars[[best_level]], nrow=1)
    
    if(trace) cat("Calculation completed in:", (proc.time() - t0)[3], "seconds\n")
  }
  
  return(list(AL = AL_out, cost = cost, Xcand = Xcand, chosen = list(level = best_level, Xnext = final_Xnext), time = unname((proc.time() - t0)[3])))
}
