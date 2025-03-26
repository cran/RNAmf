#' object to optimize the point by ALD criterion updating V1 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param fit an object of class RNAmf.
#' @return A negative V1 at Xcand.
#' @noRd
#'

obj.ALD_V1_2level <- function(Xcand, fit) { # low
  newx <- matrix(Xcand, nrow = 1)

  kernel <- fit$kernel
  constant <- fit$constant
  fit1 <- fit$fit1
  fit2 <- fit$fit2

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
      a <- Ci %*% (y2 - mu2)

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

      # mean
      predy <- mu2 + (exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))) *
                        1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                        exp(-(drop(outer(x.mu, w1.x2, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))) %*% a

      # var
      mat <- drop(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))) %o%
                    exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)]))))) * # common components
        1 / sqrt(1 + 4 * sig2 / theta[d + 1]) *
        exp(-(outer(w1.x2, w1.x2, FUN = "+") / 2 - matrix(x.mu, n, n))^2 / (theta[d + 1] / 2 + 2 * sig2)) *
        exp(-(outer(w1.x2, w1.x2, FUN = "-"))^2 / (2 * theta[d + 1]))

      VE <- - (predy - mu2)^2 + drop(t(a) %*% mat %*% a)
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
      # a <- Ci %*% (y2 + attr(y2, "scaled:center"))
      a <- Ci %*% y2

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

      # mean
      predy <- (exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))) *
                  1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                  exp(-(drop(outer(x.mu, w1.x2, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))) %*% a

      # var
      mat <- drop(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))) %o%
                    exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)]))))) * # common components
        1 / sqrt(1 + 4 * sig2 / theta[d + 1]) *
        exp(-(outer(w1.x2, w1.x2, FUN = "+") / 2 - matrix(x.mu, n, n))^2 / (theta[d + 1] / 2 + 2 * sig2)) *
        exp(-(outer(w1.x2, w1.x2, FUN = "-"))^2 / (2 * theta[d + 1]))

      VE <- - predy^2 + drop(t(a) %*% mat %*% a)
    }
  } else if (kernel == "matern1.5") {
    if (constant) {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
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
      a <- Ci %*% (y2 - mu2)

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

      # mean
      mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
      mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]

      lambda11 <- c(1, mua)
      lambda12 <- c(0, 1)
      lambda21 <- c(1, -mub)

      e1 <- cbind(matrix(1 - sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])
      e2 <- cbind(matrix(1 + sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])

      predy <- mu2 + drop(t(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                              (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                 (e1 %*% lambda11 * pnorm((mua - w1.x2) / sqrt(sig2)) +
                                    e1 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                                 exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                 (e2 %*% lambda21 * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                                    e2 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2))))) %*% a)

      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 1.5, theta = theta[d + 1])
      }
      mat <- drop(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5)) %o% t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))

      VE <- - (predy - mu2)^2 + drop(t(a) %*% mat %*% a)
    } else {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
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
      # a <- Ci %*% (y2 + attr(y2, "scaled:center"))
      a <- Ci %*% y2

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

      # mean
      mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
      mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]

      lambda11 <- c(1, mua)
      lambda12 <- c(0, 1)
      lambda21 <- c(1, -mub)

      e1 <- cbind(matrix(1 - sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])
      e2 <- cbind(matrix(1 + sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])

      predy <- drop(t(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                        (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                           (e1 %*% lambda11 * pnorm((mua - w1.x2) / sqrt(sig2)) +
                              e1 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                           exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                           (e2 %*% lambda21 * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                              e2 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2))))) %*% a)

      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 1.5, theta = theta[d + 1])
      }
      mat <- drop(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5)) %o% t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))

      VE <- - predy^2 + drop(t(a) %*% mat %*% a)
    }
  } else if (kernel == "matern2.5") {
    if (constant) {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
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
      a <- Ci %*% (y2 - mu2)

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

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

      predy <- mu2 + drop(t(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5)) *
                              (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                 (e1 %*% lambda11 * pnorm((mua - w1.x2) / sqrt(sig2)) +
                                    rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                                 exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                 (e2 %*% lambda21 * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                                    rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2))))) %*% a)

      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 2.5, theta = theta[d + 1])
      }
      mat <- drop(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5)) %o% t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))

      VE <- - (predy - mu2)^2 + drop(t(a) %*% mat %*% a)
    } else {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
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
      # a <- Ci %*% (y2 + attr(y2, "scaled:center"))
      a <- Ci %*% y2

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

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

      predy <- drop(t(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5)) *
                        (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                           (e1 %*% lambda11 * pnorm((mua - w1.x2) / sqrt(sig2)) +
                              rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                           exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                           (e2 %*% lambda21 * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                              rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2))))) %*% a)

      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 2.5, theta = theta[d + 1])
      }
      mat <- drop(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5)) %o% t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))

      VE <- - predy^2 + drop(t(a) %*% mat %*% a)
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

obj.ALD_V2_2level <- function(Xcand, fit) { # high
  newx <- matrix(Xcand, nrow = 1)

  kernel <- fit$kernel
  constant <- fit$constant
  fit1 <- fit$fit1
  fit2 <- fit$fit2

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
      a <- Ci %*% (y2 - mu2)

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

      # mean
      predy <- mu2 + (exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))) *
                        1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                        exp(-(drop(outer(x.mu, w1.x2, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))) %*% a

      # var
      mat <- drop(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))) %o%
                    exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)]))))) * # common components
        1 / sqrt(1 + 4 * sig2 / theta[d + 1]) *
        exp(-(outer(w1.x2, w1.x2, FUN = "+") / 2 - matrix(x.mu, n, n))^2 / (theta[d + 1] / 2 + 2 * sig2)) *
        exp(-(outer(w1.x2, w1.x2, FUN = "-"))^2 / (2 * theta[d + 1]))

      EV <- tau2hat - tau2hat * sum(diag(Ci %*% mat))
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
      # a <- Ci %*% (y2 + attr(y2, "scaled:center"))
      a <- Ci %*% y2

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

      # mean
      predy <- (exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))) *
                  1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                  exp(-(drop(outer(x.mu, w1.x2, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))) %*% a

      # var
      mat <- drop(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))) %o%
                    exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)]))))) * # common components
        1 / sqrt(1 + 4 * sig2 / theta[d + 1]) *
        exp(-(outer(w1.x2, w1.x2, FUN = "+") / 2 - matrix(x.mu, n, n))^2 / (theta[d + 1] / 2 + 2 * sig2)) *
        exp(-(outer(w1.x2, w1.x2, FUN = "-"))^2 / (2 * theta[d + 1]))

      EV <- tau2hat - tau2hat * sum(diag(Ci %*% mat))
    }
  } else if (kernel == "matern1.5") {
    if (constant) {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
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
      a <- Ci %*% (y2 - mu2)

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

      # mean
      mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
      mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]

      lambda11 <- c(1, mua)
      lambda12 <- c(0, 1)
      lambda21 <- c(1, -mub)

      e1 <- cbind(matrix(1 - sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])
      e2 <- cbind(matrix(1 + sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])

      predy <- mu2 + drop(t(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                              (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                 (e1 %*% lambda11 * pnorm((mua - w1.x2) / sqrt(sig2)) +
                                    e1 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                                 exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                 (e2 %*% lambda21 * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                                    e2 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2))))) %*% a)

      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 1.5, theta = theta[d + 1])
      }
      mat <- drop(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5)) %o% t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))

      EV <- tau2hat - tau2hat * sum(diag(Ci %*% mat))
    } else {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
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
      # a <- Ci %*% (y2 + attr(y2, "scaled:center"))
      a <- Ci %*% y2

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

      # mean
      mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
      mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]

      lambda11 <- c(1, mua)
      lambda12 <- c(0, 1)
      lambda21 <- c(1, -mub)

      e1 <- cbind(matrix(1 - sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])
      e2 <- cbind(matrix(1 + sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])

      predy <- drop(t(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                        (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                           (e1 %*% lambda11 * pnorm((mua - w1.x2) / sqrt(sig2)) +
                              e1 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                           exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                           (e2 %*% lambda21 * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                              e2 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2))))) %*% a)

      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 1.5, theta = theta[d + 1])
      }
      mat <- drop(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5)) %o% t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 1.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))

      EV <- tau2hat - tau2hat * sum(diag(Ci %*% mat))
    }
  } else if (kernel == "matern2.5") {
    if (constant) {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
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
      a <- Ci %*% (y2 - mu2)

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

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

      predy <- mu2 + drop(t(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5)) *
                              (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                 (e1 %*% lambda11 * pnorm((mua - w1.x2) / sqrt(sig2)) +
                                    rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                                 exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                                 (e2 %*% lambda21 * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                                    rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2))))) %*% a)

      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 2.5, theta = theta[d + 1])
      }
      mat <- drop(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5)) %o% t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))

      EV <- tau2hat - tau2hat * sum(diag(Ci %*% mat))
    } else {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.fit <- pred.matGP(fit1, newx)
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
      # a <- Ci %*% (y2 + attr(y2, "scaled:center"))
      a <- Ci %*% y2

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

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

      predy <- drop(t(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5)) *
                        (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                           (e1 %*% lambda11 * pnorm((mua - w1.x2) / sqrt(sig2)) +
                              rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2))) +
                           exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu)) / (2 * theta[d + 1]^2)) *
                           (e2 %*% lambda21 * pnorm((-mub + w1.x2) / sqrt(sig2)) +
                              rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2))))) %*% a)

      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 2.5, theta = theta[d + 1])
      }
      mat <- drop(t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5)) %o% t(cor.sep(t(newx), X2, theta[-(d + 1)], nu = 2.5))) * # constant depends on kernel
        outer(w1.x2, w1.x2, FUN = Vectorize(zeta))

      EV <- tau2hat - tau2hat * sum(diag(Ci %*% mat))
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

obj.ALD_V1_3level <- function(Xcand, fit, mc.sample, parallel = FALSE, ncore = 1) { # low

  kernel <- fit$kernel
  constant <- fit$constant
  fit.RNAmf_two_level <- fit$fit.RNAmf_two_level
  fit1 <- fit.RNAmf_two_level$fit1
  fit2 <- fit.RNAmf_two_level$fit2
  fit3 <- fit$fit3

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
    VEE.out <- foreach(i = 1:mc.sample, .combine = c) %dopar% {
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
        if (constant) {
          ### calculate the closed form ###
          X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
          w2.x3 <- fit3$X[, d + 1]
          y3 <- fit3$y
          n <- length(y3)
          theta <- fit3$theta
          tau2hat <- fit3$tau2hat
          mu3 <- fit3$mu.hat

          Ci <- fit3$Ki
          a <- Ci %*% (y3 - mu3)

          # ### scale new inputs ###
          # newx1 <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
          newx1 <- newx
          # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
          # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

          # mean
          VEE <- (exp(-distance(t(t(newx1) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)])))) *
                    1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                    exp(-(drop(outer(x.mu, w2.x3, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))) %*% a
        } else {
          ### calculate the closed form ###
          X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
          w2.x3 <- fit3$X[, d + 1]
          y3 <- fit3$y
          n <- length(y3)
          theta <- fit3$theta
          tau2hat <- fit3$tau2hat

          Ci <- fit3$Ki
          # a <- Ci %*% (y3 + attr(y3, "scaled:center"))
          a <- Ci %*% y3

          # ### scale new inputs ###
          # newx1 <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
          newx1 <- newx
          # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
          # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

          # mean
          VEE <- (exp(-distance(t(t(newx1) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)])))) *
                    1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                    exp(-(drop(outer(x.mu, w2.x3, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))) %*% a
        }
      } else if (kernel == "matern1.5") {
        if (constant) {
          ### calculate the closed form ###
          X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
          w2.x3 <- fit3$X[, d + 1]
          y3 <- fit3$y
          n <- length(y3)
          theta <- fit3$theta
          tau2hat <- fit3$tau2hat
          mu3 <- fit3$mu.hat

          Ci <- fit3$Ki
          a <- Ci %*% (y3 - mu3)

          # ### scale new inputs ###
          # newx1 <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
          newx1 <- newx
          # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
          # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

          # mean
          mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
          mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]

          lambda11 <- c(1, mua)
          lambda12 <- c(0, 1)
          lambda21 <- c(1, -mub)

          e1 <- cbind(matrix(1 - sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])
          e2 <- cbind(matrix(1 + sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])

          VEE <- drop(t(t(cor.sep(t(newx1), X3, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                          (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                e1 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                             exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                e2 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2))))) %*% a)
        } else {
          ### calculate the closed form ###
          X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
          w2.x3 <- fit3$X[, d + 1]
          y3 <- fit3$y
          n <- length(y3)
          theta <- fit3$theta
          tau2hat <- fit3$tau2hat

          Ci <- fit3$Ki
          # a <- Ci %*% (y3 + attr(y3, "scaled:center"))
          a <- Ci %*% y3

          # ### scale new inputs ###
          # newx1 <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
          newx1 <- newx
          # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
          # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

          # mean
          mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
          mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]

          lambda11 <- c(1, mua)
          lambda12 <- c(0, 1)
          lambda21 <- c(1, -mub)

          e1 <- cbind(matrix(1 - sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])
          e2 <- cbind(matrix(1 + sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])

          VEE <- drop(t(t(cor.sep(t(newx1), X3, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                          (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                e1 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                             exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                e2 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2))))) %*% a)
        }
      } else if (kernel == "matern2.5") {
        if (constant) {
          ### calculate the closed form ###
          X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
          w2.x3 <- fit3$X[, d + 1]
          y3 <- fit3$y
          n <- length(y3)
          theta <- fit3$theta
          tau2hat <- fit3$tau2hat
          mu3 <- fit3$mu.hat

          Ci <- fit3$Ki
          a <- Ci %*% (y3 - mu3)

          # ### scale new inputs ###
          # newx1 <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
          newx1 <- newx
          # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
          # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

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

          VEE <- drop(t(t(cor.sep(t(newx1), X3, theta[-(d + 1)], nu = 2.5)) *
                          (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                             exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2))))) %*% a)
        } else {
          ### calculate the closed form ###
          X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
          w2.x3 <- fit3$X[, d + 1]
          y3 <- fit3$y
          n <- length(y3)
          theta <- fit3$theta
          tau2hat <- fit3$tau2hat

          Ci <- fit3$Ki
          # a <- Ci %*% (y3 + attr(y3, "scaled:center"))
          a <- Ci %*% y3

          # ### scale new inputs ###
          # newx1 <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
          newx1 <- newx
          # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
          # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

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

          VEE <- drop(t(t(cor.sep(t(newx1), X3, theta[-(d + 1)], nu = 2.5)) *
                          (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                             exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2))))) %*% a)
        }
      }
      return(VEE) # to maximize the V1.
    }
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
        if (constant) {
          ### calculate the closed form ###
          X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
          w2.x3 <- fit3$X[, d + 1]
          y3 <- fit3$y
          n <- length(y3)
          theta <- fit3$theta
          tau2hat <- fit3$tau2hat
          mu3 <- fit3$mu.hat

          Ci <- fit3$Ki
          a <- Ci %*% (y3 - mu3)

          # ### scale new inputs ###
          # newx1 <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
          newx1 <- newx
          # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
          # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

          # mean
          VEE <- (exp(-distance(t(t(newx1) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)])))) *
                    1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                    exp(-(drop(outer(x.mu, w2.x3, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))) %*% a
        } else {
          ### calculate the closed form ###
          X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
          w2.x3 <- fit3$X[, d + 1]
          y3 <- fit3$y
          n <- length(y3)
          theta <- fit3$theta
          tau2hat <- fit3$tau2hat

          Ci <- fit3$Ki
          # a <- Ci %*% (y3 + attr(y3, "scaled:center"))
          a <- Ci %*% y3

          # ### scale new inputs ###
          # newx1 <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
          newx1 <- newx
          # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
          # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

          # mean
          VEE <- (exp(-distance(t(t(newx1) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)])))) *
                    1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                    exp(-(drop(outer(x.mu, w2.x3, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))) %*% a
        }
      } else if (kernel == "matern1.5") {
        if (constant) {
          ### calculate the closed form ###
          X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
          w2.x3 <- fit3$X[, d + 1]
          y3 <- fit3$y
          n <- length(y3)
          theta <- fit3$theta
          tau2hat <- fit3$tau2hat
          mu3 <- fit3$mu.hat

          Ci <- fit3$Ki
          a <- Ci %*% (y3 - mu3)

          # ### scale new inputs ###
          # newx1 <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
          newx1 <- newx
          # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
          # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

          # mean
          mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
          mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]

          lambda11 <- c(1, mua)
          lambda12 <- c(0, 1)
          lambda21 <- c(1, -mub)

          e1 <- cbind(matrix(1 - sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])
          e2 <- cbind(matrix(1 + sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])

          VEE <- drop(t(t(cor.sep(t(newx1), X3, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                          (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                e1 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                             exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                e2 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2))))) %*% a)
        } else {
          ### calculate the closed form ###
          X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
          w2.x3 <- fit3$X[, d + 1]
          y3 <- fit3$y
          n <- length(y3)
          theta <- fit3$theta
          tau2hat <- fit3$tau2hat

          Ci <- fit3$Ki
          # a <- Ci %*% (y3 + attr(y3, "scaled:center"))
          a <- Ci %*% y3

          # ### scale new inputs ###
          # newx1 <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
          newx1 <- newx
          # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
          # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

          # mean
          mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
          mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]

          lambda11 <- c(1, mua)
          lambda12 <- c(0, 1)
          lambda21 <- c(1, -mub)

          e1 <- cbind(matrix(1 - sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])
          e2 <- cbind(matrix(1 + sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])

          VEE <- drop(t(t(cor.sep(t(newx1), X3, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                          (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                e1 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                             exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                e2 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2))))) %*% a)
        }
      } else if (kernel == "matern2.5") {
        if (constant) {
          ### calculate the closed form ###
          X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
          w2.x3 <- fit3$X[, d + 1]
          y3 <- fit3$y
          n <- length(y3)
          theta <- fit3$theta
          tau2hat <- fit3$tau2hat
          mu3 <- fit3$mu.hat

          Ci <- fit3$Ki
          a <- Ci %*% (y3 - mu3)

          # ### scale new inputs ###
          # newx1 <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
          newx1 <- newx
          # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
          # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

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

          VEE <- drop(t(t(cor.sep(t(newx1), X3, theta[-(d + 1)], nu = 2.5)) *
                          (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                             exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2))))) %*% a)
        } else {
          ### calculate the closed form ###
          X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
          w2.x3 <- fit3$X[, d + 1]
          y3 <- fit3$y
          n <- length(y3)
          theta <- fit3$theta
          tau2hat <- fit3$tau2hat

          Ci <- fit3$Ki
          # a <- Ci %*% (y3 + attr(y3, "scaled:center"))
          a <- Ci %*% y3

          # ### scale new inputs ###
          # newx1 <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
          newx1 <- newx
          # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
          # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

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

          VEE <- drop(t(t(cor.sep(t(newx1), X3, theta[-(d + 1)], nu = 2.5)) *
                          (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                             exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                             (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2))))) %*% a)
        }
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

obj.ALD_V2_3level <- function(Xcand, fit, mc.sample, parallel = FALSE, ncore = 1) { # med

  V <- predict(fit, Xcand)$sig2
  EVE <- V +
    obj.ALD_V1_3level(Xcand, fit, mc.sample, parallel = FALSE, ncore = 1) +# -V1
    obj.ALD_V3_3level(Xcand, fit) # -V3

  return(-EVE) # to maximize the V2.
}


#' object to optimize the point by ALD criterion updating V3 with three levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param fit an object of class RNAmf.
#' @return A negative V2 at Xcand.
#' @noRd
#'

obj.ALD_V3_3level <- function(Xcand, fit) { # high
  newx <- matrix(Xcand, nrow = 1)

  kernel <- fit$kernel
  constant <- fit$constant
  fit.RNAmf_two_level <- fit$fit.RNAmf_two_level
  fit1 <- fit.RNAmf_two_level$fit1
  fit2 <- fit.RNAmf_two_level$fit2
  fit3 <- fit$fit3

  if (kernel == "sqex") {
    if (constant) {
      pred.RNAmf_two_level <- predict(fit.RNAmf_two_level, newx)
      x.mu <- pred.RNAmf_two_level$mu
      sig2 <- pred.RNAmf_two_level$sig2

      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)

      ### calculate the closed form ###
      X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
      w2.x3 <- fit3$X[, d + 1]
      y3 <- fit3$y
      n <- length(y3)
      theta <- fit3$theta
      tau2hat <- fit3$tau2hat
      mu3 <- fit3$mu.hat

      Ci <- fit3$Ki
      a <- Ci %*% (y3 - mu3)

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

      # mean
      predy <- mu3 + (exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)])))) *
                        1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                        exp(-(drop(outer(x.mu, w2.x3, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))) %*% a

      # var
      mat <- drop(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)])))) %o%
                    exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)]))))) * # common components
        1 / sqrt(1 + 4 * sig2 / theta[d + 1]) *
        exp(-(outer(w2.x3, w2.x3, FUN = "+") / 2 - matrix(x.mu, n, n))^2 / (theta[d + 1] / 2 + 2 * sig2)) *
        exp(-(outer(w2.x3, w2.x3, FUN = "-"))^2 / (2 * theta[d + 1]))

      EEV <- tau2hat - tau2hat * sum(diag(Ci %*% mat))
    } else {
      pred.RNAmf_two_level <- predict(fit.RNAmf_two_level, newx)
      x.mu <- pred.RNAmf_two_level$mu
      sig2 <- pred.RNAmf_two_level$sig2

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
      # a <- Ci %*% (y3 + attr(y3, "scaled:center"))
      a <- Ci %*% y3

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

      # mean
      predy <- (exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)])))) *
                  1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                  exp(-(drop(outer(x.mu, w2.x3, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))) %*% a

      # var
      mat <- drop(exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)])))) %o%
                    exp(-distance(t(t(newx) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)]))))) * # common components
        1 / sqrt(1 + 4 * sig2[i] / theta[d + 1]) *
        exp(-(outer(w2.x3, w2.x3, FUN = "+") / 2 - matrix(x.mu, n, n))^2 / (theta[d + 1] / 2 + 2 * sig2)) *
        exp(-(outer(w2.x3, w2.x3, FUN = "-"))^2 / (2 * theta[d + 1]))

      EEV <- tau2hat - tau2hat * sum(diag(Ci %*% mat))
    }
  } else if (kernel == "matern1.5") {
    if (constant) {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.RNAmf_two_level <- predict(fit.RNAmf_two_level, newx)
      x.mu <- pred.RNAmf_two_level$mu
      sig2 <- pred.RNAmf_two_level$sig2

      ### calculate the closed form ###
      X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
      w2.x3 <- fit3$X[, d + 1]
      y3 <- fit3$y
      n <- length(y3)
      theta <- fit3$theta
      tau2hat <- fit3$tau2hat
      mu3 <- fit3$mu.hat

      Ci <- fit3$Ki
      a <- Ci %*% (y3 - mu3)

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

      # mean
      mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
      mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]

      lambda11 <- c(1, mua)
      lambda12 <- c(0, 1)
      lambda21 <- c(1, -mub)

      e1 <- cbind(matrix(1 - sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])
      e2 <- cbind(matrix(1 + sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])

      predy <- mu3 + drop(t(t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                              (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                 (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                    e1 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                                 exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                 (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                    e2 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2))))) %*% a)

      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 1.5, theta = theta[d + 1])
      }
      mat <- drop(t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 1.5)) %o% t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 1.5))) * # constant depends on kernel
        outer(w2.x3, w2.x3, FUN = Vectorize(zeta))

      EEV <- tau2hat - tau2hat * sum(diag(Ci %*% mat))
    } else {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.RNAmf_two_level <- predict(fit.RNAmf_two_level, newx)
      x.mu <- pred.RNAmf_two_level$mu
      sig2 <- pred.RNAmf_two_level$sig2

      ### calculate the closed form ###
      X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
      w2.x3 <- fit3$X[, d + 1]
      y3 <- fit3$y
      n <- length(y3)
      theta <- fit3$theta
      tau2hat <- fit3$tau2hat

      Ci <- fit3$Ki
      # a <- Ci %*% (y3 + attr(y3, "scaled:center"))
      a <- Ci %*% y3

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

      # mean
      mua <- x.mu - sqrt(3) * sig2 / theta[d + 1]
      mub <- x.mu + sqrt(3) * sig2 / theta[d + 1]

      lambda11 <- c(1, mua)
      lambda12 <- c(0, 1)
      lambda21 <- c(1, -mub)

      e1 <- cbind(matrix(1 - sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])
      e2 <- cbind(matrix(1 + sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])

      predy <- drop(t(t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
                        (exp((3 * sig2 + 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                           (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2)) +
                              e1 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                           exp((3 * sig2 - 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                           (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                              e2 %*% lambda12 * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2))))) %*% a)

      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 1.5, theta = theta[d + 1])
      }
      mat <- drop(t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 1.5)) %o% t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 1.5))) * # constant depends on kernel
        outer(w2.x3, w2.x3, FUN = Vectorize(zeta))

      EEV <- tau2hat - tau2hat * sum(diag(Ci %*% mat))
    }
  } else if (kernel == "matern2.5") {
    if (constant) {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.RNAmf_two_level <- predict(fit.RNAmf_two_level, newx)
      x.mu <- pred.RNAmf_two_level$mu
      sig2 <- pred.RNAmf_two_level$sig2

      ### calculate the closed form ###
      X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
      w2.x3 <- fit3$X[, d + 1]
      y3 <- fit3$y
      n <- length(y3)
      theta <- fit3$theta
      tau2hat <- fit3$tau2hat
      mu3 <- fit3$mu.hat

      Ci <- fit3$Ki
      a <- Ci %*% (y3 - mu3)

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

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

      predy <- mu3 + drop(t(t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 2.5)) *
                              (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                 (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2)) +
                                    rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                                 exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                                 (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                                    rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2))))) %*% a)

      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 2.5, theta = theta[d + 1])
      }
      mat <- drop(t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 2.5)) %o% t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 2.5))) * # constant depends on kernel
        outer(w2.x3, w2.x3, FUN = Vectorize(zeta))

      EEV <- tau2hat - tau2hat * sum(diag(Ci %*% mat))
    } else {
      d <- ncol(fit1$X)
      newx <- matrix(newx, ncol = d)
      pred.RNAmf_two_level <- predict(fit.RNAmf_two_level, newx)
      x.mu <- pred.RNAmf_two_level$mu
      sig2 <- pred.RNAmf_two_level$sig2

      ### calculate the closed form ###
      X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
      w2.x3 <- fit3$X[, d + 1]
      y3 <- fit3$y
      n <- length(y3)
      theta <- fit3$theta
      tau2hat <- fit3$tau2hat

      Ci <- fit3$Ki
      # a <- Ci %*% (y3 + attr(y3, "scaled:center"))
      a <- Ci %*% y3

      # ### scale new inputs ###
      # newx <- scale_inputs(newx, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
      # x.mu <- scale_inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
      # sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

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

      predy <- drop(t(t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 2.5)) *
                        (exp((5 * sig2 + 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                           (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2)) +
                              rowSums(e1 * lambda12) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2))) +
                           exp((5 * sig2 - 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu)) / (2 * theta[d + 1]^2)) *
                           (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2)) +
                              rowSums(e2 * lambda22) * sqrt(sig2) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2))))) %*% a)

      # var
      zeta <- function(x, y) {
        zetafun(w1 = x, w2 = y, m = x.mu, s = sig2, nu = 2.5, theta = theta[d + 1])
      }
      mat <- drop(t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 2.5)) %o% t(cor.sep(t(newx), X3, theta[-(d + 1)], nu = 2.5))) * # constant depends on kernel
        outer(w2.x3, w2.x3, FUN = Vectorize(zeta))

      EEV <- tau2hat - tau2hat * sum(diag(Ci %*% mat))
    }
  }

  return(-EEV) # to maximize the V3.
}


#' @title find the next point by ALD criterion
#'
#' @description The function acquires the new point by the Active learning Decomposition (ALD) criterion.
#' It calculates the ALD criterion \eqn{\frac{V_l(\bm{x})}{\sum^l_{j=1}C_j}},
#' where \eqn{V_l(\bm{x})} is the contribution of GP emulator
#' at each fidelity level \eqn{l} and \eqn{C_j} is the simulation cost at level \eqn{j}.
#' For details, see Heo and Sung (2024, <\doi{https://doi.org/10.1080/00401706.2024.2376173}>).
#'
#' A new point is acquired on \code{Xcand}. If \code{Xcand=NULL}, a new point is acquired on unit hypercube \eqn{[0,1]^d}.
#'
#' @param Xcand vector or matrix of candidate set which could be added into the current design only used when \code{optim=FALSE}. \code{Xcand} is the set of the points where ALD criterion is evaluated. If \code{Xcand=NULL}, \eqn{100 \times d} number of points from 0 to 1 are generated by Latin hypercube design. Default is \code{NULL}.
#' @param fit object of class \code{RNAmf}.
#' @param mc.sample a number of mc samples generated for the MC approximation in 3 levels case. Default is \code{100}.
#' @param cost vector of the costs for each level of fidelity. If \code{cost=NULL}, total costs at all levels would be 1. \code{cost} is encouraged to have an ascending order of positive value. Default is \code{NULL}.
#' @param optim logical indicating whether to optimize AL criterion by \code{optim}'s gradient-based \code{L-BFGS-B} method. If \code{optim=TRUE}, \eqn{5 \times d} starting points are generated by Latin hypercube design for optimization. If \code{optim=FALSE}, AL criterion is optimized on the \code{Xcand}. Default is \code{TRUE}.
#' @param parallel logical indicating whether to compute the AL criterion in parallel or not. If \code{parallel=TRUE}, parallel computation is utilized. Default is \code{FALSE}.
#' @param ncore a number of core for parallel. It is only used if \code{parallel=TRUE}. Default is 1.
#' @param trace logical indicating whether to print the computational time for each step. If \code{trace=TRUE}, the computation time for each step is printed. Default is \code{FALSE}.
#' @return
#' \itemize{
#'   \item \code{ALD}: list of ALD criterion computed at each point of \code{Xcand} at each level if \code{optim=FALSE}. If \code{optim=TRUE}, \code{ALD} returns \code{NULL}.
#'   \item \code{cost}: a copy of \code{cost}.
#'   \item \code{Xcand}: a copy of \code{Xcand}.
#'   \item \code{chosen}: list of chosen level and point.
#'   \item \code{time}: a scalar of the time for the computation.
#' }
#' @importFrom plgp covar.sep
#' @importFrom lhs maximinLHS
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom doParallel stopImplicitCluster
#' @usage ALD_RNAmf(Xcand = NULL, fit, mc.sample = 100, cost = NULL,
#' optim = TRUE, parallel = FALSE, ncore = 1, trace=TRUE)
#' @export
#' @examples
#' \donttest{
#' library(lhs)
#' library(doParallel)
#' library(foreach)
#'
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
#' ### n1 and n2 might be changed from NestedX ###
#' ### assign n1 and n2 again ###
#' n1 <- nrow(X1)
#' n2 <- nrow(X2)
#'
#' y1 <- f1(X1)
#' y2 <- f2(X2)
#'
#' ### n=100 uniform test data ###
#' x <- seq(0, 1, length.out = 100)
#'
#' ### fit an RNAmf ###
#' fit.RNAmf <- RNAmf_two_level(X1, y1, X2, y2, kernel = "sqex")
#'
#' ### predict ###
#' predy <- predict(fit.RNAmf, x)$mu
#' predsig2 <- predict(fit.RNAmf, x)$sig2
#'
#' ### active learning with optim=TRUE ###
#' ald.RNAmf.optim <- ALD_RNAmf(
#'   Xcand = x, fit.RNAmf, cost = cost,
#'   optim = TRUE, parallel = TRUE, ncore = 2
#' )
#' print(ald.RNAmf.optim$time) # computation time of optim=TRUE
#'
#' ### active learning with optim=FALSE ###
#' ald.RNAmf <- ALD_RNAmf(
#'   Xcand = x, fit.RNAmf, cost = cost,
#'   optim = FALSE, parallel = TRUE, ncore = 2
#' )
#' print(ald.RNAmf$time) # computation time of optim=FALSE
#'
#' ### visualize ALD ###
#' oldpar <- par(mfrow = c(1, 2))
#' plot(x, ald.RNAmf$ALD$ALD1,
#'   type = "l", lty = 2,
#'   xlab = "x", ylab = "ALD criterion at the low-fidelity level",
#'   ylim = c(min(c(ald.RNAmf$ALD$ALD1, ald.RNAmf$ALD$ALD2)),
#'            max(c(ald.RNAmf$ALD$ALD1, ald.RNAmf$ALD$ALD2)))
#' )
#' points(ald.RNAmf$chosen$Xnext,
#'   ald.RNAmf$ALD$ALD1[which(x == drop(ald.RNAmf$chosen$Xnext))],
#'   pch = 16, cex = 1, col = "red"
#' )
#' plot(x, ald.RNAmf$ALD$ALD2,
#'   type = "l", lty = 2,
#'   xlab = "x", ylab = "ALD criterion at the high-fidelity level",
#'   ylim = c(min(c(ald.RNAmf$ALD$ALD1, ald.RNAmf$ALD$ALD2)),
#'            max(c(ald.RNAmf$ALD$ALD1, ald.RNAmf$ALD$ALD2)))
#' )
#' par(oldpar)}
#'
ALD_RNAmf <- function(Xcand = NULL, fit, mc.sample = 100, cost = NULL, optim = TRUE, parallel = FALSE, ncore = 1, trace=TRUE) {
  t0 <- proc.time()
  ### check the object ###
  if (!inherits(fit, "RNAmf")) {
    stop("The object is not of class \"RNAmf\" \n")
  }
  if (length(cost) != fit$level) stop("The length of cost should be the level of object")

  ### ALD ###
  if (fit$level == 2) { # level 2
    if (!is.null(cost) & cost[1] >= cost[2]) {
      warning("If the cost for high-fidelity is cheaper, acquire the high-fidelity")
    } else if (is.null(cost)) {
      cost <- c(1, 0)
    }
    if (parallel) registerDoParallel(ncore)

    fit1 <- fit$fit1
    fit2 <- fit$fit2
    constant <- fit$constant
    kernel <- fit$kernel
    g <- fit1$g

    ### Generate the candidate set ###
    if (optim){ # optim = TRUE
      Xcand <- randomLHS(5*ncol(fit1$X), ncol(fit1$X))
    }else{ # optim = FALSE
      if (is.null(Xcand)){
        Xcand <- randomLHS(100*ncol(fit1$X), ncol(fit1$X))
      }else if(is.null(dim(Xcand))){
        Xcand <- matrix(Xcand, ncol = 1)
      }
    }
    # if (ncol(Xcand) != dim(fit$fit1$X)[2]) stop("The dimension of candidate set should be equal to the dimension of the design")

    ### Calculate the contribution of GP emulator at each level ###
    t1 <- proc.time()
    if (parallel) {
      optm.mat <- foreach(i = 1:nrow(Xcand), .combine = cbind) %dopar% {
        newx <- matrix(Xcand[i, ], nrow = 1)
        return(c(
          -obj.ALD_V1_2level(newx, fit = fit),
          -obj.ALD_V2_2level(newx, fit = fit)
        ))
      }
    } else {
      optm.mat <- rbind(c(rep(0, nrow(Xcand))), c(rep(0, nrow(Xcand))))
      for (i in 1:nrow(Xcand)) {
        print(paste(i, nrow(Xcand), sep = "/"))
        newx <- matrix(Xcand[i, ], nrow = 1)

        optm.mat[1, i] <- -obj.ALD_V1_2level(newx, fit = fit)
        optm.mat[2, i] <- -obj.ALD_V2_2level(newx, fit = fit)
      }
    }
    t2 <-proc.time()
    if(trace) cat("Calculating the contribution of GP emulator at each level:", (t2 - t1)[3], "seconds\n")

    ### Find the next point ###
    if (optim) {
      X.start <- matrix(Xcand[which.max(optm.mat[1, ]), ], nrow = 1)
      optim.out <- optim(X.start, obj.ALD_V1_2level, method = "L-BFGS-B", lower = 0, upper = 1, fit = fit)
      Xnext.1 <- optim.out$par
      ALD.1 <- -optim.out$value
      t3 <-proc.time()
      if(trace) cat("Running optim for level 1:", (t3 - t2)[3], "seconds\n")

      X.start <- matrix(Xcand[which.max(optm.mat[2, ]), ], nrow = 1)
      optim.out <- optim(X.start, obj.ALD_V2_2level, method = "L-BFGS-B", lower = 0, upper = 1, fit = fit)
      Xnext.2 <- optim.out$par
      ALD.2 <- -optim.out$value
      t4 <-proc.time()
      if(trace) cat("Running optim for level 2:", (t4 - t3)[3], "seconds\n")

      ALDvalue <- c(ALD.1, ALD.2) / c(cost[1], cost[1] + cost[2])
      if (ALDvalue[2] > ALDvalue[1]) {
        level <- 2
        Xnext <- Xnext.2
      } else {
        level <- 1
        Xnext <- Xnext.1
      }
    } else {
      ALDvalue <- c(max(optm.mat[1, ]), max(optm.mat[2, ])) / c(cost[1], cost[1] + cost[2])
      if (ALDvalue[2] > ALDvalue[1]) {
        level <- 2
        Xnext <- matrix(Xcand[which.max(optm.mat[2, ]), ], nrow = 1)
      } else {
        level <- 1
        Xnext <- matrix(Xcand[which.max(optm.mat[1, ]), ], nrow = 1)
      }
    }

    chosen <- list(
      "level" = level, # next level
      "Xnext" = Xnext
    ) # next point

    ALD <- list(ALD1 = optm.mat[1, ] / cost[1], ALD2 = optm.mat[2, ] / (cost[1] + cost[2]))
  } else if (fit$level == 3) { # level 3

    if (!is.null(cost) & (cost[1] >= cost[2] | cost[2] >= cost[3])) {
      warning("If the cost for high-fidelity is cheaper, acquire the high-fidelity")
    } else if (is.null(cost)) {
      cost <- c(1, 0, 0)
    }
    if (parallel) registerDoParallel(ncore)

    fit_two_level <- fit$fit.RNAmf_two_level
    fit1 <- fit_two_level$fit1
    fit2 <- fit_two_level$fit2
    fit3 <- fit$fit3
    constant <- fit$constant
    kernel <- fit$kernel
    g <- fit1$g

    ### Generate the candidate set ###
    if (optim){ # optim = TRUE
      Xcand <- randomLHS(5*ncol(fit1$X), ncol(fit1$X))
    }else{ # optim = FALSE
      if (is.null(Xcand)){
        Xcand <- randomLHS(100*ncol(fit1$X), ncol(fit1$X))
      }else if(is.null(dim(Xcand))){
        Xcand <- matrix(Xcand, ncol = 1)
      }
    }
    # if (ncol(Xcand) != dim(fit$fit1$X)[2]) stop("The dimension of candidate set should be equal to the dimension of the design")

    ### Calculate the contribution of GP emulator at each level ###
    t1 <- proc.time()
    if (parallel) {
      optm.mat <- foreach(i = 1:nrow(Xcand), .combine = cbind) %dopar% {
        newx <- matrix(Xcand[i, ], nrow = 1)

        return(c(
          -obj.ALD_V1_3level(newx, fit = fit, mc.sample = mc.sample),
          -obj.ALD_V2_3level(newx, fit = fit, mc.sample = mc.sample),
          -obj.ALD_V3_3level(newx, fit = fit)
        ))
      }
    } else {
      optm.mat <- rbind(c(rep(0, nrow(Xcand))), c(rep(0, nrow(Xcand))), c(rep(0, nrow(Xcand))))
      for (i in 1:nrow(Xcand)) {
        print(paste(i, nrow(Xcand), sep = "/"))
        newx <- matrix(Xcand[i, ], nrow = 1)

        optm.mat[1, i] <- -obj.ALD_V1_3level(newx, fit = fit, mc.sample = mc.sample)
        optm.mat[2, i] <- -obj.ALD_V2_3level(newx, fit = fit, mc.sample = mc.sample)
        optm.mat[3, i] <- -obj.ALD_V3_3level(newx, fit = fit)
      }
    }
    t2 <- proc.time()
    if(trace) cat("Calculating the contribution of GP emulator at each level:", (t2 - t1)[3], "seconds\n")

    ### Find the next point ###
    if (optim) {
      X.start <- matrix(Xcand[which.max(optm.mat[1, ]), ], nrow = 1)
      optim.out <- optim(X.start, obj.ALD_V1_3level, method = "L-BFGS-B", lower = 0, upper = 1, fit = fit, mc.sample = mc.sample, parallel = parallel, ncore = ncore)
      Xnext.1 <- optim.out$par
      ALD.1 <- -optim.out$value
      t3 <- proc.time()
      if(trace) cat("Running optim for level 1:", (t3 - t2)[3], "seconds\n")

      X.start <- matrix(Xcand[which.max(optm.mat[2, ]), ], nrow = 1)
      optim.out <- optim(X.start, obj.ALD_V2_3level, method = "L-BFGS-B", lower = 0, upper = 1, fit = fit, mc.sample = mc.sample, parallel = parallel, ncore = ncore)
      Xnext.2 <- optim.out$par
      ALD.2 <- -optim.out$value
      t4 <- proc.time()
      if(trace) cat("Running optim for level 2:", (t4 - t3)[3], "seconds\n")

      X.start <- matrix(Xcand[which.max(optm.mat[3, ]), ], nrow = 1)
      optim.out <- optim(X.start, obj.ALD_V3_3level, method = "L-BFGS-B", lower = 0, upper = 1, fit = fit)
      Xnext.3 <- optim.out$par
      ALD.3 <- -optim.out$value
      t5 <- proc.time()
      if(trace) cat("Running optim for level 3:", (t5 - t4)[3], "seconds\n")

      ALDvalue <- c(ALD.1, ALD.2, ALD.3) / c(cost[1], cost[1] + cost[2], cost[1] + cost[2] + cost[3])
      if (ALDvalue[3] > ALDvalue[2]) {
        level <- 3
        Xnext <- Xnext.3
      } else if (ALDvalue[2] > ALDvalue[1]) {
        level <- 2
        Xnext <- Xnext.2
      } else {
        level <- 1
        Xnext <- Xnext.1
      }
    } else {
      ALDvalue <- c(max(optm.mat[1, ]), max(optm.mat[2, ]), max(optm.mat[3, ])) / c(cost[1], cost[1] + cost[2], cost[1] + cost[2] + cost[3])
      if (ALDvalue[3] > ALDvalue[2]) {
        level <- 3
        Xnext <- matrix(Xcand[which.max(optm.mat[3, ]), ], nrow = 1)
      } else if (ALDvalue[2] > ALDvalue[1]) {
        level <- 2
        Xnext <- matrix(Xcand[which.max(optm.mat[2, ]), ], nrow = 1)
      } else {
        level <- 1
        Xnext <- matrix(Xcand[which.max(optm.mat[1, ]), ], nrow = 1)
      }
    }

    chosen <- list(
      "level" = level, # next level
      "Xnext" = Xnext
    ) # next point

    ALD <- list(ALD1 = optm.mat[1, ] / cost[1], ALD2 = optm.mat[2, ] / (cost[1] + cost[2]), ALD3 = optm.mat[3, ] / (cost[1] + cost[2] + cost[3]))
  } else {
    stop("level is not 2")
  }

  if (parallel) stopImplicitCluster()
  if (optim) ALD <- NULL

  return(list(ALD = ALD, cost = cost, Xcand = Xcand, chosen = chosen, time = (proc.time() - t0)[3]))
}
