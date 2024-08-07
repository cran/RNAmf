#' fitting the model with squared exponential kernel.
#'
#' @param X vector or matrix of input locations.
#' @param y vector of response values.
#' @param g nugget parameter. Default is 1.490116e-08.
#' @param lower lower bound of theta. Default if 0.001.
#' @param upper upper bound of theta. Default if 1000.
#' @param Xscale logical indicating whether to scale X or not. Default is TRUE.
#' @param Yscale logical indicating whether to scale y or not. Only used if constant=FALSE. Default is TRUE.
#' @param constant logical indicating for constant mean (constant=TRUE) or zero mean (constant=FALSE). Default is FALSE.
#'
#' @return A list containing hyperparameters, covariance inverse matrix, X, y and logical inputs:
#' \itemize{
#'   \item \code{theta}: vector of lengthscale hyperparameter.
#'   \item \code{g}: copy of g.
#'   \item \code{Ki}: matrix of covariance inverse.
#'   \item \code{mu.hat}: optimized constant mean. If constant=FALSE, 0.
#'   \item \code{X}: copy of X. If Xscale=TRUE, scaled X.
#'   \item \code{y}: copy of y. If Yscale=TRUE, scaled y.
#'   \item \code{tau2hat}: estimated scale hyperparameter.
#'   \item \code{Xscale}: copy of Xscale.
#'   \item \code{Yscale}: copy of Yscale.
#'   \item \code{constant}: copy of constant.
#' }
#'
#' @importFrom plgp distance covar.sep
#' @importFrom stats optim quantile
#' @noRd
#' @keywords internal
#' @examples
#' \dontrun{
#' library(lhs)
#' ### synthetic function ###
#' f1 <- function(x) {
#'   sin(8 * pi * x)
#' }
#'
#' ### training data ###
#' n1 <- 15
#'
#' X1 <- maximinLHS(n1, 1)
#' y1 <- f1(X1)
#'
#' GP(X1, y1)
#' }
GP <- function(X, y, g = sqrt(.Machine$double.eps),
               Xscale = TRUE, Yscale = TRUE, constant = FALSE) {
  if (constant) {
    if (is.null(dim(X))) X <- matrix(X, ncol = 1)

    parbound <- function(XX) {
      XX <- matrix(XX, ncol = 1)

      lower <- -quantile(distance(XX)[lower.tri(distance(XX))], 0.05) / log(0.01) * (apply(XX, 2, range)[2, ] - apply(XX, 2, range)[1, ])^2
      upper <- -quantile(distance(XX)[lower.tri(distance(XX))], 0.95) / log(0.5) * (apply(XX, 2, range)[2, ] - apply(XX, 2, range)[1, ])^2

      return(c(lower, upper))
    }

    # hetGP way * 10^(d*3/4)
    XX <- (X - matrix(apply(X, 2, range)[1, ], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*% diag(1 / (apply(X, 2, range)[2, ] - apply(X, 2, range)[1, ]), ncol(X))
    lower <- pmax(10^(ncol(XX) * 3 / 4) * apply(XX, 2, parbound)[1, ], 0.2)
    upper <- 10^(ncol(XX) * 3 / 4) * apply(XX, 2, parbound)[2, ]
    init <- sqrt(lower * upper)

    # # hetGP way
    # XX <- (X - matrix(apply(X, 2, range)[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*% diag(1/(apply(X, 2, range)[2,] - apply(X, 2, range)[1,]), ncol(X))
    # lower <- apply(XX, 2, parbound)[1,]
    # upper <- apply(XX, 2, parbound)[2,]
    # init <- sqrt(lower*upper)

    if (Xscale) {
      X <- scale_inputs(X)
    } else {
      attr(X, "scaled:center") <- rep(0, ncol(X))
      attr(X, "scaled:scale") <- rep(1, ncol(X))
    }

    n <- length(y)

    nlsep <- function(par, X, Y) {
      theta <- par # lengthscale
      K <- covar.sep(X, d = theta, g = g)
      Ki <- solve(K)
      ldetK <- determinant(K, logarithm = TRUE)$modulus

      one.vec <- matrix(1, ncol = 1, nrow = n)
      mu.hat <- drop((t(one.vec) %*% Ki %*% Y) / (t(one.vec) %*% Ki %*% one.vec))

      tau2hat <- drop(t(Y - mu.hat) %*% Ki %*% (Y - mu.hat) / n)
      ll <- -(n / 2) * log(tau2hat) - (1 / 2) * ldetK
      return(drop(-ll))
    }

    gradnlsep <- function(par, X, Y) {
      theta <- par
      K <- covar.sep(X, d = theta, g = g)
      Ki <- solve(K)

      one.vec <- matrix(1, ncol = 1, nrow = n)
      mu.hat <- drop((t(one.vec) %*% Ki %*% Y) / (t(one.vec) %*% Ki %*% one.vec))

      KiY <- Ki %*% (Y - mu.hat)
      ## loop over theta components
      dlltheta <- rep(NA, length(theta))
      for (k in 1:length(dlltheta)) {
        dotK <- K * distance(X[, k]) / (theta[k]^2)
        dlltheta[k] <- (n / 2) * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dotK))
      }

      return(-c(dlltheta))
    }

    outg <- optim(init, nlsep, gradnlsep,
      method = "L-BFGS-B", lower = lower, upper = upper, X = X, Y = y
    )

    K <- covar.sep(X, d = outg$par, g = g)
    Ki <- solve(K)
    one.vec <- matrix(1, ncol = 1, nrow = n)
    mu.hat <- drop((t(one.vec) %*% Ki %*% y) / (t(one.vec) %*% Ki %*% one.vec))
    tau2hat <- drop(t(y - mu.hat) %*% Ki %*% (y - mu.hat) / nrow(X))
    theta <- outg$par
    names(theta) <- NULL

    return(list(theta = theta, g = g, Ki = Ki, mu.hat = mu.hat, X = X, y = y, tau2hat = tau2hat, Xscale = Xscale, constant = constant))
  } else {
    if (is.null(dim(X))) X <- matrix(X, ncol = 1)

    parbound <- function(XX) {
      XX <- matrix(XX, ncol = 1)

      lower <- -quantile(distance(XX)[lower.tri(distance(XX))], 0.05) / log(0.01) * (apply(XX, 2, range)[2, ] - apply(XX, 2, range)[1, ])^2
      upper <- -quantile(distance(XX)[lower.tri(distance(XX))], 0.95) / log(0.5) * (apply(XX, 2, range)[2, ] - apply(XX, 2, range)[1, ])^2

      return(c(lower, upper))
    }

    # hetGP way * 10^(d*3/4)
    XX <- (X - matrix(apply(X, 2, range)[1, ], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*% diag(1 / (apply(X, 2, range)[2, ] - apply(X, 2, range)[1, ]), ncol(X))
    lower <- pmax(10^(ncol(XX) * 3 / 4) * apply(XX, 2, parbound)[1, ], 0.2)
    upper <- 10^(ncol(XX) * 3 / 4) * apply(XX, 2, parbound)[2, ]
    init <- sqrt(lower * upper)

    # # hetGP way
    # XX <- (X - matrix(apply(X, 2, range)[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*% diag(1/(apply(X, 2, range)[2,] - apply(X, 2, range)[1,]), ncol(X))
    # lower <- - quantile(distance(XX)[lower.tri(distance(XX))], 0.05)/log(0.01) * (apply(XX, 2, range)[2,] - apply(XX, 2, range)[1,])^2
    # upper <- - quantile(distance(XX)[lower.tri(distance(XX))], 0.95)/log(0.5) * (apply(XX, 2, range)[2,] - apply(XX, 2, range)[1,])^2
    # init <- sqrt(lower*upper)

    if (Xscale) {
      X <- scale_inputs(X)
    } else {
      attr(X, "scaled:center") <- rep(0, ncol(X))
      attr(X, "scaled:scale") <- rep(1, ncol(X))
    }

    n <- length(y)
    if (Yscale) y <- scale(y, center = TRUE, scale = FALSE) # If use mean, don't scale

    nlsep <- function(par, X, Y) {
      theta <- par # lengthscale
      K <- covar.sep(X, d = theta, g = g)
      Ki <- solve(K)
      ldetK <- determinant(K, logarithm = TRUE)$modulus
      ll <- -(n / 2) * log(t(Y) %*% Ki %*% Y) - (1 / 2) * ldetK
      return(drop(-ll))
    }

    outg <- optim(init, nlsep, method = "L-BFGS-B", lower = lower, upper = upper, X = X, Y = y)

    K <- covar.sep(X, d = outg$par, g = g)
    Ki <- solve(K)
    tau2hat <- drop(t(y) %*% Ki %*% y / n)
    mu.hat <- 0
    theta <- outg$par
    names(theta) <- NULL

    return(list(theta = theta, g = g, Ki = Ki, mu.hat = mu.hat, X = X, y = y, tau2hat = tau2hat, Xscale = Xscale, Yscale = Yscale, constant = constant))
  }
}

#' predictive posterior mean and variance with squared exponential kernel.
#'
#' @param fit an object of class GP.
#' @param xnew vector or matrix of new input locations to predict.
#'
#' @return A list predictive posterior mean and variance:
#' \itemize{
#'   \item \code{mu}: vector of predictive posterior mean.
#'   \item \code{sig2}: vector of predictive posterior variance.
#' }
#'
#' @importFrom plgp covar.sep
#' @noRd
#' @keywords internal
#' @examples
#' \dontrun{
#' library(lhs)
#' ### synthetic function ###
#' f1 <- function(x) {
#'   sin(8 * pi * x)
#' }
#'
#' ### training data ###
#' n1 <- 15
#'
#' X1 <- maximinLHS(n1, 1)
#' y1 <- f1(X1)
#'
#' fit1 <- GP(X1, y1)
#'
#' ### test data ###
#' x <- seq(0, 1, 0.01)
#' pred.GP(fit1, x)
#' }
pred.GP <- function(fit, xnew) {
  constant <- fit$constant

  if (constant) {
    xnew <- as.matrix(xnew)

    Xscale <- fit$Xscale
    Ki <- fit$Ki
    theta <- fit$theta
    g <- fit$g
    X <- fit$X
    y <- fit$y
    tau2hat <- fit$tau2hat
    mu.hat <- fit$mu.hat

    if (Xscale) xnew <- scale_inputs(xnew, attr(X, "scaled:center"), attr(X, "scaled:scale"))

    KXX <- covar.sep(xnew, d = theta, g = g)
    KX <- covar.sep(xnew, X, d = theta, g = 0)

    mup2 <- mu.hat + KX %*% Ki %*% (y - mu.hat)
    Sigmap2 <- pmax(0, diag(tau2hat * (KXX - KX %*% Ki %*% t(KX))))

    return(list(mu = mup2, sig2 = Sigmap2))
  } else {
    xnew <- as.matrix(xnew)

    Xscale <- fit$Xscale
    Yscale <- fit$Yscale
    Ki <- fit$Ki
    theta <- fit$theta
    g <- fit$g
    X <- fit$X
    y <- fit$y
    tau2hat <- fit$tau2hat

    if (Xscale) xnew <- scale_inputs(xnew, attr(X, "scaled:center"), attr(X, "scaled:scale"))

    KXX <- covar.sep(xnew, d = theta, g = g)
    KX <- covar.sep(xnew, X, d = theta, g = 0)

    if (Yscale) mup2 <- KX %*% Ki %*% (y + attr(y, "scaled:center")) else mup2 <- KX %*% Ki %*% y
    Sigmap2 <- pmax(0, diag(tau2hat * (KXX - KX %*% Ki %*% t(KX))))

    return(list(mu = mup2, sig2 = Sigmap2))
  }
}
