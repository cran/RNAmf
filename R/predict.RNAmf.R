#' @title prediction of the RNAmf emulator with multiple fidelity levels.
#'
#' @description The function computes the posterior mean and variance of RNA models with multiple fidelity levels
#' by fitted model from \code{\link{RNAmf}}.
#'
#' @seealso \code{\link{RNAmf}} for model fitting.
#'
#' @details From the fitted model from \code{\link{RNAmf}},
#' the posterior mean and variance are calculated based on the closed-form expression derived by a recursive fashion.
#' The formulas depend on its kernel choices.
#' For further details, see Heo and Sung (2025, <\doi{https://doi.org/10.1080/00401706.2024.2376173}>).
#'
#' @param object An object of class \code{RNAmf} fitted by \code{\link{RNAmf}}.
#' @param x A vector or matrix of new input locations for prediction.
#' @param ... Additional arguments for compatibility with generic method \code{predict}.
#'
#' @return
#' \itemize{
#'   \item \code{mu}: A list of vectors containing the predictive posterior mean at each fidelity level.
#'   \item \code{sig2}: A list of vectors containing the predictive posterior variance at each fidelity level.
#'   \item \code{time}: A scalar indicating the computation time.
#' }
#'
#' @importFrom plgp distance
#' @importFrom stats predict
#' @rdname predict
#' @method predict RNAmf
#' @export
#' @examples
#' \donttest{
#' ### two levels example ###
#' library(lhs)
#'
#' ### Perdikaris function ###
#' f1 <- function(x) {
#'   sin(8 * pi * x)
#' }
#'
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
#' fit.RNAmf <- RNAmf(list(X1, X2), list(y1, y2), kernel = "sqex", constant=TRUE)
#'
#' ### predict ###
#' predy <- predict(fit.RNAmf, x)$mu[[2]]
#' predsig2 <- predict(fit.RNAmf, x)$sig2[[2]]
#'
#' ### RMSE ###
#' print(sqrt(mean((predy - f2(x))^2)))
#'
#' ### visualize the emulation performance ###
#' plot(x, predy,
#'   type = "l", lwd = 2, col = 3, # emulator and confidence interval
#'   ylim = c(-2, 1)
#' )
#' lines(x, predy + 1.96 * sqrt(predsig2 * length(y2) / (length(y2) - 2)), col = 3, lty = 2)
#' lines(x, predy - 1.96 * sqrt(predsig2 * length(y2) / (length(y2) - 2)), col = 3, lty = 2)
#'
#' curve(f2(x), add = TRUE, col = 1, lwd = 2, lty = 2) # high fidelity function
#'
#' points(X1, y1, pch = 1, col = "red") # low-fidelity design
#' points(X2, y2, pch = 4, col = "blue") # high-fidelity design
#'
#' ### three levels example ###
#' ### Branin function ###
#' branin <- function(xx, l){
#'   x1 <- xx[1]
#'   x2 <- xx[2]
#'   if(l == 1){
#'     10*sqrt((-1.275*(1.2*x1+0.4)^2/pi^2+5*(1.2*x1+0.4)/pi+(1.2*x2+0.4)-6)^2 +
#'     (10-5/(4*pi))*cos((1.2*x1+0.4))+ 10) + 2*(1.2*x1+1.9) - 3*(3*(1.2*x2+2.4)-1) - 1 - 3*x2 + 1
#'   }else if(l == 2){
#'     10*sqrt((-1.275*(x1+2)^2/pi^2+5*(x1+2)/pi+(x2+2)-6)^2 +
#'     (10-5/(4*pi))*cos((x1+2))+ 10) + 2*(x1-0.5) - 3*(3*x2-1) - 1
#'   }else if(l == 3){
#'     (-1.275*x1^2/pi^2+5*x1/pi+x2-6)^2 + (10-5/(4*pi))*cos(x1)+ 10
#'   }
#' }
#'
#' output.branin <- function(x, l){
#'   factor_range <- list("x1" = c(-5, 10), "x2" = c(0, 15))
#'
#'   for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
#'   branin(x[1:2], l)
#' }
#'
#' ### training data ###
#' n1 <- 20; n2 <- 15; n3 <- 10
#'
#' ### fix seed to reproduce the result ###
#' set.seed(1)
#'
#' ### generate initial nested design ###
#' X <- NestedX(c(n1, n2, n3), 2)
#' X1 <- X[[1]]
#' X2 <- X[[2]]
#' X3 <- X[[3]]
#'
#' ### n1, n2 and n3 might be changed from NestedX ###
#' ### assign n1, n2 and n3 again ###
#' n1 <- nrow(X1)
#' n2 <- nrow(X2)
#' n3 <- nrow(X3)
#'
#' y1 <- apply(X1,1,output.branin, l=1)
#' y2 <- apply(X2,1,output.branin, l=2)
#' y3 <- apply(X3,1,output.branin, l=3)
#'
#' ### n=10000 grid test data ###
#' x <- as.matrix(expand.grid(seq(0, 1, length.out = 100),seq(0, 1, length.out = 100)))
#'
#' ### fit an RNAmf ###
#' fit.RNAmf <- RNAmf(list(X1, X2, X3), list(y1, y2, y3), kernel = "sqex", constant=TRUE)
#'
#' ### predict ###
#' pred.RNAmf <- predict(fit.RNAmf, x)
#' predy <- pred.RNAmf$mu[[3]]
#' predsig2 <- pred.RNAmf$sig2[[3]]
#'
#' ### RMSE ###
#' print(sqrt(mean((predy - apply(x,1,output.branin, l=3))^2)))
#'
#' ### visualize the emulation performance ###
#' x1 <- x2 <- seq(0, 1, length.out = 100)
#' oldpar <- par(mfrow=c(1,2))
#' image(x1, x2, matrix(apply(x,1,output.branin, l=3), ncol=100),
#' zlim=c(0,310), main="Branin function")
#' image(x1, x2, matrix(predy, ncol=100),
#' zlim=c(0,310), main="RNAmf prediction")
#' par(oldpar)
#'
#' ### predictive variance ###
#' print(predsig2)}

predict.RNAmf <- function(object, x, ...) {
  t1 <- proc.time()
  if (!inherits(object, "RNAmf")) stop("The object is not of class \"RNAmf\"")

  L <- object$level
  kernel <- object$kernel
  fits <- object$fits
  d <- ncol(fits[[1]]$X)
  x <- matrix(x, ncol = d)

  mu_list  <- vector("list", L)
  sig2_list <- vector("list", L)

  # Level 1
  if (kernel == "sqex") {
    pred1 <- pred.GP(fits[[1]], x)
  } else {
    nu <- fits[[1]]$nu
    pred1 <- pred.matGP(fits[[1]], x)
  }
  mu_list[[1]]  <- pred1$mu
  sig2_list[[1]] <- pred1$sig2

  # Recursive update for k = 2..L
  for (k in 2:L) {
    fit_prev  <- fits[[k-1]]
    fit_curr  <- fits[[k]]

    X_curr <- matrix(fit_curr$X[, -(d + 1)], ncol = d)
    w_curr <- fit_curr$X[, d + 1]

    x.mu <- mu_list[[k-1]]
    sig2 <- sig2_list[[k-1]]

    y_curr <- fit_curr$y
    n <- length(y_curr)
    theta <- fit_curr$theta
    tau2hat <- fit_curr$tau2hat
    mu_curr <- fit_curr$mu.hat
    Ci <- fit_curr$Ki
    a <- Ci %*% (y_curr - mu_curr)
    a_term <- drop(a %o% a) - (Ci * tau2hat)

    # Predict at level k
    if (object$kernel == "sqex") {
      # mean
      K_x_X <- exp(-distance(t(t(x) / sqrt(theta[-(d + 1)])), t(t(X_curr) / sqrt(theta[-(d + 1)]))))
      predy <- mu_curr + (K_x_X * 1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
                            exp(-(drop(outer(x.mu, w_curr, "-")))^2 / (theta[d + 1] + 2 * sig2))) %*% a
      # var
      predsig2 <- numeric(nrow(x))
      common_term <- exp(-(outer(w_curr, w_curr, "-"))^2 / (2 * theta[d + 1]))
      for (i in 1:nrow(x)) {
        mat <- (K_x_X[i, ] %o% K_x_X[i, ]) * 1 / sqrt(1 + 4 * sig2[i] / theta[d + 1]) *
          exp(-(outer(w_curr, w_curr, "+") / 2 - x.mu[i])^2 / (theta[d + 1] / 2 + 2 * sig2[i])) *
          common_term
        predsig2[i] <- pmax(0, tau2hat - (predy[i] - mu_curr)^2 + sum(a_term * mat))
      }
    } else {
      # mean
      K_x_X <- cor.sep(x, X_curr, theta[-(d + 1)], nu=nu)
      xi_matrix <- matrix(xifun(w = rep(w_curr, each = nrow(x)),
                                m = rep(x.mu, n), s = rep(sig2, n), theta = theta[d+1], nu = nu), nrow(x), n)
      predy <- mu_curr + rowSums((K_x_X * xi_matrix) %*% a)
      # var
      predsig2 <- numeric(nrow(x))
      for(i in 1:nrow(x)) {
        zeta_mat <- matrix(zetafun(w1 = rep(w_curr, times = n), w2 = rep(w_curr, each  = n),
                                   m = x.mu[i], s = sig2[i], nu = nu, theta = theta[d+1]), n, n)
        mat <- (K_x_X[i,] %o% K_x_X[i,]) * zeta_mat
        predsig2[i] <- pmax(0, tau2hat - (predy[i] - mu_curr)^2 + sum(a_term * mat))
      }
    }
    mu_list[[k]]  <- predy
    sig2_list[[k]] <- predsig2
  }

  return(list(mu = mu_list, sig2 = sig2_list, time = (proc.time() - t1)[3]))
}
