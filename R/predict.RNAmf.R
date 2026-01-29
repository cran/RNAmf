#' @title Closed-form prediction for RNAmf model
#'
#' @description The function computes the closed-form posterior mean and variance for the RNAmf model
#' both at the fidelity levels used in model fitting using the chosen kernel.
#'
#' @param fits A fitted GP object from \code{RNAmf}.
#' @param x A vector or matrix of new input locations to predict.
#' @param kernel A character specifying the kernel type to be used. Choices are \code{"sqex"}(squared exponential kernel), \code{"matern1.5"}(Matern kernel with \eqn{\nu=1.5}), or \code{"matern2.5"}(Matern kernel with \eqn{\nu=2.5}). Default is \code{"sqex"}.
#' @param XX A list containing a pseudo-complete inputs \code{X_star}(\eqn{\left\{\mathcal{X}^*_l\right\}_{l=1}^{L}}), an original inputs \code{X_list}(\eqn{\left\{\mathcal{X}_l\right\}_{l=1}^{L}}), and a pseudo inputs \code{X_tilde}(\eqn{\left\{\widetilde{\mathcal{X}}_l\right\}_{l=1}^{L}}) for non-nested design.
#' @param pseudo_yy A list containing a pseudo-complete outputs \code{y_star}(\eqn{\left\{\mathbf{y}^*_l\right\}_{l=1}^{L}}), an original outputs \code{y_list}(\eqn{\left\{\mathbf{y}_l\right\}_{l=1}^{L}}), and a pseudo outputs \code{y_tilde}(\eqn{\left\{\widetilde{\mathbf{y}}_l\right\}_{l=1}^{L}}) imputed by \code{\link{imputer_RNA}}.
#'
#' @return A list of predictive posterior mean and variance for each level containing:
#' \itemize{
#'   \item \code{mu}: A list of predictive posterior mean at each fidelity level.
#'   \item \code{sig2}: A list of predictive posterior variance at each fidelity level.
#' }
#'
closed_form_RNA <- function(fits, x, kernel, XX=NULL, pseudo_yy=NULL) {

  if (!is.list(fits) || length(fits) < 1) stop("closed_form_RNA: 'fits' must be a valid RNA objects.")
  if (!is.numeric(x)) stop("closed_form_RNA: 'x' must be a numeric matrix or vector.")
  if (!is.null(pseudo_yy)) {
    if (is.null(XX)) stop("closed_form_RNA: 'XX' must be provided when 'pseudo_yy' is not NULL.")
    if (!is.list(XX) || is.null(XX$X_star)) stop("closed_form_RNA: 'XX' is malformed.")
    if (!is.list(pseudo_yy) || is.null(pseudo_yy$y_star)) stop("closed_form_RNA: 'pseudo_yy' is malformed.")
  }
  L <- length(fits)
  d <- ncol(fits[[1]]$X)

  mu_list   <- vector("list", L)
  sig2_list <- vector("list", L)

  if(is.null(pseudo_yy)){
    if (kernel == "sqex") {
      pred1 <- pred.GP(fits[[1]], x)
    } else {
      nu <- fits[[1]]$nu
      pred1 <- pred.matGP(fits[[1]], x)
    }
  } else {
    X1 <- XX$X_star[[1]]
    y1 <- pseudo_yy$y_star[[1]]
    n1 <- nrow(X1)

    if(kernel=="sqex"){
      K <- covar.sep(X1, d = fits[[1]]$theta, g = fits[[1]]$g)
    } else {
      nu <- fits[[1]]$nu
      K <- cor.sep(X1, theta = fits[[1]]$theta, nu = nu) + diag(fits[[1]]$g, n1)
    }

    chol_K <- chol(K)
    Ki <- chol2inv(chol_K)
    fits[[1]]$Ki <- Ki
    one.vec <- matrix(1, ncol = 1, nrow = n1)
    fits[[1]]$mu.hat <- drop(crossprod(one.vec, crossprod(Ki, y1)) / crossprod(one.vec, crossprod(Ki, one.vec)))
    fits[[1]]$tau2hat <- drop(crossprod(y1 - fits[[1]]$mu.hat, crossprod(Ki, y1 - fits[[1]]$mu.hat)) / n1)

    if(kernel == "sqex") {
      pred1 <- pred.GP(fits[[1]], x)
    } else {
      pred1 <- pred.matGP(fits[[1]], x)
    }
  }

  mu_list[[1]]   <- pred1$mu
  sig2_list[[1]] <- pred1$sig2

  for (k in 2:L) {
    fit_curr <- fits[[k]]
    X_curr  <- matrix(fit_curr$X[, -(d + 1)], ncol = d)
    theta   <- fit_curr$theta
    if(is.null(pseudo_yy)){
      w_curr  <- fit_curr$X[, d + 1]
      y_curr  <- fit_curr$y
      Ci      <- fit_curr$Ki
      mu_curr <- fit_curr$mu.hat
      tau2hat <- fit_curr$tau2hat
      n <- length(y_curr)
    } else {
      w_curr <- c(pseudo_yy$y_star[[k-1]][checkindices(XX$X_star[[k-1]], XX$X_star[[k]]), , drop = FALSE])
      y_curr <- pseudo_yy$y_star[[k]]
      g <- fit_curr$g
      n <- length(y_curr)

      X_aug <- cbind(X_curr, w_curr)
      if(kernel=="sqex"){
        K <- covar.sep(X_aug, d = theta, g = 0)
      } else {
        K <- cor.sep(X_aug, theta = theta, nu = nu)
      }
      chol_Ci <- chol(K + diag(g, n))
      Ci <- chol2inv(chol_Ci)
      one.vec <- matrix(1, ncol = 1, nrow = n)
      mu_curr <- drop(crossprod(one.vec, crossprod(Ci, y_curr)) / crossprod(one.vec, crossprod(Ci, one.vec)))
      tau2hat <- drop(crossprod(y_curr - mu_curr, crossprod(Ci, y_curr - mu_curr)) / n)
    }
    a <- crossprod(Ci, y_curr - mu_curr)

    x.mu <- mu_list[[k - 1]]
    sig2 <- sig2_list[[k - 1]]

    if (kernel == "sqex") {
      theta_d1 <- theta[d + 1]
      sqrt_theta <- sqrt(theta[-(d + 1)])
      x_s <- sweep(x, 2, sqrt_theta, "/")
      X_curr_s <- sweep(X_curr, 2, sqrt_theta, "/")
      K_x_X <- exp(-distance(x_s, X_curr_s))

      diff_sq <- drop(outer(x.mu, w_curr, "-")^2)
      denom_vec <- theta_d1 + 2 * sig2
      term_exp <- exp(-(diff_sq / denom_vec))
      term_scale <- 1 / sqrt(1 + 2 * sig2 / theta_d1)
      weights <- K_x_X * term_exp * term_scale

      predy <- mu_curr + crossprod(t(weights), a)

      predsig2 <- numeric(nrow(x))
      common_term_w <- exp(-(outer(w_curr, w_curr, "-"))^2 / (2 * theta_d1))

      for (i in 1:nrow(x)) {
        s2_i <- sig2[i]
        k_x <- K_x_X[i, ]
        mat_base <- tcrossprod(k_x)
        w_sum_half <- outer(w_curr, w_curr, "+") / 2
        exp_mat <- exp(-(w_sum_half - x.mu[i])^2 / (theta_d1 / 2 + 2 * s2_i))
        scale_val <- 1 / sqrt(1 + 4 * s2_i / theta_d1)
        mat <- mat_base * scale_val * exp_mat * common_term_w
        quad <- drop(crossprod(a, crossprod(mat, a)))
        tr_term <- sum(Ci * mat)

        predsig2[i] <- pmax(0, tau2hat - (predy[i] - mu_curr)^2 + quad - tau2hat * tr_term)
      }

    } else {
      K_x_X <- cor.sep(x, X_curr, theta[-(d + 1)], nu=nu)
      xi_vec <- xifun(w = rep(w_curr, each = nrow(x)),
                      m = rep(x.mu, n),
                      s = rep(sig2, n),
                      theta = theta[d+1], nu = nu)
      xi_matrix <- matrix(xi_vec, nrow(x), n)

      predy <- mu_curr + rowSums(crossprod(t(K_x_X * xi_matrix), a))

      predsig2 <- numeric(nrow(x))
      for(i in 1:nrow(x)) {
        zeta_mat <- matrix(zetafun(w1 = rep(w_curr, times = n),
                                   w2 = rep(w_curr, each  = n),
                                   m = x.mu[i], s = sig2[i],
                                   nu = nu, theta = theta[d+1]), n, n)
        mat <- tcrossprod(K_x_X[i, ]) * zeta_mat
        quad <- drop(crossprod(a, crossprod(mat, a)))
        tr_term <- sum(Ci * mat)

        predsig2[i] <- pmax(0, tau2hat - (predy[i] - mu_curr)^2 + quad - tau2hat * tr_term)
      }
    }
    mu_list[[k]] <- predy
    sig2_list[[k]] <- predsig2
  }
  return(list(mu = mu_list, sig2 = sig2_list))
}


#' @title prediction of the RNAmf emulator with multiple fidelity levels.
#'
#' @description The function computes the posterior mean and variance of RNA models with multiple fidelity levels
#' by fitted model from \code{\link{RNAmf}}.
#'
#' @seealso \code{\link{RNAmf}} for model fitting.
#'
#' @details The \code{predict.RNAmf} function internally calls \code{\link{closed_form_RNA}}
#' to recursively compute the closed-form posterior mean and variance at each level.
#' 
#' From the fitted model from \code{\link{RNAmf}},
#' the posterior mean and variance are calculated based on the closed-form expression derived by a recursive fashion.
#' The formulas depend on its kernel choices.
#' For further details, see Heo and Sung (2025, <\doi{https://doi.org/10.1080/00401706.2024.2376173}>).
#'
#' @param object An object of class \code{RNAmf} fitted by \code{\link{RNAmf}}.
#' @param x A vector or matrix of new input locations for prediction.
#' @param nimpute Number of imputations for non-nested designs. Default is 50.
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

predict.RNAmf <- function(object, x=NULL, nimpute=50, ...) {
  if (is.null(x)) return(fitted(object))
  t1 <- proc.time()
  if (!inherits(object, "RNAmf")) stop("The object is not of class \"RNAmf\" \n")
  if (!is.numeric(x)) stop("'x' must be numeric.")
  if (!is.numeric(nimpute) || length(nimpute) != 1 || nimpute < 1) stop("'nimpute' must be a positive integer.")

  d <- ncol(object$fits[[1]]$X)
  if (is.vector(x)) {
    if (length(x) %% d != 0) {
      stop(paste0("Input 'x' length must be a multiple of the model dimension (", d, ")."))
    }
    x <- matrix(x, ncol = d, byrow = FALSE)
  } else {
    x <- as.matrix(x)
    if (ncol(x) != d) {
      stop(paste0("Input 'x' must have ", d, " columns (matches model dimension)."))
    }
  }
  L <- object$level
  kernel <- object$kernel
  fits <- object$fits

  if(object$nested){
    pred_result <- closed_form_RNA(fits, x, kernel)
    return(list(mu = pred_result$mu, sig2 = pred_result$sig2, time = (proc.time() - t1)[3]))
  }else{
    XX <- object$XX
    yy <- object$yy
    pred1 <- object$pred1
    n <- nrow(x)

    mean_mu   <- vector("list", L); for (l in 1:L) mean_mu[[l]]   <- numeric(n)
    M2_mu     <- vector("list", L); for (l in 1:L) M2_mu[[l]]     <- numeric(n)
    mean_sig2 <- vector("list", L); for (l in 1:L) mean_sig2[[l]] <- numeric(n)

    for (m in 1:nimpute) {
      # Generate imputed dataset for the m-th imputation
      yy <- imputer_RNA(XX, yy, kernel=kernel, pred1, fits)
      pred <- closed_form_RNA(fits, x, kernel, XX, yy)

      for (l in 1:L) {
        mu_new   <- pred$mu[[l]]
        sig2_new <- pred$sig2[[l]]

        delta        <- mu_new - mean_mu[[l]]
        mean_mu[[l]] <- mean_mu[[l]] + delta / m
        M2_mu[[l]]   <- M2_mu[[l]] + delta * (mu_new - mean_mu[[l]])

        # simple running mean for sig2
        mean_sig2[[l]] <- mean_sig2[[l]] + (sig2_new - mean_sig2[[l]]) / m
      }
    }

    mu_star   <- mean_mu
    sig2_star <- vector("list", L)
    for (l in 1:L) {
      var_mu_pop   <- if (nimpute > 0) M2_mu[[l]] / nimpute else numeric(n)
      sig2_star[[l]] <- mean_sig2[[l]] + var_mu_pop
    }
    names(mu_star)   <- paste0("mu_",   1:L)
    names(sig2_star) <- paste0("sig2_", 1:L)

    return(list(mu = mu_star, sig2 = sig2_star, time = (proc.time() - t1)[3]))
  }
}
