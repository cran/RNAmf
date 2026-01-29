#' @method summary RNAmf
#' @export
summary.RNAmf <- function(object,...){
  # Create a shallow copy to avoid modifying the original object in place
  ans <- object
  # Prepend the summary class (safer than overwriting all classes)
  class(ans) <- c("summary.RNAmf", class(object))
  ans
}

#' @export
print.summary.RNAmf <- function(x, ...){
  cat("----------------------------------------------------\n")
  cat(" Recursive Non-Additive Multi-Fidelity (RNAmf) Model\n")
  cat("----------------------------------------------------\n")
  cat(sprintf("  * Levels:   %d\n", x$level))
  cat(sprintf("  * Kernel:   %s\n", x$kernel))
  cat(sprintf("  * Constant: %s\n", x$constant))
  cat(sprintf("  * nested:   %s\n", x$nested))

  # Calculate sample sizes for all levels
  sizes <- vapply(x$fits, function(f) NROW(f$X), integer(1))
  size_str <- paste(sprintf("Level %d: %d", seq_along(sizes), sizes), collapse = ", ")
  cat(sprintf("  * Samples:  %s\n", size_str))

  cat(sprintf("  * Time:     %.3f sec\n", x$time))

  cat("----------------------------------------------------\n")
  invisible(x)
}

#' @method print RNAmf
#' @export
print.RNAmf <- function(x, ...) {
  print(summary(x, ...)) # Delegate printing to the summary method
  invisible(x)           # CRITICAL: Return the original object 'x'
}

#' @method coef RNAmf
#' @export
coef.RNAmf <- function(object, ...) {
  # Return a structured list of parameters for each level
  lapply(seq_len(object$level), function(l) {
    fit <- object$fits[[l]]
    list(
      theta = fit$theta,
      mu_hat = fit$mu.hat,
      tau2_hat = fit$tau2hat
    )
  })
}

#' @method fitted RNAmf
#' @importFrom stats predict
#' @export
fitted.RNAmf <- function(object, ...) {
  # 1. Determine the input dimension 'd' from the first level (which is never augmented)
  d <- ncol(object$fits[[1]]$X)

  # 2. Iterate through each level 'l'
  lapply(seq_along(object$fits), function(l) {

    # Extract ONLY the original input columns (1:d)
    # We drop the augmented columns (previous fidelity outputs) used during training
    X_train_raw <- object$fits[[l]]$X[, 1:d, drop = FALSE]

    # Predict using the model on these training points
    # We extract 'mu[[l]]' because we want the prediction for this specific fidelity level
    predict(object, X_train_raw, ...)$mu[[l]]
  })
}

#' @method residuals RNAmf
#' @importFrom stats fitted
#' @export
residuals.RNAmf <- function(object, ...) {
  # 1. Get Fitted values (returns list of length L)
  y_hat <- fitted(object)

  # 2. Get Observed values
  # In RNAmf, the response 'y' is stored inside each fit object in the 'fits' list
  y_obs <- lapply(object$fits, function(f) f$y)

  # 3. Calculate Residuals (y_obs - y_hat)
  # Map safely subtracts corresponding elements of the two lists
  residuals_list <- Map("-", y_obs, y_hat)

  # Optional: Assign names for clarity
  names(residuals_list) <- paste0("Level", seq_along(residuals_list))

  return(residuals_list)
}

#' @method plot RNAmf
#' @importFrom graphics par lines points legend
#' @export
plot.RNAmf <- function(x, ...) {
  # 1. Determine dimensions and number of levels
  L <- x$level
  d <- ncol(x$fits[[1]]$X)

  # 2. Stop if input is not 1-dimensional
  if (d > 1) {
    stop("The plot method for 'RNAmf' currently only supports 1-dimensional inputs.")
  }

  # 3. Setup plotting layout (1 row, L columns)
  oldpar <- par(mfrow = c(1, L))
  on.exit(par(oldpar))

  # 4. Generate 1D Functional Plot
  # Create a grid of 101 points from 0 to 1
  xg <- matrix(seq(0, 1, length.out = 101), ncol = 1)

  # Get predictions for ALL levels at once
  preds <- predict(x, xg)

  for (l in 1:L) {
    # Extract training data for this level
    # Note: We take only the 1st column (original input) even if X is augmented
    X_train <- x$fits[[l]]$X[, 1]
    y_train <- x$fits[[l]]$y

    # Extract predictions for this level
    y_mu  <- preds$mu[[l]]
    y_sig <- sqrt(preds$sig2[[l]])

    # Calculate 95% Confidence Interval
    lower <- y_mu - 1.96 * y_sig
    upper <- y_mu + 1.96 * y_sig

    # Determine plot limits to ensure everything fits
    ylim <- range(c(y_train, lower, upper), na.rm = TRUE)

    # Plot Mean Curve
    plot(xg, y_mu, type = "l", lwd = 2, col = "blue", ylim = ylim,
         main = paste("Level", l, "Fit"),
         xlab = "x", ylab = "y", ...)

    # Add Confidence Interval (Dashed lines)
    lines(xg, lower, col = "blue", lty = 2)
    lines(xg, upper, col = "blue", lty = 2)

    # Add Training Points
    points(X_train, y_train, pch = 16, cex = 1, col = "black")

    # Add Legend (only on the first plot to save space)
    if (l == 1) {
      legend("topleft", legend = c("Mean", "95% CI", "Data"),
             col = c("blue", "blue", "black"),
             lty = c(1, 2, NA), pch = c(NA, NA, 16), bty = "n")
    }
  }

  invisible(x)
}

#' @title Fitting the Recursive Non-Additive model with multiple fidelity levels
#'
#' @description The function fits RNA models with designs of multiple fidelity levels.
#' The estimation method is based on MLE.
#' Available kernel choices include the squared exponential kernel, and the Matern kernel with smoothness parameter 1.5 and 2.5.
#' The function returns the fitted model at each level \eqn{1, \ldots, l}, the number of fidelity levels \eqn{l}, the kernel choice, whether constant mean or not, and the computation time.
#'
#' @seealso \code{\link{predict.RNAmf}} for prediction.
#'
#' @details Consider the model
#' \eqn{\begin{cases}
#' & f_1(\bm{x}) = W_1(\bm{x}),\\
#' & f_l(\bm{x}) = W_l(\bm{x}, f_{l-1}(\bm{x})) \quad\text{for}\quad l \geq 2,
#' \end{cases}}
#' where \eqn{f_l} is the simulation code at fidelity level \eqn{l}, and
#' \eqn{W_l(\bm{x}) \sim GP(\alpha_l, \tau_l^2 K_l(\bm{x}, \bm{x}'))} is GP model.
#' Hyperparameters \eqn{(\alpha_l, \tau_l^2, \bm{\theta_l})} are estimated by
#' maximizing the log-likelihood via an optimization algorithm "L-BFGS-B".
#' For \code{constant=FALSE}, \eqn{\alpha_l=0}.
#'
#' Covariance kernel is defined as:
#' \eqn{K_l(\bm{x}, \bm{x}')=\prod^d_{j=1}\phi(x_j,x'_j;\theta_{lj})} with
#' \eqn{\phi(x, x';\theta) = \exp \left( -\frac{ \left( x - x' \right)^2}{\theta}  \right)}
#' for squared exponential kernel; \code{kernel="sqex"},
#' \eqn{\phi(x,x';\theta) =\left( 1+\frac{\sqrt{3}|x- x'|}{\theta} \right) \exp \left( -\frac{\sqrt{3}|x- x'|}{\theta} \right)}
#' for Matern kernel with the smoothness parameter of 1.5; \code{kernel="matern1.5"} and
#' \eqn{\phi(x, x';\theta) = \left( 1+\frac{\sqrt{5}|x-x'|}{\theta} +\frac{5(x-x')^2}{3\theta^2} \right) \exp \left( -\frac{\sqrt{5}|x-x'|}{\theta} \right)}
#' for Matern kernel with the smoothness parameter of 2.5; \code{kernel="matern2.5"}.
#'
#' For further details, see Heo and Sung (2025, <\doi{https://doi.org/10.1080/00401706.2024.2376173}>).
#'
#' @param X_list A list of the matrices of input locations for all fidelity levels.
#' @param y_list A list of the vectors or matrices of response values for all fidelity levels.
#' @param kernel A character specifying the kernel type to be used. Choices are \code{"sqex"}(squared exponential), \code{"matern1.5"}, or \code{"matern2.5"}. Default is \code{"sqex"}.
#' @param constant A logical indicating for constant mean of GP (\code{constant=TRUE}) or zero mean (\code{constant=FALSE}). Default is \code{TRUE}.
#' @param init Optional vector of initial parameter values for optimization. Default is \code{NULL}.
#' @param n.iter Number of iterations for the stochastic EM algorithm for non-nested designs. Default is 50.
#' @param burn.ratio Fraction of iterations to discard as burn-in. Default is 0.75.
#' @param trace A logical indicating to print progress of iterations if \code{TRUE}, or not if \code{FALSE}. Default is \code{TRUE}.
#' @param ... Additional arguments for compatibility with \code{optim}.
#'
#' @return A list of class \code{RNAmf} with:
#' \itemize{
#'   \item \code{fits}: A list of fitted Gaussian process models \eqn{f_l(x)} at each level \eqn{1, \ldots, l}. Each element contains:
#'   \eqn{\begin{cases} & f_1 \text{ for } (X_1, y_1),\\ & f_l \text{ for } ((X_l, f_{l-1}(X_l)), y_l), \end{cases}}.
#'   \item \code{level}: The number of fidelity levels \eqn{l}.
#'   \item \code{kernel}: A copy of \code{kernel}.
#'   \item \code{constant}: A copy of \code{constant}.
#'   \item \code{nested}: A logical indicating whether the design is nested.
#'   \item \code{time}: A scalar indicating the computation time.
#' }
#' @usage RNAmf(X_list, y_list, kernel = "sqex", constant = TRUE,
#' init = NULL, n.iter = 50, burn.ratio = 0.75, trace = TRUE, ...)
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
#' ### fit an RNAmf ###
#' fit.RNAmf <- RNAmf(list(X1, X2), list(y1, y2), kernel = "sqex", constant=TRUE)
#'
#'
#' ### three levels example ###
#'
#' ### Branin function ###
#' branin <- function(xx, l){
#'   x1 <- xx[1]
#'   x2 <- xx[2]
#'   if(l == 1){
#'     10*sqrt((-1.275*(1.2*x1+0.4)^2/pi^2+5*(1.2*x1+0.4)/pi+(1.2*x2+0.4)-6)^2 +
#'               (10-5/(4*pi))*cos((1.2*x1+0.4))+ 10) +
#'               2*(1.2*x1+1.9) - 3*(3*(1.2*x2+2.4)-1) - 1 - 3*x2 + 1
#'   }else if(l == 2){
#'     10*sqrt((-1.275*(x1+2)^2/pi^2+5*(x1+2)/pi+(x2+2)-6)^2 +
#'               (10-5/(4*pi))*cos((x1+2))+ 10) + 2*(x1-0.5) - 3*(3*x2-1) - 1
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
#' ### fit an RNAmf ###
#' fit.RNAmf <- RNAmf(list(X1, X2, X3), list(y1, y2, y3), kernel = "sqex", constant=TRUE)
#'}
#'
#'

RNAmf <- function(X_list, y_list, kernel = "sqex", constant = TRUE, init = NULL, n.iter = 50, burn.ratio = 0.75, trace = TRUE, ...) {
  t0 <- proc.time()[3]
  if (!is.list(X_list)) stop("RNAmf: 'X_list' must be a list of input matrices/vectors.")
  if (!is.list(y_list)) stop("RNAmf: 'y_list' must be a list of output matrices/vectors.")
  if(kernel != "sqex") nu <- ifelse(kernel == "matern1.5", 1.5, 2.5)

  L <- length(X_list)
  d <- NCOL(X_list[[1]])
  if (length(y_list) != L) stop("RNAmf: Length of 'X_list' and 'y_list' must match.")
  if (L < 2) stop("RNAmf: Requires at least two fidelity levels (L >= 2).")
  if (!is.character(kernel) || length(kernel) != 1 || !kernel %in% c("sqex", "matern1.5", "matern2.5")) {
    stop(paste("RNAmf: 'kernel' must be one of:", paste(c("sqex", "matern1.5", "matern2.5"), collapse = ", ")))
  }
  if (!is.logical(constant)) stop("RNAmf: 'constant' must be TRUE or FALSE.")
  if (!is.numeric(n.iter) || length(n.iter) != 1 || n.iter < 1) stop("RNAmf: 'n.iter' must be a positive integer.")
  if (!is.numeric(burn.ratio) || length(burn.ratio) != 1 || burn.ratio <= 0 || burn.ratio >= 1) stop("RNAmf: 'burn.ratio' must be a numeric value strictly between 0 and 1.")
  if (!is.logical(trace)) stop("RNAmf: 'trace' must be TRUE or FALSE.")

  fits <- vector("list", L)
  fit_fun <- switch(kernel,
                    sqex = function(X, y, init=init, ...) GP(X, y, constant = constant, init=init, ...),
                    matern1.5 = function(X, y, init=init, ...) matGP(X, y, nu = nu, constant = constant, init=init, ...),
                    matern2.5 = function(X, y, init=init, ...) matGP(X, y, nu = nu, constant = constant, init=init, ...) )
  pred_fun <- switch(kernel,
                     sqex = function(fit, Xnew, cov.out=FALSE) pred.GP(fit, Xnew, cov.out=cov.out),
                     matern1.5 = function(fit, Xnew, cov.out=FALSE) pred.matGP(fit, Xnew, cov.out=cov.out),
                     matern2.5 = function(fit, Xnew, cov.out=FALSE) pred.matGP(fit, Xnew, cov.out=cov.out) )
  # Check nested
  nested <- all(mapply(checknested, X_list[-L], X_list[-1]))

  fit_core <- function(X_star, y_star, init_list) {
    fits <- vector("list", L)
    # level 1
    fits[[1]] <- fit_fun(X_star[[1]], y_star[[1]], init_list[[1]], ...)
    # levels 2..L (augmenting with prior predicted means)
    for (lvl in 2:L) {
      X_curr <- X_star[[lvl]]
      # preallocate augmented design to avoid repeated cbind
      X_aug <- X_curr
      # fill augmentation columns sequentially (no reallocation)
      for (j in 1:(lvl - 1)) {
        pj <- pred_fun(fits[[j]], X_aug[, seq_len(ncol(X_curr) + (j - 1)), drop = FALSE])
        X_aug <- cbind(X_curr, pj$mu)
      }
      fits[[lvl]] <- fit_fun(X_aug, y_star[[lvl]], init_list[[lvl]])
    }
    fits
  }

  if(nested){
    fits <- fit_core(X_star = X_list, y_star = y_list, init_list = init)
    structure(list(fits = fits, level = L, kernel = kernel, constant = constant, nested = nested, time = proc.time()[3] - t0), class = "RNAmf" )
  } else { # Non-nested
    XX <- makenested(c(X_list[1], X_list[-1]))
    y_all <- lapply(y_list, matrix)
    yy <- list(y_star  = vector("list", L), y_list  = y_all, y_tilde = vector("list", L))

    ### Draw initial y tilde ###
    ### sample initial y1_tilde from the imputer ###
    fit1  <- fit_fun(XX$X_list[[1]], y_all[[1]], init = NULL)
    pred1 <- pred_fun(fit1, XX$X_tilde[[1]], cov.out = TRUE)
    yy$y_tilde[[1]] <- pred1$mu
    yy$y_star[[1]]  <- rbind(yy$y_list[[1]], yy$y_tilde[[1]])

    ### sample initial y_2_tilde using individual GPs ###
    if (L > 2) {
      for (l in 2:(L - 1)) {
        fit_next  <- fit_fun(XX$X_list[[l]], yy$y_list[[l]], init = NULL)
        pred_next <- pred_fun(fit_next, XX$X_tilde[[l]])
        yy$y_tilde[[l]] <- pred_next$mu
        yy$y_star[[l]]  <- rbind(yy$y_list[[l]], yy$y_tilde[[l]])
      }
    }
    yy$y_star[[L]] <- yy$y_list[[L]]

    n.burnin <- ceiling(n.iter * burn.ratio) # number of burn in
    n.param <- n.iter - n.burnin # number of collected estimates
    param_mat <- rep(list(matrix(NA, nrow=n.param, ncol=d+3)), L) # theta_x, theta_y, mu.hat, tau.hat
    param_mat[[1]] <- param_mat[[1]][,-(d+3)]

    ### initial estimates
    fits <- fit_core(X_star = XX$X_star, y_star = yy$y_star, init_list = rep(list(NULL), L))

    if (trace) {
      cat(paste(vapply(seq_len(L), function(l) {
        vals <- c(fits[[l]]$theta, fits[[l]]$mu.hat, fits[[l]]$tau2hat)
        sprintf("initial (l=%d): %s", l, paste(sprintf("%.4e", vals), collapse = " "))
      }, FUN.VALUE = character(1L)), collapse = "\n"), "\n")
    }

    for (j in 1:n.iter) { # Imputation and Maximization
      # Imputation step; impute y tilde using ESS
      yy <- imputer_RNA(XX, yy, kernel=kernel, pred1, fits)

      param.init <- lapply(fits, `[[`, "theta")
      # Maximization step; optimize parameters n.iter times
      fits <- fit_core(X_star = XX$X_star, y_star = yy$y_star, init_list = param.init)

      if(j > n.burnin){
        for (l in 1:L) {
          param_mat[[l]][j-n.burnin,] <- c(fits[[l]]$theta, fits[[l]]$mu.hat, fits[[l]]$tau2hat)
        }
      }
      if (trace) {
        cat(paste(vapply(seq_len(L), function(l) {
          vals <- c(fits[[l]]$theta, fits[[l]]$mu.hat, fits[[l]]$tau2hat)
          sprintf("iter %2d/%d (l=%d): %s", j, n.iter, l, paste(sprintf("%.4e", vals), collapse = " "))
        }, FUN.VALUE = character(1L)), collapse = "\n"), "\n")
      }
    } # end of j for loop

    # average with 75% burn-in
    colnames(param_mat[[1]]) <- c(paste0("theta_", seq_len(d)), "mu_hat", "tau2_hat")
    param_mat[-1] <- lapply(param_mat[-1], function(m) { colnames(m) <- c(paste0("theta_", seq_len(d)), "theta_y", "mu_hat", "tau2_hat"); m })

    final_params <- lapply(param_mat,colMeans)

    fits[[1]]$theta <- final_params[[1]][1:d]
    fits[[1]]$mu.hat <- final_params[[1]][d+1]
    fits[[1]]$tau2hat <- final_params[[1]][d+2]
    for (l in 2:L) {
      fits[[l]]$theta <- final_params[[l]][1:(d+1)]
      fits[[l]]$mu.hat <- final_params[[l]][d+2]
      fits[[l]]$tau2hat <- final_params[[l]][d+3]
    }

    structure(list(fits = fits, level = L, kernel = kernel, constant = constant, XX = XX, yy = yy, nested = nested, pred1 = pred1, time = proc.time()[3] - t0), class = "RNAmf" )
  }
}
