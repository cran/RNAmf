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
#' @param ... Additional arguments for compatibility with \code{optim}.
#'
#' @return A list of class \code{RNAmf} with:
#' \itemize{
#'   \item \code{fits}: A list of fitted Gaussian process models \eqn{f_l(x)} at each level \eqn{1, \ldots, l}. Each element contains:
#'   \eqn{\begin{cases} & f_1 \text{ for } (X_1, y_1),\\ & f_l \text{ for } ((X_l, f_{l-1}(X_l)), y_l), \end{cases}}.
#'   \item \code{level}: The number of fidelity levels \eqn{l}.
#'   \item \code{kernel}: A copy of \code{kernel}.
#'   \item \code{constant}: A copy of \code{constant}.
#'   \item \code{time}: A scalar indicating the computation time.
#' }
#' @usage RNAmf(X_list, y_list, kernel = "sqex", constant = TRUE, ...)
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
#' ### fit an RNAmf ###
#' fit.RNAmf <- RNAmf(list(X1, X2, X3), list(y1, y2, y3), kernel = "sqex", constant=TRUE)
#'}
#'
#'

RNAmf <- function(X_list, y_list, kernel = "sqex", constant = TRUE, ...) {
  t0 <- proc.time()[3]
  L <- length(X_list)
  if (length(y_list) != L) {
    stop("Length of X_list and y_list must match.")
  }

  # Check nested
  for (lvl in 2:L) {
    if (!checknested(X_list[[lvl - 1]], X_list[[lvl]])) {
      stop(sprintf("X%s is not nested within X%s", lvl, lvl-1))
    }
  }
  if (startsWith(kernel, "matern")) nu <- as.numeric(sub("matern", "", kernel))

  fits <- vector("list", L)
  fit_fun <- switch(kernel,
                    sqex = function(X, y) GP(X, y, constant = constant),
                    matern1.5 = function(X, y) matGP(X, y, nu = nu, constant = constant),
                    matern2.5 = function(X, y) matGP(X, y, nu = nu, constant = constant) )
  pred_fun <- switch(kernel,
                     sqex = function(fit, Xnew) pred.GP(fit, Xnew)$mu,
                     matern1.5 = function(fit, Xnew) pred.matGP(fit, Xnew)$mu,
                     matern2.5 = function(fit, Xnew) pred.matGP(fit, Xnew)$mu )
  # Fit Level 1
  fits[[1]] <- fit_fun(X_list[[1]], y_list[[1]])
  # Fit level 2-L recursively
  for (lvl in 2:L) {
    X_curr <- X_list[[lvl]]
    y_curr <- y_list[[lvl]]
    X_aug <- X_curr
    for (j in 1:(lvl - 1)) {
      mu_pred <- pred_fun(fits[[j]], X_aug)
      X_aug <- cbind(X_curr, mu_pred)
    }
    fits[[lvl]] <- fit_fun(X_aug, y_curr)
  }

  structure(list(fits = fits, level = L, kernel = kernel, constant = constant, time = proc.time()[3] - t0), class = "RNAmf" )
}
