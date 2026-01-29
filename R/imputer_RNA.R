#' @title Imputation step in stochastic EM for the non-nested RNA Model
#'
#' @description The function performs the imputation step of the stochastic EM algorithm for the RNA model when the design is not nested.
#' The function generates pseudo outputs \eqn{\widetilde{\mathbf{y}}_l} at pseudo inputs \eqn{\widetilde{\mathcal{X}}_l}.
#'
#' @details The \code{imputer_RNA} function then imputes the corresponding pseudo outputs
#' \eqn{\widetilde{\mathbf{y}}_l = f_l(\widetilde{\mathcal{X}}_l)}
#' by drawing samples from the conditional normal distribution,
#' given fixed parameter estimates and previous-level outputs \eqn{Y_{l}^{*(m-1)}} for each \eqn{l},
#' at the \eqn{m}-th iteration of the EM algorithm.
#'
#' @param XX A list of design sets for all fidelity levels, containing \code{X_star}, \code{X_list}, and \code{X_tilde}.
#' @param yy A list of current observed and pseudo-responses, containing \code{y_star}, \code{y_list}, and \code{y_tilde}.
#' @param kernel A character specifying the kernel type to be used. Choices are \code{"sqex"}(squared exponential), \code{"matern1.5"}, or \code{"matern2.5"}.
#' @param pred1 Predictive results for the lowest fidelity level \eqn{f_1}. It should include \code{cov} obtained by setting \code{cov.out=TRUE}.
#' @param fits A fitted GP object from \code{RNAmf}.
#'
#' @return An updated \code{yy} list containing:
#' \itemize{
#'   \item \code{y_star}: An updated pseudo-complete outputs \eqn{\mathbf{y}^*_l}.
#'   \item \code{y_list}: An original outputs \eqn{\mathbf{y}_l}.
#'   \item \code{y_tilde}: A newly imputed pseudo outputs \eqn{\widetilde{\mathbf{y}}_l}.
#' }
#' @usage imputer_RNA(XX, yy, kernel=kernel, pred1, fits)
#' @export
#'

imputer_RNA <- function(XX, yy, kernel=kernel, pred1, fits){

  if (!is.list(XX) || is.null(XX$X_star)) stop("imputer_RNA: 'XX' must be a list containing 'X_star'.")
  if (!is.list(yy) || is.null(yy$y_list)) stop("imputer_RNA: 'yy' must be a valid output list.")
  if (!is.list(pred1) || is.null(pred1$mu) || is.null(pred1$cov)) stop("imputer_RNA: 'pred1' must contain 'mu' and 'cov' (ensure predict was with cov.out=TRUE).")
  if (!is.list(fits) || length(fits) < 1) stop("imputer_RNA: 'fits' must be a non-empty list from RNAmf.")

  L <- length(XX$X_star)
  n_list <- sapply(yy$y_list, length)

  # Sample y_1
  yy$y_tilde[[1]] <- t(mvtnorm::rmvnorm(1, mean = pred1$mu, sigma = pred1$cov))
  yy$y_star[[1]] <- rbind(yy$y_list[[1]],yy$y_tilde[[1]])

  if(L>2){
    # Sample y_2, ..., y_{L-1}
    for(l in 2:(L-1)){
      fit_l <- fits[[l]]
      if (is.null(fit_l$K) || is.null(fit_l$mu.hat)) stop(paste("imputer_RNA: Fit object at level", l, "is malformed."))

      y <- fit_l$y
      K <- fit_l$K
      g <- fit_l$g
      mu.hat <- fit_l$mu.hat
      tau2hat <- fit_l$tau2hat

      ### sampling from Y_tilde give Y_list
      list_idx <- 1:n_list[l]
      K_list <- K[list_idx, list_idx] # K corresponding to X_list
      K_tilde <- K[-list_idx, -list_idx] # K corresponding to X_tilde
      K_list_tilde <- K[list_idx, -list_idx] # K corresponding to X_list given X_tilde
      y_list <- y[list_idx]

      chol_list <- chol(K_list)
      cond_mean <- mu.hat + crossprod(K_list_tilde, backsolve(chol_list, forwardsolve(t(chol_list), y_list - mu.hat)))
      cond_var <- tau2hat * (K_tilde - crossprod(forwardsolve(t(chol_list), K_list_tilde)))
      cond_var <- (cond_var + t(cond_var)) / 2

      y_prior <- t(mvtnorm::rmvnorm(1, mean = cond_mean, sigma = cond_var))
      yy$y_tilde[[l]] <- y_prior
      yy$y_star[[l]] <- rbind(yy$y_list[[l]], yy$y_tilde[[l]])
    }
  }
  return(yy)
}
