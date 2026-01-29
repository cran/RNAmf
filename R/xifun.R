#' computing xi component in the closed form posterior mean of matern kernel.
#'
#' @param w numerical value of \eqn{y_i^{[l-1]}}.
#' @param m numerical value of \eqn{\mu^*_{l-1}(x)}.
#' @param s numerical value of \eqn{\sigma^{*2}_{l-1}(x)}.
#' @param nu numerical value of smoothness hyperparameter. It should be 1.5 or 2.5.
#' @param theta numerical value of lengthscale hyperparameter \eqn{\theta_{ly}}.
#'
#' @return calculated value of xi component.
#'
#' @importFrom stats pnorm
#' @noRd
#'

xifun <- function(w, m, s, theta, nu) {
  L <- length(w)
  sig <- sqrt(s)
  t2 <- theta^2

  if (nu == 1.5) {
    c3 <- sqrt(3)
    mua <- m - c3 * s / theta
    mub <- m + c3 * s / theta

    e1_1 <- (theta - c3 * w) / theta
    e1_2 <- c3 / theta
    e2_1 <- (theta + c3 * w) / theta

    lam11_2 <- mua
    lam21_2 <- -mub

    term1 <- exp((3*s + 2*c3*theta*(w - m)) / (2*t2)) *
      ( (e1_1 + e1_2 * lam11_2) * pnorm((mua - w)/sig) +
          e1_2 * sig / sqrt(2 * pi) * exp(-(w - mua)^2/(2*s)) )
    term2 <- exp((3*s - 2*c3*theta*(w - m)) / (2*t2)) *
      ( (e2_1 + e1_2 * lam21_2) * pnorm((w - mub)/sig) +
          e1_2 * sig / sqrt(2 * pi) * exp(-(w - mub)^2/(2*s)) )
  } else if (nu == 2.5) {
    c5 <- sqrt(5)
    mua <- m - c5 * s / theta
    mub <- m + c5 * s / theta

    e1_1 <- 1 - c5*w/theta + 5*w^2/(theta^2*3)
    e1_2 <- c5/theta - 10*w/(theta^2*3)
    e1_3 <- 5/(theta^2*3)
    e2_1 <- 1 + c5*w/theta + 5*w^2/(theta^2*3)
    e2_2 <- c5/theta + 10*w/(theta^2*3)

    lam11_2 <- mua
    lam11_3 <- mua^2 + s
    lam12_3 <- mua + w
    lam21_2 <- -mub
    lam21_3 <- mub^2 + s
    lam22_3 <- -(mub + w)

    term1 <- exp((5*s + 2*c5*theta*(w - m)) / (2*t2)) *
      ( (e1_1 + e1_2 * lam11_2 + e1_3 * lam11_3) * pnorm((mua - w)/sig) +
          (e1_2 + e1_3 * lam12_3) * sig / sqrt(2 * pi) * exp(-(w - mua)^2/(2*s)) )
    term2 <- exp((5*s - 2*c5*theta*(w - m)) / (2*t2)) *
      ( (e2_1 + e2_2 * lam21_2 + e1_3 * lam21_3) * pnorm((w - mub)/sig) +
          (e2_2 + e1_3 * lam22_3) * sig / sqrt(2 * pi) * exp(-(w - mub)^2/(2*s)) )
  }else {
    stop("xifun: 'nu' must be 1.5 or 2.5")
  }
  return(term1 + term2)
}

