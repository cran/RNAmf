#' computing zeta component in the closed form posterior variance of matern kernel.
#'
#' @param w1 numerical value of \eqn{y_i^{[l-1]}}.
#' @param w2 numerical value of \eqn{y_k^{[l-1]}}.
#' @param m numerical value of \eqn{\mu^*_{l-1}(x)}.
#' @param s numerical value of \eqn{\sigma^{*2}_{l-1}(x)}.
#' @param nu numerical value of smoothness hyperparameter. It should be 1.5 or 2.5.
#' @param theta numerical value of lengthscale hyperparameter \eqn{\theta_{ly}}.
#'
#' @return calculated value of zeta component.
#'
#' @importFrom stats pnorm
#' @noRd
#'

zetafun <- function(w1, w2, m, s, nu, theta) {
  L    <- length(w1)
  sig    <- sqrt(s)
  t2     <- theta^2

  wsmall <- pmin(w1, w2)
  wlarge <- pmax(w1, w2)

  if (nu == 1.5) {
    c3    <- sqrt(3)
    muc <- m - 2 * c3 * s/theta
    mud <- m + 2 * c3 * s/theta

    lam31_2 <- muc; lam31_3 <- muc^2+s; lam32_3 <- muc + wlarge
    lam41_2 <- m; lam41_3 <- m^2+s; lam42_3 <- m + wsmall; lam43_3 <- m + wlarge
    lam51_2 <- -mud; lam51_3 <- mud^2+s; lam52_3 <- -(mud + wsmall)

    e3_1 <- 1 + (3*wsmall*wlarge - c3*theta*(wsmall+wlarge))/t2
    e3_2 <- (2*c3*theta - 3*(wsmall+wlarge))/t2
    e3_3 <- 3/t2

    e4_1 <- 1 + ( -3*wsmall*wlarge + c3*theta*(wlarge-wsmall))/t2
    e4_2 <- 3*(wsmall+wlarge)/t2
    e4_3 <- -3/t2

    e5_1 <- 1 + (3*wsmall*wlarge + c3*theta*(wsmall+wlarge))/t2
    e5_2 <- (2*c3*theta + 3*(wsmall+wlarge))/t2
    e5_3 <- 3/t2

    term1 <- exp((6*s + c3*theta*(wsmall+wlarge - 2*m)) / t2) *
      ( (e3_1 + e3_2*lam31_2 + e3_3*lam31_3) * pnorm((muc-wlarge)/sig) +
          (e3_2 + e3_3*lam32_3) * (sig/sqrt(2*pi)) * exp(-(wlarge-muc)^2/(2*s)) )
    term2 <- exp(-(c3*(wlarge-wsmall)/theta)) *
      ( (e4_1 + e4_2*lam41_2 + e4_3*lam41_3) * (pnorm((wlarge-m)/sig) - pnorm((wsmall-m)/sig)) +
          (e4_2 + e4_3*lam42_3) * (sig/sqrt(2*pi)) * exp(-(wsmall-m)^2/(2*s)) -
          (e4_2 + e4_3*lam43_3) * (sig/sqrt(2*pi)) * exp(-(wlarge-m)^2/(2*s)) )
    term3 <- exp((6*s - c3*theta*(wsmall+wlarge - 2*m)) / t2) *
      ( (e5_1 + e5_2*lam51_2 + e5_3*lam51_3) * pnorm((wsmall-mud)/sig) +
          (e5_2 + e5_3*lam52_3) * (sig/sqrt(2*pi)) * exp(-(wsmall-mud)^2/(2*s)) )
  } else if (nu == 2.5) {
    c5  <- sqrt(5)
    muc <- m - 2*c5*s/theta
    mud <- m + 2*c5*s/theta

    lam31_2 <- muc; lam31_3 <- muc^2+s; lam31_4 <- muc^3+3*muc*s; lam31_5 <- muc^4+6*muc^2*s+3*s^2
    lam32_3 <- muc + wlarge; lam32_4 <- muc^2 + 2*s + wlarge^2 + muc*wlarge; lam32_5 <- wlarge^3 + muc^3 + muc^2*wlarge + muc*wlarge^2 + 3*s*wlarge + 5*muc*s
    lam41_2 <- m; lam41_3 <- m^2+s; lam41_4 <- m^3+3*m*s; lam41_5 <- m^4+6*m^2*s+3*s^2
    lam42_3 <- m + wsmall; lam42_4 <- m^2 + 2*s + wsmall^2 + m*wsmall; lam42_5 <- wsmall^3 + wsmall^2*m + wsmall*m^2 + m^3 + 3*s*wsmall + 5*s*m
    lam43_3 <- m + wlarge; lam43_4 <- m^2 + 2*s + wlarge^2 + m*wlarge; lam43_5 <- wlarge^3 + wlarge^2*m + wlarge*m^2 + m^3 + 3*s*wlarge + 5*s*m
    lam51_2 <- -mud; lam51_3 <- mud^2+s; lam51_4 <- -mud^3-3*mud*s; lam51_5 <- mud^4+6*mud^2*s+3*s^2
    lam52_3 <- -(mud + wsmall); lam52_4 <- mud^2 + 2*s + wsmall^2 + mud*wsmall; lam52_5 <- -mud^3 - wsmall^3 - mud^2*wsmall - mud*wsmall^2 - 3*s*wsmall - 5*mud*s

    e3_1 <- 1 + (25*wsmall^2*wlarge^2 - 3*sqrt(5)*(3*theta^3+5*theta*wsmall*wlarge)*(wsmall+wlarge) + 15*theta^2*(wsmall^2+wlarge^2+3*wsmall*wlarge)) / (9*theta^4)
    e3_2 <- (18*sqrt(5)*theta^3 + 15*sqrt(5)*theta*(wsmall^2+wlarge^2) - 75*theta^2*(wsmall+wlarge) - 50*wsmall*wlarge*(wsmall+wlarge) + 60*sqrt(5)*theta*wsmall*wlarge) / (9*theta^4)
    e3_3 <- 5*(5*wsmall^2+5*wlarge^2+15*theta^2 - 9*sqrt(5)*theta*(wsmall+wlarge)+20*wsmall*wlarge) / (9*theta^4)
    e3_4 <- 10*(3*sqrt(5)*theta - 5*(wsmall+wlarge)) / (9*theta^4)
    e3_5 <- 25/(9*theta^4)

    e4_1 <- 1 + (25*wsmall^2*wlarge^2 + 3*sqrt(5)*(3*theta^3-5*theta*wsmall*wlarge)*(wlarge-wsmall) + 15*theta^2*(wsmall^2+wlarge^2-3*wsmall*wlarge)) / (9*theta^4)
    e4_2 <- 5*(3*sqrt(5)*theta*(wlarge^2-wsmall^2) + 3*theta^2*(wsmall+wlarge) - 10*wsmall*wlarge*(wsmall+wlarge)) / (9*theta^4)
    e4_3 <- 5*(5*wsmall^2+5*wlarge^2-3*theta^2 - 3*sqrt(5)*theta*(wlarge-wsmall)+20*wsmall*wlarge) / (9*theta^4)
    e4_4 <- -50*(wsmall+wlarge) / (9*theta^4)
    e4_5 <- 25/(9*theta^4)

    e5_1 <- 1 + (25*wsmall^2*wlarge^2 + 3*sqrt(5)*(3*theta^3+5*theta*wsmall*wlarge)*(wsmall+wlarge) + 15*theta^2*(wsmall^2+wlarge^2+3*wsmall*wlarge)) / (9*theta^4)
    e5_2 <- (18*sqrt(5)*theta^3 + 15*sqrt(5)*theta*(wsmall^2+wlarge^2) + 75*theta^2*(wsmall+wlarge) + 50*wsmall*wlarge*(wsmall+wlarge) + 60*sqrt(5)*theta*wsmall*wlarge) / (9*theta^4)
    e5_3 <- 5*(5*wsmall^2+5*wlarge^2+15*theta^2 + 9*sqrt(5)*theta*(wsmall+wlarge)+20*wsmall*wlarge) / (9*theta^4)
    e5_4 <- 10*(3*sqrt(5)*theta + 5*(wsmall+wlarge)) / (9*theta^4)
    e5_5 <- 25/(9*theta^4)

    term1 <- exp((10*s + sqrt(5)*theta*(wsmall+wlarge-2*m)) / t2) *
      ( (e3_1 + e3_2*lam31_2 + e3_3*lam31_3 + e3_4*lam31_4 + e3_5*lam31_5) * pnorm((muc-wlarge)/sig) +
          (e3_2 + e3_3*lam32_3 + e3_4*lam32_4 + e3_5*lam32_5) * (sig/sqrt(2*pi)) * exp(-(wlarge-muc)^2/(2*s)) )
    term2 <- exp(-sqrt(5)*(wlarge-wsmall)/theta) *
      ( (e4_1 + e4_2*lam41_2 + e4_3*lam41_3 + e4_4*lam41_4 + e4_5*lam41_5) * (pnorm((wlarge-m)/sig) - pnorm((wsmall-m)/sig)) +
          (e4_2 + e4_3*lam42_3 + e4_4*lam42_4 + e4_5*lam42_5) * (sig/sqrt(2*pi)) * exp(-(wsmall-m)^2/(2*s)) -
          (e4_2 + e4_3*lam43_3 + e4_4*lam43_4 + e4_5*lam43_5) * (sig/sqrt(2*pi)) * exp(-(wlarge-m)^2/(2*s)) )
    term3 <- exp((10*s - sqrt(5)*theta*(wsmall+wlarge-2*m)) / t2) *
      ( (e5_1 + e5_2*lam51_2 + e5_3*lam51_3 + e5_4*lam51_4 + e5_5*lam51_5) * pnorm((wsmall-mud)/sig) +
          (e5_2 + e5_3*lam52_3 + e5_4*lam52_4 + e5_5*lam52_5) * (sig/sqrt(2*pi)) * exp(-(wsmall-mud)^2/(2*s)) )
  }else {
    stop("zetafun: 'nu' must be 1.5 or 2.5")
  }
  return(term1 + term2 + term3)
}
