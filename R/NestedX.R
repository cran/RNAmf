#' subsetting matrix to be nested
#'
#' @param X1 matrix of the design at the lower level
#' @param X2 matrix of the design at the higher level
#' @return A list containing the design and the length of its subset
#' \itemize{
#'   \item \code{X}: matrix of the design at the lower level.
#'   \item \code{le}: length of the indices of subset.
#' }
#' @noRd
#'

subsetX <- function(X1 = NULL, X2 = NULL) {
  d <- dim(as.matrix(X2))[2] # d2
  n2 <- dim(as.matrix(X2))[1] # n2
  n1 <- dim(as.matrix(X1))[1] # n1

  dist <- 0
  for (i in 1:d) {
    grid <- expand.grid(X2[, i], X1[, i])
    dist <- dist + (grid[, 1] - grid[, 2])^2
  }
  dist.mat <- matrix(dist, n2, n1)
  indice <- max.col(-(dist.mat)) # find the minimum distance of column at each row

  X1 <- matrix(X1[-indice, ], ncol=d)
  X1 <- rbind(X1, matrix(X2, ncol=d))

  return(list(X = X1, le = length(indice)))
}

#' @title Constructing nested design sets for the RNA model.
#'
#' @description The function constructs nested design sets with multiple fidelity levels
#' \eqn{\mathcal{X}_l \subseteq \cdots \subseteq \mathcal{X}_{1}} for use in \code{\link{RNAmf}}.
#'
#' @details The procedure replace the points of lower level design \eqn{\mathcal{X}_{l-1}}
#' with the closest points from higher level design \eqn{\mathcal{X}_{l}}.
#' The length of \eqn{\mathcal{X}_{l-1}} may be larger than the user specified size.
#' For details, see "\href{https://github.com/cran/MuFiCokriging}{\code{NestedDesign}}".
#'
#' @references
#' L. Le Gratiet and J. Garnier (2014). Recursive co-kriging model for design of computer experiments
#' with multiple levels of fidelity. \emph{International Journal for Uncertainty Quantification}, 4(5), 365-386;
#' \doi{doi:10.1615/Int.J.UncertaintyQuantification.2014006914}
#'
#' @param n A vector specifying the number of design points at each fidelity level \eqn{l}. Thus, the vector must have a positive value \eqn{n_1, \ldots, n_l} where \eqn{n_1 > \cdots > n_l}.
#' @param d A positive integer specifying the dimension of the design.
#' @return A list containing the nested design sets at each level, i.e., \eqn{\mathcal{X}_{1}, \ldots, \mathcal{X}_{l}}.
#' @importFrom lhs maximinLHS
#' @export
#' @examples
#' ### number of design points ###
#' n1 <- 30
#' n2 <- 15
#'
#' ### dimension of the design ###
#' d <- 2
#'
#' ### fix seed to reproduce the result ###
#' set.seed(1)
#'
#' ### generate the nested design ###
#' NX <- NestedX(c(n1, n2), d)
#'
#' ### visualize nested design ###
#' plot(NX[[1]], col="red", pch=1, xlab="x1", ylab="x2")
#' points(NX[[2]], col="blue", pch=4)
#'

NestedX <- function(n, d) { # n; vector, d; dim

  # Input validation
  if (is.unsorted(rev(n), strictly = TRUE)) {
    stop("The number of design at each level must be descending order \n")
  }
  if (!all(n > 0)) {
    stop("The number of design at each level must be positive \n")
  }
  if (length(d) != 1) {
    stop("The dimension of design at each level must be same \n")
  }

  level <- length(n)
  if (level < 2) {
    stop("The level of design should be larger than 1 \n")
  }

  # Generating initial designs
  X <- list()
  for (i in 1:level) {
    X[[i]] <- maximinLHS(n[i], d)
  }

  # Subsetting designs
  indices <- list()
  for (i in (level - 1):1) {
    SB <- subsetX(matrix(X[[i]], ncol=d), matrix(X[[i + 1]], ncol=d))
    X[[i]] <- SB$X
    n_rows <- nrow(SB$X)
    indices[[i]] <- seq(n_rows - SB$le + 1, n_rows, by = 1)
  }

  # Constructing nested designs
  X_list <- list()
  current_idx <- 1:nrow(X[[1]])
  for (k in 1:level) {
    if (k > 1) {
      current_idx <- current_idx[indices[[k-1]]]
    }
    X_list[[k]] <- matrix(X[[1]][current_idx, ], ncol = d)
  }

  return(X = X_list)
}
