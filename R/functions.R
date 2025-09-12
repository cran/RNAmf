#' scale inputs before fit the model
#'
#' @param X vector or matrix of input locations to be scaled.
#' @param back logical indicating for scale back to the original vector or matrix (back=TRUE) or scale the original vector or matrix (back=FALSE). Default is FALSE.
#'
#' @return A scaled X or original X.
#'
#' @noRd
#'

scale_inputs <- function(X, center = NULL, scale = NULL, back = FALSE) {
  if (back) {
    if (is.null(center) | is.null(scale)) stop("center and scale are required to scale back")
    X <- t(t(X) * scale + center)
  } else {
    if (is.null(center)) center <- attr(scale(X), "scaled:center")
    if (is.null(scale)) scale <- attr(scale(X), "scaled:scale")

    X <- t((t(X) - center) / scale)
    attr(X, "scaled:center") <- center
    attr(X, "scaled:scale") <- scale
  }
  return(X)
}

#' Check the design is nested
#'
#' @param XX1 vector or matrix of input locations at lower fidelity
#' @param XX2 vector or matrix of input locations at higher fidelity
#'
#' @return A logical indicating if XX2 is nested or not.
#'
#' @noRd
#'

checknested <- function(XX1, XX2) {
  checknest <- c()
  for (i in 1:nrow(XX2)) {
    checknest <- c(checknest, suppressWarnings(any(apply(XX1, 1, function(xx) {
      all.equal(XX2[i, ], xx, tolerance = sqrt(.Machine$double.eps))
    }))))
  }
  checknest[is.na(checknest)] <- FALSE
  all(checknest)
}

#' Match rows between two nested designs
#'
#' @param X_sup Matrix of the design at the lower fidelity level (the superset).
#' @param X_sub Matrix of the design at the higher fidelity level (the subset).
#' @return An integer vector of indices such that each element corresponds
#'   to the row position in \code{X_sup} matching the respective row in \code{X_sub}.
#' @details
#' This function is used to align the higher-level design points with their
#' corresponding points in the lower-level design when the designs are nested.
#' Matching is done by converting each row to a unique key string after rounding
#' to 12 decimal places, to avoid mismatches caused by floating-point precision.
#' @noRd

.match_rows <- function(X_sup, X_sub) {
  as_key <- function(M) {
    if (!is.matrix(M)) M <- as.matrix(M)
    # Optional: set a fixed rounding to avoid floating point mismatch
    M <- round(M, 12)  # fixed rounding, not user-specified
    do.call(paste, c(as.data.frame(M), sep = "\r"))
  }
  sup_key <- as_key(X_sup)
  sub_key <- as_key(X_sub)
  match(sub_key, sup_key)
}



