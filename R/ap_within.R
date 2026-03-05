# ap_within.R

#' Alternating projections within-transformation
#'
#' Residualizes columns of M w.r.t. multiple fixed effects using
#' alternating projections via collapse::fwithin.
#'
#' @param M numeric matrix (n x k)
#' @param fes list of fixed-effect vectors/factors (length >= 0)
#' @param w numeric weights (length n) or NULL
#' @param tol convergence tolerance
#' @param maxiter max AP iterations
#' @param verbose logical
#' @return numeric matrix, same dim as M
#' @keywords internal
ap_within <- function(M, fes, w = NULL, tol = 1e-8, maxiter = 30, verbose = FALSE) {
  if (!requireNamespace("collapse", quietly = TRUE)) {
    stop("Package 'collapse' is required. Please install.packages('collapse').")
  }
  M <- as.matrix(M)

  if (is.null(fes) || length(fes) == 0) {
    return(collapse::fwithin(M, g = factor(rep(1, nrow(M))), w = w))
  }

  fes_list <- lapply(fes, function(g) if (is.factor(g)) g else factor(g, exclude = NULL))

  # start with grand-mean removal
  R <- collapse::fwithin(M, g = factor(rep(1, nrow(M))), w = w)

  for (it in 1:maxiter) {
    R_old <- R
    for (g in fes_list) {
      R <- collapse::fwithin(R, g = g, w = w)
    }
    diffmax <- max(abs(R - R_old), na.rm = TRUE)
    if (isTRUE(verbose)) cat(sprintf("[AP] it=%02d  maxDelta=%.3e\n", it, diffmax))
    if (!is.finite(diffmax) || diffmax < tol) break
  }
  R
}

#' AP within-transformation (column-block version)
#'
#' @param M numeric matrix (n x k)
#' @param fes list of fixed-effect vectors/factors
#' @param w numeric weights or NULL
#' @param tol convergence tolerance
#' @param maxiter max AP iterations
#' @param verbose logical
#' @param block integer, number of columns per block
#' @return numeric matrix
#' @keywords internal
ap_within_block <- function(M, fes, w = NULL, tol = 1e-8, maxiter = 30,
                            verbose = FALSE, block = 64) {
  M <- as.matrix(M)
  k <- ncol(M)
  out <- matrix(NA_real_, nrow = nrow(M), ncol = k)
  colnames(out) <- colnames(M)
  if (k == 0) return(out)

  idxs <- split(seq_len(k), ceiling(seq_len(k) / block))
  for (ids in idxs) {
    out[, ids] <- ap_within(M[, ids, drop = FALSE], fes, w = w,
                            tol = tol, maxiter = maxiter, verbose = verbose)
  }
  out
}
