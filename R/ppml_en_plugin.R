# ppml_en_plugin.R

#' Elastic-net plugin PPML with high-dimensional fixed effects
#'
#' Plugin-lasso seeded elastic-net PPML estimator with multi-way fixed effects
#' using alternating-projections within-transformation (collapse::fwithin).
#'
#' @param data data.frame
#' @param dep character. Dependent variable name.
#' @param indep character vector or NULL.
#' @param fixed list. Fixed effects spec passed to penppml::genmodel.
#' @param cluster cluster spec passed to penppml::genmodel.
#' @param alpha elastic-net alpha in (0,1].
#' @param K integer. Number of outer EN iterations.
#' @param tol numeric. Tolerance for hdfeppml seed stage.
#' @param hdfetol numeric. HDFE tolerance for seed stage.
#' @param glmnettol numeric. Coordinate descent tolerance.
#' @param phipost logical. Passed to penppml seed stage.
#' @param verbose logical.
#' @param post logical. (kept for compatibility; not used in EN stage here)
#' @param lasso_engine character. "en" or "penppml" (kept; seed uses penppml plugin-lasso).
#' @param always_keep character vector of covariates to never penalize.
#' @param fes_mode character. currently only "ap" is implemented.
#' @param ap_tol numeric. AP convergence tolerance.
#' @param ap_maxiter integer. AP max iterations.
#'
#' @return list with coefficients, fitted mean, selection, and diagnostics.
#' @export
ppml_en_plugin <- function(
    data, dep, indep = NULL, fixed = NULL, cluster = NULL,
    alpha = 0.5, K = 15, tol = 1e-8, hdfetol = 1e-4, glmnettol = 1e-12,
    phipost = TRUE, verbose = FALSE, post = FALSE,
    lasso_engine = c("en", "penppml"),
    always_keep = NULL,
    fes_mode = c("dense", "ap"),
    ap_tol = 1e-8, ap_maxiter = 30
) {

  if (!requireNamespace("penppml", quietly = TRUE)) {
    stop("Package 'penppml' is required. Install it with install.packages('penppml').")
  }

  lasso_engine <- match.arg(lasso_engine)
  fes_mode <- match.arg(fes_mode)
  stopifnot(alpha > 0 && alpha <= 1)

  if (!requireNamespace("collapse", quietly = TRUE)) {
    stop("Please install.packages('collapse')")
  }
  if (!requireNamespace("penppml", quietly = TRUE)) {
    stop("Please install.packages('penppml')")
  }

  if (fes_mode != "ap") {
    stop("fes_mode='dense' is not implemented in this version. Use fes_mode='ap'.")
  }

  logf <- function(...) if (isTRUE(verbose)) cat(...)

  # pull internal helpers from penppml
  pen_ns <- base::asNamespace("penppml")

  genmodel         <- utils::getFromNamespace("genmodel", pen_ns)
  cluster_matrix   <- utils::getFromNamespace("cluster_matrix", pen_ns)
  hdfeppml_int_fn  <- utils::getFromNamespace("hdfeppml_int", pen_ns)
  pen_plugin_lasso <- utils::getFromNamespace("penhdfeppml_cluster_int", pen_ns)

  .hdfe_call <- function(...) {
    dots <- list(...)
    allowed <- names(formals(hdfeppml_int_fn))
    do.call(hdfeppml_int_fn, dots[intersect(names(dots), allowed)])
  }

  as_mat <- function(M) {
    M <- as.matrix(M)
    if (is.null(dim(M))) M <- matrix(M, ncol = 1)
    M
  }
  soft <- function(z, g) sign(z) * pmax(abs(z) - g, 0)

  gm <- quiet_call(genmodel(data = data, dep = dep, indep = indep, fixed = fixed, cluster = cluster))
  y <- as.numeric(gm$y)
  X_full <- as_mat(gm$x)
  fes <- gm$fes
  cl <- gm$cluster
  n <- length(y)

  keep_idx <- which(colnames(X_full) %in% (always_keep %||% character(0)))

  # ---------- SEED ----------
  logf("\n---- Seed(plugin-lasso) stage ----\n")
  seed <- tryCatch(
    quiet_call(
      pen_plugin_lasso(
        y = y, x = X_full, fes = fes, cluster = cl,
        penalty = "lasso", post = FALSE, K = 5,
        tol = tol, hdfetol = hdfetol, glmnettol = glmnettol,
        phipost = phipost, verbose = FALSE
      )
    ),
    error = function(e) NULL
  )

  if (is.null(seed)) {
    logf("[SEED] plugin-lasso failed -> fallback to theory-based seed\n")
    only_fes <- .hdfe_call(
      y = y, fes = fes, tol = tol, hdfetol = hdfetol,
      mu = NULL, saveX = TRUE, init_z = NULL,
      verbose = FALSE, maxiter = 1000, cluster = NULL, vcv = TRUE
    )
    mu <- pmax(1e-190, pmin(as.numeric(only_fes$mu), 1e190))
    z <- (y - mu) / mu + log(mu)

    z_resid <- ap_within(matrix(z), fes, w = mu, tol = ap_tol, maxiter = ap_maxiter)
    X_resid <- ap_within_block(X_full, fes, w = mu, tol = ap_tol, maxiter = ap_maxiter)

    phi <- sqrt(diag(cluster_matrix(mu * z_resid, cl, X_resid)) / n)
    k0 <- ncol(X_resid)

    lam_pkg <- 1.1 * sqrt(n) *
      stats::qnorm(1 - 0.1 / log(n + 10) / (2 * max(1, k0))) / sum(y)

    seed <- list(mu = mu, z_resid = z_resid, x_resid = X_resid, phi = phi, lambda = lam_pkg)
  }

  mu <- seed$mu
  z_resid <- seed$z_resid
  X_resid <- seed$x_resid
  phi <- seed$phi
  lam_pkg <- seed$lambda

  # normalize phi once, keep fixed in EN stage
  phi[!is.finite(phi)] <- 1
  if (length(phi) > 1) {
    q05 <- stats::quantile(phi, 0.05, na.rm = TRUE)
    q95 <- stats::quantile(phi, 0.95, na.rm = TRUE)
    phi <- pmin(pmax(phi, q05), q95)
  }
  phi <- phi / mean(phi, na.rm = TRUE)

  k_eff <- ncol(X_resid)
  sum_phi <- sum(phi, na.rm = TRUE)
  if (!is.finite(sum_phi) || sum_phi <= 0) sum_phi <- k_eff
  lambda_glmnet <- lam_pkg * sum_phi / max(1, k_eff)

  logf(sprintf("[SEED-END] k_eff=%d  lambda=%.4g  phi_bar=%.4f\n",
               k_eff, lambda_glmnet, mean(phi, na.rm = TRUE)))

  # ---------- EN STAGE ----------
  logf("\n---- EN(elastic-net) stage ----\n")

  # align X with residualized columns
  X <- X_full[, colnames(X_resid), drop = FALSE]

  penfac <- phi
  if (length(keep_idx)) penfac[keep_idx] <- 0

  b <- matrix(0, nrow = k_eff, ncol = 1, dimnames = list(colnames(X_resid), NULL))
  z <- (y - mu) / mu + log(mu)
  phi_hist <- list()

  for (it in 1:K) {
    w <- as.numeric(mu / sum(mu))
    lam_use <- lambda_glmnet
    lam1 <- alpha * lam_use
    lam2 <- (1 - alpha) * lam_use

    # coordinate descent
    for (outer in 1:200) {
      b_old <- b
      XB <- as.numeric(X_resid %*% b)

      for (j in 1:k_eff) {
        xj <- X_resid[, j]
        Hjj <- sum(w * xj * xj)

        rj <- z_resid - XB + xj * b[j]
        zj <- sum(w * xj * rj)

        if (j %in% keep_idx) {
          bj <- zj / Hjj
        } else {
          bj <- soft(zj, lam1 * penfac[j]) / (Hjj + 2 * lam2)
        }

        if (!is.finite(bj)) bj <- 0
        if (bj != b[j]) {
          XB <- XB + xj * (bj - b[j])
          b[j] <- bj
        }
      }

      if (max(abs(b - b_old)) < glmnettol) break
    }

    # IRLS update
    resid <- as.numeric(z_resid - X_resid %*% b)
    mu <- pmax(1e-190, pmin(exp(z - resid), 1e190))
    z <- (y - mu) / mu + log(mu)

    # re-residualize using AP
    z_resid <- ap_within(matrix(z), fes, w = mu, tol = ap_tol, maxiter = ap_maxiter)
    X_resid <- ap_within_block(X, fes, w = mu, tol = ap_tol, maxiter = ap_maxiter)

    # keep phi fixed
    phi_hist[[it]] <- phi

    dev_i <- mean(-2 * (y * log(y / mu) - (y - mu)), na.rm = TRUE)
    logf(sprintf("EN iter=%02d alpha=%.2f nonzero=%d lambda=%.4g phi_bar=%.4f dev=%.6f\n",
                 it, alpha, sum(b != 0), lambda_glmnet, mean(phi), dev_i))
  }

  dev <- mean(-2 * (y * log(y / mu) - (y - mu)), na.rm = TRUE)
  bic <- dev + sum(b != 0) * log(n) / n

  list(
    beta = fix_names(b),
    mu = mu,
    deviance = dev,
    bic = bic,
    phi = phi,
    lambda = lam_pkg,
    lambda_glmnet = lambda_glmnet,
    alpha = alpha,
    selected = rownames(b)[as.numeric(b) != 0],
    phi_hist = phi_hist,
    seed_snapshot = seed
  )
}
