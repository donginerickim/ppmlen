# utils.R

`%||%` <- function(x, y) if (!is.null(x)) x else y

fix_names <- function(m) {
  m <- as.matrix(m)
  rn <- rownames(m)
  if (is.null(rn)) rn <- rep("", nrow(m))
  bad <- is.na(rn) | rn == ""
  rn[bad] <- paste0("x", which(bad))
  rn <- make.unique(rn, sep = "_")
  rownames(m) <- rn
  m
}

quiet_call <- function(expr) {
  suppressWarnings(
    suppressMessages({
      .res <- NULL
      utils::capture.output({ .res <- eval.parent(substitute(expr)) })
      .res
    })
  )
}
