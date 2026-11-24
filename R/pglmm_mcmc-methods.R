#' Print method for pglmm_mcmc
#'
#' @param x An object of class \code{pglmm_mcmc}
#' @param ... Additional arguments
#'
#' @exportS3Method print pglmm_mcmc
print.pglmm_mcmc <- function(x, ...) {
  cat("--- pglmm_mcmc object containing", length(x$sig2),  " MCMC samples")
  invisible(x)
}

#' Thin draws from a single MCMC chain
#'
#' Retains every \code{thin}-th draw along the iteration dimension of a
#' \code{pglmm_mcmc} object.
#'
#' @param x A \code{pglmm_mcmc} object.
#' @param thin Positive integer thinning step (default \code{1}, keeping every
#'   draw).
#'
#' @return A thinned object of class \code{"pglmm_mcmc"}.
#' @exportS3Method thin_chain pglmm_mcmc
thin_chain <- function(x, thin = 1) {
  if (!inherits(x, "pglmm_mcmc")) {
    stop("'x' must be a 'pglmm_mcmc' object.")
  }

  if (!is.numeric(thin) || length(thin) != 1 || is.na(thin) ||
      thin < 1 || thin %% 1 != 0) {
    stop("'thin' must be a single positive integer (>= 1).")
  }
  thin <- as.integer(thin)

  # Fixed effects are required to have at least an intercept
  n_kept <- if (!is.null(dim(x$beta))) tail(dim(x$beta), 1) else NULL
  if (is.null(n_kept)) {
    stop("Could not determine the number of retained draws from 'x$beta'.")
  }

  if (thin == 1) {
    return(x)
  }

  keep_idx <- seq(from = 1, to = n_kept, by = thin)
  out <- list()

  mat_var = list("WL","WRE","muClus","pvec","gamma")
  vec_var = list("beta","Z")
  scal_var = list("sig2","zeta")

  for (var in mat_var){
    out[[var]] = x[[var]][,,keep_idx, drop = FALSE]
  }
  for (var in vec_var){
    out[[var]] = x[[var]][,keep_idx, drop = FALSE]
  }
  for (var in scal_var){
    out[[var]] = x[[var]][keep_idx,, drop = FALSE]
  }
  out$PhiClus = x[['PhiClus']][keep_idx,, drop = FALSE]
  out$names = x$names

  class(out) <- "pglmm_mcmc"
  return(out)
}

#' Combine posterior draws across multiple MCMC chains
#'
#' Thins each chain (if \code{thin > 1}) and merges retained draws from a
#' list of \code{pglmm_mcmc} objects along the iteration dimension. All
#' chains must have the same number of retained post-burn-in draws.
#'
#' @param x A list of at least two \code{pglmm_mcmc} objects.
#' @param thin Positive integer thinning step passed to \code{\link{thin_chain}}
#'   (default \code{1}).
#'
#' @return An object of class \code{"pglmm_mcmc"} containing the combined draws.
#' @exportS3Method combine_chains pglmm_mcmc
combine_chains <- function(x, thin = 1) {
  if (!is.list(x)) {
    stop("'x' must be a list.")
  }

  if (length(x) < 2) {
    stop("'x' must contain at least two 'pglmm_mcmc' objects.")
  }

  is_valid <- vapply(x, inherits, logical(1), what = "pglmm_mcmc")
  if (!all(is_valid)) {
    stop("All elements of 'x' must inherit from class 'pglmm_mcmc'.")
  }

  # Thin each chain first (also validates 'thin' and draw count)
  thinned_chains <- lapply(x, thin_chain, thin = thin)
  out = {}
  mat_var = list("WL","WRE","muClus","pvec","gamma")
  vec_var = list("beta","Z")
  scal_var = list("sig2","zeta")

  for(var in mat_var){
    out[[var]] <- do.call(abind::abind, c(lapply(thinned_chains, `[[`, var), list(along = 3)))
  }
  for(var in vec_var){
    out[[var]] <- do.call(abind::abind, c(lapply(thinned_chains, `[[`, var), list(along = 2)))
  }
  for(var in scal_var){
    out[[var]] <- do.call(abind::abind, c(lapply(thinned_chains, `[[`, var), list(along = 1)))
  }
  out$PhiClus = do.call(rbind, lapply(thinned_chains, `[[`, 'PhiClus'))


  out$names = x[[1]]$names
  class(out) <- "pglmm_mcmc"
  return(out)
}
