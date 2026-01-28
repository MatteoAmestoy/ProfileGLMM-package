#' Print method for pglmm_fit
#'
#' @param x An object of class \code{pglmm_fit}
#' @param ... Additional arguments
#'
#' @exportS3Method print pglmm_fit
print.pglmm_fit <- function(x, ...) {
  cat("--- pglmm_fit object containing the postprocessing of the MCMC samples ---\n")
  # Use a safe check (like is.null) to prevent errors during printing
  if (!is.null(x$clust$Kstar)) {
    cat("-- Representative clustering characteristics: \n")
    cat("Method used:", x$clust$mode, "\n")
    cat("Clusters found:", x$clust$Kstar, "\n")
  } else {
    cat("Representative clustering not computed.\n")
  }
  invisible(x)
}

#' Print method for pglmm_fit
#'
#' @param x An object of class \code{pglmm_fit}
#' @param ... Additional arguments
#'
#' @exportS3Method print pglmm_fit
summary.pglmm_fit <- function(x, ...) {
  cat("--- pglmm_fit summary ---\n")

  cat("-- Fixed effect estimates: \n")
  print(x$pop$betaFE)
  cat(" \n")
  if (!is.null(x$clust$Kstar)) {
    cat("-- Representative clustering estimates -- \n")
    cat("Method used:", x$clust$mode,", clusters found:", x$clust$Kstar, "\n")
    cat("Cluster centroids: \n")
    print(x$clus$cen)
    cat("Cluster interaction parameters: \n")
    print(x$clus$gamma)
  } else {
    cat("Representative clustering not computed.\n")
  }
  invisible(x)
}




#' @title Prediction of cluster memberships and outcomes
#' @description (This documentation is now for internal use only)
#' @param object An object of class \code{pglmm_fit} .
#' @param newData : A list with fields\itemize{
#' \item XFE A numeric matrix of fixed effects covariates for the prediction data.
#' \item XLat A numeric matrix of latent effect covariates.
#' \item UCont A numeric matrix or vector of continuous profile variables. Defaults to \code{NULL}.
#' \item UCat A numeric matrix or vector of categorical profile variables. Defaults to \code{NULL}.}
#' @param ... Additional arguments
#'
#' @exportS3Method predict pglmm_fit
#'
#' @examples
#' # Load post_Obj, the result of profileGLMM_postProcess()
#' data("examp")
#' post_Obj = examp$post_Obj
#'
#' # run prediction for training data
#' pred_Obj = predict(post_Obj,examp$dataProfile)
#'
#'
predict.pglmm_fit <- function(x, newData, ...) {
  return(profileGLMM_predict(x, newData$XFE, newData$XLat, newData$UCont, newData$UCat))
}


#' Print method for pglmm_fit
#'
#' @param x An object of class \code{pglmm_fit}
#' @param ... Additional arguments
#'
#' @exportS3Method print pglmm_fit
print.pglmm_fit <- function(x, ...) {
  cat("--- pglmm_fit object containing the postprocessing of the MCMC samples ---\n")
  # Use a safe check (like is.null) to prevent errors during printing
  if (!is.null(x$clust$Kstar)) {
    cat("-- Representative clustering characteristics: \n")
    cat("Method used:", x$clust$mode, "\n")
    cat("Clusters found:", x$clust$Kstar, "\n")
  } else {
    cat("Representative clustering not computed.\n")
  }
  invisible(x)
}

#' Plot method for pglmm_fit continuous covariates cluster characteristics
#'
#' @param x An object of class \code{pglmm_fit}
#' @param ... Additional arguments\itemize{
#' \item title : main title of the plot
#' \item color : palette to be used }
#'
#' @exportS3Method plot pglmm_fit
plot.pglmm_fit <- function(x, ...) {
  if (is.null(x$clust)) {
    stop("No representative clustering provided.")
  }

  # 1. Extract parameters
  args <- list(...)
  main_title <- if (!is.null(args$title)) args$title else "Cluster distributions"
  my_colors <- if (!is.null(args$color)) args$color else hcl.colors(nC, palette = "plasma")

  if(length(my_colors)!=nC){
    warning(paste0('Number of colors provided is incorrect, ',nC,' expected. reverting to default palette. '))
    my_colors <- hcl.colors(nC, palette = "plasma")
    }



  nC <- x$clust$Kstar
  cen <- x$clust$cen
  coVar <- x$clust$coVar
  d <- nrow(cen) # Number of dimensions


  # --- Internal Helper: Ellipse ---
  get_ellipse_pts <- function(mu, Sigma) {
    theta <- seq(0, 2 * pi, length.out = 100)
    circle <- rbind(cos(theta), sin(theta))
    ed <- eigen(Sigma)
    # Ensure eigenvalues aren't negative due to precision issues
    scale_matrix <- ed$vectors %*% diag(sqrt(pmax(ed$values, 0))) %*% t(ed$vectors)
    return(scale_matrix %*% circle + mu)
  }

  # --- Internal Helper: Higher Dim Plotting ---
  plot_multivariate_internal <- function(cen, coVar, colors, title_text) {
    K <- nrow(cen)
    lay_mat <- matrix(0, nrow = K, ncol = K)
    diag(lay_mat) <- 1:K

    current_id <- K + 1
    for (i in 2:K) {
      for (j in 1:(i-1)) {
        lay_mat[i, j] <- current_id
        current_id <- current_id + 1
      }
    }

    # Merge the entire upper triangle into one legend area
    legend_id <- current_id
    lay_mat[upper.tri(lay_mat)] <- legend_id

    layout(lay_mat)
    par(mar = c(3, 3, 2, 1), oma = c(1, 1, 3, 1))

    dim_limits <- lapply(1:K, function(dim_idx) {
      all_reaches <- sapply(1:nC, function(k) 3 * sqrt(coVar[dim_idx, dim_idx, k]))
      range(c(cen[dim_idx, ] - all_reaches, cen[dim_idx, ] + all_reaches))
    })

    # Diagonals
    for (i in 1:K) {
      xlims <- dim_limits[[i]]
      max_dens <- max(sapply(1:nC, function(k) dnorm(cen[i,k], cen[i,k], sqrt(coVar[i,i,k]))))
      plot(NULL, xlim = xlims, ylim = c(0, max_dens * 1.1),
           main = paste("Dim", i), yaxt = 'n', xlab="", ylab="")
      for (k in 1:nC) {
        x_seq <- seq(xlims[1], xlims[2], length.out = 100)
        lines(x_seq, dnorm(x_seq, cen[i, k], sqrt(coVar[i, i, k])), col = colors[k], lwd = 2)
      }
    }

    # Lower Triangle
    for (target in (K + 1):(legend_id - 1)) {
      pos <- which(lay_mat == target, arr.ind = TRUE)
      row_idx <- pos[1]; col_idx <- pos[2]
      plot(NULL, xlim = dim_limits[[col_idx]], ylim = dim_limits[[row_idx]], asp = 1)
      grid()
      for (k in 1:nC) {
        pts <- get_ellipse_pts(cen[c(col_idx, row_idx), k], coVar[c(col_idx, row_idx), c(col_idx, row_idx), k])
        lines(pts[1, ], pts[2, ], col = colors[k], lwd = 2)
        points(cen[col_idx, k], cen[row_idx, k], col = colors[k], pch = 19, cex = 0.5)
      }
    }

    # Legend
    plot.new()
    legend("center", legend = paste("Cluster", 1:nC),
           col = colors, lwd = 3, pch = 19, bty = "n", cex = 1.2)
    mtext(title_text, outer = TRUE, cex = 1.3, line = 1)
  }

  # --- Branching Logic ---
  if (d == 1) {
    # 1D Case
    x_min <- min(cen - 3.5 * coVar)
    x_max <- max(cen + 3.5 * coVar)
    xxx <- seq(x_min, x_max, length.out = 1000)

    plot(xxx, dnorm(xxx, mean = cen[1,1], sd = sqrt(coVar[1,1,1])), type = "n",
         ylim = c(0, max(1 / (sqrt(as.vector(coVar)) * sqrt(2 * pi)))),
         ylab = "Density", xlab = "x", main = main_title)

    for (k in 1:nC) {
      lines(xxx, dnorm(xxx, mean = cen[1, k], sd = sqrt(coVar[1, 1, k])), col = my_colors[k], lwd = 2)
    }
    legend("topright", legend = paste("Cluster", 1:nC), col = my_colors, lwd = 2, bty = "n")

  } else if (d == 2) {
    # 2D Case
    all_x <- c(); all_y <- c()
    for (k in 1:nC) {
      reach_x <- 3 * sqrt(coVar[1, 1, k])
      reach_y <- 3 * sqrt(coVar[2, 2, k])
      all_x <- c(all_x, cen[1, k] + c(-reach_x, reach_x))
      all_y <- c(all_y, cen[2, k] + c(-reach_y, reach_y))
    }

    plot(NULL, xlim = range(all_x), ylim = range(all_y), asp = 1,
         xlab = "Dim 1", ylab = "Dim 2", main = main_title)
    grid()
    for (k in 1:nC) {
      pts <- get_ellipse_pts(cen[, k], coVar[,,k])
      lines(pts[1, ], pts[2, ], col = my_colors[k], lwd = 2)
      points(cen[1, k], cen[2, k], col = my_colors[k], pch = 19, cex = 0.6)
    }
    legend("topright", legend = paste("Cluster", 1:nC), col = my_colors, lwd = 2, pch = 19, bty = "n")

  } else {
    # Higher Dimensional Case
    plot_multivariate_internal(cen, coVar, my_colors, main_title)
  }
}




