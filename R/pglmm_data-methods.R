#' Print method for pglmm_data
#'
#' @param x An object of class \code{pglmm_data}
#' @param ... Additional arguments
#'
#' @exportS3Method print pglmm_data
print.pglmm_data <- function(x, ...) {
  cat("--- pglmm Data summary ---\n")
  cat("- Number of observations : ",x$params$n,  "\n")
  cat("\n")
  cat("-- Clustering model summary --\n")
  if(x$params$catInd!=-1){
    cat("- Categorical clustering variables : '",x$d$names$UCat,"'")
    if(x$params$contInd!=-1){
      cat(" -\n- Continuous clustering variables : '",x$d$names$UCont,"'")
    }
  }else {
    cat("- Continuous clustering variables : '",x$d$names$UCont,"'")
  }
  cat("\n")
  cat("\n")
  cat("-- Outcome model summary --\n")
  if (x$regType == 0){
    cat("- Model type:  linear mixed model\n")
  }else if(x$regType == 1){
    cat("- Model type:  probit mixed model\n")
  }
  cat("- Outcome  : '",x$d$names$Y,  "' \n")
  cat("- Fixed effects  : '",x$d$names$FE,  "'\n")
  if (x$d$ZRE[1]==-1){
    cat("- No random effect \n")
  }else {
    cat("- '", x$d$names$REunit ,"' level random effects  : '", x$d$names$RE,  "'\n")
  }
  cat("- Latent clusters interacting with: '",x$d$names$Lat,  "'\n")



  invisible(x)
}
