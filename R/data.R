#' @title Simulated Data and Parameters for a exposure profile linear mixed model
#'
#' @description
#' A list containing a simulated exposure dataset (\code{df}) and the ground-truth parameters
#' (\code{theta0}) used to generate it.
#'
#' The dataset \code{df} contains \eqn{N = 4500} observations across \eqn{n_{Ind} = 1500}
#' individuals, with $n_R = 3$ repeated measures per individual.
#'
#' @format A list with 2 components:
#' \describe{
#' \item{df}{A data frame with 4,500 rows and 6 variables (the simulated data).}
#' \item{theta0}{A list of 11 elements containing the true parameters used for simulation.}
#' }
#'
#' @section \code{df} Data Variables:
#' \describe{
#' \item{X}{Continuous predictor (\eqn{\sim N(0, 1)}).}
#' \item{t}{Time-like variable (structured around 0, 1, 2).}
#' \item{indiv}{**Individual ID** (1 to 1500), the grouping factor.}
#' \item{Exp1, Exp2}{Exposure continuous predictors.}
#' \item{Y}{The **Simulated Response Variable** calculated as: \eqn{\bold{Y} = y_{Fe} + y_{Int} + y_{Re} + \epsilon}, where \eqn{\epsilon ~ N(0, 1)}.}
#' }
#'
#' @section \code{theta0} Parameters:
#' The list \code{theta0} holds the true values used to generate \code{Y}, including:
#' \itemize{
#' \item \code{Lat}: **Categorical Factor** (9 levels), defining the clusters for interaction effects.
#' \item \code{beta}: True fixed effects for the global intercept and \eqn{\bold{X}} (i.e., $(3, 2)$).
#' \item \code{alphaLat}: Vector of 18 coefficients defining the cluster-specific intercepts and slopes for \eqn{\bold{X}} within the 9 \code{Lat} categories.
#' \item \code{alphaRE}: Vector of 1500 random slopes for the time variable \eqn{\bold{t}}, drawn from $N(0, 1)$.
#' \item \code{sigma}: Residual standard deviation (1).
#' }
#'
#' @details
#' The underlying model for the response \eqn{\bold{Y}} is:
#' \deqn{\bold{Y} = \bold{X}_{Fe}\bold{\beta} + \bold{X}_{Int}\bold{\alpha}_{Lat} + \bold{X}_{Re}\bold{\alpha}_{RE} + \bold{\epsilon}}
#'
#' @source Generated synthetically by the package authors.
#' @keywords datasets
"exposure_data"

# --- Documentation for the piecewise dataset starts here ---

#' @title Simulated Data and Parameters for a Piecewise Example
#'
#' @description
#' A list containing a second simulated dataset (\code{df}) and its ground-truth
#' parameters (\code{theta0}). This dataset is generated from a **piecewise linear
#' model**, where the continuous predictor \code{x} is segmented into 6 bins, and
#' different intercept and slope coefficients are applied to each segment.
#'
#' The dataset \code{df} contains $N = 3000$ observations.
#'
#' @format A list with 2 components:
#' \describe{
#' \item{df}{A data frame with 3,000 rows and 2 variables (the simulated data).}
#' \item{theta0}{A list of 5 elements containing the true parameters used for simulation.}
#' }
#'
#' @section \code{df} Data Variables:
#' \describe{
#' \item{x}{A continuous predictor, uniformly distributed between -3 and 3.}
#' \item{Y}{The **Simulated Response Variable** defined by the piecewise linear model.}
#' }
#'
#' @section \code{theta0} Parameters:
#' The list \code{theta0} holds the true values used for simulation, including:
#' \itemize{
#' \item \code{beta}: True global intercept (i.e., (0.5)).
#' \item \code{Lat}: The categorical factor (1 to 6) derived from segmenting \code{x}.
#' \item \code{alphaLat}: Vector of $2 * 6 = 12$ coefficients defining the specific intercept and slope for \code{x} within each of the 6 segments.
#' }
#'
#' @details
#' The underlying model for the response \eqn{\bold{Y}} is:
#' \deqn{\bold{Y} = \bold{X}_{Fe}\bold{\beta} + \bold{X}_{Lat}\bold{\alpha}_{Lat} + \bold{\epsilon}}
#' where \eqn{\bold{X}_{Fe}} is the global intercept, and \eqn{\bold{X}_{Lat}}\eqn{\bold{\alpha}_{Lat}} models the piecewise relationship of \code{x} across the 6 categories defined in \code{theta0$Lat}. The error term \eqn{\bold{\epsilon} ~ N(0, 1)}.
#'
#' @source Generated synthetically by the package authors.
#' @keywords datasets
"piecewise_data"



# --- Documentation for the presaved outputs for examples  ---

#' @title List of the different outputs of the main function for examples
#'
#' @description
#' A list of the different outputs of the main function for examples
#'
#'
#' @format A list with 4 components:
#' \describe{
#' \item{dataProfile}{Output of the profileGLMM_preprocess() function example}
#' \item{MCMC_Obj}{Output of the profileGLMM_Gibbs() function example}
#' \item{post_Obj}{Output of the profileGLMM_postprocess() function example}
#' \item{pred_Obj}{Output of the profileGLMM_predict() function example}
#' }

#' @source Generated synthetically by the package authors.
#' @keywords datasets
"examp"
