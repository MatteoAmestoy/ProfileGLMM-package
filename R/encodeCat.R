#' @title One-Hot Encodes Factor Variables (FIRST Level as Reference)
#'
#' @description This function takes a dataframe, identifies all columns of class \code{factor}, and converts them into **dummy variables** using one-hot encoding via \code{stats::model.matrix}. For each factor, the function explicitly removes the first dummy variable generated, effectively making the **first level** of the factor the **reference level** (omitted category). Non-factor columns are retained as is.
#'
#' @param dataframe A \code{data.frame} containing the data to be processed, which may include factor variables.
#'
#' @returns A \code{data.frame} where:
#' \itemize{
#'   \item{All original non-factor columns are present.}
#'   \item{All original factor columns are replaced by a set of binary (0/1) dummy variables. The first level of the factor is excluded from the generated dummies, making the last level the reference.}
#' }
#' @export
#' @importFrom stats model.matrix as.formula
#'
#' @examples
#' \dontrun{
#' # Example usage (assuming 'mydata' is available):
#' # mydata$Gender = as.factor(mydata$Gender)
#' # encoded_data = encodeCat(mydata)
#'
#' }
encodeCat = function(dataframe){
  typeCol = sapply(dataframe, class)

  out_df = dataframe[, typeCol != 'factor', drop = FALSE]

  factor_cols = colnames(dataframe)[typeCol == 'factor']

  for (col_name in factor_cols) {
    # print(col_name)
    formula_str = paste0("~ ", col_name, " - 1")

    dummy_matrix = stats::model.matrix(
      as.formula(formula_str),
      data = dataframe
    )

    dummy_df = as.data.frame(dummy_matrix[,2:length(colnames(dummy_matrix))])
    # print(colnames(dummy_df))

    new_names = gsub(paste0(col_name), paste0(col_name, "."), colnames(dummy_df))
    # print(new_names)
    colnames(dummy_df) = new_names

    out_df = cbind(out_df, dummy_df)
  }

  return(out_df)
}
