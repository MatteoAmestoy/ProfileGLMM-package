library(caret) # (dummyVars) one hot encoding of factor data

#' Encodes factor variables dummies LAST is reference
#'
#' @param dataframe dataframe to convert
#'
#' @returns A dataframe with factor variable dummy encoded
#'
#' @examples
encodeCat = function(dataframe){
  ' Encodes factor variables dummies LAST is reference
-------Input:
      - dataframe: dataframe to convert
-------Output:
      - out: dataframe converted
  Note:
  '
  typeCol = sapply(dataframe, class)
  out = data.frame(dataframe[,typeCol!='factor'])
  colnames(out) = colnames(dataframe)[typeCol!='factor']
  for(col in colnames(dataframe)[typeCol=='factor']){
    dmy <- dummyVars(paste0(" ~ ",col), data = dataframe)
    dfTmp = data.frame(predict(dmy, newdata = dataframe))
    nfact = dim(dfTmp)[2]
    out[colnames(dfTmp)[2:nfact]] = dfTmp[,colnames(dfTmp)[2:nfact]]
  }
  return(out)
}
