# Function for 2D matrices returning extrema and its corresponding indices

#' @export
max2d = function(mat,epsilon=1.0,maximum=TRUE,maxRows=NA){
  # PARAMETERS:
  # epsilon       relative tolerance level in [0,1].
  #               Zero: consider only maximum values (possibly degenerated)
  #               One: consider every value
  # mat           Matrix at hand
  # maximum       Boolean value indicating, wether to search for a minimum or a maximum
  # maxRows       Upper bound for number of results
  if ( is.matrix(mat) & length(dim(mat))==2 & epsilon>=0 & epsilon<=1 ){
    dim = dim(mat)
    dimX = dim[1]
    mx = max(mat)
    mn = min(mat)
    tol = abs((mx - mn))*epsilon # convert to absolute tolerance value
    if (maximum == TRUE){
      values = mat[mat >= (mx-tol) & mat <= mx]
      indices = which(mat >= (mx-tol) & mat <= mx)
    } else {
      values = mat[mat >= mn & mat <= (mn+tol)]
      indices = which(mat >= mn & mat <= (mn+tol))
    }
    x_indices = indices %% dimX
    x_indices[x_indices == 0] = dimX        #= dimX, because otherwise values within the last row are interpreted to be in the first row
    y_indices = (indices - 0.1) %/% dimX    #-0.1, because otherwise values within the last row are interpreted to be in the next column
    y_indices = y_indices + 1               #= 1, because R indexes start with number 1
    result = data.frame("value"=values, "RowIndex"=x_indices, "ColIndex"=y_indices)
    nRow = nrow(result)
    if ( !is.na(maxRows) & nRow>maxRows ) {
      print(paste0("### Available results: ",nRow,". Only maxRows = ",maxRows," selected."))
      print("### Increase <<maxRows>> to increase number of results.")
      if(maximum == TRUE){
        result = result[order(result$value,
                              decreasing = TRUE),]
        return(result[1:maxRows,])
      } else {
        result = result[order(result$value,
                              decreasing = FALSE),]
        return(result[1:maxRows,])
      }
    }else{
      return(result)
    }
  } else {
    print("Wrong input data")
    return(data.frame("value"=NA, "RowIndex"=NA, "ColIndex"=NA))

  }
}
