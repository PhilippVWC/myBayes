#' @export
max1d = function(vec,epsilon=1.0,maximum=TRUE,maxRows=NA){
  if ( is.vector(vec) & epsilon>=0 & epsilon<=1 ){
    l = length(vec)
    if(maximum){
      vec = data.frame(sort(x = vec,
                            decreasing = TRUE,
                            index.return = TRUE))
    }else{
      vec = data.frame(sort(x = vec,
                            decreasing = FALSE,
                            index.return = TRUE))
    }
    names(x = vec) = c("value","indices")
    if (is.na(maxRows)){
      upper_index = round(epsilon*l,0)
    }else{
      upper_index = min(round(epsilon*l,0),maxRows)
    }
    candidates = vec[1:upper_index,]
    return(candidates)
  }else{
    print("Wrong input data")
    return(NULL)
  }
}
