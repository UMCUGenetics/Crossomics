convert <- function(a) {as.vector(unlist(lapply(a, function(x) {
  if (identical(unlist(x),character(0))){
    return(NA)
  } else {
    return(x)
  }
})))}
