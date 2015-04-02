
pairedPermuTest <- function(cases1, cases2, confidenceLevel=0.95) {
  # Parameters validation
  
# The confidence level must be a value between 0 and 1
  if (!is.numeric(confidenceLevel || confidenceLevel < 0 
        || confidenceLevel > 1) {
    stop("The confidence level must be a numeric value between 0 and 1.", 
      call. = FALSE)
  }
  
  call <- match.call()
  
  
  return(list(call))
}