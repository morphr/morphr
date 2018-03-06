#' Convert q to curve
#' 
#' 
#' @param q vector of coordinates of curves
#' @return numeric vector 
#' @export
q_to_curve <- function(q){
  if(is.null(ncol(q[[1]])) != TRUE){
    q = q[[1]]
  }

  n = nrow(q)
  T_col = ncol(q)
  s = pracma::linspace(0,2*pi,T_col)
  qnorm = list()
  for (i in 1:T_col){
    qnorm[i] = pracma::Norm(q[,i])
  }
  qnorm <- unlist(qnorm)
  p = list()
  for (i in 1:n){
    temp = q[i,]*qnorm
    p[[i]] <- pracma::cumtrapz(s,temp)
  }
  p <- matrix(unlist(p), ncol = length(p[[1]]), byrow = TRUE)
  return(p)
}