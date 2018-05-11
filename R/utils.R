#' Compute covariance of a set of tangent vectors (alpha_t_array) at the mean shape (qmean)
#' @param qmean matrix containing the mean shape
#' @param alpha_t_array list of tangent vectors from the mean shape
#' @return covdata list containing covariance matrix and the coefficients of eigen projections
#' @export
compute_covariance <- function(qmean, alpha_t_array){
  d = 20
  T_col = ncol(qmean)
  epsilon = 0.0001
  if(nrow(qmean) == 2){
    temp <- matrix(0,1,ncol = ncol(qmean))
    qmean <- rbind(qmean, temp)
  }
  if(nrow(alpha_t_array[[1]]) == 2){
    for(i in 1:length(alpha_t_array)){
      temp <- matrix(0,1,ncol = ncol(alpha_t_array[[i]]))
      alpha_t_array[[i]] <- rbind(alpha_t_array[[i]], temp)
    }
  }
  B = form_basis_l2_r3(d, T_col)
  B = gram_schmidt(B)
  Bnew = project_basis_tangent_space_of_c_at_q(B, qmean)
  G_O_q = form_basis_o_q(B,qmean)
  G_O_q = gram_schmidt(G_O_q)
  Bnew = gram_schmidt(Bnew)
  G = form_basis_of_tangent_space_of_s_at_q(Bnew, G_O_q)
  G = gram_schmidt(G)
  temp = project_to_basis(alpha_t_array,G)
  Aproj = temp[[2]]
  A = temp[[1]]
  C = cov(Aproj)
  temp2 = svd(C)
  S = diag(temp2$d)
  U = temp2$u
  V = temp2$v
  PC_Latent_Explained = pcacov(C)
  PC = PC_Latent_Explained[[1]]
  Latent = PC_Latent_Explained[[2]]
  Explained = PC_Latent_Explained[[3]]
  Cn = U%*%(S+epsilon*pracma::eye(ncol(S))%*%t(U))
  covdata <- list(C,Cn,U,S,V,Aproj,A,G)
  names(covdata) <- c('C', 'Cn', 'U', 'S', 'V', 'Aproj', 'A', 'G')
  
  return(covdata)
  
}





#' Convert coordinate curve to the q function
#' @param p coordinate curve
#' @return q function
#' @export
curve_to_q <- function(p){
  n <- nrow(p)
  Tcoord <- ncol(p)
  pdiff <- list()
  for( i in 1:n){
    pdiff[[i]] <- pracma::gradient(p[i,], 2*pi/(Tcoord-1))
  }
  v <- matrix(0, nrow = n, ncol = Tcoord)
  
  for(i in 1:n){
    v[i,] <- 2*pi/(Tcoord-1)*pdiff[[i]]
  }
  q <- matrix(0, nrow = n, ncol = Tcoord)
  for(i in 1:Tcoord){
    q[,i] <- v[,i]/sqrt(norm(v[,i], type = "2"))
  }
  q <- project_curve(q)
  return(q)
}

#' Compute PCA of matrix of observations
#' @param C matrix of observations
#' @return list of principal coefficients, latent (vector of eigenvalues of C), explained (percent of total variance)
#' @export
pcacov <- function(C){
  latent = svd(C)$d
  coeff = svd(C)$v
  totalvar = sum(latent)
  explained = 100*latent/totalvar
  p = nrow(coeff)
  d = ncol(coeff)
  maxind = pracma::zeros(1,ncol(coeff))
  coeff_temp = abs(coeff)
  for (i in col(coeff)){
    maxind[1,i] = which.max( coeff_temp[,i] )
  }
  maxind+pracma::linspace(0,(d-1)*p,p)
  colsign = sign(coeff[maxind+pracma::linspace(0,(d-1)*p,p)])
  colsign <- matrix(colsign,ncol(C),1)
  for (i in row(coeff)){
    coeff[i,] = coeff[i,]*colsign
  }
  return(list(coeff,latent,explained))
}




#' Numerical integration of the function y(x)
#' @param x the set of points where y is computed on
#' @param y the integrand
#' @return the result
trapz <- function(x,y,dims=1){
  if ((dims-1)>0){
    perm = c(dims:max(fdasrvf:::ndims(y),dims), 1:(dims-1))
  } else {
    perm = c(dims:max(fdasrvf:::ndims(y),dims))
  }
  
  if (fdasrvf:::ndims(y) == 0){
    m = 1
  } else {
    if (length(x) != dim(y)[dims])
      stop('Dimension Mismatch')
    y = aperm(y, perm)
    m = nrow(y)
  }
  
  if (m==1){
    M = length(y)
    out = sum(diff(x)*(y[-M]+y[-1])/2)
  } else {
    slice1 = y[as.vector(outer(1:(m-1), dim(y)[1]*( 1:prod(dim(y)[-1])-1 ), '+')) ]
    dim(slice1) = c(m-1, length(slice1)/(m-1))
    slice2 = y[as.vector(outer(2:m, dim(y)[1]*( 1:prod(dim(y)[-1])-1 ), '+'))]
    dim(slice2) = c(m-1, length(slice2)/(m-1))
    out = t(diff(x)) %*% (slice1+slice2)/2.
    siz = dim(y)
    siz[1] = 1
    out = array(out, siz)
    perm2 = rep(0, length(perm))
    perm2[perm] = 1:length(perm)
    out = aperm(out, perm2)
    ind = which(dim(out) != 1)
    out = array(out, dim(out)[ind])
  }
  
  return(out)
}





#' Compute curve length
#' @param p coordinates of curve
#' @return L length of the curve
#' @export
curve_length <- function(p){
  n = nrow(p)
  T_col = ncol(p)
  L = 0
  sq_tmp = matrix(0,n,T_col)
  del = list()
  for(j in 2:T_col){
    for(i in 1:n){
      sq_tmp[i,j-1] = (p[i,j]-p[i,j-1])^2
    }
    del[[j-1]] = (sum(sq_tmp[,j-1]))^(0.5)
    
    L = L+del[[j-1]]
  }
  return(L)
}

#' Scale, translate and rotate curve p
#' @param p curve
#' @param L scaling factor
#' @param R rotation matrix
#' @param POS translation
#' @return scaled, translated and rotated version of curve p
#' @export
repose_curve <- function(p,L,R,POS){
  n = nrow(p)
  T_col = ncol(p)
  curPos = rowMeans(p)
  p = p/curve_length(p)
  for (i in 1:n){
    p[i,] = p[i,] - curPos[i]
  }
  pnew = R%*%p
  pnew = pnew*L/curve_length(pnew)
  for (i in 1:n){
    pnew[i,] = pnew[i,] + POS[i]
  }
  return(pnew)
}

#' Save the eigen projections of the pca
#' @param covdata list of covariance matrix, eigen values returned from \code{\link{compute_covariance}}
#' @param alpha_t_array list of tangent vectors from the mean shape
#' @param eigendirs the principal eigen directions
#' @param filename exported filename
#' @return csv of eigen projections
#' @export
save_eigen_projections <- function(covdata,alpha_t_array,eigendirs,filename){
  data_to_write = list()
  name = c()
  for(i in 1:length(eigendirs)){
    data_to_write[[i]] = rep(0,length(alpha_t_array))
    name <- c(name,paste("eig",i))
  }
  names(data_to_write) <- name
  for(i in 1:length(alpha_t_array)){
    for(ctr in 1:length(eigendirs)){
      temp_ptb <- project_to_basis(alpha_t_array[i],covdata$Y)
      data_to_write[[ctr]][i] = temp_ptb[[2]]%*%covdata$U[,eigendirs[ctr]]
    }
  }
  write.csv(data_to_write, file = filename,row.names=FALSE)
  return(data_to_write)
}

#' Get the eigen projections of the pca
#' @param covdata list of covariance matrix, eigen values returned from \code{\link{compute_covariance}}
#' @param alpha_t_array list of tangent vectors from the mean shape
#' @param eigendirs the principal eigen directions
#' @return csv of eigen projections
#' @export
get_eigen_projections <- function(covdata,alpha_t_array,eigendirs){
  eigproj = list()
  name = c()
  for(i in 1:length(eigendirs)){
    eigproj[[i]] = rep(0,length(alpha_t_array))
    name <- c(name,sprintf('eig%d',i))    
  }
  names(eigproj) <- name  
  for(i in 1:length(alpha_t_array)){
    for(ctr in 1:length(eigendirs)){
      temp_ptb <- project_to_basis(alpha_t_array[i],covdata$Y)
      eigproj[[ctr]][i] = temp_ptb[[2]]%*%covdata$U[,eigendirs[ctr]]
    }
  }
  return(eigproj)
}

#' Build TPCA model from mean shape
#' @param qmean matrix containing the mean shape
#' @param alpha_t_array list of tangent vectors from the mean shape
#' @return covdata list containing covariance matrix and the coefficients of eigen projections
#' @export
build_tpca_model_from_mean <- function(qmean,alpha_t_array){
  N = nrow(alpha_t_array)
  cols = ncol(alpha_t_array)
  epsilon = 0.0001;
  Y = gram_schmidt(alpha_t_array)
  X_and_X_proj = project_to_basis(alpha_t_array,Y)
  X = X_and_X_proj[[1]]
  X_proj = X_and_X_proj[[2]]
  C = cov(X_proj)
  U = svd(C)$u
  S = pracma::zeros(length(svd(C)$d),length(svd(C)$d))
  diag(S) <- svd(C)$d
  V = svd(C)$v
  PC_Latent_Explained = pcacov(C)
  PC = PC_Latent_Explained[[1]]
  Latent = PC_Latent_Explained[[2]]
  Explained = PC_Latent_Explained[[3]]
  Cn = U%*%(S+epsilon*pracma::eye(ncol(S)))%*%t(U)
  Kinv = solve(Cn)
  covdata <- list(qmean, alpha_t_array, Y, X_proj, U, S, V, C, Cn, Kinv, PC, Latent, Explained)
  names(covdata) <- c('qmean', 'alpha_t_array', 'Y', 'X_proj', 'U', 'S', 'V', 'C', 'Cn', 'Kinv', 'PC', 'Latent', 'Explained')
  
  return(covdata)
}









