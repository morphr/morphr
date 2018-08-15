#' Form the tangent space basis of the quotient space of Q/SO(n)
#' @param B the basis for L2 (Rn)
#' @param q shape
#' @return list of basis elements of T_q(Q/SO(n))
form_basis_o_q <- function(B,q){
  d = length(B)
  n = nrow(q)
  T_col = ncol(q)
  V = form_basis_d(d,T_col)
  R1 = matrix(c(0, 1, 0, -1, 0, 0, 0, 0, 0),3,3,byrow = TRUE)
  R2 = matrix(c(0, 0, 1,0, 0, 0, -1, 0, 0),3,3,byrow = TRUE)
  R3 = matrix(c(0, 0, 0, 0, 0, 1, 0, -1, 0),3,3,byrow = TRUE)
  G = list()
  G[[1]] = R1 %*% q
  G[[2]] = R2 %*% q
  G[[3]] = R3 %*% q
  qdiff <- matrix(0,n,T_col)
  for (i in 1:n){
    qdiff[i,] = fdasrvf::gradient(q[i,],2*pi/(T_col-1))
  }
  Vdiff <- matrix(0,d,T_col)
  for (i in 1:d){
    Vdiff[i,] = fdasrvf::gradient(V[i,],2*pi/(T_col-1))
  }
  D_q <- list()
  for(i in 1:d){
    tmp1 = pracma::repmat(V[i,],n,1)
    tmp2 = pracma::repmat(Vdiff[i,],n,1)
    D_q[[i]] = qdiff*tmp1 + (1/2)*q*tmp2
  }
  O_q = c(G,D_q)
  return(O_q)
}

#' Project the basis B on the tangent space T_q(C)
#' @param B the basis for L2 (Rn)
#' @param q shape
#' @return list of basis elements of T_q(C)
project_basis_tangent_space_of_c_at_q <- function(B,q){
  Bnew = list()
  for (j in 1:length(B)){
    Bnew[[j]] = B[[j]] - innerprod_q(B[[j]],q)*q
  }
  return(Bnew)
}

#' Compute the geodesic flow at q in the direction of w for an unit time
#' @param q shape
#' @param w tangent vector at q
#' @param stp step size
#' @return list of the shapes along the geodesic flow
#' @export
geodesic_flow <- function(q,w,stp){
  n = nrow(q)
  T_col = ncol(q)
  qt = q
  lw = sqrt(innerprod_q(w,w))
  if(lw<0.001){
    alpha = list()
    alpha[[1]] = q
    return (list(qt,alpha))
  }
  alpha = list()
  alpha[[1]] = q
  for (i in 2:(stp+1)){
    qt = qt+w/stp
    qt = project_curve(qt)
    alpha[[i]] = qt
    w = project_tangent(w,qt)
    w = w*lw/sqrt(innerprod_q(w,w))
  }
  return(list(qt,alpha))
}

#' Smooth a function based on degree of smoothing
#' @param y function to be smoothed
#' @param span degree of smoothing
#' @return return a smoothed function
#' @export
smooth <- function(y, span){
  yy = y
  l = length(y)
  for (i in 1:l){
    if(i < span){
      d = i
    }else{
      d = span
    }
    w = d-1
    p2 = floor(w/2)
    if(i>(l-p2)){
      p2 = l-i
    }
    p1 = w-p2
    yy[i] = sum(y[(i-p1):(i+p2)])/d
  }
  return(yy)
}

#' Project the set of vectors alpha_t_array into a basis set Y
#' @param alpha_t_array set of tangent vectors
#' @param Y a basis set
#' @return list of projected vectors (X) and the coefficients of projections (X_proj)
#' @export
project_to_basis <- function(alpha_t_array,Y){
  n = nrow(Y[[1]])
  T_col = ncol(Y[[1]])
  X = list()
  X_proj = pracma::zeros(length(alpha_t_array),length(Y))
  for(i in 1:length(alpha_t_array)){
    X[[i]] = pracma::zeros(n,T_col)
    for(j in 1:length(Y)){
      X_proj[i,j] = innerprod_q(alpha_t_array[i],Y[[j]])
      X[[i]] = X[[i]]+X_proj[i,j]*Y[[j]]
    } 
  }
  
  return(list(X,X_proj))
}

#' Form the tangent space basis (T_q(S)) of the quotient space  S = C/(D_q x O_q)
#' @param B Orthogonalized basis for L2 (Rn)
#' @param G_O_q The basis for the tangent space T_q(O_q) (using \code{\link{form_basis_o_q}}
#' @return list of basis elements of T_q(S)
form_basis_of_tangent_space_of_s_at_q <- function(B,G_O_q){
  G = list()
  for (j in 1:length(B)){
    tmp = 0
    for (k in 1:length(G_O_q)){
      tmp = tmp+innerprod_q(B[[j]],G_O_q[[k]])*G_O_q[[k]]
    }
    G[[j]] = B[[j]] - tmp
  }
  return (G)
}

innerprod_q<-function(v1, v2){
  
  Tempt = ncol(v1)
  time=pracma::linspace(0,2*pi,Tempt)
  ip <- pracma::trapz(time,pracma::dot(v1, v2))
  return(ip)
}

#' Forms the basis for the tangent space of diffeomorphisms
#' This is actually an L2 [0,1] space
#' @param d number Fourier basis elements, d for sin(kt), cos(kt), 0<=k<=d-1 
#' @param T_col dimension of sample points
#' @return matrix of the basis
form_basis_d <- function(d,T_col){
  x = pracma::linspace(0,2*pi,T_col)
  xdarray = t(t(seq(1:d)))%*%x
  V = rbind(cos(xdarray)/sqrt(pi), sin(xdarray)/sqrt(pi))
  return(V)
}

#' Form the Schauder basis for the tangent space L2(R3)
#' @param d number Fourier basis elements, d for sin(kt), cos(kt), 0<=k<=d-1 
#' @param T_col dimension of sample points
#' @return matrix of the basis elements of length 6*d
form_basis_l2_r3 <- function(d,T_col){
  x = pracma::linspace(0,1,T_col)
  k = 0
  constB = list()
  constB[[1]] <- matrix(0,3,T_col)
  constB[[1]][1,] = sqrt(2)
  constB[[2]] <- matrix(0,3,T_col)
  constB[[2]][2,] = sqrt(2)
  constB[[3]] <- matrix(0,3,T_col)
  constB[[3]][3,] = sqrt(2)
  B = list()
  for (j in 1:d){
    B[[1+6*k]] = rbind(sqrt(2)*cos(2*pi*j*x),rep(0,T_col),rep(0,T_col))
    B[[2+6*k]] = rbind(rep(0,T_col),sqrt(2)*cos(2*pi*j*x),rep(0,T_col))
    B[[3+6*k]] = rbind(rep(0,T_col),rep(0,T_col),sqrt(2)*cos(2*pi*j*x))
    B[[4+6*k]] = rbind(sqrt(2)*sin(2*pi*j*x),rep(0,T_col),rep(0,T_col))
    B[[5+6*k]] = rbind(rep(0,T_col),sqrt(2)*sin(2*pi*j*x),rep(0,T_col))
    B[[6+6*k]] = rbind(rep(0,T_col),rep(0,T_col),sqrt(2)*sin(2*pi*j*x))
    k = k+1
  }
  
  
  return(c(constB,B))
}

#' Orthogonalized the vectors X
#' @param X list of vectors
#' @return orthogonalized vectors
#' @export
gram_schmidt <- function(X){
  epsilon = 0.000005
  rows = 1
  cols = length(X)
  i = 1
  r = 1
  Y = list()
  Y[[1]] <- X[[1]]
  while(i<=cols){
    tempvect = 0
    for(j in 1:(i-1)){
      if(i == 1){
        break
      }
      tempvect =  tempvect + innerprod_q(Y[[j]],X[[r]])*Y[[j]]
    }
    Y[[i]] = X[[r]]-tempvect
    tmp = innerprod_q(Y[[i]],Y[[i]])
    if(tmp > epsilon){
      Y[[i]] = Y[[i]]/sqrt(tmp)
      i = i+1
      r = r+1
    }else{
      if(r<i){
        r = r+1
      }else{
        break
      }
    }
  }
  return(Y)
}

#' Project the shape q to the space of closed curves, C
#' @param q shape
#' @return qnew, projected shape on C
#' @export
project_curve <- function(q){
  T1 = ncol(q)
  n = nrow(q)

  dt = 0.3

  epsilon = 1/60*2*pi/T1;
  
  e = diag(1,n)
  iter = 1
  res = rep(1,n)
  J = matrix(0,n,n)
  s = seq(0,2*pi,length.out=T1)
  qnorm = rep(0,T1)
  G = rep(0,n)
  C = rep(0,301)
  
  qnew = q
  qnew = qnew / sqrt(innerprod_q(qnew,qnew))
  
  while (fdasrvf:::pvecnorm(res) > epsilon){
    if (iter > 301){
      break
    }
    
    # Compute Jacobian
    for (i in 1:n){
      for (j in 1:n){
        J[i,j]  = 3 * pracma::trapz(s, qnew[i,]*qnew[j,])
      }
    }
    J = J + diag(1,n)
    
    for (i in 1:T1){
      qnorm[i] = fdasrvf:::pvecnorm(qnew[,i])
    }
    
    # Compute the residue
    for (i in 1:n){
      G[i] = pracma::trapz(s,qnew[i,]*qnorm)
    }
    
    res = -1*G
    

    
    
    
    C[iter] = fdasrvf:::pvecnorm(res)
    if(any(is.nan(J))||(any(is.infinite(J)))){
      cat('Projection may be inaccurate. \n')
      qnew = q
      qnew = qnew/sqrt(innerprod_q(qnew,qnew))
      break
      
    }
    if (fdasrvf:::pvecnorm(res)<epsilon)
      break
    x = solve(J,res)
    
    delG = form_basis_normal_a(qnew)
    tmp = 0
    for (i in 1:n){
      tmp = tmp + x[i]*delG[[i]]*dt
    }
    qnew = qnew + tmp
    
    iter = iter + 1
  }
  
  qnew = qnew / sqrt(innerprod_q(qnew,qnew))
  return(qnew)
}

#' Convert q to curve
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

#' Resample curve to get N points
#' @param p vector of coordinates of curves
#' @param N number of points desired
#' @return numeric vector 
#' @export
resample_curve <- function(p,N){
  n = nrow(p)
  T_col = ncol(p)
  pdiff = t(as.matrix(diff(t(p))))
  normdiff = c()
  for(i in 1:(T_col-1)){
    normdiff = c(normdiff,norm(as.matrix(pdiff[,i])))
  }
  normdiff[normdiff > 0] = 1
  pu = list()
  for(i in 1:n){
    pu[[i]] = p[i,pracma::finds(normdiff)]
  }
  pu = rbind(pu[[1]],pu[[2]])
  n = nrow(pu)
  T_col = ncol(pu)
  pgrad = matrix(rep(0,n*T_col),n,T_col)
  for(i in 1:n){
    pgrad[i,] = pracma::gradient(pu[i,])
  }
  ds = c()
  for(i in 1:T_col){
    ds = c(ds,pracma::Norm(pgrad[,i]))
  }
  S = cumsum(ds)
  newS = pracma::linspace(S[1], S[length(S)], N)
  pnew = matrix(rep(0,n,N),n,N)
  for(i in 1:n){
    pnew[i,] = pracma::interp1(S,pu[i,],newS)
  }
  
  return(pnew)
}
