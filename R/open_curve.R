#' @export
estimate_gamma_open <- function(q){
  p = morphr::q_to_curve(q)
  s = pracma::linspace(0,2*pi,ncol(q))
  pgrad = pracma::gradient(p,2*pi/ncol(q))$X
  ds = sqrt(colSums(pgrad^2))*ncol(q)
  gamma = pracma::cumtrapz(s,ds)*2*pi/max(pracma::cumtrapz(s,ds))
  return(as.vector(gamma))
}

#' @export
project_tangent_open <- function(f,q){
  return(f-q*innerprod_q(f,q))
}

#' @export
dAlpha_dt_open <- function(alpha){
  n = dim(alpha)[1]
  T_col = dim(alpha)[2]
  k = dim(alpha)[3]
  stp = k-1
  alpha_t <- rep(0, n*T_col*(stp+1))
  dim(alpha_t) <- c(n,T_col,(stp+1))
  for (tau in 2:(stp+1)){
    alpha_t[,,tau] = stp*(alpha[,,tau]-alpha[,,tau-1])
    alpha_t[,,tau] = project_tangent_open(alpha_t[,,tau], alpha[,,tau])
  }
  return(alpha_t)
}

#' @export
initialize_gamma_using_dp_open <- function(q1,q2){
  qarray_temp = list()
  qarray_temp[[1]] = q1
  qarray_temp[[2]] = q2
  save_q_shapes('DPshapedata_R.dat', qarray_temp)
  dp_match_path <- get_dp_shape_match_path()
  temp_str = paste(dp_match_path, ' DPshapedata_R.dat gamma_R.dat', sep = '')
  system(temp_str)
  gamma = load_gamma('gamma_R.dat')
  gamma = gamma/max(gamma)*2*pi
  q2n = group_action_by_gamma_cc(qarray_temp[[2]], gamma)
  return(q2n)
}

#' @export
group_action_by_gamma_cc <- function(q,gamma){
  gamma_t = pracma::gradient(as.vector(gamma),2*pi/ncol(q))
  q_composed_gamma = matrix(rep(0,nrow(q)*ncol(q)),nrow(q),ncol(q))
  for(i in 1:nrow(q)){
    q_composed_gamma[i,] = pracma::interp1(pracma::linspace(0,2*pi,ncol(q)),q[i,],gamma,method = "linear")
  }
  if(nrow(q)==2){
    sqrt_gamma_t = rbind(sqrt(gamma_t),sqrt(gamma_t))
  }else if(nrow(q)==3){
    sqrt_gamma_t = rbind(sqrt(gamma_t),sqrt(gamma_t),sqrt(gamma_t))
  }
  q2n = q_composed_gamma*sqrt_gamma_t
  return(q2n)
}

#' @export
curve_to_q_open <- function(p){
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
  q <- project_b(q)
  return(q)
}

#' @export
project_b <- function(q){
  q = q/(sqrt(pracma::trapz(pracma::linspace(0,2*pi,ncol(q)),pracma::dot(q,q))))
  return(q)
}

#' @export
compute_elastic_geodesic_open <- function(q1,q2,stp = 7,d = 5,dt = 0.1){
  q2 = project_b(q2)
  n = length(q1)
  T_col = ncol(q1)
  q2new = find_best_rotation(q1,q2)[[1]]
  temp_all = compute_geodesic_c_factor_d_open(q1,q2new,stp,d,dt)
  return(temp_all)
}

#' @export
compute_geodesic_b <- function(q1,q2,stp){
  th = acos(morphr::innerprod_q(q1,q2))
  f = q2-morphr::innerprod_q(q1,q2)*q1
  f = th*f/sqrt(morphr::innerprod_q(f,f))
  alpha = list()
  for(tau in 0:stp){
    alpha[[tau+1]] = morphr::geodesic_sphere(q1,f,tau/stp)
    alpha[[tau+1]] = project_b(alpha[[tau+1]])
  }
  alpha_m = unlist(alpha)
  dim(alpha_m) = c(dim(alpha[[1]]),stp+1)
  alpha_t = dAlpha_dt_open(alpha_m)
  E = morphr::palais_inner_product(alpha_t,alpha_t)
  L = morphr::path_length(alpha_t)
  return(list(alpha = alpha, alpha_t = alpha_t, E = E, L = L))
}

#' @export
compute_geodesic_c_factor_d_open <- function(q1,q2,stp,d,dt){
  iter = 1
  n = nrow(q1)
  T_col = ncol(q1)
  V = form_basis_d(d,T_col)
  Ednorm_sq = 100
  gamma_id = pracma::linspace(0,2*pi,T_col)
  s = pracma::linspace(0,2*pi,T_col)
  epsilon =0.000001
  EgeoC = list()
  EgeoC[1] = 100
  q1 = project_b(q1)
  q2n = project_b(q2)
  q2n = initialize_gamma_using_dp_open(q1,q2n)
  q2n = project_b(q2n)
  compute_geodesic_b_temp = compute_geodesic_b(q1,q2n,stp)
  alpha = compute_geodesic_b_temp$alpha
  alpha_t = compute_geodesic_b_temp$alpha_t
  E = compute_geodesic_b_temp$E
  L = compute_geodesic_b_temp$L
  Anormiter = 0
  geo_dist = L
  EgeoC = E
  diffEgeoC = 100
  diffL = 100

  Liter = list()
  Ednorm_sq_iter = list()
  Anormiter = list()
  if(innerprod_q(q1-q2n,q1-q2n)<1/15*epsilon){
    Anormiter = 0
    EgeoC = 0
    gamma = s
    geo_dist = 0
    alpha = list()
    for(i in 1:stp){
      alpha[[i]] = q1
    }
    alpha_t = array(0,dim=c(n,T_col,stp+1))
    return(list(alpha,alpha_t,Anormiter,EgeoC,gamma,geo_dist))
  }
  while(iter < 45 && diffL > 1/100*epsilon){
    q2n = project_b(q2n)
    compute_geodesic_b_temp = compute_geodesic_b(q1,q2n,stp)
    alpha = compute_geodesic_b_temp$alpha
    alpha_t = compute_geodesic_b_temp$alpha_t
    E = compute_geodesic_b_temp$E
    L = compute_geodesic_b_temp$L
    Liter[iter] = L
    EgeoC[iter] = E[length(E)]
    u = alpha_t[,,stp+1]
    D_q = form_basis_d_q(V,q2n)
    temp2 = project_tgt_d_q(u,D_q)
    uproj = temp2[[1]]
    a = temp2[[2]]
    a = unlist(a)
    Ednorm_sq = innerprod_q(uproj, uproj)
    Ednorm_sq_iter[iter] = Ednorm_sq
    g = a%*%V
    Anormiter[iter] = pracma::trapz(s,g^2)
    gamma_n = s-epsilon*g
    gamma_n = gamma_n - gamma_n[1]
    if(sum(gamma_n<0)>=1||sum(diff(gamma_n)<0)){
      cat("EXCEPTION: Gamma is INVALID")
      break
    }
    gamma_n = as.vector(gamma_n)
    q2n = group_action_by_gamma_cc(q2n, gamma_n)
    q2n = project_b(q2n)
    iter = iter+1
    if(iter>2){
      diffEgeoC = pracma::Norm(EgeoC[[iter-1]] - EgeoC[[iter-2]])
      diffL = pracma::Norm(Liter[[iter-1]] - Liter[[iter-2]])
    }
  }
  gamma = estimate_gamma_open(q2n)
  geo_dist = L
  return(list(alpha,alpha_t,Anormiter,EgeoC,gamma,geo_dist))
}

#' @export
geodesic_distance_all_open <- function(qarray){
  stp = 7
  dt = 0.1
  d = 5
  num_data_sets = length(qarray)
  ctr = 1
  alpha = list()
  alpha_t = list()
  Anormiter = list()
  EgeoC = list()
  gamma = list()
  geo_dist = list()
  cat('Total of ', num_data_sets,' datasets to compute.\n')
  for (i in 1:(num_data_sets-1)){
    for (j in (i+1):(num_data_sets)){
      q1 = qarray[[i]]
      q2 = qarray[[j]]
      cat("Computing elastic geodesic with", i, 'and', j, '\n')
      temp_all = compute_elastic_geodesic_open(q1,q2,stp,d,dt)
      alpha[[ctr]] = temp_all[[1]]
      alpha_t[[ctr]] = temp_all[[2]]
      Anormiter[[ctr]] = temp_all[[3]]
      EgeoC[[ctr]] = temp_all[[4]]
      gamma[[ctr]] = temp_all[[5]]
      geo_dist[[ctr]] = temp_all[[6]]
      ctr = ctr+1
    }
  }
  return(list(alpha = alpha, alpha_t = alpha_t, Anormiter = Anormiter, EgeoC = EgeoC, gamma = gamma, geo_dist = geo_dist))
}

#' @export
Resample_Curve_general_New_Preserve_Original <- function(p,N){
  n = nrow(p)
  T_col = ncol(p)
  idx = list()
  for(i in 1:T_col){
    G = p-pracma::repmat(as.matrix(p[,i]),1,T_col)
    idx[[i]] = pracma::finds(colSums(G^2)==0)
  }
  idx_rem = c()
  for(i in 1:length(idx)){
    if(length(idx[[i]])==1){
      next
    }
    idx_rem = c(idx_rem,idx[[i]][2:length(idx[[i]])])
  }
  p = p[,setdiff(c(1:T_col),idx_rem)]
  p = remove_trace_over_defects(p)
  
  
  n = nrow(p)
  T_col = ncol(p)
  pdiff = matrix(rep(0,n*T_col),n,T_col)
  for(i in 1:n){
    pdiff[i,] = pracma::gradient(p[i,])
  }
  normpdiff = c()
  for(i in 1:(T_col-1)){
    normpdiff[i] = pracma::Norm(pdiff[,i])
  }
  normpdiff[normpdiff>0] = 1
  pu = list()
  for(i in 1:n){
    pu[[i]] = p[i,pracma::finds(normpdiff)]
  }
  pu = matrix(unlist(pu),length(pu),length(pu[[1]]),byrow = TRUE)
  n = nrow(pu)
  T_col = ncol(pu)
  pgrad = matrix(rep(0,n*T_col),n,T_col)
  for(i in 1:n){
    pgrad[i,] = pracma::gradient(pu[i,])
  }
  ds = c()
  for(i in 1:T_col){
    ds[i] = pracma::Norm(pgrad[,i])
  }
  S = cumsum(ds)
  oldS = pracma::linspace(S[1],S[length(S)],length(S))
  newS = pracma::linspace(S[1],S[length(S)],N)
  newS1 = pracma::interp1(S,oldS,newS)
  pnew = matrix(rep(0,n*N),n,N)
  for(i in 1:n){
    pnew[i,] = pracma::interp1(S,pu[i,],newS1)
  }
  return(pnew)
  
  
  
}

#' @export
remove_trace_over_defects <- function(p){
  n = nrow(p)
  T_col = ncol(p)
  pdiff = matrix(rep(0,n*T_col),n,T_col)
  for(i in 1:n){
    pdiff[i,] = pracma::gradient(p[i,])
  }
  normpdiff = c()
  for(i in 1:T_col){
    normpdiff[i] = pracma::Norm(pdiff[,i])
  }
  idx = pracma::finds(normpdiff==0)
  if(length(idx) == 0){
    status = 0
    pnew = p
    return(p)
  }else{
    print("Check")
    idx_to_remove = c()
    for(ctr in 1:length(idx)){
      idx_to_remove = c(idx_to_remove,idx[ctr])
      N = min(idx[ctr]-1,T_col-idx[ctr]-1)
      for(i in 1:N){
        if(normpdiff[idx[ctr]+i] == normpdiff[idx[ctr]-i]){
          idx_to_remove = c(idx_to_remove,idx[ctr]+i,idx[ctr]-i)
        }
      }
    }
    idx_to_remove = sort(idx_to_remove)
    idx_to_keep = setdiff(c(1:Tcol),idx_to_remove)
    pnew = p[,idx_to_keep]
    pnew = remove_duplicate_points(pnew)
    return(pnew)
  }
}

#' @export
remove_duplicate_points <- function(pnew){
  n = nrow(p)
  T_col = ncol(p)
  idx = list()
  for(i in 1:T_col){
    G = p-pracma::repmat(as.matrix(p[,i]),1,T_col)
    idx[[i]] = pracma::finds(colSums(G^2)==0)
  }
  idx_rem = c()
  for(i in 1:length(idx)){
    if(length(idx[[i]])==1){
      next
    }
    idx_rem = c(idx_rem,idx[[i]][2:length(idx[[i]])])
  }
  p = p[,setdiff(c(1:T_col),idx_rem)]
  idx = setdiff(c(1:T_col),idx_rem)
  return(p)
}

