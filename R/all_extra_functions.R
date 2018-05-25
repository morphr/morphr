back_parallel_transport_c <- function(wfinal, alpha){
    n = dim(alpha)[1]
    T_col = dim(alpha)[2]
    k = dim(alpha)[3]
    stp = k-1
    wtilde <- rep(0, n*T_col*(stp+1))
    dim(wtilde) <- c(n,T_col,(stp+1))
    wtilde[,,stp+1] = wfinal
    twtilde = wtilde
    for (tau in  seq(stp,1,-1)){
        wtilde[,,tau] = parallel_transport_c(wtilde[,,tau+1], alpha[,,tau])
        twtilde[,,tau] = (tau-1)/stp*wtilde[,,tau]
    }
    return(list(wtilde, twtilde))
}

compute_elastic_geodesic <- function(q1,q2,stp = 7,d = 5,dt = 0.1){
    q2 = project_curve(q2)
    n = length(q1)
    T_col = ncol(q1)
    q2 = regroup(q1,q2)
    q2new = find_rotation_and_seed_unique(q1,q2)
    temp_all = compute_geodesic_c_factor_d(q1,q2new,stp,d,dt)
    return(temp_all)
}

compute_geodesic_c_factor_d <- function(q1,q2,stp,d,dt){
    iter = 1
    n = nrow(q1)
    T_col = ncol(q1)
    V = form_basis_d(d,T_col)
    Ednorm_sq = 100
    gamma_id = pracma::linspace(0,2*pi,T_col)
    s = pracma::linspace(0,2*pi,T_col)
    epsilon =0.1
    EgeoC = list()
    EgeoC[1] = 100
    q2n = regroup(q1,q2)
    q2n = project_curve(q2n)
    diffEgeoC = 100
    diffL = 100
    q2n = initialize_gamma_using_dp(q1,q2)
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
    while(iter < 45 && diffL > 1/50*epsilon){
        q2n = regroup(q1,q2n)
        q2n = project_curve(q2n)
        temp = compute_geodesic_c(q1,q2n,stp,dt)
        alpha = temp[[1]]
        alpha_t = temp[[2]]
        E = unlist(temp[[3]])
        L = unlist(temp[[4]])
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
        gamma_n=gamma_n/max(gamma_n)*2*pi
        if(sum(gamma_n<0)>=1||sum(diff(gamma_n)<0)){
            cat("EXCEPTION: Gamma is INVALID")
            break
        }
        gamma_n = as.vector(gamma_n)
        q2n = group_action_by_gamma(q2n, gamma_n)
        q2n = project_curve(q2n)
        iter = iter+1
        if(iter>2){
            diffEgeoC = pracma::Norm(EgeoC[[iter-1]] - EgeoC[[iter-2]])
            diffL = pracma::Norm(Liter[[iter-1]] - Liter[[iter-2]])
        }
    }
    gamma = estimate_gamma(q2n)
    geo_dist = L
    return(list(alpha,alpha_t,Anormiter,EgeoC,gamma,geo_dist))
}

compute_geodesic_c <- function(q0,q1,stp,dt){
    ##Delete this
    #q0 = q1
    #q1 = q2n
    ##Delete this
    q1 = regroup(q0,q1)
    q1 = project_curve(q1)
    alpha = geodesic_q(q0,q1,stp)
    iter = 1
    E = list()
    E[1] = 1000
    vnorm = 100
    delE = 10
    while(iter<10 && vnorm > 10^-3){
        alpha_t <- dAlpha_dt(alpha)
        L = path_length(alpha_t)
        E[iter] = palais_inner_product(alpha_t, alpha_t)
        w = cov_int_alpha_t(alpha,alpha_t)
        wtilde = back_parallel_transport_c(w[,,stp+1],alpha)[[1]]
        twtilde = back_parallel_transport_c(w[,,stp+1],alpha)[[2]]
        v = w-twtilde
        vnorm = palais_inner_product(v,v)
        alpha = path_update(alpha, -v, 0.5)
        alpha[,,dim(alpha)[3]] = q1
        iter = iter+1
    }
    return(list(alpha, alpha_t, E, L))
}

cov_int_alpha_t <- function(alpha,alpha_t){
    n = dim(alpha_t)[1]
    T_col = dim(alpha_t)[2]
    k = dim(alpha_t)[3]
    stp = k-1
    w <- rep(0, n*T_col*(stp+1))
    dim(w) <- c(n,T_col,(stp+1))
    for(tau in 2:(stp+1)){
        wprev = parallel_transport_c(w[,,tau-1], alpha[,,tau])
        w[,,tau] = wprev+alpha_t[,,tau]/stp
        w[,,tau] = project_tangent(w[,,tau], alpha[,,tau])
    }
    return(w)
}

dAlpha_dt <- function(alpha){
    n = dim(alpha)[1]
    T_col = dim(alpha)[2]
    k = dim(alpha)[3]
    stp = k-1
    alpha_t <- rep(0, n*T_col*(stp+1))
    dim(alpha_t) <- c(n,T_col,(stp+1))
    for (tau in 2:(stp+1)){
        alpha_t[,,tau] = stp*(alpha[,,tau]-alpha[,,tau-1])
        alpha_t[,,tau] = project_tangent(alpha_t[,,tau], alpha[,,tau])
    }
    return(alpha_t)
}

estimate_gamma <- function(q){
    p = q_to_curve(q)
    n = nrow(p)
    T_col = ncol(p)
    pdiff = t(apply( p , 1 , diff ))
    ds = sqrt(apply( pdiff^2 , 2 , sum ))*T_col
    gamma = cumsum(ds)*2*pi/max(cumsum(ds))
    return(gamma)
}

find_rotation_and_seed_unique <- function(q1,q2){
    N = 30
    n = nrow(q1)
    T_col = ncol(q1)
    Ltwo = list()
    Rlist = list()
    for (ctr in 0:T_col){
        q2n = shift_f(q2,ctr)
        q2new = find_best_rotation(q1,q2n)[[1]]
        R = find_best_rotation(q1,q2n)[[2]]
        Ltwo[ctr+1] = innerprod_q(q1-q2new, q1-q2new)
        Rlist[[ctr+1]] = R
    }
    Ltwo_idx = which.min(Ltwo)
    q2new = shift_f(q2, Ltwo_idx-1)
    q2new = Rlist[[Ltwo_idx]]%*%q2new
    return(q2new)
}

find_best_rotation <- function(q1,q2){
    n = nrow(q1)
    T_col = ncol(q1)
    A = q1%*%t(q2)
    U = svd(A)$u
    S = diag(svd(A)$d)
    V = svd(A)$v
    eps = 2.2204*10^(-16)
    if(abs(det(U)*det(V) - 1) < 10*eps){
        S = pracma::eye(n)
    }else{
        S = pracma::eye(n)
        S = -S
    }
    R = U%*%S%*%t(V)
    q2new = R%*%q2
    return(list(q2new, R))
}

find_mean_shape <- function(qarray){
    n = nrow(qarray[[1]])
    T_col = ncol(qarray[[1]])
    N = length(qarray)
    stp = 7
    dt = 0.1
    d = 10
    qmean = pracma::zeros(n,T_col)
    for (i in 1:N){
        qmean = qmean+qarray[[i]]
    }
    qmean = qmean/N
    qmean = project_curve(qmean)
    qshapes = list()
    norm_alpha_t_mean = list()
    sum_sq_dist_iter = list()
    for (iter in 1:6){
        cat('Iteration ', iter, '\n')
        alpha_t_mean = pracma::zeros(n,T_col)
        sum_sq_dist = 0
        qshapes[[1]] = qmean
        for(i in 1:N){
            cat(i, '\n')
            qshapes[[2]] = qarray[[i]]
            temp_all_2 = geodesic_distance_all(qshapes)
            alpha = temp_all_2[[1]]
            alpha_t = temp_all_2[[2]]
            Anormiter = temp_all_2[[3]]
            Egeo = temp_all_2[[4]]
            gamma = temp_all_2[[5]]
            geo_dist = temp_all_2[[6]]
            alpha_t_mean = alpha_t_mean+alpha_t[[1]][,,2]/N
            sum_sq_dist = sum_sq_dist+geo_dist[[1]]^2
        }
        norm_alpha_t_mean[iter] = sqrt(innerprod_q(alpha_t_mean, alpha_t_mean))
        sum_sq_dist_iter[iter] = sum_sq_dist
        qmean = geodesic_flow(qmean, alpha_t_mean, stp)[[1]]
    }
    qshapes[[1]] = qmean
    gamma_array = list()
    alpha_t_array = list()
    alpha_array = list()
    Anorm_array = list()
    Egeo_array = list()
    geo_dist_array = list()
    for(i in 1:N){
        qshapes[[2]] = qarray[[i]]
        temp_all_3 = geodesic_distance_all(qshapes)
        alpha = temp_all_3[[1]]
        alpha_t = temp_all_3[[2]]
        Anormiter = temp_all_3[[3]]
        Egeo = temp_all_3[[4]]
        gamma_array[[i]] = temp_all_3[[5]]
        geo_dist = temp_all_3[[6]]
        alpha_t_array[[i]] = alpha_t[[1]][,,2]
        alpha_array[[i]] = alpha[[1]]
        Anorm_array[[i]] = Anormiter[[1]]
        Egeo_array[[i]] = Egeo[[1]]
        geo_dist_array[[i]] = geo_dist[[1]]
    }
    sum_sq_dist_iter = unlist(sum_sq_dist_iter)
    sum_sq_dist_iter = sum_sq_dist_iter/N
    return(list(qmean=qmean,alpha_array=alpha_array, alpha_t_array=alpha_t_array,norm_alpha_t_mean=norm_alpha_t_mean,gamma_array=gamma_array,sum_sq_dist_iter=sum_sq_dist_iter,Egeo_array=Egeo_array,geo_dist_array=geo_dist_array))
}

form_basis_d_q <- function(V,q){
    d = dim(V)[1]
    n = dim(q)[1]
    T_col = dim(q)[2]
    D_q <- rep(0, n*T_col*d)
    dim(D_q) <- c(n,T_col,d)
    qdiff <- pracma::zeros(n,T_col)
    Vdiff <- pracma::zeros(d,T_col)
    for(i in 1:n){
        qdiff[i,] = pracma::gradient(q[i,],2*pi/(T_col-1))
    }
    for(i in 1:d){
        Vdiff[i,] = pracma::gradient(V[i,],2*pi/(T_col-1))
    }
    for(i in 1:d){
        tmp1 = pracma::repmat(V[i,],n,1)
        tmp2 = pracma::repmat(Vdiff[i,],n,1)
        D_q[,,i] = qdiff*tmp1+1/2*q*tmp2
    }
    return(D_q)
}

form_basis_d2 <- function(d,T_col){
    x = pracma::linspace(0,2*pi, T_col)
    temp = 1:d
    dim(temp) <- c(d,1)
    xdarray = temp%*%x
    V = rbind(cos(xdarray)/sqrt(pi), sin(xdarray)/sqrt(pi))
    return(V)
}

form_basis_normal_a <- function(q){
    n = nrow(q)
    T_col = ncol(q)
    e = pracma::eye(n)
    Ev <- rep(0, n*T_col*n)
    dim(Ev) <- c(n,T_col,n)
    for (i in 1:n){
        Ev[,,i] = pracma::repmat(e[,i],1,T_col)
    }
    qnorm = list()
    for(i in 1:T_col){
        qnorm[i] = pracma::Norm(q[,i])
    }
    qnorm = unlist(qnorm)
    delG = list()
    for(i in 1:n){
        tmp1 = pracma::repmat(q[i,]/qnorm,n,1)
        tmp2 = pracma::repmat(qnorm, n, 1)
        delG[[i]] = tmp1*q+tmp2*Ev[,,i]
    }
    return(delG)
}

geodesic_distance_all <- function(qarray){
    stp = 6
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
        for (j in (i+1):num_data_sets){
            q1 = qarray[[i]]
            q2 = qarray[[j]]
            cat("Computing elastic geodesic with", i, 'and', j, '\n')
            temp_all = compute_elastic_geodesic(q1,q2,stp,d,dt)
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

geodesic_q <- function(q1,q2,stp){
    th = acos(innerprod_q(q1,q2))
    alpha_new <- rep(0, nrow(q1)*ncol(q2)*(stp+1))
    dim(alpha_new) <- c(nrow(q1),ncol(q2),(stp+1))
    for( tau in 0:stp){
        t = tau/stp
        if (th<0.0001){
            alpha_new[,,tau+1] = q1
        }else{
            alpha_new[,,tau+1] = (sin(th - t*th)*q1 + sin(t*th)*q2)/sin(th)
            alpha_new[,,tau+1] = project_curve(alpha_new[,,tau+1])
        }
    }
    return(alpha_new)
}

geodesic_sphere <- function(x_init,g,dt){
    gnorm = sqrt(innerprod_q(g,g))
    X = cos(dt*gnorm)*x_init+sin(dt*gnorm)*g/gnorm
    return(X)
}

group_action_by_gamma <- function(q,gamma){
    n = nrow(q)
    T_col = ncol(q)
    gamma_t = pracma::gradient(gamma,2*pi/(T_col-1))
    q_composed_gamma = matrix(rep(0,n*T_col),nrow=n,ncol=T_col)
    for(i in 1:n){
        q_composed_gamma[i,] = pracma::interp1(pracma::linspace(0,2*pi,T_col), q[i,], gamma, 'nearest')
    }
    sqrt_gamma_t = pracma::repmat(sqrt(gamma_t),n,1)
    qn = q_composed_gamma*sqrt_gamma_t
    return(qn)
}


initialize_gamma_using_dp <- function(q1,q2){
  qarray_temp = list()
  qarray_temp[[1]] = q1
  qarray_temp[[2]] = q2
  save_q_shapes('DPshapedata_R.dat', qarray_temp)
  dp_match_path <- get_dp_shape_match_path()
  temp_str = paste(dp_match_path, ' DPshapedata_R.dat gamma_R.dat', sep = '')
  system(temp_str)
  gamma = load_gamma('gamma_R.dat')
  gamma = gamma/max(gamma)*2*pi
  q2n = group_action_by_gamma(qarray_temp[[2]], gamma)
  return(q2n)
}


save_q_shapes <- function(v_szfilename, qarray_temp){
  # setwd("~/srp99/Research/dinosaur_code")
  N = length(qarray_temp)
  n = as.integer(nrow(qarray_temp[[1]]))
  T_col = as.integer(ncol(qarray_temp[[2]]))
  x <- as.integer(c(N,n,T_col))
  write.filename = file(v_szfilename, "wb")
  writeBin(x, write.filename)
  for(i in 1:N){
    for(j in 1:n){
      writeBin(qarray_temp[[i]][j,], write.filename, size = 4)
    }
  }
  close(write.filename)
  success = 0
  return(success)
}

load_gamma <- function(filename){
  read.filename <- file('gamma_R.dat', 'rb')
  N = readBin(read.filename,"integer",n = 1)
  gamma = readBin(read.filename,"double",n = 100, size = 4)
  close(read.filename)
  return(gamma)
}

palais_inner_product <- function(v1,v2){
    stp = dim(v1)[3]
    s = pracma::linspace(0,1,stp)
    v_inner = list()
    for(i in 1:stp){
        v_inner[i] = innerprod_q(v1[,,i], v2[,,i])
    }
    v_inner = unlist(v_inner)
    val = pracma::trapz(s,v_inner)
    return(val)
}

parallel_transport_c <- function(w,q3){
    lw = sqrt(innerprod_q(w,w))
    if(lw<0.0001){
        w_new = w
    }else{
        w_new = project_tangent(w,q3)
        w_new = w_new*lw/sqrt(innerprod_q(w_new,w_new))
    }
    return(w_new)
}

path_length <- function(alpha_t){
    stp = dim(alpha_t)[3]
    s = pracma::linspace(0,1,stp)
    v_sqrt_inner = list()
    for (i in 1:stp){
        v_sqrt_inner[i] = sqrt(innerprod_q(alpha_t[,,i],alpha_t[,,i]))
    }
    v_sqrt_inner = unlist(v_sqrt_inner)
    L = pracma::trapz(s,v_sqrt_inner)
    return(L)
}

path_update <- function(alpha, v, dt){
    k = dim(alpha)[3]
    n = dim(alpha)[1]
    T_col = dim(alpha)[2]
    stp = k-1
    alpha_new <- rep(0, n*T_col*(stp+1))
    dim(alpha_new) <- c(n,T_col,(stp+1))
    for(tau in 1:(stp+1)){
        if(sqrt(innerprod_q(v[,,tau], v[,,tau]))<0.0001){
            alpha_new[,,tau] = alpha[,,tau]
        }else{
            alpha_new[,,tau] = geodesic_sphere(alpha[,,tau], v[,,tau], dt)
            alpha_new[,,tau] = project_curve(alpha_new[,,tau])
        }
    }
    return(alpha_new)
}

project_tangent <- function(f,q){
    n = nrow(q)
    T_col = ncol(q)
    w = f-innerprod_q(f,q)*q
    e = pracma::eye(n)
    g = form_basis_normal_a(q)
    Evorth = gram_schmidt(g)
    Ev <- rep(0, nrow(Evorth[[1]])*ncol(Evorth[[1]])*n)
    dim(Ev) <- c(nrow(Evorth[[1]]),ncol(Evorth[[1]]),n)
    for(i in 1:n){
        Ev[,,i] = Evorth[[i]]
    }
    sum = 0
    for(i in 1:n){
        sum = sum+innerprod_q(w,Ev[,,i])*Ev[,,i]
    }
    fnew = w-sum
    return(fnew)
}

project_tgt_d_q <- function(u,D_q){
    n = dim(D_q)[1]
    T_col = dim(D_q)[2]
    d = dim(D_q)[3]
    uproj = 0
    a = list()
    for(i in 1:d){
        a[i] = innerprod_q(u,D_q[,,i])
        uproj = uproj+a[[i]]*D_q[,,i]
    }
    return(list(uproj,a))
}

regroup <- function(q1,q2){
    n = length(q1)
    T_col = ncol(q1)
    E = list()
    for (tau in 0:(T_col-1)){
        q2n = shift_f(q2,tau)
        E[tau+1] = innerprod_q(q1-q2n, q1-q2n)
    }
    E = unlist(E)
    idx = which.min(E)
    q2n = shift_f(q2, idx-1)
    return(q2n)
}

shift_f <- function(p,tau){
    n = nrow(p)
    T_col = ncol(p)
    if(tau == 0){
        pn = p
        return(pn)
    }
    if(tau>0){
        if(tau+1>T_col){
            pn = p
            return(p)
        }
        pn = p[,(tau+1):T_col]
        pn = cbind(pn, p[,1:tau])
    }else{
        t = abs(tau)+1
        pn = p[,t:T_col]
        pn = cbind(pn, p[,1:(t-1)])
        return(pn)
    }
}
grad_gamma <- function(gamma,c){
  g = rep(0,100)
  h = c*seq(1:100)
  n = 100
  g[1] = (gamma[2]-gamma[1])/(h[2]-h[1])
  g[100] = (gamma[100]-gamma[99])/(h[100]-h[99])
  g[2:99] = gamma[3:100]-gamma[1:98]/(h[3:100]-h[1:98])
  return(g)
}

