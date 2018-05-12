#' Plot the curve p with color value
#' @param p coordinates of curve
#' @param colorstr color value
#' @param l Default value is False, which means returning a plot
#' @export
plot_curve <- function(p, colorstr, l = FALSE, filename = ''){
  colorval <- NULL
  if(colorstr == 'k'){
    colorval = "black"
  }else if(colorstr == 'b'){
    colorval = "blue"
  }else if(colorstr == 'r'){
    colorval = "red"
  }else if(colorstr == 'm'){
    colorval = "magenta"
  }else if(colorstr == 'g'){
    colorval = "green"
  }else{
    colorval = "black"
  }
  styleval <- '-'
  n = nrow(p)
  T_col = ncol(p)
  if(n == 2){
    if (l == TRUE){
      lines(p[1,],p[2,],type = "l",col = colorval,lwd=1)
    }else{
      plot(p[1,],p[2,],type = "l",col = colorval,lwd=1,axes=FALSE, xlab = '', ylab = '', main = filename )
    }
  }
  if(n == 3){
    if (l == TRUE){
      lines(p[1,],p[2,],type = "l",col = colorval,lwd=1)
    }else{
      plot(p[1,],p[2,],type = "l",col = colorval,lwd=1,axes=FALSE, xlab = '', ylab = '', main = filename )
    }
    for (i in 1:T_col){
      text(p[1,i],p[2,i],toString(i),cex = 1)
    }
    text(p[1,1],p[2,1],toString(0))
  }
  
}

#' Utility function for shiny app
#' @param maxplots maximum number of plots
#' @param input_n number of iterations
#' @param X shape
#' @return plots
get_plot_output_list <- function(max_plots, input_n, X, filenames) {
  # Insert plot output objects the list
  
  
  
  plot_output_list <- lapply(1:input_n, function(i) {
    plotname <- filenames[i]
    plot_output_object <- plotOutput(plotname, height = 280, width = 250)
    plot_output_object <- renderPlot({
      plot(X[[i]][1,],X[[i]][2,], type = "l",axes=FALSE,xlab = '', ylab = '', main = plotname)
      
    },height = 200, width = 300)
  })
  
  do.call(tagList, plot_output_list) # needed to display properly.
  
  return(plot_output_list)
}


plot_eigenshapes <- function(qmean, alpha_t_array, covdata, eigdir, dt = 0.15) {
  
   X <- get_eigenshapes(qmean, alpha_t_array, covdata, eigdir, dt)
   plot_curve_list(X$X, color_code = X$color_code, titletxt = paste("Eigen shape variation along ",eigdir)) 
}

get_eigenshapes <- function(qmean, alpha_t_array, covdata, eigdir, dt = 0.15) {

  U = covdata$U
  S = covdata$S
  C = covdata$Cn
  Y = covdata$Y
  N = nrow(U)
  zerovect_new <- matrix(0,1,40)
  gphi=zerovect_new
  gth=zerovect_new
  n = nrow(qmean)
  T_col = ncol(qmean)
  stp = 8
  t = eigdir
  Unew = U[,t]
  sigma = diag(S)
  theta = pi/2
  if(n==2){
    R = pracma::eye(2)
  }else if(n==3){
    R = pracma::eye(2)
  }
  
  
  dx = 100.05
  dy = 100.05
  stp = 8
  ctr = 1
  N = 10
  delta = 2
  epsilon_array = list()
  epsilon_array[1] = 0
  Xtemp = list()
  for(epsilon in seq(-1*delta,delta,2*delta/N)){
    epsilon_array[ctr] = epsilon
    gphi = epsilon*sqrt(sigma[eigdir])*t(U[,eigdir])
    final_alpha_t = matrix(0,n,T_col)
    for (i in 1:length(Y)){
      final_alpha_t = final_alpha_t +gphi[i] * Y[[i]]
    }
    temp_temp = geodesic_flow(qmean,final_alpha_t,stp)
    qfinal = temp_temp[1]
    qfinal = qfinal[[1]]
    alpha_final = temp_temp[2]
    if(is.null(ncol(qfinal))){
      
      Xtemp[[ctr]] = matrix(0,2,100)
      ctr = ctr+1
    }else{
      Xtemp[[ctr]] = q_to_curve(R%*%qfinal)
      ctr = ctr+1
    }
    
  }

  X <- list()
  color_code <- c()
  for(i in 1:length(Xtemp)){
    X[[i]] <- repose_curve(Xtemp[[i]],1,R,c(0,i*dt,i*dt))
    color_code[i] <- "blue"
    if(epsilon_array[[i]] == 0){
      color_code[i] <- "red"
    }
  }
  return(list(X=X, color_code=color_code))
}




display_eigen_and_original_shapes_along_extremes <- function(X,qmean,alpha_t_array,covdata,eigdir,fignum){
  U = covdata$U
  S = covdata$S
  C = covdata$Cn
  Y = covdata$Y
  N = nrow(U)
  zerovect_new <- matrix(0,1,40)
  gphi=zerovect_new
  gth=zerovect_new
  n = nrow(qmean)
  T_col = ncol(qmean)
  stp = 8
  t = eigdir
  Unew = U[,t]
  sigma = diag(S)
  theta = pi/2
  if(n==2){
    R = pracma::eye(2)
  }else if(n==3){
    R = pracma::eye(2)
  }
  
  
  dx = 100.05
  dy = 100.05
  stp = 8
  ctr = 1
  N = 10
  delta = 2
  epsilon_array = list()
  epsilon_array[1] = 0
  Xtemp = list()
  for(epsilon in seq(-1*delta,delta,2*delta/N)){
    epsilon_array[ctr] = epsilon
    gphi = epsilon*sqrt(sigma[eigdir])*t(U[,eigdir])
    final_alpha_t = matrix(0,n,T_col)
    for (i in 1:length(Y)){
      final_alpha_t = final_alpha_t +gphi[i] * Y[[i]]
    }
    temp_temp = geodesic_flow(qmean,final_alpha_t,stp)
    qfinal = temp_temp[1]
    qfinal = qfinal[[1]]
    alpha_final = temp_temp[2]
    if(is.null(ncol(qfinal))){
      
      Xtemp[[ctr]] = matrix(0,2,100)
      ctr = ctr+1
    }else{
      Xtemp[[ctr]] = q_to_curve(R%*%qfinal)
      ctr = ctr+1
    }
    
  }
  dt = 0.15
  plot(c(-1,1),c(0,2), type="n", main = paste("Shape variation along Eigen-axis",eigdir)) 
  for(i in 1:length(Xtemp)){
    pfinal = repose_curve(Xtemp[[i]],1,R,c(0,i*dt,i*dt))
    if(epsilon_array[[i]] == 0){
      plot_curve(pfinal,'r',l = TRUE)
    }else{
      plot_curve(pfinal,'b',l = TRUE)
    }
  }
  
  
}


#' Display variation of pca
#' @param pathname location of the file
#' @param alpha_t_array list of tangent vectors from the mean shape
#' @export
plot_pca_variation <- function(alpha_t_array, qmean, qarray, eigdir, dt = 0.15){
  cat('Building TPCA model ...')
  covdata = build_tpca_model_from_mean(qmean,alpha_t_array)
  cat('done.\n')
  plot_eigenshapes(qmean, alpha_t_array, covdata, eigdir, dt = dt)
  # display_eigen_and_original_shapes_along_extremes(X,qmean,alpha_t_array,covdata,1,8)
}

plot_curve_list <- function(X,  color_code = c(), vertical = T, titletxt="") {

  n <- nrow(X[[1]])
  N <- length(X)
  if (length(color_code) == 0){
    color_code = rep('blue', N)
  }
  
  if(n==2){
    R = pracma::eye(2)
  }else if(n==3){
    R = pracma::eye(2)
  }
  
  ii <- 1
  df <- data.frame(x = X[[ii]][1, ], y = X[[ii]][2, ])
  p <- ggplot2::ggplot(df, ggplot2::aes(x , y) ) + ggplot2::geom_path(size = 1, color = color_code[ii]) 
  for (ii in 2:N) {
    df <- data.frame(x = X[[ii]][1, ], y = X[[ii]][2, ])
    p <- p + ggplot2::geom_path(data = df, color = color_code[ii], size=1)
  }
  p <- p + ggplot2::coord_fixed() + ggplot2::scale_y_reverse()+  ggplot_axis_off() + 
    ggplot2::labs(title=titletxt) + ggplot2::theme(title = ggplot2::element_text(size = ggplot2::rel(1))) 
  
  return(p)
}


mdsplot <- function(alpha_t_array, geo_dis, X, color_code_path_name, titletxt ="MDS scatter plot"){
  
  
  Taxoncolorcodes <- readr::read_csv(color_code_path_name,col_names = TRUE)
  filenames <- Taxoncolorcodes$specimen_ID
  color_code <- Taxoncolorcodes$color
  geo_dis <- as.vector(unlist(geo_dis))
  distance_matrix <- pracma::squareform(geo_dis)
  example_NMDS=vegan::metaMDS(distance_matrix,k=2,trymax=200)
  Y = example_NMDS$points
  theta = pi
  R = pracma::eye(2)
  class_color_code <- data.frame(label=levels(factor(Taxoncolorcodes$class)), 
                                 colors = RColorBrewer::brewer.pal(nlevels(factor(Taxoncolorcodes$class)), name = "Set1"))
  class_color_map <- match(Taxoncolorcodes$class, class_color_code$label)
  L = 0.08
  i = 1
  xpos = Y[i,1]
  ypos = -Y[i,2]
  zpos = 0
  pfinal = repose_curve(q_to_curve(curve_to_q(X[[i]])),0.12,R,c(xpos,ypos,0))
  df <- data.frame(x = pfinal[1, ], y = pfinal[2,])
  df2 <- data.frame(labelname = filenames[i], x = xpos, y = ypos)
  p <- ggplot2::ggplot(df, ggplot2::aes(x , y) ) + ggplot2::geom_path(data = df, color = as.character(class_color_code$colors[class_color_map[i]]), size=1) + 
    ggplot2::annotate("text", x = xpos, y = ypos, label = filenames[i], size = 2, fontface = 2)
  for(i in 2:ncol(distance_matrix)){
    xpos = Y[i,1]
    ypos = -Y[i,2]
    zpos = 0
    pfinal = repose_curve(q_to_curve(curve_to_q(X[[i]])),0.12,R,c(xpos,ypos,0))
    df <- data.frame(x = pfinal[1, ], y = pfinal[2,])
    df2 <- data.frame(labelname = filenames[i], x = xpos, y = ypos)
    p <- p+ggplot2::geom_path(data = df, color = as.character(class_color_code$colors[class_color_map[i]]), size=1) + 
      ggplot2::annotate("text", x = xpos, y = ypos, label = filenames[i], size = 2, fontface = 2)
  }
  p <- p + ggplot2::coord_fixed() + ggplot2::scale_y_reverse() +
    ggplot2::labs(x = 'X', y = 'Y') + 
    ggplot2::ggtitle(titletxt) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot_set_axis()
  return(p)
}

#Modification
deformation_field_per_class <- function(pathname){
  data <- R.matlab::readMat(pathname)
  for (ii in 1:4) {
    X1 <-data$pmean.array[[ii]][[1]]
    # alpha_t <- data$alpha.t.array[[1]][[1]]
    # alpha_t_norm <- apply(alpha_t, 2, function(x) sqrt(sum(x^2)))
    
    alpha_t_mag_smooth = as.vector(data$alpha.t.mag.smooth.array[[ii]][[1]])
    # alpha_t_mag = as.vector(data$alpha.t.mag)
    dat = data.frame(x = X1[1,], y = X1[2,], alpha_t_mag_smooth = alpha_t_mag_smooth)
    
    jet_color <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    
    class_id <- as.character(data$class.id[[ii]][[1]])
    p <- ggplot2::ggplot(dat, ggplot2::aes(x , y, color = alpha_t_mag_smooth)) + ggplot2::geom_path(size = 7) + 
      ggplot2::scale_y_reverse() + ggplot2::scale_color_gradientn(colors =jet_color(256)) + ggplot2::theme(legend.text=ggplot2::element_text(size=14)) +
      ggplot2::labs(color="Deformation") + ggplot2::labs(title="Humeri_Including_Indeterminate2") + ggplot2::theme(title = ggplot2::element_text(size = ggplot2::rel(1.8))) +
      ggplot2::xlab('') + ggplot2::ylab('') + ggplot2::coord_fixed()
    
    
    #ggsave(sprintf('protoceratops_class_%s_deform.pdf', class_id), width = 11, height = 8.5)
    print(p)
    
  }
  
  
  
  
  
}

#' Display shape orientation
#' @param pathname location of the ucf file containing hape
#' @export
verify_shapes <- function(pathname){
  fid = read.table(pathname,stringsAsFactors = FALSE)[[1]]
  i = 1
  for (i in 1:length(fid)){
    fname = fid[i]
    X = read_ucf_multiple_levels(fname)[[1]]
    X = t(X)
    plot_curve(X,'r',filename = fname)

  }
}


#' Display Rshiny
#' @export
run_shiny <- function() {
  appDir <- system.file("shinyapps", "myapp", package = "morphr")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `morphr`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}




deformation_field_all <- function(alpha_t_array, pmean, qmean,X){
  alpha_t_mean = 0
  for (jj in 1:length(alpha_t_array)){
    alpha_t_mean = alpha_t_mean + alpha_t_array[[jj]]
  }
  alpha_t_mean = alpha_t_mean/ncol(alpha_t_mean)
  alpha_t_mag = c()
  for (i in 1:ncol(alpha_t_mean)){
    alpha_t_mag[i] = pracma::Norm(alpha_t_mean[,i])
  }
  alpha_t_mag = (alpha_t_mag - min(alpha_t_mag))/(max(alpha_t_mag) - min(alpha_t_mag))
  alpha_t_mag_smooth = smooth(alpha_t_mag, 3)
  alpha_t_mag_smooth = (alpha_t_mag_smooth - min(alpha_t_mag_smooth))/(max(alpha_t_mag_smooth) - min(alpha_t_mag_smooth))
  plot(c(0,100),c(0,1), type="n")
  x_coord = c(1:100)
  lines(x_coord,alpha_t_mag_smooth,type = "l",col = "blue")
  lines(x_coord,alpha_t_mag,type = "l",col = "red")
  pmean = q_to_curve(qmean)
  x = pmean[1,]
  y = pmean[2,]
  z = rep(0,length(x))
  X1 <- X[[1]]
  X1 <- pmean
  #alpha_t_norm <- apply(alpha_t_array, 2, function(x) sqrt(sum(x^2)))
  alpha_t_mag_smooth = as.vector(alpha_t_mag_smooth)
  alpha_t_mag = as.vector(alpha_t_mag)
  dat = data.frame(x = X1[1,], y = X1[2,], alpha_t_mag = alpha_t_mag, alpha_t_mag_smooth = alpha_t_mag_smooth)
  
  jet_color <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  
  p <- ggplot2::ggplot(dat, ggplot2::aes(x , y, color = alpha_t_mag_smooth)) + ggplot2::geom_path(size = 7) + 
    ggplot2::scale_y_reverse() + ggplot2::scale_color_gradientn(colors =jet_color(256)) + ggplot2::theme(legend.text=ggplot2::element_text(size=14)) +
    ggplot2::labs(color="Deformation") + ggplot2::labs(title="Humeri_Including_Indeterminate2") + ggplot2::theme(title = ggplot2::element_text(size = ggplot2::rel(1.8))) +
    ggplot2::xlab('') + ggplot2::ylab('') + ggplot2::coord_fixed()
  
  print(p)
  
}


plot_dendrogram <- function(geo_dist, color_code_path_name){
  geo_dist = unlist(geo_dist)
  Taxoncolorcodes <- readr::read_csv(color_code_path_name,col_names = TRUE)
  filenames <- Taxoncolorcodes$specimen_ID
  color_code <- Taxoncolorcodes$color
  D.clust <- hclust(as.dist(pracma::squareform(geo_dist)),method="average")
  D.clust$labels <- filenames
  plot(D.clust, xlab = '', main='Dwarf Morphometric Placement', sub='', ylab='')
}

plotpca <- function(qmean, alpha_t_array, qarray, colorpath, eigdir = c(1, 2), titletxt = "") {
  Taxoncolorcodes <- readr::read_csv(colorpath,col_names = TRUE)
  filenames <- Taxoncolorcodes$specimen_ID
  color_code <- Taxoncolorcodes$color
  covdata = build_tpca_model_from_mean(qmean,alpha_t_array)
  eigproj_landmarks <- get_eigen_projections(covdata,alpha_t_array,1:5)
  R = pracma::eye(2)

  class_color_code <- data.frame(label=levels(factor(Taxoncolorcodes$class)), 
                                 colors = RColorBrewer::brewer.pal(nlevels(factor(Taxoncolorcodes$class)), name = "Set1"))
  class_color_map <- match(Taxoncolorcodes$class, class_color_code$label)
  
  ii <- 1
  # color_code[ii] <- "red"
  xpos = eigproj_landmarks[[ eigdir[1] ]][ii]
  ypos = eigproj_landmarks[[ eigdir[2] ]][ii]
  pfinal = repose_curve(q_to_curve(qarray[[ii]]),0.12,R,c(xpos,ypos,0))
  df <- data.frame(x = pfinal[1, ], y = pfinal[2,])
  df2 <- data.frame(labelname = filenames[ii], x = xpos, y = ypos)
  p <- ggplot2::ggplot(df, ggplot2::aes(x , y) ) + ggplot2::geom_path(size = 1, color = as.character(class_color_code$colors[class_color_map[ii]])) +
    ggplot2::annotate("text", x = xpos, y = ypos, label = filenames[ii], size = 3, fontface = 2)

  for (ii in 2:length(eigproj_landmarks$eig1)) {
    xpos = eigproj_landmarks[[ eigdir[1] ]][ii]
    ypos = eigproj_landmarks[[ eigdir[2] ]][ii]
    # color_code[ii] <- "red"
    pfinal = repose_curve(q_to_curve(qarray[[ii]]),0.12,R,c(xpos,ypos,0))
    df <- data.frame(x = pfinal[1, ], y = pfinal[2,])
    df2 <- data.frame(labelname = filenames[ii], x = xpos, y = ypos)
    p <- p + ggplot2::geom_path(data = df, color = as.character(class_color_code$colors[class_color_map[ii]]), size=1) + 
      ggplot2::annotate("text", x = xpos, y = ypos, label = filenames[ii], size = 2, fontface = 2)
  }
  p <- p + ggplot2::coord_fixed() + ggplot2::scale_y_reverse() +
        ggplot2::labs(x = sprintf('Eigen axis %d', eigdir[1]), y = sprintf('Eigen axis %d', eigdir[2])) + 
        ggplot2::ggtitle(titletxt) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot_set_axis()
  return(p)  
  
}


PCA_plot <- function(alpha_t_array, qmean, qarray, colorpath){
  Taxoncolorcodes <- readr::read_csv(colorpath,col_names = TRUE)
  filenames <- Taxoncolorcodes$specimen_ID
  color_code <- Taxoncolorcodes$color
  covdata = build_tpca_model_from_mean(qmean,alpha_t_array)
  eigproj_landmarks <- get_eigen_projections(covdata,alpha_t_array,1:5)
  R = pracma::eye(2)
  plot(c(min(eigproj_landmarks[[1]]),max(eigproj_landmarks[[1]])),c(min(eigproj_landmarks[[2]]),max(eigproj_landmarks[[2]])), type="n", xlab="Eigen Axis 1",ylab="Eigen Axis 2", main = "PCA scatter plot") 
  for(i in 1:length(eigproj_landmarks$`eig1`)){
    xpos = eigproj_landmarks[[1]][i]
    ypos = eigproj_landmarks[[2]][i]
    pfinal = repose_curve(q_to_curve(qarray[[i]]),0.12,R,c(xpos,ypos,0))
    plot_curve(pfinal,color_code[i],l = TRUE)
    #print(c(xpos,ypos))
    #print(filenames[i+1])#depending on heading
    text(xpos, ypos, labels = filenames[i],cex=0.3)
  }
  return(eigproj_landmarks)
}

#' Plot the geodesic path
#' @param alpha list of shapes along the geodesic 
#' @param color color name
#' @param L1 length(size) of first curve
#' @param L2 length(size) of second curve
#' @return figure handle
#' @export
plot_geodesic_path <- function(alpha, color = "cycle", L1 = 1, L2 = 1) {
  dim_alpha <- dim(alpha)
  n <- dim_alpha[1]
  T1 <- dim_alpha[2]
  k <- dim_alpha[3]
  
  if (n == 3)
    stop('Ploting for 3D curves not supported.')
  curve_path <- plyr::aaply(alpha, 3, q_to_curve)
  
  # Center all curves to origin 0 and shift them succesively by delta
  # Linearly interpolate the scale from L1 to L2 and rescale all the curves
  L <- seq(L1, L2, length=k)
  dt <- 0.25*(L1 + L2)/2
  shifts <- matrix(c(1:k*dt, rep(0, k)), nrow = 2, ncol = k, byrow = T)
  for (ii in 1:k) {
    curve_path[ii,,] <- sweep(curve_path[ii,,], 1, rowMeans(curve_path[ii,,]), FUN="-")
    curve_path[ii,,] <- L[ii]*curve_path[ii,,]
    curve_path[ii,,] <- sweep(curve_path[ii,,], 1, shifts[, ii], FUN="+")
  }
  
  if (color == "cycle") {
    jet_color <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    colors <- jet_color(k)
  }
  else
    colors <- rep("red", k)
  
  ii <- 1
  df <- data.frame(x = curve_path[ii,1, ], y = curve_path[ii,2,])
  p <- ggplot2::ggplot(df, ggplot2::aes(x , y) ) + ggplot2::geom_path(size = 1, color = colors[ii]) 
  for (ii in 2:k) {
    df <- data.frame(x = curve_path[ii,1, ], y = curve_path[ii,2,])
    p <- p + ggplot2::geom_path(data = df, color = colors[ii], size=1)
  }
  p <- p + ggplot2::coord_fixed() + ggplot2::scale_y_reverse()+  ggplot_axis_off()
  print(p)
  return(p)
}

#' Plot magnitude of the deformation field
#' @param alpha_t tangent vector for the deformation field
#' @param beta coordinate curve that will overlay the deformation field
#' @param titletext caption for the figure
#' @return figure handle
#' @export
plot_deformation_field <- function(alpha_t, beta, titletext="") {
  alpha_t_mag = c()  
  alpha_t <- alpha_t[,,2]
  for (i in 1:ncol(alpha_t)){
    alpha_t_mag[i] = pracma::Norm(alpha_t[,i])
  }
  alpha_t_mag = (alpha_t_mag - min(alpha_t_mag))/(max(alpha_t_mag) - min(alpha_t_mag))
  alpha_t_mag_smooth = smooth(alpha_t_mag, 3)
  alpha_t_mag_smooth = (alpha_t_mag_smooth - min(alpha_t_mag_smooth))/(max(alpha_t_mag_smooth) - min(alpha_t_mag_smooth))
  x = beta[1,]
  y = beta[2,]
  z = rep(0,length(x))
  X1 <- beta
  alpha_t_mag_smooth = as.vector(alpha_t_mag_smooth)
  alpha_t_mag = as.vector(alpha_t_mag)
  dat = data.frame(x = X1[1,], y = X1[2,], alpha_t_mag = alpha_t_mag, alpha_t_mag_smooth = alpha_t_mag_smooth)
  
  jet_color <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  
  p <- ggplot2::ggplot(dat, ggplot2::aes(x , y, color = alpha_t_mag_smooth)) + ggplot2::geom_path(size = 4) + 
    ggplot2::scale_y_reverse() + ggplot2::scale_color_gradientn(colors =jet_color(256)) + ggplot2::theme(legend.text=ggplot2::element_text(size=14)) +
    ggplot2::labs(color="Deformation") + ggplot2::labs(title=titletext) + ggplot2::theme(title = ggplot2::element_text(size = ggplot2::rel(1.8))) +
    ggplot2::xlab('') + ggplot2::ylab('') + ggplot2::coord_fixed() + ggplot_axis_off()
  
  print(p)
  return(p)
  
}
  
  
ggplot_axis_off <- function() {
  ggplot2::theme(
    axis.text.x=ggplot2::element_blank(),
    axis.text.y=ggplot2::element_blank(),
    axis.ticks=ggplot2::element_blank(),
    axis.title.x=ggplot2::element_blank(),
    axis.title.y=ggplot2::element_blank(),
    plot.margin=grid::unit(c(0,0,0,0), "mm")
  )
}

ggplot_set_axis <- function() {
  return(
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = ggplot2::rel(2)),
      axis.text.x = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(2)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(2))
      )
  )

}
