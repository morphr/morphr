#' Plot the curve p with color value
#' @param p coordinates of curve
#' @param colorstr color value
#' @param l Default value is False, which means returning a plot
#' @export
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
  }
  styleval <- '-'
  n = nrow(p)
  T_col = ncol(p)
  if(n == 2){
    if (l == TRUE){
      lines(p[1,],p[2,],type = "l",col = colorval)
    }else{
      plot(p[1,],p[2,],type = "l",col = colorval)
    }
  }
  if(n == 3){
    if (l == TRUE){
      lines(p[1,],p[2,],type = "l",col = colorval)
    }else{
      plot(p[1,],p[2,],type = "l",col = colorval)
    }
    for (i in 1:T_col){
      text(p[1,i],p[2,i],toString(i),cex = 0.7)
    }
    text(p[1,1],p[2,1],toString(0))
  }
  
}

#' Utility function for shiny app
#' @param maxplots maximum number of plots
#' @param input_n number of iterations
#' @param X shape
#' @return plots
get_plot_output_list <- function(max_plots, input_n, X) {
  # Insert plot output objects the list
  
  
  
  plot_output_list <- lapply(1:input_n, function(i) {
    plotname <- paste("Original Shape", i, sep="")
    plot_output_object <- plotOutput(plotname, height = 280, width = 250)
    plot_output_object <- renderPlot({
      plot(X[[i]][1,],X[[i]][2,], type = "l",axes=FALSE,xlab = '', ylab = '', main = plotname)
      
    })
  })
  
  do.call(tagList, plot_output_list) # needed to display properly.
  
  return(plot_output_list)
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
plot_pca_variation <- function(pathname){
  data_ppa <- R.matlab::readMat(pathname)
  alpha_t_array <- list()
  for(i in 1:length(data_ppa$alpha.array)) {
    alpha_t_array[[i]]<- data_ppa$alpha.t.array[[i]][[1]]
  }
  qmean <- data_ppa$qmean
  qarray <- data_ppa$qarray
  cat('Building TPCA model ...')
  covdata = build_tpca_model_from_mean(qmean,alpha_t_array)
  cat('done.\n')
  display_eigen_and_original_shapes_along_extremes(X,qmean,alpha_t_array,covdata,1,8)
}


# Needs modification
#mat_file_path_name = '/Users/xlzzxl/Desktop/files/drive-download-20180305T024235Z-001/all_mean.mat'
#csv_datafile = '/Users/xlzzxl/Desktop/files/drive-download-20180305T024235Z-001/Taxoncolorcodes.csv'
mdsplot <- function(mat_file_path_name, color_code_path_name){
  data_mds <- R.matlab::readMat(mat_file_path_name)
  alpha_t_array <- list()
  for(i in 1:length(data_mds$alpha.array)) {
    alpha_t_array[[i]]<- data_mds$alpha.t.array[[i]][[1]]
  }
  #qmean <- data_mds$qmean
  #qarray <- data_mds$qarray
  geo_dis <- data_mds$geo.dist
  Taxoncolorcodes <- read.csv(color_code_path_name, header = F, stringsAsFactors = FALSE)
  filenames <- Taxoncolorcodes[,1]
  color_code <- Taxoncolorcodes[,3]
  geo_dis <- as.vector(unlist(geo_dis))
  distance_matrix <- pracma::squareform(geo_dis)
  example_NMDS=vegan::metaMDS(distance_matrix,k=2,trymax=200)
  Y = example_NMDS$points
  theta = pi
  R = pracma::eye(2)
  L = 0.08
  labels = list()
  plot(c(-0.2,0.17),c(-0.15,0.15), type="n", xlab="X",ylab="Y", main = "MDS scatter plot including Dwarf") 
  for(i in 1:ncol(distance_matrix)){
    xpos = (-1)*Y[i,1]
    ypos = (-1)*Y[i,2]
    zpos = 0
    labels[[i]] = filenames[[i]]
    pfinal = repose_curve(q_to_curve(curve_to_q(data_mds$X[[i]][[1]])),0.12,R,c(xpos,ypos,0))
    plot_curve(pfinal,color_code[i],l = TRUE)
    text(xpos, ypos, labels = filenames[i],cex=0.3)
  }
  
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
    plot(c(min(X[1,]),max(X[1,])),c(min(X[2,]),max(X[2,])), type="n") 
    plot_curve(X,'r')
  }
}


#' Display Rshiny
#' @export
run_shiny <- function() {
  appDir <- system.file("shiny-examples", "myapp", package = "morphr")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `morphr`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}



#need modification
deformation_field_all <- function(pathname){
  data_df <- R.matlab::readMat(pathname)
  alpha_t_mean = data_df$alpha.t.mean
  alpha_t_mean = 0
  alpha_t_array <- list()
  for(i in 1:length(data_df$alpha.array)) {
    alpha_t_array[[i]]<- data_df$alpha.t.array[[i]][[1]]
  }
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
  pmean = data_df$pmean
  qmean = data_df$qmean
  pmean = q_to_curve(qmean)
  x = pmean[1,]
  y = pmean[2,]
  z = rep(0,length(x))
  X1 <- data_df$X[[1]][[1]]
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

#Modification
plot_dendrogram <- function(mat_file_path_name, color_code_path_name){
  data <- R.matlab::readMat("/Users/xlzzxl/srp99/ontogeny_Hypacrosaurus_stebingeri/all_mean.mat")
  geo_dist <- as.vector(unlist(data$geo.dist))
  taxon_codes <- read.csv('/Users/xlzzxl/srp99/ontogeny_Hypacrosaurus_stebingeri/Taxoncolorcodes.csv', header = F)
  
  D.clust <- hclust(as.dist(squareform(geo_dist)),method="average")
  D.clust$labels <- taxon_codes$V1
  plot(D.clust, xlab = '', main='Dwarf Morphometric Placement', sub='', ylab='')
}
