server <- function(input, output) {
  
  # output$value <- renderPrint({ input$text })
  
  
  output$contents <- renderText({
    req(input$file1)
    
    
    fid = read.table(input$file1$datapath, stringsAsFactors = FALSE)[[1]]
    
    return("File successfully uploaded")
    
    
  })
  
  
  
  
  
  
  output$mean_shape <- renderPlot({
    req(input$Mean_shape_plotButton)
    req(input$file1)
    req(input$var)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making plot", value = 0)
    X = main_closed(input$file1$datapath)
    ar = array(unlist(X),c(2,100,length(X)))
    if(input$var == "Open"){
      progress$inc(1)
      b = fdasrvf::curve_karcher_mean(ar,mode = "O",maxit = 6)
      
    }else{
      progress$inc(1)
      b = fdasrvf::curve_karcher_mean(ar,mode = "C",maxit = 6)
      
    }
    plot(b$betamean[1,],b$betamean[2,],type = "l",axes=FALSE,xlab = '', ylab = '', main = "Mean Shape")
    
  })
  
  
  
  
  
  output$geo_dist <- renderTable({
    req(input$goButton)
    req(input$file1)
    req(input$var)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating distances", value = 0)
    X = main_closed(input$file1$datapath)
    output_result <- matrix(rep(0,length(X)^2),nrow=length(X))
    n = (length(X)-1)*length(X)
    for(i in 1:(length(X)-1)){
      for (j in (i+1):length(X)){
        temp = NULL
        progress$inc(1, detail = paste("Calculating distance between", i, "and", j))
        
        if(input$var == "Open"){
          temp = fdasrvf::calc_shape_dist(X[[i]],X[[j]], mode = "O")
        }else{
          temp = fdasrvf::calc_shape_dist(X[[i]],X[[j]], mode = "C")
        }
        output_result[i,j] = temp
        output_result[j,i] = temp
        
        
      }
    }
    rc_name <- c()
    for(i in 1:length(X)){
      a <- paste("shape",i, sep = "")
      rc_name <- c(rc_name,a)
    }
    rownames(output_result) <- rc_name
    colnames(output_result) <- rc_name
    return(output_result)
  })
  
  
  observe({
    
    req(input$file1)
    #req(input$var == "Original Curves")
    X = main_closed(input$file1$datapath)
    output$plots <- renderUI({ get_plot_output_list(length(X),length(X),X) })
  })
  
  
  
  
  
  
  
  output$PCA_plot <- renderPlot({
    req(input$PCA_plotButton)
    req(input$file_ready)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making PCA plot", value = 1)
    data_pca <- R.matlab::readMat(input$mat_path)
    alpha_t_array <- list()
    for(i in 1:length(data_pca$alpha.array)) {
      alpha_t_array[[i]]<- data_pca$alpha.t.array[[i]][[1]]
    }
    qmean <- data_pca$qmean
    qarray <- data_pca$qarray
    Taxoncolorcodes <- readr::read_csv(input$color_path,col_names = FALSE)
    filenames <- Taxoncolorcodes$X1
    color_code <- Taxoncolorcodes$X3
    covdata = build_tpca_model_from_mean(qmean,alpha_t_array)
    eigproj_landmarks = save_eigen_projections(covdata,alpha_t_array,1:5,'eig_proto.csv')
    R = pracma::eye(2)
    plot(c(-0.15,0.15),c(-0.2,0.15), type="n", xlab="Eigen Axis 1",ylab="Eigen Axis 2", main = "PCA scatter plot of Humeri including indeterminate") 
    for(i in 1:length(eigproj_landmarks$`eig 1`)){
      xpos = eigproj_landmarks[[1]][i]
      ypos = eigproj_landmarks[[2]][i]
      pfinal = repose_curve(q_to_curve(qarray[[i]]),0.12,R,c(xpos,ypos,0))
      plot_curve(pfinal,color_code[i],l = TRUE)
      text(xpos, ypos, labels = filenames[i],cex=0.3)
    }
    
    
    
    
    
  })
  
  
  
  
  output$mds_plot <- renderPlot({
    req(input$mds_plotButton)
    req(input$file_ready)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making MDS plot", value = 1)
    data_mds <- R.matlab::readMat(input$mat_path)
    alpha_t_array <- list()
    for(i in 1:length(data_mds$alpha.array)) {
      alpha_t_array[[i]]<- data_mds$alpha.t.array[[i]][[1]]
    }
    qmean <- data_mds$qmean
    qarray <- data_mds$qarray
    geo_dis <- data_mds$geo.dist
    Taxoncolorcodes <- readr::read_csv(input$color_path,col_names = FALSE)
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
      plot_curve(pfinal,color_code[i], l = TRUE)
      text(xpos, ypos, labels = filenames[i],cex=0.3)
    }
    
    
    
    
    
    
  })
  
  
  output$dendogram_plot <- renderPlot({
    req(input$dendogram_plotButton)
    req(input$file_ready)
    data_ddg <- R.matlab::readMat(input$mat_path)
    Taxoncolorcodes <- readr::read_csv(input$color_path,col_names = FALSE)
    geo_dist <- as.vector(unlist(data_ddg$geo.dist))
    D.clust <- hclust(as.dist(pracma::squareform(geo_dist)),method="average")
    D.clust$labels <- Taxoncolorcodes$X1
    plot(D.clust, xlab = '', main='Dwarf Morphometric Placement', sub='', ylab='')
  })
  
  
  
  output$dfa_plot <- renderPlot({
    req(input$dfa)
    req(input$file_ready)
    data_df <- R.matlab::readMat(input$mat_path)
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
    
    
    ggplot2::ggplot(dat, ggplot2::aes(x , y, color = alpha_t_mag_smooth)) + ggplot2::geom_path(size = 7) + 
      ggplot2::scale_y_reverse() + ggplot2::scale_color_gradientn(colors =jet_color(256)) + ggplot2::theme(legend.text=ggplot2::element_text(size=14)) +
      ggplot2::labs(color="Deformation") + ggplot2::labs(title="Humeri_Including_Indeterminate2") + ggplot2::theme(title = ggplot2::element_text(size = ggplot2::rel(1.8))) +
      ggplot2::xlab('') + ggplot2::ylab('') + ggplot2::coord_fixed()
    
    
    
  })
  
  
  
  
  
  
  
}