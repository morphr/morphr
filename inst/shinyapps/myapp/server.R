server <- function(input, output) {
  
  # output$value <- renderPrint({ input$text })
  
  
  output$contents <- renderText({
    req(input$file1)
    req(input$Color_File)
    
    
    
    return("Files successfully uploaded")
    
    
  })
  
  
  
  
  
  
  output$mean_shape <- renderPlot({
    req(input$Mean_shape_plotButton)
    req(input$file1)
    req(input$Color_File)
    req(input$var)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making plot", value =1)
    X = main(input$file1$datapath)
    qarray = list()
    for(i in 1:length(X)){
      qarray[[i]] = curve_to_q(X[[i]])
    }
    if(input$var == "Open"){
      #Modification
      progress$inc(1)
      b = fdasrvf::curve_karcher_mean(ar,mode = "O",maxit = 6)
      
    }else{
      progress$inc(1)
      temp_all_mean_shape = find_mean_shape(qarray)
      saveRDS(temp_all_mean_shape, file = "all_mean_shape.rds")
      qmean = temp_all_mean_shape[[1]]
      alpha_array = temp_all_mean_shape[[2]]
      alpha_t_array = temp_all_mean_shape[[3]]
      norm_alpha_t_mean = temp_all_mean_shape[[4]]
      gamma_array = temp_all_mean_shape[[5]]
      sum_sq_dist = temp_all_mean_shape[[6]]
      qmean_new = project_curve(qmean)
      pmean = q_to_curve(qmean_new)
      pmean_new = rbind(pmean[1,],-pmean[2,])
      plot_curve(pmean,'blue', l = FALSE, filename = "Mean Shape")
    }
    
  })
  
  
  
  
  
  output$geo_dist <- plotly::renderPlotly({
    req(input$file1)
    req(input$Color_File)
    req(input$var)
    req(input$goButton)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating distances", value = 1)
    X = main(input$file1$datapath)
    qarray = list()
    for(i in 1:length(X)){
      qarray[[i]] = curve_to_q(X[[i]])
    }
    output_result <- matrix(rep(0,length(X)^2),nrow=length(X))
    n = (length(X)-1)*length(X)
    all_geodesic = geodesic_distance_all(qarray)
    saveRDS(all_geodesic, file = "all_geodsic.rds")
    geo_dist = all_geodesic[[6]]
    output_result[lower.tri(output_result, diag=FALSE)] <- unlist(geo_dist)
    output_result <- t(output_result)
    for(i in 1:(length(X)-1)){
      for (j in (i+1):length(X)){
        
        output_result[j,i] = output_result[i,j]
        
        
      }
    }
    rc_name <- c()
    Taxoncolorcodes <- readr::read_csv(input$Color_File$datapath,col_names = FALSE)
    filenames <- Taxoncolorcodes$X1
    for(i in 1:length(X)){
      rc_name <- c(rc_name,filenames[i])
    }
    rownames(output_result) <- rc_name
    colnames(output_result) <- rc_name
    return(plotly::plot_ly(x = rc_name, y = rc_name ,z = output_result, type = "heatmap"))
  })
  
  
  observe({
    
    req(input$file1)
    req(input$Color_File)
    #req(input$var == "Original Curves")
    X = main(input$file1$datapath)
    Taxoncolorcodes <- readr::read_csv(input$Color_File$datapath,col_names = TRUE)
    filenames <- Taxoncolorcodes$specimen_ID
    color_code <- Taxoncolorcodes$color
    class_color_code <- data.frame(label=levels(factor(Taxoncolorcodes$class)), 
                                   colors = colorRampPalette(RColorBrewer::brewer.pal(9, name = "Set1"))(nlevels(factor(Taxoncolorcodes$class))))
    
    class_color_map <- match(Taxoncolorcodes$class, class_color_code$label)
    output$plots <- renderPlot({
      if(length(X)<=8){
        par(mfrow=c(1,length(X)))
        par(mar=c(1,1,1,1))
      }else{
        rr <- ceiling(length(X)/7)
        par(mfrow=c(rr,7))
        par(mar=c(1,1,1,1))
      }
      for(i in 1:length(X)){
        plot_curve(X[[i]],as.character(class_color_code$colors[class_color_map[i]]), filename = filenames[i])
      }
    }
    )
  })
  
  
  
  
  
  
  
  output$PCA_plot <- renderPlot({
    req(input$Mean_shape_plotButton)
    req(input$PCA_plotButton)
    req(input$Color_File)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making PCA plot", value = 1)
    X = main(input$file1$datapath)
    qarray = list()
    for(i in 1:length(X)){
      qarray[[i]] = curve_to_q(X[[i]])
    }
    temp_all_mean_shape = readRDS(file = "all_mean_shape.rds")
    plotpca(temp_all_mean_shape$qmean, temp_all_mean_shape$alpha_t_array, qarray,
            input$Color_File$datapath, eigdir = c(1, 2), titletxt = "PCA Plot")
    
    
    
    
    
  })
  
  
  
  
  output$mds_plot <- renderPlot({
    req(input$Mean_shape_plotButton)
    req(input$goButton)
    req(input$mds_plotButton)
    req(input$Color_File)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making MDS plot", value = 1)
    X = main(input$file1$datapath)
    qarray = list()
    for(i in 1:length(X)){
      qarray[[i]] = curve_to_q(X[[i]])
    }
    temp_all = readRDS(file = "all_geodsic.rds")
    temp_all_mean_shape = readRDS(file = "all_mean_shape.rds")
    mdsplot(temp_all$alpha_t, temp_all$geo_dist,X,input$Color_File$datapath )

    
    
    
    
    
    
  })
  
  
  output$dendogram_plot <- renderPlot({
    req(input$dendogram_plotButton)
    req(input$goButton)
    req(input$Color_File)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making Dendogram", value = 1)
    X = main(input$file1$datapath)
    qarray = list()
    for(i in 1:length(X)){
      qarray[[i]] = curve_to_q(X[[i]])
    }
    temp_all = readRDS(file = "all_geodsic.rds")
    geo_dist = temp_all[[6]]
    plot_dendrogram(geo_dist, input$Color_File$datapath)
  })
  
  
  
  output$dfa_plot <- renderPlot({
    req(input$Mean_shape_plotButton)
    req(input$dfa)
    req(input$Color_File)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making Deformation Field Plot", value = 1)
    X = main(input$file1$datapath)
    qarray = list()
    for(i in 1:length(X)){
      qarray[[i]] = curve_to_q(X[[i]])
    }
    temp_all_mean_shape = readRDS(file = "all_mean_shape.rds")
    qmean = temp_all_mean_shape[[1]]
    alpha_t_array = temp_all_mean_shape[[3]]
    pmean = q_to_curve(project_curve(qmean))
    deformation_field_all(alpha_t_array, pmean, qmean,X)
    
    
    
  })
  
  output$PCA_variation_plot <- renderPlot({
    req(input$Mean_shape_plotButton)
    req(input$PCA_variation_plotButton)
    req(input$Color_File)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Making PCA Variation Plot", value = 1)
    X = main(input$file1$datapath)
    qarray = list()
    for(i in 1:length(X)){
      qarray[[i]] = curve_to_q(X[[i]])
    }
    temp_all_mean_shape = readRDS(file = "all_mean_shape.rds")
    plot_pca_variation(temp_all_mean_shape$alpha_t_array, temp_all_mean_shape$qmean, eigdir = 1,qarray, dt = 0)
    
    
    
    
  })
  
  
  
  
  
  
  
}
