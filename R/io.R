#' Read shape coordinates from svg files
#' @param file_path location of svg files
#' @param N number of points (default value 100)
#' @param r number of iterations (default value 3)
#' @return coordinates of shape
#' @export
read_svg_file <- function(file_path, N = 100, r = 10){
  #Function require XML library. Read in the command in the SVG file
  file_type = "O"
  doc <- XML::htmlParse(file_path)
  p <- XML::xpathSApply(doc, "//path", XML::xmlGetAttr, "d")
  if(stringr::str_detect(p, "z")){
    file_type = "C"
  }
  if(file_type == "C"){
    element <- strsplit(p, "\\, |\\,| ")[[1]]
    
    #Create the starting point and make sure we create the starting point correctly
    segment = c()
    curr_pos = c(0,0)
    pos = c(as.numeric(element[2]),as.numeric(element[3]))
    curr_pos = curr_pos + pos
    start_pos = curr_pos
    control1 = c(as.numeric(element[5]), as.numeric(element[6]))
    control2 = c(as.numeric(element[7]), as.numeric(element[8]))
    end = c(as.numeric(element[9]), as.numeric(element[10]))
    control1 = control1+curr_pos
    control2 = control2+curr_pos
    end = end + curr_pos
    segment_part = matrix(c(start_pos[1],control1[1],control2[1],end[1],start_pos[2],control1[2],control2[2],end[2]),4,2)
    segment = list()
    segment[[1]] = segment_part
    curr_pos = end
    
    #All the rest of the points are dependent on the starting point. Calulate all of them in a loop and append to segment as a 4*2 matrix
    j=2
    for(i in seq(from = 11, to = length(element)-1, by = 6)){
      start_pos = curr_pos
      control1 = c(as.numeric(element[i]), as.numeric(element[i+1]))
      control2 = c(as.numeric(element[i+2]), as.numeric(element[i+3]))
      end = c(as.numeric(element[i+4]), as.numeric(element[i+5]))
      control1 = control1+curr_pos
      control2 = control2+curr_pos
      end = end + curr_pos
      curr_pos = end
      segment_part = matrix(c(start_pos[1],control1[1],control2[1],end[1],start_pos[2],control1[2],control2[2],end[2]),4,2)
      segment[[j]] = segment_part
      j = j+1
    }
    
    #All the points we generate if we choose to generate r point for each segment of Bezier curve. Need knotR package
    point <- c()
    for(i in 1:length(segment)){
      xy <- knotR::bezier(segment[[i]], n = r) #This is the problem! We should always specify n = (a number)
      point <- rbind(point, xy)
    }
    
    #Use fdasrvf package to create equal distance point based on specified N and plot those points
    point_resample <- resample_curve(t(point), N = N)
    return(point_resample)
    
  }else{
    p_modified <- stringr::str_extract_all(p,"\\-?[0-9.cC]+")[[1]]
    element <- unlist(strsplit(p_modified, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE))
    
    
    
    segment = c()
    curr_pos = c(0,0)
    pos = c(as.numeric(element[1]),as.numeric(element[2]))
    curr_pos = curr_pos + pos
    start_pos = curr_pos
    control1 = c(as.numeric(element[4]), as.numeric(element[5]))
    control2 = c(as.numeric(element[6]), as.numeric(element[7]))
    end = c(as.numeric(element[8]), as.numeric(element[9]))
    control1 = control1+curr_pos
    control2 = control2+curr_pos
    end = end + curr_pos
    segment_part = matrix(c(start_pos[1],control1[1],control2[1],end[1],start_pos[2],control1[2],control2[2],end[2]),4,2)
    segment = list()
    segment[[1]] = segment_part
    curr_pos = end
    j=2
    
    for(i in seq(from = 10, to = length(element), by = 7)){
      
      start_pos = curr_pos
      control1 = c(as.numeric(element[i+1]), as.numeric(element[i+2]))
      control2 = c(as.numeric(element[i+3]), as.numeric(element[i+4]))
      end = c(as.numeric(element[i+5]), as.numeric(element[i+6]))
      if(element[i] != "C"){
        control1 = control1+curr_pos
        control2 = control2+curr_pos
        end = end + curr_pos
      }
      curr_pos = end
      segment_part = matrix(c(start_pos[1],control1[1],control2[1],end[1],start_pos[2],control1[2],control2[2],end[2]),4,2)
      segment[[j]] = segment_part
      j = j+1
    }
    
    
    point <- c()
    for(i in 1:length(segment)){
      xy <- knotR::bezier(segment[[i]], n=4)
      point <- rbind(point, xy)
    }
    #plot(point)
    
    point_resample <- resample_curve(t(point), N = N)
    return(point_resample)
    
  }
  
}

#' Read shape coordinates from ucf files
#' @param file_path location of ucf files
#' @param ndim number of dimensions to read. A ucf file can contain multidimensional curves
#' @return list of coordinates of shape
#' @export
read_ucf_multiple_levels <- function(filepath, ndim = 2){
  fid <- file(filepath, "r")
  file_content <- readLines(fid)
  close(fid)  
  j = 1
  num_levels = NULL
  while(j<=length(file_content)){
    if(file_content[j] == "<levels>"){
      num_levels = as.numeric(file_content[j+1])
    }
    j = j+1
  }
  N = NULL
  Xtemp = list()
  X = list()
  for(i in 1:num_levels){
    k = 1
    while(k<=length(file_content)){
      if(file_content[k] == "<point_num=>"){
        N = as.numeric(file_content[k+1])
        for(l in 1:N){
          temp <- as.numeric(strsplit(file_content[k+2+l], " ")[[1]])
          Xtemp[[l]] = temp
          
        }
        break
      }
      
      k = k+1
    }
    
    Xtemp <- t(matrix(unlist(Xtemp),3,N))
    Xtemp <- Xtemp[, 1:ndim]
    X[[i]] = t(Xtemp)
  }
  if (num_levels == 1) {
    X <- X[[i]]
  }

  return (X)
  
}

main_closed <- function(path){
  fid = read.table(path, stringsAsFactors = FALSE)[[1]]
  i = 1
  X = NULL
  while(i<=length(fid)){
    fname <- fid[i]
    if(stringr::str_detect(fname,"\\.ucf")){
      Xtemp = read_ucf_multiple_levels(fname)
      X[[i]] = t(Xtemp[[1]][,1:2])
      i = i+1
    }else{
      Xtemp = read_svg_file(fname)
      X[[i]] <- t(Xtemp[,1:2])
      i = i+1
    }
    
    
  }
  
  qarray = list()
  for(i in 1:length(X)){
    qarray[[i]] = curve_to_q(X[[i]])
  }
  
  n = nrow(qarray[[1]])
  T_col = ncol(qarray[[1]])
  return(X)
  
}

#' @export
read_curve <- function(curve_file, ndim=2) {
  if (! file.exists(curve_file) ) {
    stop(sprintf('Curve file %s does not exist.', curve_file), call. = FALSE)
  }
  
  ext <- tools::file_ext(curve_file)
  switch(ext,
         "svg" = {return( read_svg_file(curve_file) )},
         "ucf" = {return( read_ucf_multiple_levels(curve_file, ndim = ndim) )}
         )
}




