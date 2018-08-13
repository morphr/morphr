p_adjust <- function (pvalues, method = "fdr") {
  valid_methods <- c(p.adjust.methods, "perm")
  if (!(method %in% valid_methods)) {
    warning(sprintf("%s is not a valid multiple comparisons method. Using no correction.", 
                    method), call. = FALSE)
    return(pvalues)
  }
  if (method %in% p.adjust.methods && method != "perm") {
    return(p.adjust(abs(pvalues), method))
  }
}

sign_tvalues <- function (tvalues){
  tvalues_sign <- sign(tvalues)
  tvalues_sign[tvalues_sign == 0] <- 1
  return(tvalues_sign)
}


shape_ttest <- function(X1, X2, paired=FALSE) {
  
  n1 <- dim(X1)[1]
  n2 <- dim(X2)[1]
  
  if (paired == TRUE) {
    D = X1 - X2
    D_mean <- colMeans(D)
    D_dev <- sweep(D, 2, D_mean)
    s1 <- sqrt(colSums(D_dev^2)/(n1-1))
    tvalues <- D_mean/(s1/sqrt(n1))
    pvalues <- 2*pt(abs(tvalues),n1-1, lower.tail = FALSE)
    pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  }
  else {
    X1_mean <- colMeans(X1)
    X2_mean <- colMeans(X2)
    
    X1_dev <- sweep(X1, 2, X1_mean)
    X2_dev <- sweep(X2, 2, X2_mean)
    
    SS1 <- colSums(X1_dev^2) # sum of squares of group 1
    SS2 <- colSums(X2_dev^2) # sum of squares of group 2
    
    s1_sq <- SS1/(n1-1)
    s2_sq <- SS2/(n2-1)
    
    se_diff <- sqrt(s1_sq/n1 + s2_sq/n2)
    tvalues <- (X1_mean - X2_mean)/(se_diff + .Machine$double.eps)
    # Calculate the degrees of freedom using the Welch–Satterthwaite approximation
    deg <- (s1_sq/n1 + s2_sq/n2)^2/(s1_sq^2/(n1^2*(n1-1)) + s2_sq^2/(n2^2*(n2-1)))
    pvalues <- 2*pt(abs(tvalues), deg, lower.tail = FALSE)
    pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  }
  return(list("tvalues"=tvalues, "pvalues"=pvalues))
}

shape_corr <- function(X, Y) {
  
  N <- length(Y)
  X_dev <- sweep(X, 2, colMeans(X))
  Y_dev <- Y - mean(Y)
  corr_coeff <- as.numeric((Y_dev %*% X_dev)/sqrt(colSums(X_dev^2)*sum(Y_dev^2)))
  tvalues <- corr_coeff * sqrt((N-2)/(1-corr_coeff^2 + .Machine$double.eps))
  pvalues <- 2*pt(abs(tvalues), N-2, lower.tail = FALSE)
  pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  return(list("tvalues"=tvalues, "pvalues"=pvalues, "corr_coeff"=corr_coeff))
}

lm_vec <- function(main_effect="", covariates, demographic,y) {
  data_array = as.matrix(demographic[,y])
  N <- length(data_array)
  Np <- length(c(main_effect,covariates))
  temp = paste(y,"~ ",sep = " ")
  if(main_effect==""){
    X <- model.matrix(as.formula(paste(temp, paste(covariates, collapse= "+"))),demographic)
    Npfull = length(covariates)
    Npnull = length(covariates)
    unique = ""
    fullmodel = paste(covariates, collapse= "+")
    X_design_full = X
    X_design_null = X
    fullvars = c("",covariates)
  }else{
    all_variable = c(main_effect,covariates)
    X <- model.matrix(as.formula(paste(temp, paste(all_variable, collapse= "+"))),demographic)
    Npfull = length(all_variable)
    Npnull = length(covariates)
    unique = main_effect
    fullmodel = paste(all_variable, collapse= "+")
    X_design_full = X
    X_design_null = X[,c(1,3:ncol(X))]
    fullvars = c(main_effect,covariates)
  }
  X_hat <- solve(t(X) %*% X) %*% t(X) # pre hat matrix
  beta_coeff <- X_hat %*% data_array  # beta coefficients
  Y <- X %*% beta_coeff  # predicted response
  rss <- colSums((data_array - Y)^2) # residual sum of squares
  
  if (main_effect == "")
    main_effect = "(Intercept)" # If main_effect is empty, return the parameters of the Intercept
  if(main_effect%in%sqrt(diag(solve(t(X) %*% X)))){
    se <- sqrt(diag(solve(t(X) %*% X)))[[main_effect]] * sqrt(rss / (N-Np-1))
    tvalues <- as.numeric(beta_coeff[main_effect, ]/(se + .Machine$double.eps)) # tvalue
  }else{
    se <- sqrt(diag(solve(t(X) %*% X)))[[2]] * sqrt(rss / (N-Np-1))
    tvalues <- as.numeric(beta_coeff[2, ]/(se + .Machine$double.eps)) # tvalue
  }
  pvalues <- 2*pt(abs(tvalues), N-Np-1, lower.tail = FALSE) # pvalu
  residuals <- data_array - Y
  pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  return(list(pvalues = pvalues, tvalues = tvalues, beta_coeff = beta_coeff, 
              rss = rss, residuals = residuals,Npfull = Npfull, Npnull = Npnull,
              fullmodel = fullmodel, unique = unique,fullvars = fullvars,
              X_design_full = X_design_full, X_design_null = X_design_null,se = se))
}

anova_vec <- function(lm_full, lm_null, demographics,data_array) {
  
  N <- length(data_array)
  Fstat <- (lm_null$rss - lm_full$rss)/lm_full$rss * (N - lm_full$Npfull - 1)/(lm_full$Npfull - lm_null$Npnull)  # F statistic
  
  model_unique_idx <- which(lm_null$unique %in% lm_null$fullvars) + 1    # Add 1, because the first column in the design matrix is the intercept
  
  se_full_unique <- sqrt(diag(solve(t(lm_full$X_design_full) %*% lm_full$X_design_full)))[model_unique_idx] *
    sqrt(lm_full$rss / (N - lm_full$Npfull - 1))
  
  tvalues <- lm_full$beta_coeff[model_unique_idx, ]/(se_full_unique + .Machine$double.eps)
  pvalues <- 1 - pf(Fstat, lm_full$Npfull - lm_null$Npnull, N - lm_full$Npfull - 1)
  pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  tvalues_sign <- sign_tvalues(tvalues)
  return(list(pvalues = pvalues, tvalues = tvalues, tvalues_sign = tvalues_sign, 
              se = se_full_unique, Fstat = Fstat))
}

shape_anova <- function(main_effect="", covariates="", demographic, y) {
  
  
  message('Running the statistical model. This may take a while...', appendLF = FALSE)
  lm_full <- lm_vec(main_effect = main_effect, covariates = covariates,demographic = demographic,y = y)
  lm_null <- lm_vec(main_effect = "", covariates = covariates,demographic = demographic,y = y)
  anova_model <- anova_vec(lm_full,lm_null,demographic,as.matrix(demographic[, main_effect]))
  anova_model$pvalues[is.nan(anova_model$pvalues)] <- 1
  anova_model$pvalues <- anova_model$pvalues*anova_model$tvalues_sign
  anova_model$tvalues[abs(anova_model$pvalues) >= 0.05] <- 0
  anova_model$pvalues_adjusted <- p_adjust(anova_model$pvalues)
  message('Done.')
  return(anova_model)
}

calculate_p_value <- function(main_effect="", covariates="", loaded_data,demographic_path,size = 150){
  original_demographic = readr::read_csv(demographic_path)
  original_demographic = as.data.frame(original_demographic)
  p_value = c()
  for(i in 1:size){
    p_value[i] = shape_anova(main_effect,covariates,loaded_data,paste("V",i,sep = ''))$pvalues_adjusted[[1]]
  }
  return(p_value)
}
