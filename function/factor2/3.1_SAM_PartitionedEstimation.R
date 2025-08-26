##==============================================================================
##
## Script name:             00_Function_LSAM_PartitionedEstimation.R
##
## Purpose of script:       Compute loadings and residual variances for each
##                          measurement block separately and combine estimates
##                          into unified matrices
##
##==============================================================================

partitioned.estimation <- function(data = NULL, method = NULL) {
  
  # Determine estimator to be used in the cfa() function
  cfa.estimator <- switch(method,
                          "SAM_MLB" = "ML",
                          "SAM_ULSB" = "ULS",
                          "SAM_MultGroup" = "guttman",
                          "SAM_FABIN2_ULS" = "fabin2",
                          "SAM_Bentler_ULS" = "bentler",
                          "SAM_JS_aggregated" = "JSA")
  
  # Determine additional arguments for cfa() function
  est.args <- switch(method,
                     "SAM_MLB" = list(),
                     "SAM_ULSB" = list(),
                     "SAM_MultGroup" = list(),
                     "SAM_FABIN2_ULS" = list(thetapsi.method = "ULS"),
                     "SAM_Bentler_ULS" = list(GLS = FALSE),
                     "SAM_JS_aggregated" = list())
  
  # Count groups of indicators to define measurement blocks
  grp <- stringr::str_extract(colnames(data), "[A-Za-z]+")
  grp2 <- factor(grp, levels = c("x","z","y"))
  ind.names.split <- split(colnames(data), grp2)
  
  grps <- length(unique(grp))
  fac.names <- paste0("f", unique(grp))
  
  # Create LAMBDA and THETA matrix (all measurement blocks combined)
  LAMBDA <- matrix(0, nrow = ncol(data), ncol = grps)
  dimnames(LAMBDA) <- list(colnames(data), fac.names)
  THETA <- matrix(0, nrow = ncol(data), ncol = ncol(data))
  dimnames(THETA) <- list(colnames(data), colnames(data))
  
  # Initialize convergence issues
  con <- 0L
  
  # Estimate loadings and residual variances per block separately
  for(i in 1:grps) {
    
    data.sub <- data[ , ind.names.split[[i]]]
    
    if (ncol(as.data.frame(data.sub)) == 1) {
      data.sub <- as.data.frame(data.sub)
      names(data.sub) <- "y"
    }
    
    model <- paste0("f", unique(grp)[i], " =~ ", 
                    paste0(ind.names.split[[i]], collapse = " + "))
    
    fit <- lavaan::cfa(model = model, 
                   data = as.data.frame(data.sub),
                   estimator = cfa.estimator,
                   bounds = "standard",
                   estimator.args = est.args)
    
    if(!inherits(fit, "try-error")){
      
      if(lavaan::lavInspect(fit, "converged")) { con <- con + 1L }
      
      param_table <- lavaan::parameterTable(fit)
      united_table <- tidyr::unite(param_table, col = "param", lhs:rhs, sep = "")
      selected_table <- dplyr::select(united_table, param, est)
      est <- tibble::deframe(selected_table)
      
      lambda.sub <- est[grep("=~", names(est), value = TRUE)]
      theta.sub <- est[grep("~~f|=~", names(est), value = TRUE, invert = TRUE)]
      
      # Assign to appropriate elements in LAMBDA and THETA
      LAMBDA[ind.names.split[[i]] , fac.names[i]] <- lambda.sub
      if (ncol(data.sub) != 1) {
        diag(THETA[ind.names.split[[i]], ind.names.split[[i]]]) <- theta.sub
      }
      
    } else {
      
      LAMBDA[ind.names.split[[i]] , fac.names[i]] <- rep(0, 1)
      #diag(THETA[ind.names.split[[i]], ind.names.split[[i]]]) <- matrix(0,1,1)
      
    }
  }
  
  # Compute percentage of converged solutions
  conv <- con / grps
  
  return(list(LAMBDA = LAMBDA, THETA = THETA, convergence.pct = conv))
  
}
