##==============================================================================
##
## Script name:             00_Function_LSAM_ComputeSummaryStats.R
##
## Purpose of script:       Compute summary statistics E(eta) and var(eta) via
##                          explicit formulae Rosseel & Loh (2022), using three
##                          different mapping matrices (ML, GLS, ULS) and using
##                          the Fuller (1987) approach to ensure var(eta) is 
##                          positive definite
##
##==============================================================================

LSAM_SumStats <- function(S = NULL, 
                          sample.nobs = NULL,
                          LAMBDA = NULL,
                          THETA = NULL,
                          mapping = NULL) {
  
  if(mapping == "ML") {
    
    # do not allow for zero values on the diagonal (cannot be inverted)
    zero.theta.idx <- which(abs(diag(THETA)) < 1e-06)
    if (length(zero.theta.idx) > 0L) {
      tmp <- THETA[-zero.theta.idx, -zero.theta.idx, drop = FALSE]
      tmp.inv <- solve(tmp)
      THETA.inv <- matrix(0, nrow = nrow(THETA), ncol = ncol(THETA))
      THETA.inv[-zero.theta.idx, -zero.theta.idx] <- tmp.inv
      diag(THETA.inv)[zero.theta.idx] <- 1
    } else {
      THETA.inv <- solve(THETA)
    }
    
    # compute M-matrix
    M.hat <- solve(t(LAMBDA) %*% THETA.inv %*% LAMBDA) %*%
      t(LAMBDA) %*% THETA.inv
    
  } else if(mapping == "GLS") {
    
    M.hat <- solve(t(LAMBDA) %*% solve(S) %*% LAMBDA) %*%
      t(LAMBDA) %*% solve(S)
    
  } else if(mapping == "ULS") {
    
    M.hat <- solve(t(LAMBDA) %*% LAMBDA) %*% t(LAMBDA)
    
  } else {
    
    warning("choose 'ML', 'GLS', or 'ULS' for the mapping argument")
    
  }
  
  # Rosseel & Loh (2021), p. 18 eq. 21
  # Apply correction is lambda.corr is smaller than the threshold
  S.star <- M.hat %*% S %*% t(M.hat)
  THETA.star <- M.hat %*% THETA %*% t(M.hat)
  lambda.hat <- lavaan:::lav_matrix_symmetric_diff_smallest_root(S.star, THETA.star)
  cutoff <- 1 + 1/(sample.nobs-1)
  
  if(lambda.hat < cutoff) {
    
    lambda.corr <- lambda.hat - 1/(sample.nobs-1)
    # cat("Applying Fuller correction \n")
    var_eta <- S.star - (lambda.corr * THETA.star)
    
  } else {
    
    # cat("No Fuller correction \n")
    var_eta <- S.star - THETA.star
    
  }
  
  if(!is.null(dimnames(LAMBDA))) {
    dimnames(var_eta) <- list(colnames(LAMBDA), colnames(LAMBDA))
  }
  
  return(var_eta)
  
}
