PL <- function(T_loadings, E_loadings){
  cells <- length(T_loadings)
  
  zero_true <- T_loadings == 0
  zero_est <- E_loadings == 0
  correct_zero <- sum(zero_est & zero_true)
  
  nonzero_true <- T_loadings != 0
  nonzero_est <- E_loadings != 0
  correct_nonzero <- sum(nonzero_est & nonzero_true)
  
  rate <- (correct_zero + correct_nonzero) / cells
  
  return(rate)
}


