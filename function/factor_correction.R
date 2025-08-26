#### function for transform the loading matrix
check_lambda <- function(loadings){
  # For each observed variable, "max_abs_indices" corresponds to the factor with the largest absolute loading
  max_abs_indices <- apply(loadings, 1, function(x) which.max(abs(x)))
  # Create a same size matrix as "loadings" which is filled with zeros
  transformed_loadings <- matrix(0, nrow = nrow(loadings), ncol = ncol(loadings))
  # Assign +1 or -1 to the largest absolute loading in each observed variable
  for (u in 1:nrow(loadings)) {
    max_index <- max_abs_indices[u]
    max_value <- loadings[u, max_index]
    transformed_loadings[u, max_index] <- ifelse(max_value > 0, 1, -1)
  }
  return(transformed_loadings)
}

#### correct directions
check_direction <- function(check_load, train_score, true_load, factor_n){
  perm <- gtools::permutations(factor_n, factor_n)
  absdiff <- c()
  for (p in 1:nrow(perm)) {
    corsign <- sign(diag(cor(true_load, check_load[, perm[p,]])))
    L_res <- (check_load[, perm[p,]]) %*% diag(corsign)
    absdiff[p] <- sum(rowSums(abs(true_load - L_res)))
  }
  bestperm <- which.min(absdiff)
  loadings <- check_load[, perm[bestperm,]]
  corsign <- sign(diag(cor(true_load, loadings)))
  changed_loadings <- loadings %*% diag(corsign)
  
  if(sum(is.na(train_score)) == 1){
    train_score_T <- NA
  }else{
    train_score_T <- (train_score[, perm[bestperm,]]) %*% diag(corsign)
  }
  
  return(list(loadings = changed_loadings, 
              score = train_score_T))
}
