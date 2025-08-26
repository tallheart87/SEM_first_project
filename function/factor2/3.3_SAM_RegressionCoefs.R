##==============================================================================
##
## Script name:             00_Function_RegressionCoefs.R
##
## Purpose of script:       Compute regression coefficients from cov matrix
##
##==============================================================================

LSAM_regcoef <- function(model = NULL, sumstats = NULL) {
  
  # extract regressions from model syntax
  regressions <- gsub(" ", "", grep(pattern = " ~ f", 
                      strsplit(model, "\n")[[1]], invert = FALSE, value = TRUE))
  
  # extract dependent variables
  DVs <- sub("~.*", "", regressions) 
  IVs <- strsplit(sub(".*~", "", regressions), split = "\\+")
  
  # create empty vector
  betas <- numeric(length = sum(lengths(IVs)))
  names(betas) <- paste0(paste0(rep(DVs, times = lengths(IVs)), sep = "~"), 
                         unlist(IVs))
  
  # compute betas
  for(i in 1:length(DVs)) {
    
    betas[(1 + sum(lengths(IVs[1:i-1]))):sum(lengths(IVs[1:i]))] <- 
      sumstats[DVs[i], IVs[[i]]] %*% solve(sumstats[IVs[[i]], IVs[[i]]])
    
  }

  # return betas
  return(betas)
  
}
