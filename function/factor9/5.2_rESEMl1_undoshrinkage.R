#2. REGULARIZED ESEM WITH THE LASSO
rESEMl1_undoshrinkage <- function(DATA, R, P, MaxIter, eps){
  lambda <- 0 ###undo shrinkage
  I <- dim(DATA)[1]
  J <- dim(DATA)[2]
  convAO <- 0
  iter <- 1
  Lossc <- 1
  Lossvec <- Lossc   #Used to compute convergence criterium
  svd1 <- svd(DATA, R, R)
  P2 <- svd1$v %*% diag(svd1$d[1:R])/sqrt(I)
  P1 <- matrix(rnorm(J*R), nrow = J, ncol = R) #initialize P 
  loadings <- .3*P1 + .7*P2
  scores <- svd1$u 
  diffT <- 0
  diffP <- 0
  while (convAO == 0) {
    iter0 <- 1
    Losst <- 1
    Lossvec0 <- Losst
    convT0 <- 0
    Lossvec1 <- 1
    #1. Update component scores 
    while(convT0 == 0){
      Lossu1old <- LOSS(DATA,scores,loadings,lambda)
      E <- DATA - scores%*%t(loadings) 
      for (r in 1:R){
        Er <- E + scores[,r]%*%t(loadings[,r])
        num <- Er%*%loadings[,r]
        scores[,r] <- sqrt(I)*num/sqrt(sum(num^2))
        Lossu1 <- LOSS(DATA,scores,loadings,lambda)
        diffT <- c(diffT,Lossu1old-Lossu1)
        Lossu1old <- Lossu1
        Lossvec1 <- c(Lossvec1, Lossu1)
      }
      #t(scores)%*%scores
      #Calculate loss
      Lossu0 <- LOSS(DATA,scores,loadings,lambda)
      Lossvec0 <- c(Lossvec0,Lossu0)
      # check convergence
      if (iter0 > MaxIter) {
        convT0 <- 1
      }
      if (abs(Losst-Lossu0) < eps){
        convT0 <- 1
      }
      iter0 <- iter0 + 1
      Losst <- Lossu0
    }
    Loss <- LOSS(DATA, scores, loadings, lambda)
    
    #2. Update loadings
    Lossu1old <- LOSS(DATA,scores,loadings,lambda)
    E <- DATA - scores%*%t(loadings) 
    for (r in 1:R){
      Er <- E+scores[,r]%*%t(loadings[,r])
      crosstEr <- t(Er)%*%scores[,r]
      loadings[,r]<-sign(crosstEr)*apply(cbind(abs(crosstEr)-lambda/2,0),1,max)/I
    }
    loadings[P==0]<-0###fixed zero
    
    #Calculate loss
    Lossu <- LOSS(DATA,scores,loadings,lambda)
    Lossvec <- c(Lossvec,Lossu)
    if (iter > MaxIter) {
      convAO <- 1
    }
    if (abs(Lossc-Lossu) < eps) {
      convAO <- 1
    }
    iter <- iter + 1
    Lossc <- Lossu
  }
  return_rlslv <- list()
  return_rlslv$loadings <- loadings
  return_rlslv$scores <- scores
  return_rlslv$Loss <- Loss
  return_rlslv$Lossvec <- Lossvec
  
  return(return_rlslv)
}
