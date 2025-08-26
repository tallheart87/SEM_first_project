INITLOADINGS <- function(DATA, R, CARD){
  J <- dim(DATA)[2]
  I <- dim(DATA)[1]
  svd1 <- svd(DATA, R, R)
  P1 <- matrix(rnorm(J*R), ncol = R, nrow = J)
  P2 <- svd1$v %*% diag(svd1$d[1:R], nrow = R, ncol = R)/sqrt(I)
  P <- P1*.3 + P2*.7
  if (length(CARD)>1){
    CARDc <- J-CARD
    for (r in 1:R){
      ind <- sort(abs(P[,r]), index.return = TRUE)
      P[ind$ix[1:CARDc[r]],r] <- 0
    } 
  } else {
    CARDc <- J*R-CARD
    ind <- sort(abs(P), index.return = TRUE)
    P[ind$ix[1:CARDc]] <- 0
  }
  loadings <- P
  return(loadings)
}

RESIDUAL <- function(DATA, SCORES, LOADINGS){
  XHAT <- SCORES%*%t(LOADINGS)
  res <- sum(rowSums((XHAT-DATA)^2))
  
  return(res)
}

CCrESEM <- function(DATA, R, CARD, MaxIter, eps) {
  J <- dim(DATA)[2]
  I <- dim(DATA)[1]
  convAO <- 0
  iter <- 1
  Lossc <- 1
  Lossvec <- Lossc        #Used to compute convergence criterium
  loadings <- INITLOADINGS(DATA, R, CARD)
  A <- loadings
  scores <- matrix(nrow = I, ncol = R)
  if (sum(length(CARD)) > 1){
    CARDc <- J-CARD #Number of loadings -per PC- that have to be zero
    while (convAO == 0) {
      #1. Update component scores
      XP <- DATA %*% loadings/ (sqrt(I))
      scores <- sqrt(I) * svd(XP, R, R)$u %*% t(svd(XP, R, R)$v)
      A <- t(DATA) %*% scores/ I
      #Calculate loss
      Loss <- RESIDUAL(DATA,scores,loadings)
      Lossvec <- c(Lossvec,Loss)
      #2. Update loadings
      loadings <- A
      for (r in 1:R){
        ind <- sort(abs(A[,r]), index.return = TRUE)
        loadings[ind$ix[1:CARDc[r]],r] <- 0
      } 
      #Calculate loss
      Lossu <- RESIDUAL(DATA,scores,loadings)
      Lossvec <- c(Lossvec,Lossu)
      if (iter > MaxIter) {
        convAO <- 1
      }
      if (abs(Lossc-Lossu) < eps) {
        convAO <- 1
      }
      iter <- iter + 1
      Lossc <- Lossu
    }} else {
      CARDc <- J*R-CARD #Number of loadings -over all PCs- that have to be zero
      while (convAO == 0) {
        #1. Update component scores
        XP <- DATA %*% loadings / (sqrt(I))
        scores <- sqrt(I) * svd(XP, R, R)$u %*% t(svd(XP, R, R)$v)
        A <- t(DATA) %*% scores / I
        #Calculate loss
        Loss <- RESIDUAL(DATA,scores,loadings)
        Lossvec <- c(Lossvec,Loss)
        #2. Update loadings
        loadings <- A
        ind <- sort(abs(A), index.return = TRUE)
        loadings[ind$ix[1:CARDc]] <- 0
        #Calculate loss
        Lossu <- RESIDUAL(DATA,scores,loadings)
        Lossvec <- c(Lossvec,Lossu)
        if (iter > MaxIter ) {
          convAO <- 1
        }
        if (abs(Lossc-Lossu) < eps) {
          convAO <- 1
        }
        iter <- iter + 1
        Lossc <- Lossu
      }
    }
  uslpca <- list('scores' = scores, 'loadings' = loadings, 'Lossvec' = Lossvec, 'Residual' = Loss)
  return(uslpca)
}