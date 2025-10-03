## First data generation (de Rooij) (correlated factor) ## 

# Factor correction
source("function/factor_correction.R")
#Criteria
source("function/PL_criteria.R")
# SAM
source("function/factor9/3.1_SAM_PartitionedEstimation.R")
source("function/factor9/3.2_SAM_ComputeSummaryStats.R")
source("function/factor9/3.3_SAM_RegressionCoefs.R")
# rESEM
source("function/factor9/5_rESEM_function.R")
source("function/factor9/5.2_rESEMl1_undoshrinkage.R")
source("function/factor9/5.1_rESEM_l1.R")
# PLS_prediction
source("function/factor9/7.1_simplePLS.R")
source("function/factor9/7.2_PLSpredict.R")
source("function/factor9/7.3_get_intercepts.R")
#-------------------------------------------------------------------------------
# Data generation
#-------------------------------------------------------------------------------
rep = 1
ss = 4
m = 2
s = 2

measurement = c("M:weak", "M:strong")
structural = c("S:weak", "S:strong")
skewtheta = FALSE
direct = FALSE
samplesize <- c(100,200,500,1000)
repetitions <- 30
repi <- 1:30
factor_N <- 9
item_N <- 4 
item_total <- 36

criteria <- lapply(repi, function(r) {
  lapply(measurement, function(m) {
    lapply(structural, function(s) {
      setNames(vector("list", length(samplesize)), as.character(samplesize))
    }) |> setNames(structural)
  }) |> setNames(measurement)
}) |> setNames(paste0("rep", repi))

# variance covariance matrix of latent variables from CERQ study
vcov.theta = matrix(c(
  1.000, 0.467, 0.425, 0.153, 0.465, 0.455, 0.445, 0.283, 0.101,
  0.467, 1.000, 0.471, 0.356, 0.493, 0.492, 0.513, 0.355, 0.256,
  0.425, 0.471, 1.000, 0.109, 0.652, 0.418, 0.142, 0.531, 0.327,
  0.153, 0.356, 0.109, 1.000, 0.324, 0.416, 0.564, 0.022, 0.063,
  0.465, 0.493, 0.652, 0.324, 1.000, 0.811, 0.430, 0.104, 0.114,
  0.455, 0.492, 0.418, 0.416, 0.811, 1.000, 0.685, -0.145, -0.070,
  0.445, 0.513, 0.142, 0.564, 0.430, 0.685, 1.000, -0.070, 0.060,
  0.283, 0.355, 0.531, 0.022, 0.104, -0.145, -0.070, 1.000, 0.664,
  0.101, 0.256, 0.327, 0.063, 0.114, -0.070, 0.060, 0.664, 1.000
), 9, 9)


signbetas = sign(matrix(c( 2.712,
                           -0.136,
                           0.048,
                           0.513,
                           0.525,
                           -2.671,
                           0.042,
                           5.324,
                           -0.461), 9, 1))


for(rep in 1:repetitions){
  for(ss in 1:4){
    for(m in 1:2){
      for(s in 1:2){
        cat("This is repetition:", rep, "from", repetitions, "\n")
        
        n2 = samplesize[ss] # training sample (ss = 100, 200, 500, 1000)
        n1 = 1000 # test sample (1000)
        meas = measurement[m]
        struc = structural[s]
        
        n = n1 + n2 #(ss + 1000)
        idx1 = 1:n1 # test data index
        idx2 = (n1 + 1): n # training sample index
        
        #-------------------------------------------------------------------------------
        ## First data generation (de Rooij) (correlated factor) ## 
        set.seed(0922)
        vareta = 1.1
        
        while(vareta > 1){
          # latent variables
          theta = mvtnorm::rmvnorm(n, rep(0,factor_N), vcov.theta)
          if(skewtheta){
            theta = scale(exp(theta)) # skew distribution of latent variables
          }
          
          #indicators 
          X = matrix(NA, n, item_total)
          #Y = matrix(NA, n, 16)
          
          # can adjust factor loadings for strong/wide/weak measurement model
          if(meas == "M:strong"){
            lambda = runif(item_total, min = 0.5, max = 0.8)
          }
          else if(meas == "M:weak"){
            lambda = runif(item_total, min = 0.2, max = 0.5)
          }
          
          # get different format loadings (with 0)
          lambda1 <- matrix(lambda, nrow = item_N, ncol = factor_N)
          eta <- lambda1
          p <- nrow(eta)
          k <- ncol(eta)
          P_loadings <- sapply(1:k, function(j) {
            c(rep(0, (j-1)*p), eta[, j], rep(0, (k-j)*p))
          })
          P_loadings <- matrix(P_loadings, nrow = p*k, ncol = k,
                               dimnames = list(paste0("row", 1:(p*k)), colnames(eta)))
          
          # can adjust regression coefficients
          if(!direct){
            if(struc == "S:strong"){
              beta = signbetas * matrix(runif(factor_N, min = 0.25, max = 0.40), factor_N, 1)
            }
            else if(struc == "S:weak"){
              beta = signbetas * matrix(runif(factor_N, min = 0.15, max = 0.25), factor_N, 1)
            }
          }
          else{
            beta = signbetas * matrix(runif(factor_N, min = 0.15, max = 0.25), factor_N, 1)
          }
          
          # means
          mx = runif(item_total, min = 1.50, max = 3.0)
          #my = runif(16, min = 1.15, max = 2.15)
          
          j = 0
          for(f in 1:factor_N){
            for(r in 1:item_N){
              j = j + 1
              X[ , j] = mx[j] + lambda[j] * theta[, f] + rnorm(n, mean = 0, sd = sqrt(1 - lambda[j]^2))
            }
          }
          
          if(!direct){
            eta = theta %*% beta
          }
          # with direct effects the differences in the two S conditions are about the strength of the direct effects
          else{
            dbetas = matrix(0, ncol(X), 1)
            if(struc == "S:strong"){
              dbetas[seq(3, 36, by = 4), 1] = (runif(9, min = 0.15, max = 0.25) * sample(c(-1,1), 9, replace = TRUE))
            }
            else if(struc == "S:weak"){
              dbetas[seq(3, 36, by = 4), 1] = (runif(9, min = 0.10, max = 0.15) * sample(c(-1,1), 9, replace = TRUE))
            }
            eta = theta %*% beta + X %*% dbetas
          }
          vareta = var(eta)
          #cat("the variance is", vareta, "\n")
        }
        
        Y = eta + rnorm(n, mean = 0, sd = sqrt(1 - vareta))
        
        #lambday = runif(16, min = 0.2, max = 0.7)
        
        #for(r in 1:16){
        #  Y[ , r] = my[r] + lambday[r] * eta + rnorm(n, mean = 0, sd = sqrt(1 - lambday[r]^2))
        #}
        
        #-------------------------------------------------------------------------------
       
        data = data.frame(cbind(X,Y))
        colnames(data) = c(paste("x", c(1:item_total), sep = ""), "y")
        # Split data
        n <- nrow(data)
        ## set as data.frame
        X <- as.data.frame(X)
        y <- as.data.frame(Y)
        ## Training data
        X_train <- X[idx2, ] 
        Y_train <- Y[idx2, ] 
        ## Test data
        X_test <- X[idx1, ]
        Y_test <- Y[idx1, ]
        ## Combine
        df_train <- cbind(X_train, Y_train) |> setNames(c(paste0("x", 1:item_total), "y"))
        df_test <- cbind(X_test, Y_test) |> setNames(c(paste0("x", 1:item_total), "y"))
        
        #-------------------------------------------------------------------------------
        # Criteria
        #-------------------------------------------------------------------------------
        
        methods <- c("SumScore_E", 
                     "SEM_Reg", 
                     "SAM",
                     "SGCCA","rESEM", 
                     "SEM_BASED", "PLS", 
                     "GLM", "elastic")
        # Prediction statistics: MAE, RMSE, out-of-sample R2, correlation
        MAE <- matrix(NA, nrow = 1, ncol = length(methods))
        colnames(MAE) <- methods
        RMSE <- matrix(NA, nrow = 1, ncol = length(methods))
        colnames(RMSE) <- methods
        OFS <- matrix(NA, nrow = 1, ncol = length(methods))
        colnames(OFS) <- methods
        
        # Measurement model part
        congruence <- matrix(NA, nrow = ncol(P_loadings), ncol = length(methods))
        colnames(congruence) <- methods
        PL_rate <- matrix(NA, nrow = 1, ncol = length(methods))
        colnames(PL_rate) <- methods
        
        # Structure model part
        bias_beta0 <- matrix(NA, nrow = 1, ncol = length(methods))
        colnames(bias_beta0) <- methods
        bias_beta <- matrix(NA, nrow = length(beta), ncol = length(methods))
        colnames(bias_beta) <- methods
        rownames(bias_beta) <- paste0("b",c(1:length(beta)))
        bias_betaALL <- matrix(NA, nrow = 1, ncol = length(methods))
        colnames(bias_betaALL) <- methods
        
        #-------------------------------------------------------------------------------
        # Classical approach
        #-------------------------------------------------------------------------------
        lambda <- lambda1
        sample_size <- n1*2
        beta0 <- 0
        #1. Sum score
        
        ##(1) Empirical method
        ### Step 1. Exploratory factor analysis
        #### lambda for X
        result <- psych::fa(X_train, nfactors = ncol(lambda), rotate='varimax', fm = "minres")# defult is minimum residual factoring
        loadings_EFA <- as.matrix(result$loadings[,1:ncol(lambda)])
        
        ### Step 2. Set the unit weight for the indicators
        #### function for transform the loading matrix
        weight_Emp <- check_lambda(loadings_EFA)
        #### correct directions
        weight_Emp <- check_direction(check_load = weight_Emp, train_score = NA, true_load = P_loadings, factor_n = ncol(lambda))$loadings
        
        ### Step 3: Calculate sum score
        sum_score <- as.matrix(X_train) %*% weight_Emp
        ### Step 4. Estimate the regression coefficients
        reg_mod <- lm(Y_train ~ sum_score[,1] + sum_score[,2] + sum_score[,3] + sum_score[,4] + sum_score[,5] + sum_score[,6] + sum_score[,7] + sum_score[,8] + sum_score[,9])
        reg_coef <- unname(reg_mod$coefficients)
        reg_coef[is.na(reg_coef)] <- 0
        ### Step 5. Calculate the predicted variable based on the test data
        #### test score
        test_score <- as.matrix(X_test) %*% weight_Emp
        Y_hat <- rep(reg_coef[1], length(Y_test)) +
          test_score[,1] * reg_coef[2] +
          test_score[,2] * reg_coef[3] +
          test_score[,3] * reg_coef[4] +
          test_score[,4] * reg_coef[5] +
          test_score[,5] * reg_coef[6] +
          test_score[,6] * reg_coef[7] +
          test_score[,7] * reg_coef[8] +
          test_score[,8] * reg_coef[9] +
          test_score[,9] * reg_coef[10]
        ### Step 6. Criteria
        #### Prediction
        MAE[1, "SumScore_E"] <- sum(abs(Y_test-Y_hat)) / (sample_size/2)
        RMSE[1, "SumScore_E"] <- sqrt(sum((Y_test-Y_hat)^2) / (sample_size/2))
        OFS[1, "SumScore_E"] <- (sum((Y_test-Y_hat)^2) / sum((Y_test)^2))
        
        #### Measurement model part
        congruence[1:ncol(P_loadings), "SumScore_E"] <- diag(abs(psych::factor.congruence(weight_Emp, P_loadings)))
        PL_rate[1, "SumScore_E"] <- PL(P_loadings, weight_Emp)
        
        #### Regression coefficient
        bias_beta0[1, "SumScore_E"] <- abs(reg_coef[1]-beta0)
        bias_beta[, "SumScore_E"] <- abs(reg_coef[-1]-beta)
        bias_betaALL[1, "SumScore_E"] <- mean(c(bias_beta0[1, "SumScore_E"], bias_beta[, "SumScore_E"]))
        
        #-------------------------------------------------------------------------------
        
        #2. CB-SEM, Regression method
        
        colnames(df_train) <- c(paste0("x", 1:36), "y")
        colnames(df_test) <- c(paste0("x", 1:36), "y")
        
        ## Step 1. Set the model
        SEM_Model <- '
          # Measurement model
          l1 =~ x1 + x2 + x3 + x4
          l2 =~ x5 + x6 + x7 + x8
          l3 =~ x9 + x10 + x11 + x12
          l4 =~ x13 + x14 + x15 + x16
          l5 =~ x17 + x18 + x19 + x20
          l6 =~ x21 + x22 + x23 + x24
          l7 =~ x25 + x26 + x27 + x28
          l8 =~ x29 + x30 + x31 + x32
          l9 =~ x33 + x34 + x35 + x36

          # Structural model
          y ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 + l8 + l9

           # Fix variance to 1
          l1 ~~ 1 * l1
          l2 ~~ 1 * l2
          l3 ~~ 1 * l3
          l4 ~~ 1 * l4
          l5 ~~ 1 * l5
          l6 ~~ 1 * l6
          l7 ~~ 1 * l7
          l8 ~~ 1 * l8
          l9 ~~ 1 * l9

        '
        ## Step 2. Estimation of CB-SEM
        fit <- lavaan::sem(model= SEM_Model, data = df_train, meanstructure = T, warn = FALSE)
        
        if (lavaan::lavInspect(fit, "converged") == F){
          MAE[1, "SEM_Reg"] <- NA
          RMSE[1, "SEM_Reg"] <- NA
          OFS[1, "SEM_Reg"] <- NA
          
          #### Measurement model part
          congruence[1:ncol(P_loadings), "SEM_Reg"] <- NA
          PL_rate[1, "SEM_Reg"] <- NA
          
          #### Regression coefficient
          bias_beta0[1, "SEM_Reg"] <- NA
          bias_beta[, "SEM_Reg"] <- NA
          bias_betaALL[1, "SEM_Reg"] <- NA
        }else{
          reg_coef <- c(lavaan::coef(fit)["y~1"], lavaan::coef(fit)[paste0("y~l",c(1:9))]) # regression coefficient
          
          ## Step 3. Get the factor score of test data
          test_score1 <- lavaan::lavPredict(fit, newdata = df_test, method = "regression", transform = TRUE)
          test_score1 <- lavaan::lavPredict(fit, newdata = df_test, method = "Bartlett", transform = TRUE)
          res_loading <- lavaan::inspect(fit, "est")$lambda
          res_loading <- res_loading[1:36,1:9]
          #cov_matrix_factors <- cov(test_score1)
          
          ### Step 4. Calculate the predicted variable
          Y_hat <- reg_coef[1] + test_score1 %*% matrix(reg_coef[-1], factor_N, 1)
          ### Step 5. Criteria
          #### Prediction
          MAE[1, "SEM_Reg"] <- sum(abs(Y_test-Y_hat)) / (sample_size/2)
          RMSE[1, "SEM_Reg"] <- sqrt(sum((Y_test-Y_hat)^2) / (sample_size/2))
          OFS[1, "SEM_Reg"] <- (sum((Y_test-Y_hat)^2) / sum((Y_test)^2))
          
          #### Measurement model part
          congruence[1:ncol(P_loadings), "SEM_Reg"] <- diag(abs(psych::factor.congruence(res_loading, P_loadings)))
          PL_rate[1, "SEM_Reg"] <- PL(P_loadings, res_loading)
          
          #### Regression coefficient
          bias_beta0[1, "SEM_Reg"] <- abs(reg_coef[1]-beta0)
          bias_beta[, "SEM_Reg"] <- abs(reg_coef[-1]-beta)
          bias_betaALL[1, "SEM_Reg"] <- mean(c(bias_beta0[1, "SEM_Reg"], bias_beta[, "SEM_Reg"]))
        }
        
        
        #-------------------------------------------------------------------------------
        # SEM Based (Two-stage)
        #-------------------------------------------------------------------------------
        
        #3. Structural After Measurement (SAM)
        ## step 0. data name
        colnames(df_train) <- c(unlist(lapply(letters[1:9], function(l) paste0(l, 1:4))), "y")
        
        ## Step 1. Estimate measurement model with CFA
        SAM_Model <- '
  # Measurement model
  Fa =~ a1 + a2 + a3 + a4
  Fb =~ b1 + b2 + b3 + b4
  Fc =~ c1 + c2 + c3 + c4
  Fd =~ d1 + d2 + d3 + d4
  Fe =~ e1 + e2 + e3 + e4
  Ff =~ f1 + f2 + f3 + f4
  Fg =~ g1 + g2 + g3 + g4
  Fh =~ h1 + h2 + h3 + h4
  Fi =~ i1 + i2 + i3 + i4

  # Structural model
  Fy ~ Fa + Fb + Fc + Fd + Fe + Ff + Fg + Fh + Fi

  # Fix latent variances to 1
  Fa ~~ 1*Fa
  Fb ~~ 1*Fb
  Fc ~~ 1*Fc
  Fd ~~ 1*Fd
  Fe ~~ 1*Fe
  Ff ~~ 1*Ff
  Fg ~~ 1*Fg
  Fh ~~ 1*Fh
  Fi ~~ 1*Fi
'
        ### Compute loadings and residual variances for each measurement block separately
        # sub.est <- partitioned.estimation(data = df_train, method = "SAM_Bentler_ULS")
        mm_est <- partitioned.estimation(data = df_train, method = "SAM_MLB")
        #### loadings
        res_loading <- mm_est$LAMBDA[1:36, 1:9]
        res_loading <- check_direction(check_load = res_loading, train_score = NA, true_load = P_loadings, factor_n = ncol(P_loadings))$loadings
        #### residuals
        res_residual <- mm_est$THETA[1:36,1:36]
        if (all(colSums(res_residual != 0) > 0) == F){
          MAE[1, "SEM_Reg"] <- NA
          RMSE[1, "SEM_Reg"] <- NA
          OFS[1, "SEM_Reg"] <- NA
          
          #### Measurement model part
          congruence[1:ncol(P_loadings), "SEM_Reg"] <- NA
          PL_rate[1, "SEM_Reg"] <- NA
          
          #### Regression coefficient
          bias_beta0[1, "SEM_Reg"] <- NA
          bias_beta[, "SEM_Reg"] <- NA
          bias_betaALL[1, "SEM_Reg"] <- NA
        }else{
          ### Obtain estimated variance–covariance matrix
          covariance_matrix <- cov(df_train) * (nrow(df_train) - 1) / nrow(df_train)
          ### Compute summary statistics E(eta) and var(eta)
          summ.stats <- LSAM_SumStats(S = covariance_matrix,
                                      sample.nobs = sample_size,
                                      LAMBDA = mm_est$LAMBDA,
                                      THETA = mm_est$THETA,
                                      mapping = "ML")
          ## Step 2. Estimate structural model (based on the variance-covariance matrix of the factor scores)
          reg_coef <- LSAM_regcoef(model = SAM_Model, sumstats = summ.stats)
          ### calculate beta0
          # mu_Score <- colMeans(t(theta))
          # mu_Y <- colMeans(Y_train)
          # reg_coef <-  c(mu_Y - sum(mu_Score * reg_coef), reg_coef)
          reg_coef <- c(0, reg_coef)
          ## Step 3. Get the factor score (Bartlett method) of test data
          test_score <- t(solve(res_residual) %*% res_loading %*% solve(t(res_loading) %*% solve(res_residual) %*% res_loading)) %*% t(X_test)
          test_score <- t(test_score)
          
          ### Step 4. Calculate the predicted variable
          Y_hat <- reg_coef[1] + test_score %*% matrix(reg_coef[-1], factor_N, 1)
          ### Step 5. Criteria
          #### Prediction
          MAE[1, "SAM"] <- sum(abs(Y_test - Y_hat)) / (sample_size / 2)
          RMSE[1, "SAM"] <- sqrt(sum((Y_test - Y_hat)^2) / (sample_size / 2))
          OFS[1, "SAM"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)
          
          #### Measurement model part
          congruence[1:ncol(P_loadings), "SAM"] <- diag(abs(psych::factor.congruence(res_loading, P_loadings)))
          PL_rate[1, "SAM"] <- PL(P_loadings, res_loading)
          
          #### Regression coefficient
          bias_beta0[1, "SAM"] <- abs(reg_coef[1]-beta0)
          bias_beta[, "SAM"] <- abs(reg_coef[-1]-beta)
          bias_betaALL[1, "SAM"] <- mean(c(bias_beta0[1, "SAM"], bias_beta[, "SAM"]))
        }
        
        
        #-------------------------------------------------------------------------------
        
        #4. Sparse generalized canonical correlation analysis (SGCCA)
        ## Step 0. Transform data
        colnames(df_train) <- c(paste0("x", 1:36), "y")
        colnames(df_test) <- c(paste0("x", 1:36), "y")
        
        x_train <- df_train[, 1:4]
        z_train<- df_train[, 5:8]
        a_train <- df_train[, 9:12] 
        b_train <- df_train[, 13:16] 
        c_train <- df_train[, 17:20] 
        d_train <- df_train[, 21:24] 
        e_train <- df_train[, 25:28] 
        f_train <- df_train[, 29:32] 
        g_train <- df_train[, 33:36] 
        y_train <- df_train[, 37]
        
        x_test <- df_test[, 1:4]
        z_test <- df_test[, 5:8]
        a_test <- df_test[, 9:12] 
        b_test <- df_test[, 13:16] 
        c_test <- df_test[, 17:20] 
        d_test <- df_test[, 21:24] 
        e_test <- df_test[, 25:28] 
        f_test <- df_test[, 29:32] 
        g_test <- df_test[, 33:36] 
        y_test <- df_test[, 37]
        
        
        ## Step 1. Estimate measurement model
        ### Settings
        blocks <- list(
          bl1  = x_train,
          bl2  = z_train,
          bl3  = a_train,
          bl4  = b_train,
          bl5  = c_train,
          bl6  = d_train,
          bl7  = e_train,
          bl8  = f_train,
          bl9  = g_train,
          bl10 = y_train
        )
        
        connection <- matrix(0, 10, 10)
        connection[1:9, 10] <- 1
        connection[10, 1:9] <- 1
        diag(connection) <- 0
        
        cv <- RGCCA::rgcca_cv(
          blocks            = blocks,
          connection        = connection,
          method            = "sgcca",
          response          = 10,              
          par_type          = "sparsity",      
          par_length        = 10,              # 10 random candidates/block in [0.5, 1]
          ncomp             = 1,               # 1 component/block 
          scheme            = "centroid",
          scale             = TRUE,
          scale_block       = TRUE,
          validation        = "kfold",         
          k                 = 5,               # 5-fold CV
          prediction_model  = "glmnet",        
          metric            = "MAE")     
        
        sparsity <- c(cv$best_params)
        ### Run SGCCA
        fit_sgcca <- RGCCA::rgcca(blocks = blocks,
                                  connection = connection,
                                  method = "sgcca",
                                  sparsity = sparsity,
                                  # each block remains separate rather than combining into a meta-block
                                  superblock = FALSE, 
                                  # Extracts one component per block
                                  ncomp = rep(1,10), 
                                  scheme = "factorial", 
                                  # extracted components are orthogonal (uncorrelated)
                                  comp_orth = TRUE, 
                                  verbose = F)
        ### Get the weights
        w_x <- fit_sgcca$astar$bl1
        w_z <- fit_sgcca$astar$bl2
        w_a <- fit_sgcca$astar$bl3
        w_b <- fit_sgcca$astar$bl4
        w_c <- fit_sgcca$astar$bl5
        w_d <- fit_sgcca$astar$bl6
        w_e <- fit_sgcca$astar$bl7
        w_f <- fit_sgcca$astar$bl8
        w_g <- fit_sgcca$astar$bl9
        ### get the full weight matrix
        w_list <- list(w_x, w_z, w_a, w_b, w_c, w_d, w_e, w_f, w_g)
        W <- matrix(0, nrow = nrow(P_loadings), ncol = ncol(P_loadings),
                    dimnames = dimnames(P_loadings))
        for (j in seq_len(ncol(P_loadings))) {
          idx <- which(P_loadings[, j] != 0)
          W[idx, j] <- w_list[[j]]
        }
        
        ## Step 2. Calculate the factor score for training data
        train_score_x <- as.matrix(x_train) %*% w_x
        train_score_z <- as.matrix(z_train) %*% w_z
        train_score_a <- as.matrix(a_train) %*% w_a
        train_score_b <- as.matrix(b_train) %*% w_b
        train_score_c <- as.matrix(c_train) %*% w_c
        train_score_d <- as.matrix(d_train) %*% w_d
        train_score_e <- as.matrix(e_train) %*% w_e
        train_score_f <- as.matrix(f_train) %*% w_f
        train_score_g <- as.matrix(g_train) %*% w_g
        ## Step 3. Estimate structural model
        train_df <- data.frame(
          Y_train,
          train_score_x, train_score_z,
          train_score_a, train_score_b, train_score_c,
          train_score_d, train_score_e, train_score_f, train_score_g
        )
        res_reg <- lm(Y_train ~ ., data = train_df)
        reg_coef <- unname(res_reg$coefficients)
        ## Step 4. Calculate the factor score for test data
        test_score_x <- as.matrix(x_test) %*% w_x
        test_score_z <- as.matrix(z_test) %*% w_z
        test_score_a <- as.matrix(a_test) %*% w_a
        test_score_b <- as.matrix(b_test) %*% w_b
        test_score_c <- as.matrix(c_test) %*% w_c
        test_score_d <- as.matrix(d_test) %*% w_d
        test_score_e <- as.matrix(e_test) %*% w_e
        test_score_f <- as.matrix(f_test) %*% w_f
        test_score_g <- as.matrix(g_test) %*% w_g
        ## Step 5. Calculate the predicted variable
        Y_hat <- rep(reg_coef[1], length(Y_test)) +
          test_score_x * reg_coef[2] +
          test_score_z * reg_coef[3] +
          test_score_a * reg_coef[4] +
          test_score_b * reg_coef[5] +
          test_score_c * reg_coef[6] +
          test_score_d * reg_coef[7] +
          test_score_e * reg_coef[8] +
          test_score_f * reg_coef[9] +
          test_score_g * reg_coef[10]
        
        ## Step 6. Criteria
        ### Prediction
        MAE[1, "SGCCA"] <- sum(abs(Y_test - Y_hat)) / (sample_size / 2)
        RMSE[1, "SGCCA"] <- sqrt(sum((Y_test - Y_hat)^2) / (sample_size / 2))
        OFS[1, "SGCCA"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)
        
        #### Measurement model part
        congruence[1:ncol(P_loadings), "SGCCA"] <- diag(abs(psych::factor.congruence(W, P_loadings)))
        PL_rate[1, "SGCCA"] <- PL(P_loadings, W)
        
        #### Regression coefficient
        bias_beta0[1, "SGCCA"] <- abs(reg_coef[1]-beta0)
        bias_beta[, "SGCCA"] <- abs(reg_coef[-1]-beta)
        bias_betaALL[1, "SGCCA"] <- mean(c(bias_beta0[1, "SGCCA"], bias_beta[, "SGCCA"]))
        
        
        #-------------------------------------------------------------------------------
        
        #5. Regularized Exploraty Structural Equation Modeling (rESEM) (cardinality constrains)
        X_train_s <- scale(X_train)
        
        ## Step 1. Estimate the measurement part
        # res_method <- MULTISTART_CCrESEM(DATA = as.matrix(X_train_s), R = 9, CARD = rep(4,9), MaxIter = 100, eps = 10^-4, nstarts = 20)
        
        lasso_start <- c(39, 64, 120, 250)
        lasso <- findLasso(dat = as.matrix(X_train_s), Ptrue = P_loadings, maxItr = 20, lassou = lasso_start[ss])
        res_method <- MULTISTART_rESEMl1(DATA = as.matrix(X_train_s), R = 9, MaxIter = 20, eps = 10^-6, nstarts = 50, lambda = lasso$lasso)
        res_method <- rESEMl1_undoshrinkage(DATA = as.matrix(X_train_s), R = 9, P = res_method$loadings, MaxIter = 100, eps = 10^-6)
        #res_method <- MULTISTART_rESEMl1(DATA = as.matrix(X_train_s), R = 9, MaxIter = 20, eps = 10^-6, nstarts = 50, lambda = 245)
        #sum(round(res_method$loadings,3) != 0)
        ### get the scores and loadings
        train_score <- res_method$scores
        res_loading <- res_method$loadings
        ### correct directions
        res_direction <- check_direction(check_load = res_loading, train_score = train_score, true_load = P_loadings, factor_n = ncol(P_loadings))
        res_loading <- res_direction$loadings
        train_score_T <- res_direction$score
        
        ## Step 2. Estimate the structural part
        res_reg <- lm(Y_train ~  train_score_T[,1] + train_score_T[,2] + 
                        train_score_T[,3] + train_score_T[,4] + 
                        train_score_T[,5] + train_score_T[,6] +
                        train_score_T[,7] + train_score_T[,8] +
                        train_score_T[,9] )
        reg_coef <- unname(res_reg$coefficients)
        
        ## Step 3. Calculate the factor scores
        
        test_score <- as.matrix(scale(X_test)) %*% res_loading %*% solve(t(res_loading) %*% res_loading)
        ## Step 4. Calculate the predicted variable
        Y_hat <- reg_coef[1] + test_score %*% matrix(reg_coef[-1], factor_N, 1)
        
        ## Step 5. Criteria
        ### Prediction
        MAE[1, "rESEM"] <- sum(abs(Y_test - Y_hat)) / (sample_size / 2)
        RMSE[1, "rESEM"] <- sqrt(sum((Y_test - Y_hat)^2) / (sample_size / 2))
        OFS[1, "rESEM"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)
        
        #### Measurement model part
        congruence[1:ncol(P_loadings), "rESEM"] <- diag(abs(psych::factor.congruence(res_loading, P_loadings)))
        PL_rate[1, "rESEM"] <- PL(P_loadings, res_loading)
        
        #### Regression coefficient
        bias_beta0[1, "rESEM"] <- abs(reg_coef[1]-beta0)
        bias_beta[, "rESEM"] <- abs(reg_coef[-1]-beta)
        bias_betaALL[1, "rESEM"] <- mean(c(bias_beta0[1, "rESEM"], bias_beta[, "rESEM"]))
        
        
        #-------------------------------------------------------------------------------
        # SEM Based (One-stage)
        #-------------------------------------------------------------------------------
        
        #6. SEM-Based Prediction Rule
        colnames(df_train) <- c(paste0("x", 1:36), "y")
        colnames(df_test) <- c(paste0("x", 1:36), "y")
        ## Step 1. Build the model
        SEMRule_Model <- '
  # Measurement model
  l1 =~ x1 + x2 + x3 + x4
  l2 =~ x5 + x6 + x7 + x8
  l3 =~ x9 + x10 + x11 + x12
  l4 =~ x13 + x14 + x15 + x16
  l5 =~ x17 + x18 + x19 + x20
  l6 =~ x21 + x22 + x23 + x24
  l7 =~ x25 + x26 + x27 + x28
  l8 =~ x29 + x30 + x31 + x32
  l9 =~ x33 + x34 + x35 + x36
  
  # Structural model
  y ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 + l8 + l9
  
  # Fix variance to 1
  l1 ~~ 1 * l1
  l2 ~~ 1 * l2
  l3 ~~ 1 * l3
  l4 ~~ 1 * l4
  l5 ~~ 1 * l5
  l6 ~~ 1 * l6
  l7 ~~ 1 * l7
  l8 ~~ 1 * l8
  l9 ~~ 1 * l9
'
        ## Step 2. Estimation of CB-SEM
        fit <- lavaan::sem(SEMRule_Model , data = df_train , meanstructure = TRUE, warn = FALSE)
        reg_coef <- c(lavaan::coef(fit)["y~1"], lavaan::coef(fit)[paste0("y~l",c(1:9))]) # regression coefficient
        res_loading <- lavaan::inspect(fit, "est")$lambda
        res_loading <- res_loading[1:36,1:9]
        ## Step 3. Calculate the predicted variable
        Y_hat <- lavaan::lavPredictY(fit, newdata = df_test, ynames = "y", xnames = c(paste("x", c(1:36), sep = "")))
        ### Step 4. Criteria
        #### Prediction
        MAE[1, "SEM_BASED"] <- sum(abs(Y_test - Y_hat)) / (sample_size / 2)
        RMSE[1, "SEM_BASED"] <- sqrt(sum((Y_test - Y_hat)^2) / (sample_size / 2))
        OFS[1, "SEM_BASED"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)
        
        #### Measurement model part
        congruence[1:ncol(P_loadings), "SEM_BASED"] <- diag(abs(psych::factor.congruence(res_loading, P_loadings)))
        PL_rate[1, "SEM_BASED"] <- PL(P_loadings, res_loading)
        
        #### Regression coefficient
        bias_beta0[1, "SEM_BASED"] <- abs(reg_coef[1]-beta0)
        bias_beta[, "SEM_BASED"] <- abs(reg_coef[-1]-beta)
        bias_betaALL[1, "SEM_BASED"] <- mean(c(bias_beta0[1, "SEM_BASED"], bias_beta[, "SEM_BASED"]))
        
        
        #-------------------------------------------------------------------------------
        
        #7. PLS-SEM
        colnames(df_train) <- c(paste0("x", 1:36), "y")
        colnames(df_test) <- c(paste0("x", 1:36), "y")
        ## Step 1. Create matrix of the measurement and structural Model
        smMatrix <- matrix(c("X","Y",
                             "Z","Y",
                             "A","Y",
                             "B","Y",
                             "C","Y",
                             "D","Y",
                             "E","Y",
                             "F","Y",
                             "G","Y"),
                           nrow = 9, ncol = 2, byrow = T, 
                           dimnames = list(1:9,c("source", "target")))
        mmMatrix <- matrix(c("X","x1","R",
                             "X","x2","R",
                             "X","x3","R",
                             "X","x4","R",
                             "Z","z1","R",
                             "Z","z2","R",
                             "Z","z3","R",
                             "Z","z4","R",
                             "A","a1","R",
                             "A","a2","R",
                             "A","a3","R",
                             "A","a4","R",
                             "B","b1","R",
                             "B","b2","R",
                             "B","b3","R",
                             "B","b4","R",
                             "C","c1","R",
                             "C","c2","R",
                             "C","c3","R",
                             "C","c4","R",
                             "D","d1","R",
                             "D","d2","R",
                             "D","d3","R",
                             "D","d4","R",
                             "E","e1","R",
                             "E","e2","R",
                             "E","e3","R",
                             "E","e4","R",
                             "F","f1","R",
                             "F","f2","R",
                             "F","f3","R",
                             "F","f4","R",
                             "G","g1","R",
                             "G","g2","R",
                             "G","g3","R",
                             "G","g4","R",
                             "Y","y","R"),
                           nrow = 37, ncol = 3, byrow = T,
                           dimnames = list(1:37,c("latent","measurement","type")))
        ## Step 2. Prediction on the model
        NameVariables <- mmMatrix[,"measurement"]
        colnames(df_train) <- NameVariables
        colnames(df_test) <- NameVariables
        predTrain <- PLSpredict(df_train, df_test, smMatrix, mmMatrix, maxIt = 300, stopCriterion = 9)
        Y_hat <- predTrain$predictedMeasurements
        # get loadings and regression coefficients
        PLS_sim <- simplePLS(df_train, smMatrix, mmMatrix, maxIt = 300, stopCriterion = 9)
        PLS_loadings <- PLS_sim$outer_loadings[-37,-10]
        # get not standardized regression coefficient
        reg_coef <- get_intercepts(PLS_sim)
        ## Step 3. Criteria
        ### Prediction
        MAE[1, "PLS"] <- sum(abs(Y_test - Y_hat)) / (sample_size / 2)
        RMSE[1, "PLS"] <- sqrt(sum((Y_test - Y_hat)^2) / (sample_size / 2))
        OFS[1, "PLS"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)
        
        #### Measurement model part
        congruence[1:ncol(P_loadings), "PLS"] <- diag(psych::factor.congruence(PLS_loadings[1:nrow(P_loadings), 1:9], P_loadings))
        PL_rate[1, "PLS"] <- PL(P_loadings, PLS_loadings)
        
        ### Regression coefficient
        bias_beta0[1, "PLS"] <- abs(reg_coef[1]-beta0)
        bias_beta[, "PLS"] <- abs(reg_coef[-1]-beta)
        bias_betaALL[1, "PLS"] <- mean(c(bias_beta0[1, "PLS"], bias_beta[, "PLS"]))
        
        
        #-------------------------------------------------------------------------------
        # Machine learning approaches
        #-------------------------------------------------------------------------------
        
        #8. Linear regression
        colnames(df_train) <- c(paste0("x", 1:36), "y")
        colnames(df_test) <- c(paste0("x", 1:36), "y")
        ## Step 1. Estimate linear model based on training data
        Lmodel <- lm(y ~., data = df_train)
        ### Get the regression coefficients
        reg_coef <- Lmodel$coefficients
        ## Step 2. Calculate the predicted variable
        X_test_int <- cbind(1, as.matrix(X_test))
        Y_hat <- X_test_int %*% reg_coef
        # Y_hat <- predict(Lmodel, newdata = df_test)
        
        ## Step 3. Criteria
        ### Prediction
        MAE[1, "GLM"] <- sum(abs(Y_test - Y_hat)) / (sample_size / 2)
        RMSE[1, "GLM"] <- sqrt(sum((Y_test - Y_hat)^2) / (sample_size / 2))
        OFS[1, "GLM"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)
        
        #-------------------------------------------------------------------------------
        
        #9. Elastic net 
        ## Step 1. Data transform
        X_train <- as.matrix(X_train)
        Y_train <- as.matrix(Y_train)
        X_test <- as.matrix(X_test)
        ## Step 2. Training data
        ### Choose the alpha through cross-validation
        alphas <- seq(0, 1, by = 0.1)
        cvs <- rep(NA, 11)
        for (i in 1:length(alphas)) {
          cvs[i] <- min(
            glmnet::cv.glmnet(X_train, Y_train, alpha  = alphas[i])$cvm)
        }
        aa <- which.min(cvs) # the alpha having the smallest output
        ### Choose the lambda through cross-validation
        cv_fit <- glmnet::cv.glmnet(X_train, Y_train, alpha = alphas[aa])
        # Extract lambda
        best_lambda <- cv_fit$lambda.min # the best lambda having the smallest result
        ## Step 3. Fit elestic net
        ### Fit final model on entire training set with that lambda and alpha
        final_model <- glmnet::glmnet(X_train, Y_train, alpha = alphas[aa], lambda = best_lambda)
        ### the regression coefficient 
        reg_coef <- coef(final_model)
        ## Step 4. Calculate the predicted variable
        Y_hat <- predict(final_model, newx = X_test, s = best_lambda)
        
        ## Step 5. Criteria
        ### Prediction
        MAE[1, "elastic"] <- sum(abs(Y_test - Y_hat)) / (sample_size / 2)
        RMSE[1, "elastic"] <- sqrt(sum((Y_test - Y_hat)^2) / (sample_size / 2))
        OFS[1, "elastic"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)
        
        #-------------------------------------------------------------------------------
        # The end
        #-------------------------------------------------------------------------------
        
        criteria[[rep]][[measurement[m]]][[structural[s]]][[as.character(samplesize[ss])]] <- list(
          MAE         = MAE,
          RMSE        = RMSE,
          OFS         = OFS,
          congruence  = congruence,
          PL_rate     = PL_rate,
          beta_mean   = bias_betaALL)
      }
    }
  }
}


rep_ids <- paste0("rep", 1:10)

sample_sizes <- c("100", "200", "500", "1000")

methods <- c("SumScore_E", 
             "SEM_Reg", 
             "SAM",
             "SGCCA","rESEM", 
             "SEM_BASED", "PLS", 
             "GLM", "elastic")

## MAE, RMSE, beta
mean_MAE_tbl <- sapply(sample_sizes, function(n) {
  mat <- sapply(rep_ids, function(r) {
    # this should return a named numeric vector (methods as names)
    criteria[[r]][["M:Strong"]][["S:Strong"]][[n]][["MAE"]]
  })
  rowMeans(as.matrix(mat), na.rm = TRUE)
})
rownames(mean_MAE_tbl) <- c("SumScore_E", 
                            "SEM_Reg",
                            "SAM",
                            "SGCCA","rESEM", 
                            "SEM_BASED", "PLS", 
                            "GLM", "elastic")

## congruence
mean_congruence <- list()
for (ss in sample_sizes) {
  # collect colMeans for each replication
  cong_mat <- sapply(rep_ids, function(r) {
    colMeans(criteria[[r]][["M:weak"]][["S:weak"]][[ss]][["congruence"]],na.rm = T)
  })
  
  # average across reps
  mean_congruence[[ss]] <- rowMeans(cong_mat)
  mean_congruence <- sapply(sample_sizes, function(n) { cbind(mean_congruence[[n]])})
}
rownames(mean_congruence) <- c("SumScore_E", 
                               #"SEM_Reg", 
                               "SGCCA","rESEM", 
                               "SEM_BASED", "PLS", 
                               "GLM", "elastic")


which.max(as.matrix(mean_MAE_tbl), na.rm = TRUE)
write.csv(mean_congruence,"0_congruence_studyWW.csv")
#saveRDS(criteria, "criteria_rep30_0909.rds")
criteria <- readRDS("criteria_rep30_0909.rds")
saveRDS(criteria,"criteria_rep10_0923_deRooij.rds")
