#-------------------------------------------------------------------------------
# Load package
#-------------------------------------------------------------------------------
# "cbsem" package
#url <- "https://cran.r-project.org/src/contrib/Archive/cbsem/cbsem_1.0.0.tar.gz"
#install.packages(url, repos = NULL, type = "source")

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
#install.packages("foreach")
#install.packages("doParallel")
library(foreach)
library(doParallel)

# factor, item settings and sample size
## x facotrs
factor_N <- 9
item_N <- 4
item_total <- factor_N * item_N
## y factor
Y_n <- 1
XY_n <- factor_N + Y_n

# do Parallel
cl <- makeCluster(max(1, parallel::detectCores() - 1))
registerDoParallel(cl)
# replications
rep <- 40
# conditions
measurement = c("M:weak", "M:strong")
structural = c("S:weak", "S:strong")
sample_condition <- c(100,200,500,1000)
# criteria
criteria <- lapply(measurement, function(m) {
  lapply(structural, function(s) {
    setNames(vector("list", length(sample_condition)), as.character(sample_condition))
    }) |> setNames(structural)
  }) |> setNames(measurement)

start <- Sys.time()
# start
result <- foreach(i = 1:rep) %dopar% {
  
  for(ss in 1:4){
    for(m in 1:2){
      for(s in 1:2){
        cat("SampleSize:",ss,"measurement:",m,"structural:",s,"\n")
        # measurement and strucrtural conditions
        meas = measurement[m]
        struc = structural[s]
        # sample size
        train_n <- sample_condition[ss]
        test_n <- 1000
        sample_size <- train_n + test_n
        
        # loadings and betas
        if (meas == "M:weak"){
          P_loadings <- readRDS("loadingsW_0922.rds")
          lambda <- readRDS("lambda1W_0922.rds")
        }else{
          P_loadings <- readRDS("loadings_0922.rds")
          lambda <- readRDS("lambda1_0922.rds")
        }
        if (struc == "S:weak"){
          beta <- readRDS("betaW_0922.rds")
        }else{
          beta <- readRDS("beta_0922.rds")
        }
        beta0 <- 0
        allCoef <- rbind(P_loadings, matrix(beta,1,9))
        cov_matrix <- readRDS("covarianceSS.rds")
        
        # data generation for PLS-SEM based model
        ## Structural model among 10 composites
        B <- matrix(0, nrow = XY_n, ncol = XY_n)
        B[XY_n, 1:factor_N] <- beta
        ## Which composite each indicator belongs to
        indicatorx <- rep(1:factor_N, each = item_N)    # 9*4 X-indicators 
        indicatory <- 1          # 1 Y-indicator (for the single endogenous comp)
        ## Reflective loadings
        lambdax <- c(matrix(unlist(lambda), ncol = 1))
        lambday <- 1
        ## Covariance of exogenous composites (eta1, eta2)
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
        Sxixi <- vcov.theta
        ## covariance matrix of a GSC model
        out <- cbsem::gscmcov(
          B = B,
          indicatorx = indicatorx,
          indicatory = indicatory,
          lambdax = lambdax,
          lambday = lambday,
          Sxixi = Sxixi
        )
        Sigma <- out$S
        ## sample data
        data <- mvtnorm::rmvnorm(sample_size, rep(0, nrow(Sigma)), Sigma)
        colnames(data) = c(paste("x", c(1:item_total), sep = ""), "y")
        X <- data[,1:nrow(P_loadings)]
        Y <- data[,item_total+1]
        
        # split data
        idx_train <- 1:train_n # train index
        idx_test <- (train_n + 1): sample_size # test index
        ## Training data
        X_train <- X[idx_train, ] 
        Y_train <- Y[idx_train] 
        ## Test data
        X_test <- X[idx_test, ]
        Y_test <- Y[idx_test]
        ## Combine
        df_train <- cbind(X_train, Y_train) 
        df_test <- cbind(X_test, Y_test)
        colnames(df_train) <- c(paste0("x", 1:36), "y")
        colnames(df_test) <- c(paste0("x", 1:36), "y")
        
        #-------------------------------------------------------------------------------
        # Criteria
        #-------------------------------------------------------------------------------
        
        methods <- c("SumScore", 
                     "SEM_Reg", "SEM_Bar",
                     "SAM","SAM_Reg",
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
        #1. Sum score
        
        ## Step 1. Exploratory factor analysis
        result <- psych::fa(X_train, 
                            nfactors = factor_N, 
                            rotate='varimax', 
                            fm = "minres") # defult is minimum residual factoring
        loadings_EFA <- as.matrix(result$loadings[,1:factor_N])
        
        ## Step 2. Set the unit weight for the indicators
        ### check_lambda() function for transform the loading matrix
        weight_Emp <- check_lambda(loadings_EFA)
        ### check_direction() correct directions
        weight_Emp <- check_direction(check_load = weight_Emp, 
                                      train_score = NA, 
                                      true_load = P_loadings, 
                                      factor_n = factor_N
        )$loadings
        
        ## Step 3: Calculate sum score
        sum_score <- as.matrix(X_train) %*% weight_Emp
        
        ## Step 4. Estimate the regression coefficients
        reg_mod <- lm(Y_train ~ ., data = as.data.frame(sum_score))
        reg_coef <- unname(reg_mod$coefficients)
        reg_coef[is.na(reg_coef)] <- 0
        
        ## Step 5. Calculate the predicted variable based on the test data
        ### test score
        test_score <- as.matrix(X_test) %*% weight_Emp
        Y_hat <- reg_coef[1] + test_score %*% reg_coef[-1]
        
        ## Step 6. Criteria
        ### Prediction
        MAE[1, "SumScore"] <- sum(abs(Y_test-Y_hat)) / test_n
        RMSE[1, "SumScore"] <- sqrt(sum((Y_test-Y_hat)^2) / test_n)
        OFS[1, "SumScore"] <- (sum((Y_test-Y_hat)^2) / sum((Y_test)^2))
        
        ### Measurement model part
        congruence[1:ncol(P_loadings), "SumScore"] <- diag(abs(psych::factor.congruence(weight_Emp, P_loadings)))
        PL_rate[1, "SumScore"] <- PL(P_loadings, weight_Emp)
        
        ### Regression coefficient
        bias_beta0[1, "SumScore"] <- abs(reg_coef[1]-beta0)
        bias_beta[, "SumScore"] <- abs(reg_coef[-1]-beta)
        bias_betaALL[1, "SumScore"] <- mean(c(bias_beta0[1, "SumScore"], bias_beta[, "SumScore"]))
        
        #-------------------------------------------------------------------------------
        #2. CB-SEM, Regression method
        
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
        fit <- lavaan::sem(model= SEM_Model, 
                           data = df_train, 
                           meanstructure = T, 
                           warn = FALSE)
        reg_coef <- c(lavaan::coef(fit)["y~1"], 
                      lavaan::coef(fit)[paste0("y~l",c(1:factor_N))]) # regression coefficient
        
        ## Step 3. Get the factor score of test data
        res_loading <- lavaan::inspect(fit, "est")$lambda
        res_loading <- res_loading[1:item_total,1:factor_N]
        ### (1) regression
        test_score <- lavaan::lavPredict(fit, newdata = df_test, method = "regression", transform = TRUE)
        ### (2) Bartlett
        test_score1 <- lavaan::lavPredict(fit, newdata = df_test, method = "Bartlett", transform = TRUE)
        
        ## Step 4. Calculate the predicted variable
        ### (1) regression
        Y_hat <- reg_coef[1] + test_score %*% matrix(reg_coef[-1], factor_N, 1)
        ### (2) Bartlett
        Y_hat1 <- reg_coef[1] + test_score1 %*% matrix(reg_coef[-1], factor_N, 1)
        
        ## Step 5. Criteria
        ### Prediction
        #### (1) regression
        MAE[1, "SEM_Reg"] <- sum(abs(Y_test-Y_hat)) / test_n
        RMSE[1, "SEM_Reg"] <- sqrt(sum((Y_test-Y_hat)^2) / test_n)
        OFS[1, "SEM_Reg"] <- (sum((Y_test-Y_hat)^2) / sum((Y_test)^2))
        #### (2) Bartlett
        MAE[1, "SEM_Bar"] <- sum(abs(Y_test-Y_hat1)) / test_n
        RMSE[1, "SEM_Bar"] <- sqrt(sum((Y_test-Y_hat1)^2) / test_n)
        OFS[1, "SEM_Bar"] <- (sum((Y_test-Y_hat1)^2) / sum((Y_test)^2))
        
        ### Measurement model part
        congruence[1:ncol(P_loadings), "SEM_Reg"] <- diag(abs(psych::factor.congruence(res_loading, P_loadings)))
        PL_rate[1, "SEM_Reg"] <- PL(P_loadings, res_loading)
        
        congruence[1:ncol(P_loadings), "SEM_Bar"] <- diag(abs(psych::factor.congruence(res_loading, P_loadings)))
        PL_rate[1, "SEM_Bar"] <- PL(P_loadings, res_loading)
        
        ### Regression coefficient
        bias_beta0[1, "SEM_Reg"] <- abs(reg_coef[1]-beta0)
        bias_beta[, "SEM_Reg"] <- abs(reg_coef[-1]-beta)
        bias_betaALL[1, "SEM_Reg"] <- mean(c(bias_beta0[1, "SEM_Reg"], bias_beta[, "SEM_Reg"]))
        
        bias_beta0[1, "SEM_Bar"] <- abs(reg_coef[1]-beta0)
        bias_beta[, "SEM_Bar"] <- abs(reg_coef[-1]-beta)
        bias_betaALL[1, "SEM_Bar"] <- mean(c(bias_beta0[1, "SEM_Bar"], bias_beta[, "SEM_Bar"]))
        
        #-------------------------------------------------------------------------------
        # SEM Based (Two-stage)
        #-------------------------------------------------------------------------------
        #3. Structural After Measurement (SAM)
        
        ## step 0. data name
        df_trainSAM <- df_train
        colnames(df_trainSAM) <- c(unlist(lapply(letters[1:factor_N], function(l) paste0(l, 1:4))), "y")
        
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
'
        
        ### Compute loadings and residual variances for each measurement block separately
        # sub.est <- partitioned.estimation(data = df_train, method = "SAM_Bentler_ULS")
        mm_est <- partitioned.estimation(data = df_trainSAM, 
                                         method = "SAM_MLB")
        #### loadings
        res_loading <- mm_est$LAMBDA[1:item_total, 1:factor_N]
        res_loading <- check_direction(check_load = res_loading, 
                                       train_score = NA, 
                                       true_load = P_loadings, 
                                       factor_n = ncol(P_loadings)
        )$loadings
        #### residuals
        res_residual <- mm_est$THETA[1:item_total,1:item_total]
        if (all(colSums(res_residual != 0) > 0) == F){
          ## Last step. Criteria
          ### prediction
          MAE[1, "SEM_Reg"] <- NA
          RMSE[1, "SEM_Reg"] <- NA
          OFS[1, "SEM_Reg"] <- NA
          
          ### Measurement model part
          congruence[1:ncol(P_loadings), "SEM_Reg"] <- NA
          PL_rate[1, "SEM_Reg"] <- NA
          
          ### Regression coefficient
          bias_beta0[1, "SEM_Reg"] <- NA
          bias_beta[, "SEM_Reg"] <- NA
          bias_betaALL[1, "SEM_Reg"] <- NA
        }else{
          ### Obtain sample estimated variance–covariance matrix
          covariance_matrix <- cov(df_trainSAM) * (train_n - 1) / train_n
          ### Compute summary statistics E(eta) and var(eta)
          summ.stats <- LSAM_SumStats(S = covariance_matrix,
                                      sample.nobs = train_n,
                                      LAMBDA = mm_est$LAMBDA,
                                      THETA = mm_est$THETA,
                                      mapping = "ML")
          
          ## Step 2. Estimate structural model (based on the variance-covariance matrix of the factor scores)
          reg_coef <- LSAM_regcoef(model = SAM_Model, 
                                   sumstats = summ.stats) # unstandardized beta
          
          ## Step 3. Get the factor score (Bartlett method) of test data
          test_score <- t(solve(res_residual) %*% res_loading %*% solve(t(res_loading) %*% solve(res_residual) %*% res_loading)) %*% t(X_test)
          test_score <- t(test_score)
          ## regression method factor score
          test_scoreReg <- summ.stats[1:9,1:9] %*% t(res_loading) %*% solve((res_loading %*% summ.stats[1:9,1:9] %*% t(res_loading) + res_residual)) %*% t(X_test)
          test_scoreReg <- t(test_scoreReg)
          
          ### get intercept
          #### Bartlett 
          train_score <- t(solve(res_residual) %*% res_loading %*% solve(t(res_loading) %*% solve(res_residual) %*% res_loading)) %*% t(X_train)
          train_score <- t(train_score)
          MeanY <- 0 # Y loading = 1 and residual is 0
          MeanX <- apply(train_score, 2, mean)
          intercept <- MeanY - sum(reg_coef*MeanX)
          reg_coefBar <- c(intercept, reg_coef)
          #### Regression
          train_scoreReg <- summ.stats[1:9,1:9] %*% t(res_loading) %*% solve((res_loading %*% summ.stats[1:9,1:9] %*% t(res_loading) + res_residual)) %*% t(X_train)
          train_scoreReg <- t(train_scoreReg)
          MeanY <- 0 # Y loading = 1 and residual is 0
          MeanXReg <- apply(train_scoreReg, 2, mean)
          interceptReg <- MeanY - sum(reg_coef*MeanXReg)
          reg_coefReg <- c(interceptReg, reg_coef)
          
          ### Step 4. Calculate the predicted variable
          Y_hat <- reg_coefBar[1] + test_score %*% matrix(reg_coefBar[-1], factor_N, 1)
          #### Regression
          Y_hatReg <- reg_coefReg[1] + test_score %*% matrix(reg_coefReg[-1], factor_N, 1)
          
          
          ### Step 5. Criteria
          #### Bartlett
          #### Prediction
          MAE[1, "SAM"] <- sum(abs(Y_test - Y_hat)) / test_n
          RMSE[1, "SAM"] <- sqrt(sum((Y_test - Y_hat)^2) / test_n)
          OFS[1, "SAM"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)
          
          #### Measurement model part
          congruence[1:ncol(P_loadings), "SAM"] <- diag(abs(psych::factor.congruence(res_loading, P_loadings)))
          PL_rate[1, "SAM"] <- PL(P_loadings, res_loading)
          
          #### Regression coefficient
          bias_beta0[1, "SAM"] <- abs(reg_coefBar[1]-beta0)
          bias_beta[, "SAM"] <- abs(reg_coefBar[-1]-beta)
          bias_betaALL[1, "SAM"] <- mean(c(bias_beta0[1, "SAM"], bias_beta[, "SAM"]))
          
          #### Regression
          #### Prediction
          MAE[1, "SAM_Reg"] <- sum(abs(Y_test - Y_hatReg)) / test_n
          RMSE[1, "SAM_Reg"] <- sqrt(sum((Y_test - Y_hatReg)^2) / test_n)
          OFS[1, "SAM_Reg"] <- sum((Y_test - Y_hatReg)^2) / sum(Y_test^2)
          
          #### Measurement model part
          congruence[1:ncol(P_loadings), "SAM_Reg"] <- diag(abs(psych::factor.congruence(res_loading, P_loadings)))
          PL_rate[1, "SAM_Reg"] <- PL(P_loadings, res_loading)
          
          #### Regression coefficient
          bias_beta0[1, "SAM_Reg"] <- abs(reg_coefReg[1]-beta0)
          bias_beta[, "SAM_Reg"] <- abs(reg_coefReg[-1]-beta)
          bias_betaALL[1, "SAM_Reg"] <- mean(c(bias_beta0[1, "SAM"], bias_beta[, "SAM"]))
          
        }
        
        #-------------------------------------------------------------------------------
        #4. Sparse generalized canonical correlation analysis (SGCCA)
        
        ## Step 0. Set matrix
        ### set blocks
        #### a~i are predicting factors with each 4 items; j is the predicted factor
        idx <- split(1:item_total, rep(letters[1:factor_N], each = item_N))
        blocks <- c(
          setNames(lapply(idx, \(i) as.matrix(df_train[, i, drop = FALSE])), paste0("bl", 1:factor_N)),
          list(bl10 = df_train[,item_total+1])
        )
        ### set connection
        connection <- matrix(0, XY_n, XY_n)
        connection[1:factor_N, XY_n] <- 1
        connection[XY_n, 1:factor_N] <- 1
        diag(connection) <- 0
        
        ## Step 1. Get sparsity
        ### use train data to run cross validation to get the sparsity
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
        
        ## Step 2. SGCCA
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
        ### extract weights
        w_list <- setNames(lapply(paste0("bl", 1:factor_N), function(b) fit_sgcca$astar[[b]]), paste0("bl", 1:factor_N))
        ### turn into wight matrix
        W <- matrix(0, nrow = nrow(P_loadings), ncol = ncol(P_loadings),
                    dimnames = dimnames(P_loadings))
        for (j in seq_len(ncol(P_loadings))) {
          idx2 <- which(P_loadings[, j] != 0)
          W[idx2, j] <- w_list[[j]]
        }
        
        # Step 3. Calculate the factor score for training data
        train_score <- lapply(paste0("bl", 1:factor_N), function(b) blocks[[b]] %*% w_list[[b]])
        train_score <- do.call(cbind, train_score)
        train_df <- data.frame(train_score, Y_train)
        
        ## Step 4. Estimate structural model
        res_reg <- lm(Y_train ~ ., data = train_df)
        reg_coef <- unname(res_reg$coefficients)
        
        ## Step 5. Calculate the factor score for test data
        test_data <- c(
          setNames(lapply(idx, \(i) as.matrix(df_test[, i, drop = FALSE])), paste0("bl", 1:factor_N)),
          list(bl10 = df_test[,item_total+1])
        )
        test_score <- lapply(paste0("bl", 1:factor_N), function(b) test_data[[b]] %*% w_list[[b]])
        test_score <- do.call(cbind, test_score)
        
        ## Step 6. Calculate the predicted variable
        Y_hat <- reg_coef[1] + test_score %*% matrix(reg_coef[-1], factor_N, 1)
        
        ## Step 7. Criteria
        ### Prediction
        MAE[1, "SGCCA"] <- sum(abs(Y_test - Y_hat)) / test_n
        RMSE[1, "SGCCA"] <- sqrt(sum((Y_test - Y_hat)^2) / test_n)
        OFS[1, "SGCCA"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)
        
        #### Measurement model part
        congruence[1:ncol(P_loadings), "SGCCA"] <- diag(abs(psych::factor.congruence(W, P_loadings)))
        PL_rate[1, "SGCCA"] <- PL(P_loadings, W)
        
        #### Regression coefficient
        bias_beta0[1, "SGCCA"] <- abs(reg_coef[1]-beta0)
        bias_beta[, "SGCCA"] <- abs(reg_coef[-1]-beta)
        bias_betaALL[1, "SGCCA"] <- mean(c(bias_beta0[1, "SGCCA"], bias_beta[, "SGCCA"]))
        
        #-------------------------------------------------------------------------------
        #5. Regularized Exploratory Structural Equation Modeling (rESEM) (lasso)
        
        ## Step 0. Scale x data
        X_train_s <- scale(X_train)
        
        ## Step 1. Estimate the measurement part
        lasso_start <- c(39, 64, 120, 250)
        lasso <- findLasso(dat = as.matrix(X_train_s), 
                           Ptrue = P_loadings, 
                           maxItr = 20, 
                           lassou = lasso_start[ss])
        res_method <- MULTISTART_rESEMl1(DATA = as.matrix(X_train_s), 
                                         R = factor_N, 
                                         MaxIter = 20, 
                                         eps = 10^-6, 
                                         nstarts = 50, 
                                         lambda = lasso$lasso)
        res_method <- rESEMl1_undoshrinkage(DATA = as.matrix(X_train_s), 
                                            R = factor_N, 
                                            P = res_method$loadings, 
                                            MaxIter = 100, 
                                            eps = 10^-6)
        ### get the scores and loadings
        train_score <- res_method$scores
        res_loading <- res_method$loadings
        ### correct directions
        res_direction <- check_direction(check_load = res_loading, train_score = train_score, true_load = P_loadings, factor_n = ncol(P_loadings))
        #### get the correct direction loadings and score
        res_loading <- res_direction$loadings
        train_score <- res_direction$score
        train_df <- data.frame(train_score, Y_train)
        
        ## Step 2. Estimate the structural part
        res_reg <- lm(Y_train ~ ., data = train_df)
        reg_coef <- unname(res_reg$coefficients)
        
        ## Step 3. Calculate the factor scores
        test_score <- as.matrix(scale(X_test)) %*% res_loading %*% solve(t(res_loading) %*% res_loading)
        
        ## Step 4. Calculate the predicted variable
        Y_hat <- reg_coef[1] + test_score %*% matrix(reg_coef[-1], factor_N, 1)
        
        ## Step 5. Criteria
        ### Prediction
        MAE[1, "rESEM"] <- sum(abs(Y_test - Y_hat)) / test_n
        RMSE[1, "rESEM"] <- sqrt(sum((Y_test - Y_hat)^2) / test_n)
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
        reg_coef <- c(lavaan::coef(fit)["y~1"], lavaan::coef(fit)[paste0("y~l",c(1:factor_N))]) # regression coefficient
        res_loading <- lavaan::inspect(fit, "est")$lambda
        res_loading <- res_loading[1:item_total,1:factor_N]
        ## Step 3. Calculate the predicted variable
        Y_hat <- lavaan::lavPredictY(fit, newdata = df_test, ynames = "y", xnames = c(paste("x", c(1:item_total), sep = "")))
        ### Step 4. Criteria
        #### Prediction
        MAE[1, "SEM_BASED"] <- sum(abs(Y_test - Y_hat)) / test_n
        RMSE[1, "SEM_BASED"] <- sqrt(sum((Y_test - Y_hat)^2) / test_n)
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
        df_trainPLS <- df_train
        df_testPLS <- df_test
        colnames(df_trainPLS) <- NameVariables
        colnames(df_testPLS) <- NameVariables
        predTrain <- PLSpredict(df_trainPLS, df_testPLS, 
                                smMatrix, mmMatrix, 
                                maxIt = 300, stopCriterion = 9)
        Y_hat <- predTrain$predictedMeasurements
        ### get loadings and regression coefficients
        PLS_sim <- simplePLS(df_trainPLS, 
                             smMatrix, mmMatrix, 
                             maxIt = 300, stopCriterion = 9)
        PLS_loadings <- PLS_sim$outer_loadings[-(item_total+1),-(factor_N+1)]
        ### get regression coefficient
        reg_coef <- get_intercepts(PLS_sim)
        
        ## Step 3. Criteria
        ### Prediction
        MAE[1, "PLS"] <- sum(abs(Y_test - Y_hat)) / test_n
        RMSE[1, "PLS"] <- sqrt(sum((Y_test - Y_hat)^2) / test_n)
        OFS[1, "PLS"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)
        
        ### Measurement model part
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
        
        ## Step 1. Estimate linear model based on training data
        Lmodel <- lm(y ~., data = as.data.frame(df_train))
        ### Get the regression coefficients
        reg_coef <- Lmodel$coefficients
        
        ## Step 2. Calculate the predicted variable
        X_test_int <- cbind(1, as.matrix(X_test))
        Y_hat <- X_test_int %*% reg_coef
        
        ## Step 3. Criteria
        ### Prediction
        MAE[1, "GLM"] <- sum(abs(Y_test - Y_hat)) / test_n
        RMSE[1, "GLM"] <- sqrt(sum((Y_test - Y_hat)^2) / test_n)
        OFS[1, "GLM"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)
        
        #-------------------------------------------------------------------------------
        #9. Elastic net 
        
        ## Step 1. Training data
        ### Choose the alpha through cross-validation
        alphas <- seq(0, 1, by = 0.1)
        cvs <- rep(NA, length(alphas))
        for (i in 1:length(alphas)) {
          cvs[i] <- min(
            glmnet::cv.glmnet(X_train, Y_train, alpha  = alphas[i])$cvm)
        }
        aa <- which.min(cvs) # the alpha having the smallest output
        ### Choose the lambda through cross-validation
        cv_fit <- glmnet::cv.glmnet(X_train, Y_train, alpha = alphas[aa])
        # Extract lambda
        best_lambda <- cv_fit$lambda.min # the best lambda having the smallest result
        
        ## Step 2. Fit elastic net
        ### Fit final model on entire training set with that lambda and alpha
        final_model <- glmnet::glmnet(X_train, Y_train, alpha = alphas[aa], lambda = best_lambda)
        ### the regression coefficient 
        reg_coef <- coef(final_model)
        
        ## Step 3. Calculate the predicted variable
        Y_hat <- predict(final_model, newx = X_test, s = best_lambda)
        
        ## Step 4. Criteria
        ### Prediction
        MAE[1, "elastic"] <- sum(abs(Y_test - Y_hat)) / test_n
        RMSE[1, "elastic"] <- sqrt(sum((Y_test - Y_hat)^2) / test_n)
        OFS[1, "elastic"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)
        
        
        
        #-------------------------------------------------------------------------------
        # The end
        #-------------------------------------------------------------------------------
        
        criteria[[measurement[m]]][[structural[s]]][[as.character(sample_condition[ss])]] <- list(
          MAE         = MAE,
          RMSE        = RMSE,
          OFS         = OFS,
          congruence  = congruence,
          PL_rate     = PL_rate,
          beta_mean   = bias_betaALL)
        
      }
    }
  }
  saveRDS(criteria, paste0("criteria_", format(Sys.time(), "%m-%d-%H-%M-%S"), ".rds"))
}
end_time <- Sys.time()
time <- (end_time - start); time
stopCluster(cl)

criteria <- readRDS("Replication_1110-03-12-22-23.rds")


