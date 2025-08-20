# Factor correction
source("function/factor_correction.R")
# SAM
source("function/factor2/3.1_SAM_PartitionedEstimation.R")
source("function/factor2/3.2_SAM_ComputeSummaryStats.R")
source("function/factor2/3.3_SAM_RegressionCoefs.R")
# rESEM
source("function/factor2/5_rESEM_function.R")
# PLS_prediction
source("function/factor2/7_simplePLS.R")
source("function/factor2/7_PLSpredict.R")

## settings
### sample size
sample_size <- 1000
### number of indicators for a factor
itemNum <- 10
### measurement part
#### first loadings
loading1 <- c(rep(sqrt(.6), itemNum), rep(0, 10))
loading2 <- c(rep(0, 10),rep(sqrt(.6), itemNum))
P_loadings <- matrix(cbind(loading1, loading2), nrow = length(loading1), ncol = 2)
### structural part
beta0 <- 0
beta1 <- 1
beta2 <- 1

# ## generate data1
Uniq <- diag(1 - rowSums(P_loadings ^ 2))
SIGMA <- P_loadings %*% t(P_loadings) + Uniq
X <- MASS::mvrnorm(n = sample_size,
                   mu = rep(0, nrow(P_loadings)),
                   Sigma = SIGMA,
                   empirical = TRUE
)
Score <- t(solve(Uniq) %*% P_loadings %*% solve(t(P_loadings) %*% solve(Uniq) %*% P_loadings)) %*% t(X)
Y <- rep(beta0, sample_size) + beta1 * Score[1, ] + beta2 * Score[2, ] + rnorm(n = sample_size, mean = 0, sd = 0.01)
### The dataset
data <- cbind(X,Y)
colnames(data) <- c(paste("x", c(1:itemNum), sep = ""), paste("z", c(1:itemNum), sep = ""),"y")

## generate data2
residVar <- 0.01^2
Uniq <- diag(1 - rowSums(P_loadings ^ 2))
allUniq <- diag(c(1 - rowSums(P_loadings ^ 2), residVar))
allCoef <- rbind(P_loadings, c(beta1, beta2))
Phi <- matrix(c(1,0,0,1),2,2)
allSIGMA <- allCoef %*% Phi %*% t(allCoef) + allUniq
data <- MASS::mvrnorm(n = sample_size,
                         mu = rep(0, nrow(allCoef)),
                         Sigma = allSIGMA,
                         empirical = TRUE)
colnames(data) <- c(paste("x", c(1:itemNum), sep = ""), paste("z", c(1:itemNum), sep = ""),"y")
cov(data)

X <- data[,1:20]
Y <- data[,21]

## generate data3
### factor score
Score <- cbind(
  F1 = rnorm(sample_size, mean = 0, sd = 1),
  F2 = rnorm(sample_size, mean = 0, sd = 1)
)
### Construct the uniqueness (measurement error) covariance
Uniq <- diag(1 - rowSums(P_loadings ^ 2))
### error 
errors <- MASS::mvrnorm(
  n    = sample_size,
  mu   = rep(0, nrow(P_loadings)),
  Sigma = Uniq,
  empirical = TRUE
)
### observation x
X <- Score %*% t(P_loadings) + errors
### predicted value Y
Y <- beta0 +
     beta1 * Score[, "F1"] +
     beta2 * Score[, "F2"] +
     rnorm(sample_size, mean = 0, sd = 0.01)
### The dataset
data <- cbind(X,Y)
colnames(data) <- c(
  paste0("x", 1:itemNum),  # x1–x10 load on F1
  paste0("z", 1:itemNum),  # z1–z10 load on F2
  "y"
)

## Split data
n <- nrow(data)
### training set index
train_idx <- sample(n, round(n/2)) 
test_idx  <- setdiff(1:sample_size, train_idx)
### set as data.frame
X <- as.data.frame(X)
Y <- as.data.frame(Y)
### Training data
X_train <- X[train_idx, ] 
Y_train <- Y[train_idx, ] 
### Test data
X_test <- X[test_idx, ]
Y_test <- Y[test_idx, ]
## Combine
n_ind <- length(loading1)
df_train <- cbind(X_train, Y_train) |> setNames(c(paste("x", c(1:itemNum), sep = ""), paste("z", c(1:itemNum), sep = ""),"y"))
df_test <- cbind(X_test, Y_test) |> setNames(c(paste("x", c(1:itemNum), sep = ""), paste("z", c(1:itemNum), sep = ""),"y"))

#-------------------------------------------------------------------------------
# Criteria
#-------------------------------------------------------------------------------

methods <- c("SumScore_E", "SumScore_T", 
             "SEM_Reg", "SEM_Bar", 
             "SAM", "SGCCA","rESEM", 
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

# Structure model part
bias_beta0 <- matrix(NA, nrow = 1, ncol = length(methods))
colnames(bias_beta0) <- methods
bias_beta1 <- matrix(NA, nrow = 1, ncol = length(methods))
colnames(bias_beta1) <- methods
bias_beta2 <- matrix(NA, nrow = 1, ncol = length(methods))
colnames(bias_beta2) <- methods
bias_betaALL <- matrix(NA, nrow = 1, ncol = length(methods))
colnames(bias_betaALL) <- methods

#-------------------------------------------------------------------------------
# Classical approach
#-------------------------------------------------------------------------------

#1. Sum score

##(1) Empirical method
### Step 1. Exploratory factor analysis
result <- psych::fa(X_train, nfactors = ncol(P_loadings), rotate='varimax', fm = "minres")# defult is minimum residual factoring
loadings_EFA <- as.matrix(result$loadings)

### Step 2. Set the unit weight for the indicators
#### function for transform the loading matrix
weight_Emp <- check_lambda(loadings_EFA)
#### correct directions
weight_Emp <- check_direction(check_load = weight_Emp, train_score = NA, true_load = P_loadings, factor_n = ncol(P_loadings))$loadings

### Step 3: Calculate sum score
sum_score <- as.matrix(X_train) %*% weight_Emp
### Step 4. Estimate the regression coefficients
reg_mod <- lm(Y_train ~ sum_score[,1] + sum_score[,2])
reg_coef <- unname(reg_mod$coefficients)
### Step 5. Calculate the predicted variable based on the test data
#### test score
test_score <- as.matrix(X_test) %*% weight_Emp
Y_hat <- rep(reg_coef[1],length(Y_test)) + test_score[,1]*reg_coef[2] + test_score[,2]*reg_coef[3]
### Step 6. Criteria
#### Prediction
MAE[1, "SumScore_E"] <- sum(abs(Y_test-Y_hat)) / (sample_size/2)
RMSE[1, "SumScore_E"] <- sqrt(sum((Y_test-Y_hat)^2) / (sample_size/2))
OFS[1, "SumScore_E"] <- (sum((Y_test-Y_hat)^2) / sum((Y_test)^2))

#### Measurement model part
congruence[1:ncol(P_loadings), "SumScore_E"] <- diag(psych::factor.congruence(weight_Emp, P_loadings))

#### Regression coefficient
bias_beta0[1, "SumScore_E"] <- abs(reg_coef[1]-beta0)
bias_beta1[1, "SumScore_E"] <- abs(reg_coef[2]-beta1)
bias_beta2[1, "SumScore_E"] <- abs(reg_coef[3]-beta2)
bias_betaALL[1, "SumScore_E"] <- abs(sum(bias_beta0[1, "SumScore_E"], bias_beta1[1, "SumScore_E"], bias_beta2[1, "SumScore_E"]))

#-------------------------------------------------------------------------------

##(2) Theory-based method
### Step1: Set the unit weight for the indicators
#### function for transform the loading matrix
weight_theory <- check_lambda(loadings_EFA)
#### correct directions
weight_theory <- check_direction(check_load = weight_Emp, train_score = NA, true_load = P_loadings, factor_n = ncol(P_loadings))$loadings
### Step2: Calculate sum score
sum_score <- as.matrix(X_train) %*% weight_theory
### Step 3. Estimate the regression coefficients
reg_mod <- lm(Y_train ~ sum_score[,1] + sum_score[,2])
reg_coef <- unname(reg_mod$coefficients)
### Step 4. Calculate the predicted variable based on the test data
#### test score
test_score <-  as.matrix(X_test) %*% weight_theory
Y_hat <- rep(reg_coef[1],length(Y_test)) + test_score[,1]*reg_coef[2] + test_score[,2]*reg_coef[3]
### Step 5. Criteria
#### Prediction
MAE[1, "SumScore_T"] <- sum(abs(Y_test-Y_hat)) / (sample_size/2)
RMSE[1, "SumScore_T"] <- sqrt(sum((Y_test-Y_hat)^2) / (sample_size/2))
OFS[1, "SumScore_T"] <- (sum((Y_test-Y_hat)^2) / sum((Y_test)^2))

#### Measurement model part
congruence[1:ncol(P_loadings), "SumScore_T"] <- diag(psych::factor.congruence(weight_theory, P_loadings))

#### Regression coefficient
bias_beta0[1, "SumScore_T"] <- abs(reg_coef[1]-beta0)
bias_beta1[1, "SumScore_T"] <- abs(reg_coef[2]-beta1)
bias_beta2[1, "SumScore_T"] <- abs(reg_coef[3]-beta2)
bias_betaALL[1, "SumScore_T"] <- abs(sum(bias_beta0[1, "SumScore_T"], bias_beta1[1, "SumScore_T"], bias_beta2[1, "SumScore_T"]))

#-------------------------------------------------------------------------------

#2. Traditional SEM, Regression method

## Step 2. Set the model
SEM_Model <- '
  #Structural model
  y ~ factor1 + factor2
  
  # Measurement model for factor1
  factor1 =~ NA*x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
  factor2 =~ NA*z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10
  
  # Fix variance of factor1 to 1
  factor1 ~~ 1 * factor1
  factor2 ~~ 1 * factor2
'
## Step 3. Estimation of CB-SEM
fit <- lavaan::sem(model= SEM_Model, data = df_train, meanstructure = T)
reg_coef <- cbind(lavaan::coef(fit)["y~1"], lavaan::coef(fit)["y~factor1"], lavaan::coef(fit)["y~factor2"]) # regression coefficient
### Step 4. Get the factor score of test data
test_score <- lavaan::lavPredict(fit, newdata = df_test, method = "regression")
###
res_loading <- lavaan::inspect(fit, "est")$lambda
res_loading <- res_loading[1:20,1:2]
res_residual <- diag(lavaan::inspect(fit, "est")$theta)[-21]
res_residual <- diag(res_residual)
test_score1 <- t(solve(res_loading %*% t(res_loading) + res_residual) %*% res_loading) %*% t(X_test)
test_score1 <- t(test_score1)
###

### Step 5. Calculate the predicted variable
Y_hat <- rep(reg_coef[1],length(Y_test)) + test_score1[,1]*reg_coef[2] + test_score1[,2]*reg_coef[3]
### Step 6. Criteria
#### Prediction
MAE[1, "SEM_Reg"] <- sum(abs(Y_test-Y_hat)) / (sample_size/2)
RMSE[1, "SEM_Reg"] <- sqrt(sum((Y_test-Y_hat)^2) / (sample_size/2))
OFS[1, "SEM_Reg"] <- (sum((Y_test-Y_hat)^2) / sum((Y_test)^2))

#### Measurement model part
congruence[1:ncol(P_loadings), "SEM_Reg"] <- diag(psych::factor.congruence(res_loading, P_loadings))

#### Regression coefficient
bias_beta0[1, "SEM_Reg"] <- abs(reg_coef[1]-beta0)
bias_beta1[1, "SEM_Reg"] <- abs(reg_coef[2]-beta1)
bias_beta2[1, "SEM_Reg"] <- abs(reg_coef[3]-beta2)
bias_betaALL[1, "SEM_Reg"] <- abs(sum(bias_beta0[1, "SEM_Reg"], bias_beta1[1, "SEM_Reg"], bias_beta2[1, "SEM_Reg"]))


#2.1 Traditional SEM, Bartlett method
### Step 4. Get the factor score of test data
test_score <- lavaan::lavPredict(fit, df_test, method = "Bartlett")
test_score1 <- t(solve(res_residual) %*% res_loading %*% solve(t(res_loading) %*% solve(res_residual) %*% res_loading)) %*% t(X_test)
test_score1 <- t(test_score1)
### Step 5. Calculate the predicted variable
Y_hat <- rep(reg_coef[1],length(Y_test)) + test_score1[,1]*reg_coef[2] + test_score1[,2]*reg_coef[3]
### Step 6. Criteria
#### Prediction
MAE[1, "SEM_Bar"] <- sum(abs(Y_test - Y_hat)) / (sample_size / 2)
RMSE[1, "SEM_Bar"] <- sqrt(sum((Y_test - Y_hat)^2) / (sample_size / 2))
OFS[1, "SEM_Bar"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)

#### Measurement model part
congruence[1:ncol(P_loadings), "SEM_Bar"] <- diag(psych::factor.congruence(res_loading, P_loadings))

#### Regression coefficient
bias_beta0[1, "SEM_Bar"] <- abs(reg_coef[1] - beta0)
bias_beta1[1, "SEM_Bar"] <- abs(reg_coef[2] - beta1)
bias_beta2[1, "SEM_Bar"] <- abs(reg_coef[3] - beta2)
bias_betaALL[1, "SEM_Bar"] <- abs(
  sum(
    bias_beta0[1, "SEM_Bar"],
    bias_beta1[1, "SEM_Bar"],
    bias_beta2[1, "SEM_Bar"]
  )
)


#-------------------------------------------------------------------------------
# SEM Based (Two-stage)
#-------------------------------------------------------------------------------

#3. Structural After Measurement (SAM) 

## Step 1. Estimate measurement model with CFA
SAM_Model <- '
  #Structural model
  fy ~ fx + fz
  
  # Measurement model for factor1
  fx =~ NA*x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
  fz =~ NA*z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10
  
  # Fix variance of factor1 to 1
  fx ~~ 1 * fx
  fz ~~ 1 * fz
'
### Compute loadings and residual variances for each measurement block separately
# sub.est <- partitioned.estimation(data = df_train, method = "SAM_Bentler_ULS")
mm_est <- partitioned.estimation(data = df_train, method = "SAM_MLB")
#### loadings
res_loading <- mm_est$LAMBDA[1:20, 1:2]
res_loading <- check_direction(check_load = res_loading, train_score = NA, true_load = P_loadings, factor_n = ncol(P_loadings))$loadings
#### residuals
res_residual <- mm_est$THETA[1:20,1:20]
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
mu_Score <- colMeans(t(Score))
mu_Y <- colMeans(Y)
beta_vec <- c(reg_coef[1], reg_coef[2])
reg_coef <- c(mu_Y - sum(mu_Score * beta_vec), reg_coef[1], reg_coef[2])
## Step 3. Get the factor score (Bartlett method) of test data
test_score <- t(solve(res_residual) %*% res_loading %*% solve(t(res_loading) %*% solve(res_residual) %*% res_loading)) %*% t(X_test)
test_score <- t(test_score)
### Step 4. Calculate the predicted variable
Y_hat <- rep(reg_coef[1],length(Y_test)) + test_score[,1] * reg_coef[2] + test_score[,2] * reg_coef[3]
### Step 5. Criteria
#### Prediction
MAE[1, "SAM"] <- sum(abs(Y_test - Y_hat)) / (sample_size / 2)
RMSE[1, "SAM"] <- sqrt(sum((Y_test - Y_hat)^2) / (sample_size / 2))
OFS[1, "SAM"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)

#### Measurement model part
congruence[1:ncol(P_loadings), "SAM"] <- diag(psych::factor.congruence(res_loading, P_loadings))

#### Regression coefficient
#bias_beta0[1, "SAM"] <- abs(reg_coef[1] - beta0)
bias_beta1[1, "SAM"] <- abs(reg_coef[1] - beta1)
bias_beta2[1, "SAM"] <- abs(reg_coef[2] - beta2)
bias_betaALL[1, "SAM"] <- abs(
  sum(
    bias_beta0[1, "SAM"],
    bias_beta1[1, "SAM"],
    bias_beta2[1, "SAM"]
  )
)

#-------------------------------------------------------------------------------

#4. Sparse generalized canonical correlation analysis (SGCCA)
## Step 0. Transform data
x_train <- df_train[, 1:10]
z_train<- df_train[, 11:20]
y_train <- df_train[, 21]

x_test <- df_test[, 1:10]
z_test<- df_test[, 11:20]
y_test <- df_test[, 21]

## Step 1. Estimate measurement model
### Settings
blocks <- list(bl1 = x_train, bl2 = z_train, bl3 = y_train)
connection <- matrix(c(0, 0, 1, 0, 1, 0, 1, 0, 0), nrow = 3, byrow = TRUE)
sparsity <- c(1,1,1)
### Run SGCCA
fit_sgcca <- RGCCA::rgcca(blocks = blocks,
                          connection = connection,
                          method = "sgcca",
                          sparsity = sparsity,
                          # each block remains separate rather than combining into a meta-block
                          superblock = FALSE, 
                          # Extracts one component per block
                          ncomp = c(1,1,1), 
                          scheme = "factorial", 
                          # extracted components are orthogonal (uncorrelated)
                          comp_orth = TRUE, 
                          verbose = F)
### Get the weights
w_x <- fit_sgcca$astar$bl1
w_z <- fit_sgcca$astar$bl2
weight_SGCCA <- cbind(c(w_x, rep(0,itemNum)), c(rep(0,itemNum) ,w_z))
## Step 2. Calculate the factor score for training data
train_score_x <- as.matrix(x_train) %*% w_x
train_score_z <- as.matrix(z_train) %*% w_z
## Step 3. Estimate structural model
res_reg <- lm(Y_train ~ train_score_x + train_score_z)
reg_coef <- unname(res_reg$coefficients)
## Step 4. Calculate the factor score for test data
test_score_x <- as.matrix(x_test)%*%w_x
test_score_z <- as.matrix(z_test)%*%w_z
## Step 5. Calculate the predicted variable
Y_hat <- rep(reg_coef[1],length(Y_test)) + test_score_x*reg_coef[2] + test_score_z*reg_coef[3]
## Step 6. Criteria
### Prediction
MAE[1, "SGCCA"] <- sum(abs(Y_test - Y_hat)) / (sample_size / 2)
RMSE[1, "SGCCA"] <- sqrt(sum((Y_test - Y_hat)^2) / (sample_size / 2))
OFS[1, "SGCCA"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)

#### Measurement model part
congruence[1:ncol(P_loadings), "SGCCA"] <- diag(psych::factor.congruence(weight_SGCCA, P_loadings))

#### Regression coefficient
bias_beta0[1, "SGCCA"] <- abs(reg_coef[1] - beta0)
bias_beta1[1, "SGCCA"] <- abs(reg_coef[2] - beta1)
bias_beta2[1, "SGCCA"] <- abs(reg_coef[3] - beta2)
bias_betaALL[1, "SGCCA"] <- abs(
  sum(
    bias_beta0[1, "SGCCA"],
    bias_beta1[1, "SGCCA"],
    bias_beta2[1, "SGCCA"]
  )
)

#-------------------------------------------------------------------------------

#5. Regularized Exploraty Structural Equation Modeling (rESEM) (cardinality constrains)

## Step 1. Estimate the measurement part
res_method <- MULTISTART_CCrESEM(DATA = as.matrix(X_train), R = 2, CARD = 20, MaxIter = 200, eps = 10^-6)
### get the scores and loadings
train_score <- res_method$scores
res_loading <- res_method$loadings
### correct directions
res_direction <- check_direction(check_load = res_loading, train_score = train_score, true_load = P_loadings, factor_n = ncol(P_loadings))
res_loading <- res_direction$loadings
train_score_T <- res_direction$score

## Step 2. Estimate the structural part
res_reg <- lm(Y_train ~ train_score_T[,1] + train_score_T[,2])
reg_coef <- unname(res_reg$coefficients)

## Step 3. Calculate the factor scores
test_score <- as.matrix(X_test) %*% res_loading %*% solve(t(res_loading) %*% res_loading)
## Step 4. Calculate the predicted variable
Y_hat <- rep(reg_coef[1],length(Y_test)) + test_score[,1]*reg_coef[2] + test_score[,2]*reg_coef[3]
## Step 5. Criteria
### Prediction
MAE[1, "rESEM"] <- sum(abs(Y_test - Y_hat)) / (sample_size / 2)
RMSE[1, "rESEM"] <- sqrt(sum((Y_test - Y_hat)^2) / (sample_size / 2))
OFS[1, "rESEM"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)

#### Measurement model part
congruence[1:ncol(P_loadings), "rESEM"] <- diag(psych::factor.congruence(res_loading, P_loadings))

#### Regression coefficient
bias_beta0[1, "rESEM"] <- abs(reg_coef[1] - beta0)
bias_beta1[1, "rESEM"] <- abs(reg_coef[2] - beta1)
bias_beta2[1, "rESEM"] <- abs(reg_coef[3] - beta2)
bias_betaALL[1, "rESEM"] <- abs(
  sum(
    bias_beta0[1, "rESEM"],
    bias_beta1[1, "rESEM"],
    bias_beta2[1, "rESEM"]
  )
)

#-------------------------------------------------------------------------------
# SEM Based (One-stage)
#-------------------------------------------------------------------------------

#6. SEM-Based Prediction Rule
## Step 1. Build the model
SEMRule_Model <- '
  # Structural model
  y ~ factor1 + factor2
  
  # Measurement model for factor1
  factor1 =~ NA*x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
  factor2 =~ NA*z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10
  
  # Fix variance of factor1 to 1
  factor1 ~~ 1 * factor1
  factor2 ~~ 1 * factor2
'
## Step 2. Estimation of CB-SEM
fit <- lavaan::sem(SEMRule_Model , data = df_train , meanstructure = TRUE)
reg_coef <- cbind(lavaan::coef(fit)["y~1"], lavaan::coef(fit)["y~factor1"], lavaan::coef(fit)["y~factor2"]) # regression coefficient
## Step 3. Calculate the predicted variable
Y_hat <- lavaan::lavPredictY(fit, newdata = df_test, xnames = c(paste("x", c(1:itemNum), sep = ""), paste("z", c(1:itemNum), sep = "")), ynames = c("y"))
### Step 4. Criteria
#### Prediction
MAE[1, "SEM_BASED"] <- sum(abs(Y_test - Y_hat)) / (sample_size / 2)
RMSE[1, "SEM_BASED"] <- sqrt(sum((Y_test - Y_hat)^2) / (sample_size / 2))
OFS[1, "SEM_BASED"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)

#### Measurement model part
congruence[1:ncol(P_loadings), "SEM_BASED"] <- diag(psych::factor.congruence(weight_Emp, P_loadings))

#### Regression coefficient
bias_beta0[1, "SEM_BASED"] <- abs(reg_coef[1] - beta0)
bias_beta1[1, "SEM_BASED"] <- abs(reg_coef[2] - beta1)
bias_beta2[1, "SEM_BASED"] <- abs(reg_coef[3] - beta2)
bias_betaALL[1, "SEM_BASED"] <- abs(
  sum(
    bias_beta0[1, "SEM_BASED"],
    bias_beta1[1, "SEM_BASED"],
    bias_beta2[1, "SEM_BASED"]
  )
)

#-------------------------------------------------------------------------------

#7. PLS-SEM
## Step 1. Create matrix of the measurement and structural Model
smMatrix <- matrix(c("X","Y","Z","Y"),
                   nrow = 2, ncol = 2, byrow = T, 
                   dimnames = list(1:2,c("source", "target")))
mmMatrix <- matrix(c("X","x1","R",
                     "X","x2","R",
                     "X","x3","R",
                     "X","x4","R",
                     "X","x5","R",
                     "X","x6","R",
                     "X","x7","R",
                     "X","x8","R",
                     "X","x9","R",
                     "X","x10","R",
                     "Z","z1","R",
                     "Z","z2","R",
                     "Z","z3","R",
                     "Z","z4","R",
                     "Z","z5","R",
                     "Z","z6","R",
                     "Z","z7","R",
                     "Z","z8","R",
                     "Z","z9","R",
                     "Z","z10","R",
                     "Y","y","R"),
                   nrow = 21, ncol = 3, byrow = T,
                   dimnames = list(1:21,c("latent","measurement","type")))
## Step 2. Prediction on the model
predTrain <- PLSpredict(df_train, df_test, smMatrix, mmMatrix, 300,9)
Y_hat <- predTrain$predictedMeasurements

## Step 3. Criteria
### Prediction
MAE[1, "PLS"] <- sum(abs(Y_test - Y_hat)) / (sample_size / 2)
RMSE[1, "PLS"] <- sqrt(sum((Y_test - Y_hat)^2) / (sample_size / 2))
OFS[1, "PLS"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)

#### Measurement model part
congruence[1:ncol(P_loadings), "PLS"] <- diag(psych::factor.congruence(PLS_loadings[1:nrow(P_loadings), c("X","Z")], P_loadings))

#### Regression coefficient
bias_beta0[1, "PLS"] <- abs(reg_coef["Y"] - beta0)
bias_beta1[1, "PLS"] <- abs(reg_coef["X"] - beta1)
bias_beta2[1, "PLS"] <- abs(reg_coef["Z"] - beta2)
bias_betaALL[1, "PLS"] <- abs(
  sum(
    bias_beta0[1, "PLS"],
    bias_beta1[1, "PLS"],
    bias_beta2[1, "PLS"]
  )
)

#-------------------------------------------------------------------------------
# Machine learning approaches
#-------------------------------------------------------------------------------

#8. Linear regression
## Step 1. Estimate linear model based on test data
Lmodel <- lm(y ~., data = df_train)
### Get the regression coefficients
reg_coef <- Lmodel$coefficients
## Step 2. Calculate the predicted variable
X_test_int <- cbind(1, as.matrix(X_test))
Y_hat <- X_test_int %*% reg_coef

## Step 3. Criteria
### Prediction
MAE[1, "GLM"] <- sum(abs(Y_test - Y_hat)) / (sample_size / 2)
RMSE[1, "GLM"] <- sqrt(sum((Y_test - Y_hat)^2) / (sample_size / 2))
OFS[1, "GLM"] <- sum((Y_test - Y_hat)^2) / sum(Y_test^2)

#### Measurement model part
# congruence[1:ncol(P_loadings), "GLM"] <- diag(psych::factor.congruence(weight_Emp, P_loadings))

#### Regression coefficient
# bias_beta0[1, "GLM"] <- abs(reg_coef[1] - beta0)
# bias_beta1[1, "GLM"] <- abs(reg_coef[2] - beta1)
# bias_beta2[1, "GLM"] <- abs(reg_coef[3] - beta2)
# bias_betaALL[1, "GLM"] <- abs(
#   sum(
#     bias_beta0[1, "GLM"],
#     bias_beta1[1, "GLM"],
#     bias_beta2[1, "GLM"]
#   )
# )

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

#### Measurement model part
congruence[1:ncol(P_loadings), "elastic"] <- diag(psych::factor.congruence(weight_Emp, P_loadings))

#### Regression coefficient
# bias_beta0[1, "elastic"] <- abs(reg_coef[1] - beta0)
# bias_beta1[1, "elastic"] <- abs(reg_coef[2] - beta1)
# bias_beta2[1, "elastic"] <- abs(reg_coef[3] - beta2)
# bias_betaALL[1, "elastic"] <- abs(
#   sum(
#     bias_beta0[1, "elastic"],
#     bias_beta1[1, "elastic"],
#     bias_beta2[1, "elastic"]
#   )
# )
