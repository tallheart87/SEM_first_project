# Factor correction
source("function/factor_correction.R")
# SAM
source("function/factor9/3.1_SAM_PartitionedEstimation.R")
source("function/factor9/3.2_SAM_ComputeSummaryStats.R")
source("function/factor9/3.3_SAM_RegressionCoefs.R")
# rESEM
source("function/factor9/5_rESEM_function.R")
# PLS_prediction
source("function/factor9/7_simplePLS.R")
source("function/factor9/7_PLSpredict.R")
#-------------------------------------------------------------------------------
# Data generation
#-------------------------------------------------------------------------------

# sample size
SS <- c(200,400,1000,2000)
sample_size <- SS[4]

# measurement model
lambda <- matrix(runif(36, min = 0.5, max = 0.8), nrow = 4, ncol = 9)
colnames(lambda) <- c(paste0("eta", c(1:9)))
## get different format loadings (with 0)
eta <- as.matrix(lambda)
p <- nrow(eta)
k <- ncol(eta)
P_loadings <- sapply(1:k, function(j) {
  c(rep(0, (j-1)*p), eta[, j], rep(0, (k-j)*p))
})
P_loadings <- matrix(P_loadings, nrow = p*k, ncol = k,
                     dimnames = list(paste0("row", 1:(p*k)), colnames(eta)))
# structural model
beta0 <- 0
beta <- runif(9, min = 0.25, max = 0.4)

# data generate
Uniq <- diag(1 - rowSums(P_loadings ^ 2))
SIGMA <- P_loadings %*% t(P_loadings) + Uniq
X <- MASS::mvrnorm(n = sample_size,
                   mu = rep(0, nrow(P_loadings)),
                   Sigma = SIGMA,
                   empirical = TRUE
)
Score <- t(solve(Uniq) %*% P_loadings %*% solve(t(P_loadings) %*% solve(Uniq) %*% P_loadings)) %*% t(X)
Y <- rep(beta0, sample_size) + beta[1] * Score[1, ] + 
                               beta[2] * Score[2, ] + 
                               beta[3] * Score[3, ] +
                               beta[4] * Score[4, ] + 
                               beta[5] * Score[5, ] + 
                               beta[6] * Score[6, ] +
                               beta[7] * Score[7, ] + 
                               beta[8] * Score[8, ] + 
                               beta[9] * Score[9, ] + rnorm(n = sample_size, mean = 0, sd = 0.01)


# The dataset
data <- cbind(X,Y)
colnames(data) <- c(
  paste0("x", 1:36),  # x1–x10 load on F1
  "y"  # z1–z10 load on F2
)

# Split data
n <- nrow(data)
## training set index
train_idx <- sample(n, round(n/2)) 
test_idx  <- setdiff(1:n, train_idx)
## set as data.frame
X <- as.data.frame(X)
Y <- as.data.frame(Y)
## Training data
X_train <- X[train_idx, ] 
Y_train <- Y[train_idx, ] 
## Test data
X_test <- X[test_idx, ]
Y_test <- Y[test_idx, ]
## Combine
df_train <- cbind(X_train, Y_train) |> setNames(c(paste0("x", 1:36), "y"))
df_test <- cbind(X_test, Y_test) |> setNames(c(paste0("x", 1:36), "y"))

#-------------------------------------------------------------------------------
# Criteria
#-------------------------------------------------------------------------------

methods <- c("SumScore_E", 
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
congruence[1:ncol(P_loadings), "SumScore_E"] <- diag(psych::factor.congruence(weight_Emp, P_loadings))

#### Regression coefficient
bias_beta0[1, "SumScore_E"] <- abs(reg_coef[1]-beta0)
bias_beta1[1, "SumScore_E"] <- abs(reg_coef[2]-beta1)
bias_beta2[1, "SumScore_E"] <- abs(reg_coef[3]-beta2)
bias_betaALL[1, "SumScore_E"] <- abs(sum(bias_beta0[1, "SumScore_E"], bias_beta1[1, "SumScore_E"], bias_beta2[1, "SumScore_E"]))

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
fit <- lavaan::sem(model= SEM_Model, data = df_train, meanstructure = T)
reg_coef <- c(lavaan::coef(fit)["y~1"], lavaan::coef(fit)[paste0("y~l",c(1:9))]) # regression coefficient

## Step 3. Get the factor score of test data
test_score <- lavaan::lavPredict(fit, newdata = df_test, method = "regression")
###
res_loading <- lavaan::inspect(fit, "est")$lambda
res_loading <- res_loading[1:36,1:9]
res_residual <- diag(lavaan::inspect(fit, "est")$theta)[-37]
res_residual <- diag(res_residual)
test_score1 <- t(solve(res_loading %*% t(res_loading) + res_residual) %*% res_loading) %*% t(X_test)
test_score1 <- t(test_score1)
###

### Step 4. Calculate the predicted variable
Y_hat <- rep(reg_coef[1],length(Y_test)) + 
  test_score1[,1]*reg_coef[2] + 
  test_score1[,2]*reg_coef[3] +
  test_score1[,3]*reg_coef[4] + 
  test_score1[,4]*reg_coef[5] +
  test_score1[,5]*reg_coef[6] + 
  test_score1[,6]*reg_coef[7] +
  test_score1[,7]*reg_coef[8] + 
  test_score1[,8]*reg_coef[9] +
  test_score1[,9]*reg_coef[10] 
### Step 5. Criteria
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
Y_hat <- rep(reg_coef[1],length(Y_test)) + 
             test_score1[,1]*reg_coef[2] + 
             test_score1[,2]*reg_coef[3] +
             test_score1[,3]*reg_coef[4] + 
             test_score1[,4]*reg_coef[5] +
             test_score1[,5]*reg_coef[6] + 
             test_score1[,6]*reg_coef[7] +
             test_score1[,7]*reg_coef[8] + 
             test_score1[,8]*reg_coef[9] +
             test_score1[,9]*reg_coef[10] 
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
reg_coef <- c(mu_Y - sum(mu_Score * reg_coef), reg_coef)
## Step 3. Get the factor score (Bartlett method) of test data
test_score <- t(solve(res_residual) %*% res_loading %*% solve(t(res_loading) %*% solve(res_residual) %*% res_loading)) %*% t(X_test)
test_score <- t(test_score)
### Step 4. Calculate the predicted variable
Y_hat <- rep(reg_coef[1],length(Y_test)) + 
  test_score[,1]*reg_coef[2] + 
  test_score[,2]*reg_coef[3] +
  test_score[,3]*reg_coef[4] + 
  test_score[,4]*reg_coef[5] +
  test_score[,5]*reg_coef[6] + 
  test_score[,6]*reg_coef[7] +
  test_score[,7]*reg_coef[8] + 
  test_score[,8]*reg_coef[9] +
  test_score[,9]*reg_coef[10] 
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
df_train <- cbind(X_train, Y_train) |> setNames(c(paste0("x", 1:36), "y"))
df_test <- cbind(X_test, Y_test) |> setNames(c(paste0("x", 1:36), "y"))

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
connection <- diag(10)[, 10:1]
sparsity <- rep(1,10)
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
res_method <- MULTISTART_CCrESEM(DATA = as.matrix(X_train), R = 9, CARD = rep(4,9), MaxIter = 200, eps = 10^-6)
### get the scores and loadings
train_score <- res_method$scores
res_loading <- res_method$loadings
### correct directions
res_direction <- check_direction(check_load = res_loading, train_score = train_score, true_load = P_loadings, factor_n = ncol(P_loadings))
res_loading <- res_direction$loadings
train_score_T <- res_direction$score

## Step 2. Estimate the structural part
res_reg <- lm(Y_train ~ train_score_T[,1] + train_score_T[,2] + 
                        train_score_T[,3] + train_score_T[,4] + 
                        train_score_T[,5] + train_score_T[,6] +
                        train_score_T[,7] + train_score_T[,8] +
                        train_score_T[,9] )
reg_coef <- unname(res_reg$coefficients)

## Step 3. Calculate the factor scores
test_score <- as.matrix(X_test) %*% res_loading %*% solve(t(res_loading) %*% res_loading)
## Step 4. Calculate the predicted variable
Y_hat <- rep(reg_coef[1],length(Y_test)) + 
             test_score[,1]*reg_coef[2] + 
             test_score[,2]*reg_coef[3] +
             test_score[,3]*reg_coef[4] + 
             test_score[,4]*reg_coef[5] +
             test_score[,5]*reg_coef[6] + 
             test_score[,6]*reg_coef[7] +
             test_score[,7]*reg_coef[8] + 
             test_score[,8]*reg_coef[9] +
             test_score[,9]*reg_coef[10] 

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
fit <- lavaan::sem(SEMRule_Model , data = df_train , meanstructure = TRUE)
reg_coef <- c(lavaan::coef(fit)["y~1"], lavaan::coef(fit)[paste0("y~l",c(1:9))]) # regression coefficient
## Step 3. Calculate the predicted variable
Y_hat <- lavaan::lavPredictY(fit, newdata = df_test, xnames = c(paste("x", c(1:36), sep = ""), "y"))
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

