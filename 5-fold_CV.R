# 5-fold cv --------


library(glmnet)
library(selectiveInference)
library(ggplot2)
library(MASS)
rm(list = ls())

n = 60
p = 45

# Equicorrelated Covariance
Sigma <- matrix(0.8, n, n)
diag(Sigma) <- 1
mean_vec <- rep(0, n)

# Building a design matrix X 
set.seed(123) 
X <- mvrnorm(p, mu = mean_vec, Sigma = Sigma)
X <- t(X)
X = scale(X, TRUE, TRUE)

# Generating betas 
set.seed(123)
num_non_zero_betas = 10
beta = c(rep(2, num_non_zero_betas), rep(0, p - num_non_zero_betas)) 
mu = X %*% beta

# Error variance 
sigma = 1

SEED = 111



## 5-fold CV ----------------
B = 1e4
coverage_prob_cv = numeric(B) 
confidence_len_cv = numeric(B)
coverage_prob_cv_act = numeric(B) 
confidence_len_cv_act = numeric(B) 
coverage_prob_cv_inact = numeric(B) 
confidence_len_cv_inact = numeric(B) 

f1_scores = numeric(B)
lambdas = numeric(B)
beta_M_CV = matrix(NA, nrow = B, ncol = 5)  

pb = txtProgressBar(min = 0, max = 100, style = 3, width = 50, char = "=") # For progress bar 
set.seed(SEED)

for(b in 1:B) {
  y = mu + sigma * rnorm(n)
  lasso = cv.glmnet(X, y, standardize = FALSE, nfolds = 5)
  lambda = lasso$lambda.min 
  lambdas[b] = lambda
  beta_lasso = coef(lasso, x = X, y = y, s = lambda, exact = TRUE)[-1]
  
  # Define the selection model 
  M = which(beta_lasso != 0)  
  X_M = X[, M]
  
  beta_M = solve(t(X_M) %*% X_M) %*% t(X_M) %*% mu 
  active_M = beta_M[1:num_non_zero_betas]
  inactive_M = beta_M[-c(1:num_non_zero_betas)]
  
  selInf = fixedLassoInf(X, y, beta_lasso, lambda * n, sigma = sigma) 
  CIs = selInf$ci
  active_CIs = CIs[1:num_non_zero_betas,]
  inactive_CIs = CIs[-c(1:num_non_zero_betas),]
  
  coverage_prob_cv[b] = sum((beta_M > CIs[,1])*(beta_M < CIs[,2]))/length(M) 
  confidence_len_cv[b] = sum((CIs[,2] - CIs[,1]))/length(M)
  coverage_prob_cv_act[b] = sum((active_M > active_CIs[,1])*(active_M < active_CIs[,2]))/length(active_M) 
  confidence_len_cv_act[b] = sum((active_CIs[,2] - active_CIs[,1]))/length(active_M)
  
  if (length(inactive_M) == 1) {
    coverage_prob_cv_inact[b] = sum((inactive_M > inactive_CIs[1])*(inactive_M < inactive_CIs[2]))/length(inactive_M) 
    confidence_len_cv_inact[b] = sum((inactive_CIs[2] - inactive_CIs[1]))/length(inactive_M)
  } else if (length(inactive_M) == 0) {
    coverage_prob_cv_inact[b] = 0
    confidence_len_cv_inact[b] = 0
  } else {
    coverage_prob_cv_inact[b] = sum((inactive_M > inactive_CIs[,1])*(inactive_M < inactive_CIs[,2]))/length(inactive_M) 
    confidence_len_cv_inact[b] = sum((inactive_CIs[,2] - inactive_CIs[,1]))/length(inactive_M)
  }
  
  
  # Compute F1 score
  true_positive = sum(beta != 0 & beta_lasso != 0)
  false_positive = sum(beta == 0 & beta_lasso != 0)
  false_negative = sum(beta != 0 & beta_lasso == 0)
  
  f1_scores[b] = 2 * true_positive / (2 * true_positive + false_positive + false_negative)
  
  setTxtProgressBar(pb, 100 * b / B) # For progress bar 
}

# Plotting lambdas
hist_lambdas <- ggplot(data.frame(lambdas), aes(lambdas)) + 
  geom_histogram(binwidth=0.0005, fill="#69b3a2", color="#e9ecef", alpha=0.9) + 
  labs(x=bquote(lambda~"s"), y = "Count") +
  theme_classic()

hist_lambdas
ggsave("hist_lambdas_4.pdf", plot = hist_lambdas, width = 6, height = 4)

mean(f1_scores)
fixed_lambda_f1 = lambdas[which.min(abs(mean(f1_scores) - f1_scores))]
fixed_lambda_f1

fixed_lambda_avg = mean(lambdas)
fixed_lambda_avg

mean(coverage_prob_cv)
median(confidence_len_cv)

Finite_CI_prob_cv = sum(is.infinite(confidence_len_cv) == FALSE) / B 
Finite_CI_prob_cv_act = sum(is.infinite(confidence_len_cv_act) == FALSE) / B 
Finite_CI_prob_cv_inact = sum(is.infinite(confidence_len_cv_inact) == FALSE) / B 



## Fixed lambda (having the closest F1 score) ------------- 


B = 1e4
coverage_prob_fixed_f1 = numeric(B) 
confidence_len_fixed_f1 = numeric(B)
coverage_prob_fixed_act_f1 = numeric(B) 
confidence_len_fixed_act_f1 = numeric(B) 
coverage_prob_fixed_inact_f1 = numeric(B) 
confidence_len_fixed_inact_f1 = numeric(B) 

f1_scores_fixed_f1 = numeric(B)
beta_M_fixed_f1 = matrix(NA, nrow = B, ncol = 5)

pb = txtProgressBar(min = 0, max = 100, style = 3, width = 50, char = "=") # For progress bar 
set.seed(SEED)

for(b in 1:B){
  
  B = 1e4
  y = mu + sigma*rnorm(n)
  lasso = glmnet(X, y, standardize = FALSE) 
  lambda = fixed_lambda_f1
  beta_lasso = coef(lasso, x = X, y = y, s = lambda, exact = TRUE)[-1]
  
  # Define the selection model 
  M = which(beta_lasso != 0)  
  X_M = X[ , M]
  
  beta_M = solve(t(X_M)%*%X_M)%*%t(X_M)%*% mu
  active_M = beta_M[1:num_non_zero_betas]
  inactive_M = beta_M[-c(1:num_non_zero_betas)]
  
  selInf = fixedLassoInf(X, y, beta_lasso, lambda*n, sigma = sigma) 
  CIs = selInf$ci
  active_CIs = CIs[1:num_non_zero_betas,]
  inactive_CIs = CIs[-c(1:num_non_zero_betas),]
  
  coverage_prob_fixed_f1[b] = sum((beta_M > CIs[,1])*(beta_M < CIs[,2]))/length(M) 
  confidence_len_fixed_f1[b] = sum((CIs[,2] - CIs[,1]))/length(M)
  coverage_prob_fixed_act_f1[b] = sum((active_M > active_CIs[,1])*(active_M < active_CIs[,2]))/length(active_M) 
  confidence_len_fixed_act_f1[b] = sum((active_CIs[,2] - active_CIs[,1]))/length(active_M)
  
  if (length(inactive_M) == 1) {
    coverage_prob_fixed_inact_f1[b] = sum((inactive_M > inactive_CIs[1])*(inactive_M < inactive_CIs[2]))/length(inactive_M) 
    confidence_len_fixed_inact_f1[b] = sum((inactive_CIs[2] - inactive_CIs[1]))/length(inactive_M)
  } else if (length(inactive_M) == 0) {
    coverage_prob_fixed_inact_f1[b] = 0
    confidence_len_fixed_inact_f1[b] = 0
  } else {
    coverage_prob_fixed_inact_f1[b] = sum((inactive_M > inactive_CIs[,1])*(inactive_M < inactive_CIs[,2]))/length(inactive_M) 
    confidence_len_fixed_inact_f1[b] = sum((inactive_CIs[,2] - inactive_CIs[,1]))/length(inactive_M)
  }
  
  # Compute F1 score
  true_positive = sum(beta != 0 & beta_lasso != 0)
  false_positive = sum(beta == 0 & beta_lasso != 0)
  false_negative = sum(beta != 0 & beta_lasso == 0)
  
  f1_scores_fixed_f1[b] = true_positive / (true_positive + 1/2 * (false_positive + false_negative)) 
  
  setTxtProgressBar(pb, 100*b/B) # For progress bar 
}

mean(f1_scores_fixed_f1)
mean(f1_scores)

mean(coverage_prob_fixed_f1)
median(confidence_len_fixed_f1)


Finite_CI_prob_fixed_f1 = sum(is.infinite(confidence_len_fixed_f1) == FALSE) / B 
Finite_CI_prob_fixed_act_f1 = sum(is.infinite(confidence_len_fixed_act_f1) == FALSE) / B 
Finite_CI_prob_fixed_inact_f1 = sum(is.infinite(confidence_len_fixed_inact_f1) == FALSE) / B 


## Fixed lambda (Average) -------------


B = 1e4
coverage_prob_fixed_avg = numeric(B) 
confidence_len_fixed_avg = numeric(B)
coverage_prob_fixed_act_avg = numeric(B) 
confidence_len_fixed_act_avg = numeric(B) 
coverage_prob_fixed_inact_avg = numeric(B) 
confidence_len_fixed_inact_avg = numeric(B) 

f1_scores_fixed_avg = numeric(B)
beta_M_fixed_avg = matrix(NA, nrow = B, ncol = 5)

pb = txtProgressBar(min = 0, max = 100, style = 3, width = 50, char = "=") # For progress bar 
set.seed(SEED)

for(b in 1:B){
  
  B = 1e4
  y = mu + sigma*rnorm(n)
  lasso = glmnet(X, y, standardize = FALSE) 
  lambda = fixed_lambda_avg
  beta_lasso = coef(lasso, x = X, y = y, s = lambda, exact = TRUE)[-1]
  
  # Define the selection model 
  M = which(beta_lasso != 0)  
  X_M = X[ , M]
  
  beta_M = solve(t(X_M)%*%X_M)%*%t(X_M)%*% mu
  active_M = beta_M[1:num_non_zero_betas]
  inactive_M = beta_M[-c(1:num_non_zero_betas)]
  
  selInf = fixedLassoInf(X, y, beta_lasso, lambda*n, sigma = sigma) 
  CIs = selInf$ci
  active_CIs = CIs[1:num_non_zero_betas,]
  inactive_CIs = CIs[-c(1:num_non_zero_betas),]
  
  coverage_prob_fixed_avg[b] = sum((beta_M > CIs[,1])*(beta_M < CIs[,2]))/length(M) 
  confidence_len_fixed_avg[b] = sum((CIs[,2] - CIs[,1]))/length(M)
  coverage_prob_fixed_act_avg[b] = sum((active_M > active_CIs[,1])*(active_M < active_CIs[,2]))/length(active_M) 
  confidence_len_fixed_act_avg[b] = sum((active_CIs[,2] - active_CIs[,1]))/length(active_M)
  
  if (length(inactive_M) == 1) {
    coverage_prob_fixed_inact_avg[b] = sum((inactive_M > inactive_CIs[1])*(inactive_M < inactive_CIs[2]))/length(inactive_M) 
    confidence_len_fixed_inact_avg[b] = sum((inactive_CIs[2] - inactive_CIs[1]))/length(inactive_M)
  } else if (length(inactive_M) == 0) {
    coverage_prob_fixed_inact_avg[b] = 0
    confidence_len_fixed_inact_avg[b] = 0
  } else {
    coverage_prob_fixed_inact_avg[b] = sum((inactive_M > inactive_CIs[,1])*(inactive_M < inactive_CIs[,2]))/length(inactive_M) 
    confidence_len_fixed_inact_avg[b] = sum((inactive_CIs[,2] - inactive_CIs[,1]))/length(inactive_M)
  }
  
  # Compute F1 score
  true_positive = sum(beta != 0 & beta_lasso != 0)
  false_positive = sum(beta == 0 & beta_lasso != 0)
  false_negative = sum(beta != 0 & beta_lasso == 0)
  
  f1_scores_fixed_avg[b] = true_positive / (true_positive + 1/2 * (false_positive + false_negative)) 
  
  setTxtProgressBar(pb, 100*b/B) # For progress bar 
}

mean(f1_scores_fixed_avg)
mean(f1_scores)

mean(coverage_prob_fixed_avg)
median(confidence_len_fixed_avg)


Finite_CI_prob_fixed_avg = sum(is.infinite(confidence_len_fixed_avg) == FALSE) / B 
Finite_CI_prob_fixed_act_avg = sum(is.infinite(confidence_len_fixed_act_avg) == FALSE) / B 
Finite_CI_prob_fixed_inact_avg = sum(is.infinite(confidence_len_fixed_inact_avg) == FALSE) / B 


## Visualisation ---------------

# Comparing F1 score 

F1_score_df <- data.frame(
  Iteration = 1:B,
  CV = f1_scores,
  Fixed_by_F1 = f1_scores_fixed_f1,
  Fixed_by_avg = f1_scores_fixed_avg
)

F1_score_melt <- melt(F1_score_df, id.vars = "Iteration")


box_plot_F1 <- ggplot(F1_score_melt, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  labs(x = expression(paste(lambda,"s")),
       y = "F1 scores",
       fill = expression(paste(lambda,"s"))) +
  theme_minimal()

box_plot_F1

ggsave("box_plot_F1_4.pdf", plot = box_plot_F1, width = 6, height = 4)


# Comparing CI lengths 

## Total 
CI_length_Total_df <- data.frame(
  Iteration = 1:B,
  CV_Total = confidence_len_cv,
  Fixed_F1_Total = confidence_len_fixed_f1,
  Fixed_Avg_Total = confidence_len_fixed_avg
)

CI_length_Total_melt <- melt(CI_length_Total_df, id.vars = "Iteration")

box_plot_CI_Total <- ggplot(CI_length_Total_melt, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  labs(
    x = expression(paste(lambda, "s")),
    y = "CI lengths",
    fill = expression(paste(lambda, "s"))
  ) +
  theme_minimal()

print(box_plot_CI_Total)

ggsave("box_plot_CI_Total_4.pdf", plot = box_plot_CI_Total, width = 6, height = 4)


## Active coefficients 

CI_length_Act_df <- data.frame(
  Iteration = 1:B,
  CV_Act = confidence_len_cv_act,
  Fixed_F1_Act = confidence_len_fixed_act_f1,
  Fixed_Avg_Act = confidence_len_fixed_act_avg
)

CI_length_Act_melt <- melt(CI_length_Act_df, id.vars = "Iteration")

box_plot_CI_Act <- ggplot(CI_length_Act_melt, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  labs(
    x = expression(paste(lambda, "s")),
    y = "CI lengths",
    fill = expression(paste(lambda, "s"))
  ) +
  theme_minimal()

print(box_plot_CI_Act)

ggsave("box_plot_CI_Act_4.pdf", plot = box_plot_CI_Act, width = 6, height = 4)


## Inactive coefficients 

CI_length_Inact_df <- data.frame(
  Iteration = 1:B,
  CV_Inact = confidence_len_cv_inact,
  Fixed_F1_Inact = confidence_len_fixed_inact_f1,
  Fixed_Avg_Inact = confidence_len_fixed_inact_avg
)

CI_length_Inact_melt <- melt(CI_length_Inact_df, id.vars = "Iteration")

box_plot_CI_Inact <- ggplot(CI_length_Inact_melt, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  labs(
    x = expression(paste(lambda, "s")),
    y = "CI lengths",
    fill = expression(paste(lambda, "s"))
  ) +
  theme_minimal()

print(box_plot_CI_Inact)

ggsave("box_plot_CI_Inact_4.pdf", plot = box_plot_CI_Inact, width = 6, height = 4)



# --- Table 

# CV 
cv_result <- data.frame(Coverage_average = c(mean(coverage_prob_cv), mean(coverage_prob_cv_act), mean(coverage_prob_cv_inact)),
                        CI_len_median = c(median(confidence_len_cv), median(confidence_len_cv_act), median(confidence_len_cv_inact)),
                        Finite_CI_prob = c(Finite_CI_prob_cv, Finite_CI_prob_cv_act, Finite_CI_prob_cv_inact)
)
rownames(cv_result) <- c("Total", "Active (M)", "Inactive (-M)")
round(cv_result, 3)

# Fixed_F1 

fixed_f1_result <- data.frame(Coverage_average = c(mean(coverage_prob_fixed_f1), mean(coverage_prob_fixed_act_f1), mean(coverage_prob_fixed_inact_f1)),
                              CI_len_median = c(median(confidence_len_fixed_f1), median(confidence_len_fixed_act_f1), median(confidence_len_fixed_inact_f1)),
                              Finite_CI_prob = c(Finite_CI_prob_fixed_f1, Finite_CI_prob_fixed_act_f1, Finite_CI_prob_fixed_inact_f1)
)
rownames(fixed_f1_result) <- c("Total", "Active (M)", "Inactive (-M)")
round(fixed_f1_result, 3)

# Fixed_avg 

fixed_avg_result <- data.frame(Coverage_average = c(mean(coverage_prob_fixed_avg), mean(coverage_prob_fixed_act_avg), mean(coverage_prob_fixed_inact_avg)),
                               CI_len_median = c(median(confidence_len_fixed_avg), median(confidence_len_fixed_act_avg), median(confidence_len_fixed_inact_avg)),
                               Finite_CI_prob = c(Finite_CI_prob_fixed_avg, Finite_CI_prob_fixed_act_avg, Finite_CI_prob_fixed_inact_avg)
)
rownames(fixed_avg_result) <- c("Total", "Active (M)", "Inactive (-M)")
round(fixed_avg_result, 3)