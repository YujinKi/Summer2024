# p exceeds n with all non-zero betas --------- 


library(glmnet)
library(selectiveInference)
library(ggplot2)
library(MASS)
library(reshape2)

n = 60
p = 80

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
num_non_zero_betas = p
beta = rep(2, num_non_zero_betas)
mu = X %*% beta

# Error variance 
sigma = 1

SEED = 111



## 10-fold CV ---------------
B = 1e4
coverage_prob_cv = numeric(B) 
confidence_len_cv = numeric(B)

f1_scores = numeric(B)
lambdas = numeric(B)
seleted_models = numeric(B)
beta_M_CV = matrix(NA, nrow = B, ncol = 5)  

pb = txtProgressBar(min = 0, max = 100, style = 3, width = 50, char = "=") # For progress bar 
set.seed(SEED)

for(b in 1:B) {
  y = mu + sigma * rnorm(n)
  lasso = cv.glmnet(X, y, standardize = FALSE)
  lambda = lasso$lambda.min 
  lambdas[b] = lambda
  beta_lasso = coef(lasso, x = X, y = y, s = lambda, exact = TRUE)[-1]
  
  # Define the selection model 
  M = which(beta_lasso != 0)  
  seleted_models[b] = length(M)
  
  if(length(M) == 0) {
    coverage_prob_cv[b] = NA
    confidence_len_cv[b] = NA
    next
  }
  
  X_M = X[, M]
  
  
  # Compute F1 score
  true_positive = sum(beta != 0 & beta_lasso != 0)
  false_positive = sum(beta == 0 & beta_lasso != 0)
  false_negative = sum(beta != 0 & beta_lasso == 0)
  
  f1_scores[b] = 2 * true_positive / (2 * true_positive + false_positive + false_negative)
  
  setTxtProgressBar(pb, 100 * b / B) # For progress bar 
}

# Plotting lambdas
hist_lambdas <- ggplot(data.frame(lambdas), aes(lambdas)) +
  geom_histogram(binwidth=0.05, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  labs(title = bquote("Histogram of cross-validation" ~ lambda),x=bquote(lambda~"s"), y = "Count") +
  theme_classic()

hist_lambdas
ggsave("hist_lambdas_5.pdf", plot = hist_lambdas, width = 6, height = 4)


mean(f1_scores)
fixed_lambda_f1 = lambdas[which.min(abs(mean(f1_scores) - f1_scores))]
fixed_lambda_f1

fixed_lambda_avg = mean(lambdas)
fixed_lambda_avg

incld_prob = sum(seleted_models < n)/B
incld_prob

hist_selected_models_cv <- ggplot(data.frame(seleted_models), aes(seleted_models)) +
  geom_histogram(binwidth=1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_vline(xintercept=60, color="red", linetype = 2) +
  annotate("text", x=58.5, y=2200, label="n = 60", size=5, colour = "red", angle = 90) +
  labs(x="The number of selected model", y = "Count") +
  theme_classic()

hist_selected_models_cv
ggsave("hist_selected_models_cv.pdf", plot = hist_selected_models_cv, width = 6, height = 4)




## Fixed lambda (having the closest F1 score) ----------- 


B = 1e4
coverage_prob_fixed_f1 = numeric(B) 
confidence_len_fixed_f1 = numeric(B)
coverage_prob_fixed_act_f1 = numeric(B) 
confidence_len_fixed_act_f1 = numeric(B) 
coverage_prob_fixed_inact_f1 = numeric(B) 
confidence_len_fixed_inact_f1 = numeric(B) 

f1_scores_fixed_f1 = numeric(B)
seleted_models_f1 = numeric(B)
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
  seleted_models_f1[b] = length(M)
  
  if(length(M) == 0) {
    coverage_prob_cv[b] = NA
    confidence_len_cv[b] = NA
    next
  }
  
  X_M = X[, M]
  
  
  # Compute F1 score
  true_positive = sum(beta != 0 & beta_lasso != 0)
  false_positive = sum(beta == 0 & beta_lasso != 0)
  false_negative = sum(beta != 0 & beta_lasso == 0)
  
  f1_scores_fixed_f1[b] = true_positive / (true_positive + 1/2 * (false_positive + false_negative)) 
  
  setTxtProgressBar(pb, 100*b/B) # For progress bar 
}

mean(f1_scores_fixed_f1)
mean(f1_scores)

selected_prob_f1 = sum(seleted_models_f1 < n)/B
selected_prob_f1

hist_selected_models_f1 <- ggplot(data.frame(seleted_models_f1), aes(seleted_models_f1)) +
  geom_histogram(binwidth=1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_vline(xintercept=60, color="red", linetype = 2) +
  annotate("text", x=58.5, y=2200, label="n = 60", size=5, colour = "red", angle = 90) +
  labs(x="The number of selected model", y = "Count") +
  theme_classic()

hist_selected_models_f1
ggsave("hist_selected_models_f1.pdf", plot = hist_selected_models_f1, width = 6, height = 4)



## Fixed lambda (Average) -----------


B = 1e4
coverage_prob_fixed_avg = numeric(B) 
confidence_len_fixed_avg = numeric(B)


f1_scores_fixed_avg = numeric(B)
seleted_models_avg = numeric(B)
beta_M_fixed_avg = matrix(NA, nrow = B, ncol = 5)

pb = txtProgressBar(min = 0, max = 100, style = 3, width = 50, char = "=") # For progress bar 
set.seed(SEED)

for(b in 1:B){
  
  B = 1e4
  y = mu + sigma*rnorm(n)
  lasso = glmnet(X, y, standardize = FALSE, thresh = 1e-10) 
  lambda = fixed_lambda_avg
  beta_lasso = coef(lasso, x = X, y = y, s = lambda, exact = TRUE)[-1]
  
  # Define the selection model 
  M = which(beta_lasso != 0)  
  seleted_models_avg[b] = length(M)
  
  if(length(M) == 0) {
    coverage_prob_cv[b] = NA
    confidence_len_cv[b] = NA
    next
  }
  
  X_M = X[, M]
  
 
  # Compute F1 score
  true_positive = sum(beta != 0 & beta_lasso != 0)
  false_positive = sum(beta == 0 & beta_lasso != 0)
  false_negative = sum(beta != 0 & beta_lasso == 0)
  
  f1_scores_fixed_avg[b] = true_positive / (true_positive + 1/2 * (false_positive + false_negative)) 
  
  setTxtProgressBar(pb, 100*b/B) # For progress bar 
}

mean(f1_scores_fixed_avg)
mean(f1_scores)

selected_prob_avg = sum(seleted_models_avg < n)/B
selected_prob_avg

hist_selected_models_avg <- ggplot(data.frame(seleted_models_avg), aes(seleted_models_avg)) +
  geom_histogram(binwidth=1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_vline(xintercept=60, color="red", linetype = 2) +
  annotate("text", x=58.5, y=2200, label="n = 60", size=5, colour = "red", angle = 90) +
  labs(x="The number of selected model", y = "Count") +
  theme_classic()

hist_selected_models_avg
ggsave("hist_selected_models_avg.pdf", plot = hist_selected_models_avg, width = 6, height = 4)



## Visualisation -------------

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
  labs(title = expression("Box plot of F1 scores"),
       x = expression(paste(lambda,"s")),
       y = "F1 scores",
       fill = expression(paste(lambda,"s"))) +
  theme_minimal()

box_plot_F1

ggsave("box_plot_F1_5.pdf", plot = box_plot_F1, width = 6, height = 4)


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
    title = "Box Plots of Confidence Interval Lengths (Total coefficietns)",
    x = expression(paste(lambda, "s")),
    y = "CI lengths",
    fill = expression(paste(lambda, "s"))
  ) +
  theme_minimal()

print(box_plot_CI_Total)

ggsave("box_plot_CI_Total_5.pdf", plot = box_plot_CI_Total, width = 6, height = 4)



# --- Table 

# CV 
cv_result <- data.frame(Coverage_average = c(mean(coverage_prob_cv), mean(coverage_prob_fixed_f1), mean(coverage_prob_fixed_avg)),
                        CI_len_median = c(median(confidence_len_cv), median(confidence_len_fixed_f1), median(confidence_len_fixed_avg)),
                        Finite_CI_prob = c(Finite_CI_prob_cv, Finite_CI_prob_fixed_f1, Finite_CI_prob_fixed_avg)
)
rownames(cv_result) <- c("10-fold CV", "Fixed by F1 score", "Fixed by Average")
round(cv_result, 3)




