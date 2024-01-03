library(tidyverse)
theme_set(theme_bw())

EM_sim_CIR <- function(total_time, step_length, param, X_0){
  beta <- param[1]
  mu <- param[2]
  sigma <- param[3]
  
  # Calculate the number of steps correctly
  N <- as.integer(total_time / step_length) 
  if (total_time %% step_length != 0) {
    N <- N + 1
  }
  
  # Initialize the process
  X <- numeric(N + 1) # N + 1 to include the initial point
  X[1] <- X_0
  dW <- rnorm(N, mean = 0, sd = sqrt(step_length))
  
  # Euler-Maruyama method to simulate the process
  for (t in 2:(N + 1)) {
    X[t] <- X[t - 1] - beta * (X[t - 1] - mu) * step_length + sigma * sqrt(X[t - 1]) * dW[(t - 1)]
  }
  
  # Euler-Maruyama method to simulate -log(X_t) directly, but with same wiener increments
  neg_logX <- numeric(N + 1)
  neg_logX[1] <- -log(X_0)
  
  for(t in 2:(N + 1)) {
    neg_logX[t] <- neg_logX[t - 1] +
      (beta * (1 - mu*exp(neg_logX[t - 1])) + 1/2 * sigma^2 * exp(neg_logX[t - 1])) * step_length -  
      sigma * exp(neg_logX[t - 1]/2) * dW[(t - 1)]
  }
  
  lamperti_X <- numeric(N + 1)
  lamperti_X[1] <- 2 * sqrt(X_0)
  
  for(t in 2:(N + 1)) {
    lamperti_X[t] <- lamperti_X[t - 1] +
      1/2 * (- beta * lamperti_X[t - 1] + (4 * beta * mu - sigma^2)/lamperti_X[t - 1]) * step_length +
      sigma * dW[(t - 1)]
  }
  
  neg_log_lamperti_X <- numeric(N + 1)
  neg_log_lamperti_X[1] <- -2 * exp(-neg_logX[1]/2)
  
  for(t in 2:(N + 1)) {
    neg_log_lamperti_X[t] <- neg_log_lamperti_X[t - 1] +
      1 / 2 * (-beta * neg_log_lamperti_X[t - 1] + 1/neg_log_lamperti_X[t - 1] * (4*mu - sigma^2)) * step_length -
      sigma * dW[(t - 1)]
  }
  
  # Create time vector
  t_vec <- seq(0, length.out = N + 1, by = step_length)
  
  # Create tibble
  tibble(t = t_vec, X = X,
         manual_negative_log = -log(X), ito_negative_log = neg_logX,
         manual_lamperti = 2*sqrt(X), lamperti_X,
         manual_negative_log_lamperti = -2*exp(-neg_logX/2), neg_log_lamperti_X,
         X_closed_form = (4 * beta * X) / (sigma^2 * (1 - exp(-beta * step_length))))
}

#CIR models

# Baseline closed form solution implementation. Note the data needs to be multiplied by 4 * beta/(sigma^2 * (1 - exp(-beta * t)))
closedform_likelihood_CIR <- function(param, data, dt){
  # Initialize parameters
  beta <- param[1]
  mu <- param[2]
  sigma <- param[3]
  
  # Prepare data
  N <- length(data)
  lower_x <- data[1:(N - 1)]
  upper_x <- data[2:N]
  
  likelihood_value <- -sum(dchisq(upper_x, df = 4 * mu * beta / sigma^2, ncp = lower_x * exp(-beta * dt), log = TRUE))
  likelihood_value
}

closedform_likelihood_CIR(initial_param, sim_data$X_closed_form, actual_dt)


# Direct Euler Maruyama
euler_maruyama_likelihood_CIR <- function(param, data, dt){
  # Initialize parameters
  beta <- param[1]
  mu <- param[2]
  sigma <- param[3]
  
  # Prepare data
  N <- length(data)
  lower_x <- data[1:(N - 1)]
  upper_x <- data[2:N]
  
  likelihood_value <- 1/2 * sum(log(2*pi*sigma^2 * lower_x * dt) + 
                                    (upper_x - lower_x + beta * (lower_x - mu) * dt)^2 /
                                    (sigma^2 * lower_x * dt))
  likelihood_value
}

# Lamperti-transformed CIR models

# Euler-maruyama
euler_maruyama_likelihood_CIR_lamperti <- function(param, data, dt){
  # Initialize parameters
  beta <- param[1]
  mu <- param[2]
  sigma <- param[3]
  
  # Prepare data
  N <- length(data)
  lower_x <- data[1:(N - 1)]
  upper_x <- data[2:N]
  
  likelihood_value <- 1/2 * sum(log(2 * pi * dt * sigma^2) +
                                  (upper_x - lower_x * (1 - beta * dt / 2) + ((-4 * beta * mu - sigma^2) * dt)/(2 * lower_x))^2 /
                                  (sigma^2 * dt))
                        
  likelihood_value
}

# Lie Trotter

# runge_kutta4_solver <- function(df, x0, y0, h){
#   k1 <- df(x0, y0)
#   k2 <- df(x0 + h/2, y0 + h*k1/2)
#   k3 <- df(x0 + h/2, y0 + h*k2/2)
#   k4 <- df(x0 + h, y0 + h*k3)
#   
#   yn <- y0 + h/6 * (k1 + 2*k2+ 2*k3+ k4)
#   yn
# }

LT_likelihood_CIR_lamperti <- function(param, data, dt){
  beta <- param[1]
  mu <- param[2]
  sigma <- param[3]
  
  B <- 2 * beta * mu - sigma^2/2
  
  N <- length(data)
  lower_x <- data[1:(N - 1)]
  upper_x <- data[2:N]
  
  f_delta_t <- sqrt(2 * B * dt + lower_x^2)
  
  N * (log(sigma^2/beta) + log(1 - exp(-beta * dt))) +
    beta * sum((upper_x - exp(-beta * dt / 2) * f_delta_t)^2) /
    (sigma^2 * (1 - exp(-beta * dt)))
}


# Strang
strang_likelihood_CIR_lamperti <- function(param, data, dt){
  beta <- param[1]
  mu <- param[2]
  sigma <- param[3]
  
  B <- 2 * beta * mu - sigma^2/2
  
  N <- length(data)
  lower_x <- data[1:(N - 1)]
  upper_x <- data[2:N]
  
  f_delta_t <- sqrt(B * dt + lower_x^2)
  
  f_inverse_delta_t <- sqrt(upper_x^2 - B * dt)
  
  f_inverse_delta_diff <- upper_x/f_inverse_delta_t
  
  
  N * (log(sigma^2/beta) + log(1 - exp(-beta * dt))) +
    beta * sum((f_inverse_delta_t - exp(-beta * dt / 2) * f_delta_t)^2) / 
    (sigma^2 * (1 - exp(-beta * dt))) -
    2 * sum(log(abs(f_inverse_delta_diff)))
}

# Negative log models

euler_maruyama_likelihood_CIR_neglog <- function(param, data, dt){
  beta <- param[1]
  mu <- param[2]
  sigma <- param[3]
  
  N <- length(data)
  lower_x <- data[1:(N - 1)]
  upper_x <- data[2:N]
  
  likelihood_value <- 1/2 * sum(
    log(2 * pi * sigma^2 * dt) + lower_x +
      (upper_x - lower_x + (- beta * (1 - mu * exp(lower_x)) + 1/2 * sigma^2 * exp(lower_x)) * dt)^2 / 
      (sigma^2 * exp(lower_x) * dt)
  )
  likelihood_value
}

# Simulation and optimization of likelihoods
beta <- 5
mu <- 20
sigma <- 4
X_0 <- mu

# Check boundary
2*beta*mu>=sigma^2


actual_dt <- 0.001

sim_data <-  EM_sim_CIR(109, actual_dt, c(beta, mu, sigma), X_0)

#sim_data %>% ggplot(aes(x = t, y = X)) + geom_line()

initial_param <- c(beta, mu, sigma)

# CIR models

# Closed form solution
result_closedform_likelihood_CIR <- optim(par = initial_param, fn = closedform_likelihood_CIR,
                                          data = sim_data$X_closed_form, dt = actual_dt)


# Direct maruyama

result_EM_CIR <- optim(par = initial_param, fn = euler_maruyama_likelihood_CIR,
      data = sim_data$X, dt = actual_dt)

result_EM_CIR$par

# Lamperti transform

result_EM_CIR_lamperti <- optim(par = initial_param, fn = euler_maruyama_likelihood_CIR_lamperti,
                       data = sim_data$lamperti_X, dt = actual_dt)

result_EM_CIR_lamperti$par


# -log of CIR models

# Direct maruyama


result_EM_CIR_neglog <- optim(par = initial_param, fn = euler_maruyama_likelihood_CIR_neglog,
                       data = sim_data$ito_negative_log, dt = actual_dt)

result_EM_CIR_neglog$par

result_LT_likelihood_CIR_lamperti <- optim(par = initial_param, fn = LT_likelihood_CIR_lamperti,
                              data = sim_data$lamperti_X, dt = actual_dt)

result_LT_likelihood_CIR_lamperti$par

result_strang_likelihood_CIR_lamperti <- optim(par = initial_param, fn = strang_likelihood_CIR_lamperti,
                                           data = sim_data$lamperti_X, dt = actual_dt)

result_strang_likelihood_CIR_lamperti$par



# 
# ### Experiments
# 
# 
# # Experiment with ARE
# dt_values <- c(0.0001, 0.00001, 0.00001)
# results <- tibble(dt = double(), beta_ARE = double(), mu_ARE = double(), sigma_ARE = double())
# num_runs <- 10
# results <- tibble(dt = double(), run = integer(), beta_ARE = double(), mu_ARE = double(), sigma_ARE = double())
# 
# # Loop over dt values
# for (dt in dt_values) {
#   for (run in 1:num_runs) {
#     # Simulate data
#     sim_data <- EM_sim_CIR(total_time = 10, step_length = dt, param = c(beta, mu, sigma), X_0 = X_0)
#     
#     # Estimate parameters
#     result_EM <- optim(par = initial_param, fn = euler_maruyama_likelihood_CIR,
#                        method = "L-BFGS-B",
#                        lower = c(-Inf, 0, 0),
#                        upper = c(Inf, Inf, Inf), 
#                        data = sim_data$X, dt = dt)
#     
#     # Calculate ARE
#     ARE <- abs(result_EM$par - initial_param) / initial_param
#     
#     # Store results
#     results <- bind_rows(results, tibble(dt = dt, run = run, beta_ARE = ARE[1], mu_ARE = ARE[2], sigma_ARE = ARE[3]))
#   }
# }
# 
# # Calculate statistics
# stats <- results %>%
#   group_by(dt) %>%
#   summarize(
#     beta_median = median(beta_ARE),
#     beta_q025 = quantile(beta_ARE, probs = 0.025),
#     beta_q975 = quantile(beta_ARE, probs = 0.975),
#     mu_median = median(mu_ARE),
#     mu_q025 = quantile(mu_ARE, probs = 0.025),
#     mu_q975 = quantile(mu_ARE, probs = 0.975),
#     sigma_median = median(sigma_ARE),
#     sigma_q025 = quantile(sigma_ARE, probs = 0.025),
#     sigma_q975 = quantile(sigma_ARE, probs = 0.975)
#   )
# 
# # Plotting
# ggplot(stats, aes(x = 1/dt)) +
#   geom_line(aes(y = beta_median, color = "Beta Median")) +
#   geom_ribbon(aes(ymin = beta_q025, ymax = beta_q975, fill = "Beta"), alpha = 0.3) +
#   geom_line(aes(y = mu_median, color = "Mu Median")) +
#   geom_ribbon(aes(ymin = mu_q025, ymax = mu_q975, fill = "Mu"), alpha = 0.3) +
#   geom_line(aes(y = sigma_median, color = "Sigma Median")) +
#   geom_ribbon(aes(ymin = sigma_q025, ymax = sigma_q975, fill = "Sigma"), alpha = 0.3) +
#   scale_color_manual(values = c("Beta Median" = "blue", "Mu Median" = "green", "Sigma Median" = "red")) +
#   scale_fill_manual(values = c("Beta" = "blue", "Mu" = "green", "Sigma" = "red")) +
#   labs(title = "ARE vs Input Size", x = "Input Size (1/dt)", y = "ARE", color = "Parameter", fill = "Parameter") +
#   theme_minimal()
# 
