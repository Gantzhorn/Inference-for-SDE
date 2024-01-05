library(tidyverse)
theme_set(theme_bw())

sim_doublewell <- function(total_time, step_length, param, X_0){
  
  beta1 <- param[1]
  beta2 <- param[2]
  beta4 <- param[3]
  sigma <- param[4]
  
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
    X[t] <- X[t - 1] - (beta4 * X[t - 1]^3 - beta2 * X[t - 1] - beta1) * step_length + sigma * dW[t - 1]
  }
  
  X_ito <- numeric(N + 1)
  X_ito[1] <- exp(-X_0)
  
  for (t in 2:(N + 1)) {
    X_ito[t] <- X_ito[t - 1] + X_ito[t - 1] * (beta2 * log(X_ito[t - 1]) - beta4 * log(X_ito[t - 1])^3 - beta1 + 1/2 * sigma^2) * step_length - 
      sigma * X_ito[t - 1] * dW[t - 1]
  }
  
  t_vec <- seq(0, length.out = N + 1, by = step_length)
  
  # Create tibble
  tibble(t = t_vec, X = X, X_ito_manual = exp(-X), X_ito)
}

beta1 <- 1
beta2 <- 5
beta4 <- 1
sigma <- 3
x0 <- 1
actual_dt <- 0.0001
trueparam <- c(beta1, beta2, beta4, sigma)


sample_doublewell <- sim_doublewell(10, actual_dt, trueparam, x0)

sample_doublewell %>% ggplot(aes(x = t, y = X)) + geom_step()

sample_doublewell %>% ggplot(aes(x = t, y = X_ito)) + geom_step()

euler_maruyama_likelihood_doublewell <- function(param, data, dt){
  # Initialize parameters
  beta1 <- param[1]
  beta2 <- param[2]
  beta4 <- param[3]
  sigma <- exp(param[4])
  
  # Prepare data
  N <- length(data)
  lower_x <- data[1:(N - 1)]
  upper_x <- data[2:N]
  
  mu <- lower_x - (beta4 * lower_x^3 - beta2 * lower_x - beta1) * dt
  sd <- sigma * sqrt(dt)
  # likelihood_value <- 1/2 * sum(log(2*pi*sigma^2 * dt) + 
  #                                 (upper_x - lower_x - (beta4 * lower_x^3 - beta2 * lower_x - beta1) * dt)^2 /
  #                                 (sigma^2 * dt))
  -sum(dnorm(upper_x, mean = mu, sd = sd, log = TRUE))
  # likelihood_value
}

LT_likelihood_double_well <- function(param, data, dt){
  
  beta1 <- param[1]
  beta2 <- param[2]
  beta4 <- param[3]
  sigma <- param[4]
  
  
  N <- length(data)
  lower_x <- data[1:(N - 1)]
  upper_x <- data[2:N]
  
  f_delta_t <- 1 / sqrt(2 * beta4 * dt + lower_x^(-2))
  
  mu <- exp(beta2 * dt) * (f_delta_t + beta1 / beta2) - beta1 / beta2
  
  omega  <- sigma^2 / (2 * beta2) * (exp(2 * beta2 * dt) - 1)
  
  -sum(dnorm(upper_x, mean = mu, sd = sqrt(omega), log = TRUE))
}


optim(par = trueparam, fn = LT_likelihood_double_well,
      data = sample_doublewell$X, dt = actual_dt)

Strang_likelihood_double_well <- function(param, data, dt){
  
  beta1 <- param[1]
  beta2 <- param[2]
  beta4 <- param[3]
  sigma <- param[4]
  
  
  N <- length(data)
  lower_x <- data[1:(N - 1)]
  upper_x <- data[2:N]
  
  f_delta_t <- 1 / sqrt(beta4 * dt + lower_x^(-2))
  
  f_inverse_delta_t <- 1 / sqrt(1/upper_x^2 - beta4 * dt)
  
  f_prime_inverse_delta_t <- 1 / (upper_x^3 * (1 / upper_x^2 - beta4 * dt)^(3/2)) 
  
  mu <- exp(beta2 * dt) * (f_delta_t + beta1 / beta2) - beta1 / beta2
  
  omega  <- sigma^2 / (2 * beta2) * (exp(2 * beta2 * dt) - 1)
  
  -sum(dnorm(f_inverse_delta_t, mean = mu, sd = sqrt(omega), log = TRUE)) - sum(log(abs(f_prime_inverse_delta_t)))
}

optim(par = trueparam, fn = Strang_likelihood_double_well,
      data = sample_doublewell$X, dt = actual_dt)

euler_maruyama_likelihood_transformed_doublewell <- function(param, data, dt){
  # Initialize parameters
  beta1 <- param[1]
  beta2 <- param[2]
  beta4 <- param[3]
  sigma <- param[4]
  
  # Prepare data
  N <- length(data)
  lower_x <- data[1:(N - 1)]
  upper_x <- data[2:N]
  
likelihood_value <- 1/2 * sum(log(2*pi*sigma^2 * lower_x^2 * dt) + 
                                  (upper_x - lower_x * (1 + (beta2 * log(lower_x) - beta4 * log(lower_x)^3 - beta1 + sigma^2/2) * dt))^2 /
                                  (sigma^2 * lower_x^2 * dt))
  likelihood_value
}

## Optimization

result_euler_maruyama_likelihood_doublewell <-  optim(par = trueparam, fn = euler_maruyama_likelihood_doublewell,
                                                      data = sample_doublewell$X, dt = actual_dt)

result_euler_maruyama_likelihood_doublewell$par


result_euler_maruyama_likelihood_transformed_doublewell <-  optim(par = trueparam, fn = euler_maruyama_likelihood_transformed_doublewell,
                                                      data = sample_doublewell$X_ito, dt = actual_dt)

result_euler_maruyama_likelihood_transformed_doublewell$par
