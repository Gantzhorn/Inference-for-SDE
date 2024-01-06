DW_functions <- function() {
  list(
    f = function(x, par) par[1] + par[2] * x - par[3] * x^3,
    df = function(x, par) par[2] - 3 * par[3] * x^2,
    dtf = function(x, par) 0,
    d2f = function(x, par) -6 * par[3] * x,
    sigma = function(x, par) exp(par[4]),
    dsigma = function(x, par) 0,
    # SPLITTING SCHEME
    split_A = function(par) par[2],
    split_b = function(par) -par[1] / par[2],
    split_f = function(x, par) -par[3] * x^3,
    split_df = function(x, par) -3 * par[3] * x^2,
    # MISC
    transform = function(par) c(par[1], par[2], par[3], exp(par[4])),
    inverse_transform = function(par) c(par[1:3], log(par[4])),
    nparms = 4
  )
}

DW2_functions <- function() {
  list(
    f = function(x, par) par[1] + par[2] * x - par[3] * x^3,
    df = function(x, par) par[2] - 3 * par[3] * x^2,
    dtf = function(x, par) 0,
    d2f = function(x, par) -6 * par[3] * x,
    sigma = function(x, par) exp(par[4]),
    dsigma = function(x, par) 0,
    # SPLITTING SCHEME
    split_A = function(par) par[2],
    split_b = function(par) 0,
    split_f = function(x, par) -par[3] * x^3 + par[1],
    split_df = function(x, par) -3 * par[3] * x^2,
    # MISC
    transform = function(par) c(par[1], par[2], par[3], exp(par[4])),
    inverse_transform = function(par) c(par[1:3], log(par[4])),
    nparms = 4
  )
}

exp_DW_functions <- function(c) {
  list(
    f = function(x, par) exp(c) * x * (-par[3] * (log(x) + c)^3 + par[2] * (log(x) + c) - par[1] + exp(2 * par[4]) / 2),
    df = function(x, par) exp(c) * (-par[3] * (log(x) + c)^3 + par[2] * (log(x) + c) - par[1] + exp(2 * par[4]) / 2 - 3 * par[3] * (log(x) + c)^2 + par[2]),
    dtf = function(x, par) 0,
    d2f = function(x, par) exp(c) * (-3 * par[3] * (log(x) + c)^2 + par[2] - 6 * par[3] * (log(x) + c)) / x,
    sigma = function(x, par) exp(par[4]) * exp(c) * x,
    dsigma = function(x, par) exp(par[4]) * exp(c),
    # NO SPLITTING SCHEME
    # MISC
    transform = function(par) c(par[1], par[2], par[3], exp(par[4])),
    inverse_transform = function(par) c(par[1:3], log(par[4])),
    nparms = 4
  )
}


# Er håbløst ustabil. Kun dynamiske effekter er interssante 
Log_CIR_functions <- function(c) {
  list(
    f = function(x, par) -exp(par[1]) * (1 - par[2] * exp(-x - c)) - 0.5 * exp(2 * par[3]) * exp(-x - c),
    df = function(x, par) -exp(par[1]) * par[2] * exp(-x - c) + 0.5 * exp(2 * par[3]) * exp(-x - c),
    dtf = function(x, par) 0,
    d2f = function(x, par) exp(par[1]) * par[2] * exp(-x - c) - 0.5 * exp(2 * par[3]) * exp(-x - c),
    sigma = function(x, par) exp(par[3]) * exp(-0.5 * (x + c)),
    dsigma = function(x, par) - 0.5 * exp(par[3]) * exp(-0.5 * (x + c)),
    # NO SPLITTING SCHEME
    # MISC
    transform = function(par) c(exp(par[1]), par[2], exp(par[3])),
    inverse_transform = function(par) c(log(par[1]), par[2], log(par[3])),
    nparms = 3
  )
}

CIR_functions <- function() {
  list(
    f = function(x, par) -exp(par[1]) * (x - par[2]),
    df = function(x, par) -exp(par[1]),
    dtf = function(x, par) 0,
    d2f = function(x, par) 0,
    sigma = function(x, par) exp(par[3]) * sqrt(x),
    dsigma = function(x, par) exp(par[3]) / (2 * sqrt(x)),
    # NO SPLITTING SCHEME
    # MISC
    transform = function(par) c(exp(par[1]), par[2], exp(par[3])),
    inverse_transform = function(par) c(log(par[1]), par[2], log(par[3])),
    nparms = 3
  )
}

Lamperti_CIR_functions <- function() {
  list(
    f = function(x, par) -exp(par[1]) * (x / 2 - 2 * par[2] / x) - exp(2 * par[3]) / (2 * x),
    df = function(x, par) -exp(par[1]) * (0.5 + 2 * par[2] / x^2) + exp(2 * par[3]) / (2 * x^2),
    dtf = function(x, par) 0,
    d2f = function(x, par) (4 * par[2] * exp(par[1]) - exp(2 * par[3])) / x^3,
    sigma = function(x, par) exp(par[3]),
    dsigma = function(x, par) 0,
    # SPLITTING SCHEME
    split_A = function(par) - 0.5 * exp(par[1]),
    split_b = function(par) 0,
    split_f = function(x, par) (4 * exp(par[1]) * par[2] - exp(2 * par[3])) / (2 * x),
    split_df = function(x, par) -(4 * exp(par[1]) * par[2] - exp(2 * par[3])) / (2 * x^2),
    # MISC
    transform = function(par) c(exp(par[1]), par[2], exp(par[3])),
    inverse_transform = function(par) c(log(par[1]), par[2], log(par[3])),
    nparms = 3
  )
}

local_linearization <- function(par, x, fs, dt) {
  x0 <- x[1:(length(x) - 1)]
  x1 <- x[2:length(x)]
  f <- fs$f(x0, par)
  df <- fs$df(x0, par)
  d2f <- fs$d2f(x0, par)
  dtf <- fs$dtf(x0, par)
  sigma <- fs$sigma(x0, par)
  
  r0 <- expm1(df * dt) / df
  LL_mean <- x0 + r0 * f + (r0 - dt / df) * (dtf + 0.5 * sigma^2 * d2f)
  LL_sd <- sigma * sqrt(expm1(2 * df * dt) / (2 * df))
  
  -sum(dnorm(x1, mean = LL_mean, sd = LL_sd, log = TRUE))
}

rlocal_linearization <- function(par, x0, t0, t1, fs, dt) {
  N <- (t1 - t0) / dt
  ff <- fs$f
  fdf <- fs$df
  fd2f <- fs$d2f
  fdtf <- fs$dtf
  fsigma <- fs$sigma
  
  x <- numeric(N + 1)
  x[1] <- x0
  xi <- rnorm(N, mean = 0, sd = 1)
  
  for(i in 2:(N + 1)) {
    f <- ff(x[i - 1], par)
    df <- fdf(x[i - 1], par)
    d2f <- fd2f(x[i - 1], par)
    dtf <- fdtf(x[i - 1], par)
    sigma <- fsigma(x[i - 1], par)
    
    r0 <- expm1(df * dt) / df
    LL_mean <- x[i - 1] + r0 * f + (r0 - dt / df) * (dtf + 0.5 * sigma^2 * d2f)
    LL_sd <- sigma * sqrt(expm1(2 * df * dt) / (2 * df))
    
    x[i] <- LL_mean + LL_sd * xi[i - 1]
  }
  x
}

resid_local_linearization <- function(par, x, fs, dt) {
  x0 <- x[1:(length(x) - 1)]
  x1 <- x[2:length(x)]
  f <- fs$f(x0, par)
  df <- fs$df(x0, par)
  d2f <- fs$d2f(x0, par)
  dtf <- fs$dtf(x0, par)
  sigma <- fs$sigma(x0, par)
  
  r0 <- expm1(df * dt) / df
  LL_mean <- x0 + r0 * f + (r0 - dt / df) * (dtf + 0.5 * sigma^2 * d2f)
  LL_sd <- sigma * sqrt(expm1(2 * df * dt) / (2 * df))
  
  qqnorm(pnorm(x1, mean = LL_mean, sd = LL_sd))
}

euler_maruyama <- function(par, x, fs, dt) {
  x0 <- x[1:(length(x) - 1)]
  x1 <- x[2:length(x)]
  f <- fs$f(x0, par)
  sigma <- fs$sigma(x0, par)
  EM_mean <- x0 + f * dt
  EM_sd <- sigma * sqrt(dt)
  
  -sum(dnorm(x1, mean = EM_mean, sd = EM_sd, log = TRUE))
}

reuler_maruyama <- function(par, x0, t0, t1, fs, dt) {
  N <- (t1 - t0) / dt
  f <- fs$f
  sigma <- fs$sigma
  x <- numeric(N + 1)
  x[1] <- x0
  xi <- rnorm(N, mean = 0, sd = sqrt(dt))
  for(i in 2:(N + 1)) {
    x[i] <- x[i - 1] + f(x[i - 1], par) * dt + sigma(x[i - 1], par) * xi[i - 1]
  }
  x
}

resid_euler_maruyama <- function(par, x, fs, dt) {
  x0 <- x[1:(length(x) - 1)]
  x1 <- x[2:length(x)]
  f <- fs$f(x0, par)
  sigma <- fs$sigma(x0, par)
  EM_mean <- x0 + f * dt
  EM_sd <- sigma * sqrt(dt)
  
  qqnorm(pnorm(x1, mean = EM_mean, sd = EM_sd))
}

rmilstein <- function(par, x0, t0, t1, fs, dt) {
  N <- (t1 - t0) / dt
  f <- fs$f
  sigma <- fs$sigma
  dsigma <- fs$dsigma
  x <- numeric(N + 1)
  x[1] <- x0
  xi <- rnorm(N, mean = 0, sd = sqrt(dt))
  for(i in 2:(N + 1)) {
    fsigma <- sigma(x[i - 1], par)
    fdsigma <- dsigma(x[i - 1], par)
    x[i] <- x[i - 1] + f(x[i - 1], par) * dt + fsigma * xi[i - 1] +
      0.5 * fdsigma * fsigma * (xi[i - 1]^2 - dt)
  }
  x
}
# Kun til CIR-modellen
martingale <- function(x, dt) {
  N <- length(x)
  M <- length(x) - 1
  X0 <- x[1:M]
  X1 <- x[2:N]
  S1 <- sum(X1)
  S2 <- sum(1 / X0)
  S3 <- sum(X1 / X0)
  
  ebdel <- (M * S3 - S1 * S2) / (M^2 - S1 * S2)
  beta <- -log(ebdel) / dt
  
  mu <- (S3 - M * ebdel) / ((1 - ebdel) * S2)
  
  S4 <- sum((X1 - mu + ebdel * (mu - X0))^2 / X0)
  S5 <- sum((X1 * (ebdel - ebdel^2) / beta + mu * (1 - ebdel)^2 / (2 * beta)) / X1)
  
  sugma <- S4 / S5
  
  c(beta = beta, mu = mu, sigma = sqrt(sugma))
}

# Only for CIR
closedform_likelihood_CIR <- function(par, x, t, fs, dt){
  # Initialize parameters
  beta <- exp(par[1])
  mu <- par[2]
  sigma <- exp(par[3])
  
  # Prepare data
  N <- length(x)
  x <- 4 * beta / (sigma^2 * (1 - exp(-beta * t))) * x
  lower_x <- x[1:(N - 1)]
  upper_x <- x[2:N]
  
  -sum(dchisq(upper_x, df = 4 * mu * beta / sigma^2, ncp = lower_x * exp(-beta * dt), log = TRUE))
}

runge_kutta_nystrom <- function(y0, dy0, h, f, n = 1) {
  h <- h / n
  
  for(i in 1:n) {
    k1 <- f(0, dy0, y0)
    
    dy1 <- dy0 + 0.5 * h * k1
    y1 <- y0 + 0.5 * h * (dy0 + dy1) / 2
    k2 <- f(0.5 * h, dy1, y1)
    
    dy2 <- dy0 + 0.5 * h * k2
    y2 <- y0 + 0.5 * h * (dy0 + dy2) / 2
    k3 <- f(0.5 * h, dy2, y2)
    
    dy3 <- dy0 + k3 * h
    y3 <- y0 + h * (dy0 + dy3) / 2
    k4 <- f(h, dy3, y3)
    
    dy0 <- dy0 + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    y0 <- y0 + h * (dy1 + 2 * dy2 + 2 * dy3 + dy0) / 6
  }
  
  list(y = y0, dy = dy0)
}

runge_kutta <- function(y0, h, f, n = 1) {
  h <- h / n
  
  for(i in 1:n) {
    k1 <- f(0, y0)
    k2 <- f(0.5 * h, y0 + 0.5 * h * k1)
    k3 <- f(0.5 * h, y0 + 0.5 * h * k2)
    k4 <- f(h, y0 + h * k3)
    
    y0 <- y0 + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
  }
  y0
}

strang <- function(par, x, fs, dt) {
  x0 <- x[1:(length(x) - 1)]
  x1 <- x[2:length(x)]
  A <- fs$split_A(par)
  b <- fs$split_b(par)
  
  sigma <- fs$sigma(x0, par)
  
  diff_f <- function(t, y) fs$split_f(y, par)
  
  inv_f <- runge_kutta(x1, -dt / 2, diff_f)
  inv_f2 <- runge_kutta(x1 + 0.01, -dt / 2, diff_f)
  inv_f3 <- runge_kutta(x1 - 0.01, -dt / 2, diff_f)
  # inv_f4 <- runge_kutta(x1 + 2 * 0.01, -dt / 2, diff_f)
  # inv_f5 <- runge_kutta(x1 - 2 * 0.01, -dt / 2, diff_f)
  f <- runge_kutta(x0, dt / 2, diff_f)
  df <- (inv_f2 - inv_f3) / (2 * 0.01) # Richardson Extrapolation
  # df <- (inv_f5 - 8 * inv_f3 + 8 * inv_f2 - inv_f4) / (12 * 0.01) # Richardson Extrapolation 2
  # df <- (inv_f2 - inv_f) / 0.1
  
  mu <- exp(A * dt) * (f - b) + b
  omega <- sigma * sqrt(expm1(2 * A * dt) / (2 * A))
  
  -sum(dnorm(inv_f, mean = mu, sd = omega, log = TRUE)) - sum(log(abs(df)))
}

rstrang <- function(par, x0, t0, t1, fs, dt) {
  N <- (t1 - t0) / dt
  sigma <- fs$sigma
  A <- fs$split_A(par)
  b <- fs$split_b(par)
  sigma <- fs$sigma
  diff_f <- function(t, y) fs$split_f(y, par)
  
  x <- numeric(N + 1)
  x[1] <- x0
  xi <- rnorm(N, mean = 0, sd = 1)
  
  for(i in 2:(N + 1)) {
    fsigma <- sigma(x[i - 1], par)
    f <- runge_kutta(x[i - 1], dt / 2, diff_f)
    mu <- exp(A * dt) * (f - b) + b
    omega <- fsigma * sqrt(expm1(2 * A * dt) / (2 * A))
    x[i] <- runge_kutta(mu + omega * xi[i - 1], dt / 2, diff_f)
  }
  x
}

resid_strang <- function(par, x, fs, dt) {
  x0 <- x[1:(length(x) - 1)]
  x1 <- x[2:length(x)]
  A <- fs$split_A(par)
  b <- fs$split_b(par)
  
  sigma <- fs$sigma(x0, par)
  
  diff_f <- function(t, y) fs$split_f(y, par)
  
  inv_f <- runge_kutta(x1, -dt / 2, diff_f)
  inv_f2 <- runge_kutta(x1 + 0.01, -dt / 2, diff_f)
  inv_f3 <- runge_kutta(x1 - 0.01, -dt / 2, diff_f)
  # inv_f4 <- runge_kutta(x1 + 2 * 0.01, -dt / 2, diff_f)
  # inv_f5 <- runge_kutta(x1 - 2 * 0.01, -dt / 2, diff_f)
  f <- runge_kutta(x0, dt / 2, diff_f)
  df <- (inv_f2 - inv_f3) / (2 * 0.01) # Richardson Extrapolation
  # df <- (inv_f5 - 8 * inv_f3 + 8 * inv_f2 - inv_f4) / (12 * 0.01) # Richardson Extrapolation 2
  # df <- (inv_f2 - inv_f) / 0.1
  
  mu <- exp(A * dt) * (f - b) + b
  omega <- sigma * sqrt(expm1(2 * A * dt) / (2 * A))
  
  qqnorm(pnorm(inv_f, mean = mu, sd = omega))
}

lie_trotter <- function(par, x, fs, dt) {
  x0 <- x[1:(length(x) - 1)]
  x1 <- x[2:length(x)]
  A <- fs$split_A(par)
  b <- fs$split_b(par)
  
  sigma <- fs$sigma(x0, par)
  
  diff_f <- function(t, y) fs$split_f(y, par)
  f <- runge_kutta(x0, dt, diff_f)
  
  mu <- exp(A * dt) * (f - b) + b
  omega <- sigma * sqrt(expm1(2 * A * dt) / (2 * A))
  
  -sum(dnorm(x1, mean = mu, sd = omega, log = TRUE))
}

rlie_trotter <-  function(par, x0, t0, t1, fs, dt) {
  N <- (t1 - t0) / dt
  sigma <- fs$sigma
  A <- fs$split_A(par)
  b <- fs$split_b(par)
  sigma <- fs$sigma
  
  diff_f <- function(t, y) fs$split_f(y, par)
  
  x <- numeric(N + 1)
  x[1] <- x0
  xi <- rnorm(N, mean = 0, sd = 1)
  
  for(i in 2:(N + 1)) {
    fsigma <- sigma(x[i - 1], par)
    f <- runge_kutta(x[i - 1], dt, diff_f)
    mu <- exp(A * dt) * (f - b) + b
    omega <- fsigma * sqrt(expm1(2 * A * dt) / (2 * A))
    x[i] <- mu + omega * xi[i - 1]
  }
  x
}

resid_lie_trotter <- function(par, x, fs, dt) {
  x0 <- x[1:(length(x) - 1)]
  x1 <- x[2:length(x)]
  A <- fs$split_A(par)
  b <- fs$split_b(par)
  
  sigma <- fs$sigma(x0, par)
  
  diff_f <- function(t, y) fs$split_f(y, par)
  f <- runge_kutta(x0, dt, diff_f)
  
  mu <- exp(A * dt) * (f - b) + b
  omega <- sigma * sqrt(expm1(2 * A * dt) / (2 * A))
  
  qqnorm(pnorm(x1, mean = mu, sd = omega))
}