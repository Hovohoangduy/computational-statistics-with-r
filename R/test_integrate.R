source(file="R/integ_trap.R")

x <- seq(0, pi, length.out = 101)
f1 <- function(x) x * sin(x)
integ_trap(x = x, f = f1)

integrate(f = f1, lower = 0, upper = pi)

## monte carlo
x_uniform <- runif(n = 10000, 0, pi)
f1_uni <- pi*f1(x_uniform)
mean(f1_uni)

f2 <- function(x, mu, sigma) {
  t <- 1 - x
  res <- exp(-0.5 * (1/sigma^2) * (log(x/t) - mu)^2)/
    (sqrt(2*pi)*sigma * x * t)
  return(res)
}

a <- 0
x2 <- seq(10^(-9), 1/(1+exp(-a)), length.out = 101)
integ_trap(x = x2, f = f2, mu = 2, sigma = 1.2)
pnorm(0, mean = 2, sd = 1.2)

f2_1 <- f2(x = x2, mu = 2, sigma = 1.2)

x_uni_2 <- runif(n = 1000, 0, 1/(1+exp(-a)))
f2_uni <- (1/(1+exp(-a))) * f2(x=x_uni_2, mu = 2, sigma = 1.2)
mean(f2_uni)

## pi estimate
pi_est <- function(n, seed = 1105) {
  set.seed(seed)
  x1 <- runif(n, -1, 1)
  x2 <- runif(n, -1, 1)
  res <- 4*mean(x1^2 + x2^2 <= 1)
  return(res)
}
pi1 <- pi_est(n=400)
pi1

pi_est_n <- sapply(c(100, 200, 500, 1000), function(x) pi_est(x))
pi_est_n

n_sim <- 50
u_sim <- runif(n_sim)
x_2 <- -log(1 - u_sim*(1 - exp(-1)))
I2_est <- (1 - exp(-1)) * mean(1 / (1 + x_2^2))
I2_est
