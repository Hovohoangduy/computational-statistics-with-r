set.seed(123)

n <- 300
y <- c(
  rnorm(0.3 * n, mean = -1, sd=sqrt(0.1)),
  rnorm(0.7 * n, mean = 1, sd = 0.1)
)

kernel_gaussian <- function(u) {
  dnorm(u)
}

kernel_epanechinikov <- function(u) {
  0.75 * (1 - u^2) * (abs(u) <= 1)
}

kernel_manual <- function(x, data, h, kernel = kernel_gaussian) {
  n_x <- length(x)
  f_hat <- numeric(n_x)
  for (i in seq_len(n_x)) {
    x0 <- x[i]
    f_hat[i] <- mean(kernel(x0 - data) / h) / h
  }
  return(f_hat)
}

x_grid <- seq(-4, 4, length.out = 500)
h <- 0.3

f_hat <- kernel_manual(x_grid, y, h)
print(f_hat)

plot(x_grid, f_hat, type = "l", lwd = 2, col = "blue", xlab = "x", ylab = "Density", main = "Kernel Density Estimation")
hist(y, prob = TRUE, col = "gray90", border = "white", add = TRUE)