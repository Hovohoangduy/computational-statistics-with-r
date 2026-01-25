kernel_gauss <- function(u) {
  kern <- exp(-0.5 * u^2) / sqrt(2 * pi)
  return(kern)
}

lp_reg <- function(x, y, x_eval, h, p, kernel = "gauss") {
  n <- length(x)
  Z <- matrix(1, n, p + 1)
  for (j in 1:p) {
    Z[, j + 1] <- (x - x_eval)^j
  }
  w <- kernel_gauss((x - x_eval) / h)
  W = diag(w)
  M <- solve(t(Z) %*% W %*% Z) %*% t(Z) %*% W
  beta_hat <- M %*% y
  return(list(beta = beta_hat, S_row = M[1,]))
}

cv_lp <- function(x, y, h, p) {
  n <- length(y)
  y_hat <- numeric(n)
  S_diag <- numeric(n)
  for (i in 1:n) {
    fit <- lp_reg(x, y, x_eval = x[i], h = h, p = p)
    y_hat[i] <- fit$beta[1]
    S_diag[i] <- fit$S_row[i]
  }
  mean(((y - y_hat) / (1 - S_diag))^2)
}

cv_lp <- Vectorize(FUN = cv_lp, vectorize.args = "h")

# GCV
gcv_lp <- function(x, y, h, p) {
  n <- length(x)
  y_hat <- numeric(n)
  trace_S <- 0
  for (i in 1:n) {
    fit <- lp_reg(x, y, x_eval = x[i], h = h, p = p)
    y_hat[i] <- fit$beta[1]
    trace_S <- trace_S + fit$S_row[i]
  }
  sum((y - y_hat)^2) / (n - trace_S)^2
}

gcv_lp = Vectorize(FUN = gcv_lp, vectorize.args = "h")

poly_reg <- function(x, y, x_eval, h = NULL, p, method = c("CV", "GCV")) {
  method <- match.arg(method)
  
  if (is.null(h)) {
    h_grid <- seq(sd(x) / 10, sd(x), length.out = 30)
    if (method == "CV") {
      cv_vals <- sapply(h_grid, cv_lp, x = x, y = y, p = p)
      h <- h_grid[which.min(cv_vals)]
    }
    
    if (method == "GCV") {
      gcv_vals <- sapply(h_grid, cv_lp, x = x, y = y, p = p)
      h <- h_grid(which.min(cv_vals))
    }
  }
  m_hat <- lp_reg(x, y, x_eval, h, p)
  return(list(m_hat = m_hat, h = h, p = p))
}



### -----------
x <- sort(runif(100, 0, 10))
y <- sin(x) * rnorm(100, 0, 0.2)
lp_reg(x = x, y = y, x_eval = 5, h = 0.5, p = 1)

poly_reg(x = x, y = y, x_eval = 5, h = NULL, p = 1)


