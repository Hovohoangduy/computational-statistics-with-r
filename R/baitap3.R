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
  if (sum(w > 1e-6) < (p + 1)) return(NULL)
  W = diag(w)
  if (rcond(A) < 1e-12) return(NULL)
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
    if (is.null(fit)) return(Inf)
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
    h_min <- diff(range(x)) / 5
    h_max <- diff(range(x))
    h_grid <- seq(h_min, h_max, length.out = 30)
    if (method == "CV") {
      cv_vals <- cv_lp(x, y, h_grid, p)
      h <- h_grid[which.min(cv_vals)]
    }
    
    if (method == "GCV") {
      gcv_vals <- gcv_lp(x, y, h_grid, p)
      h <- h_grid[which.min(gcv_vals)]
    }
  }
  m_hat <- lp_reg(x, y, x_eval, h, p)
  return(list(m_hat = m_hat, h = h, p = p))
}



### -----------
x <- sort(runif(100, 0, 10))
y <- sin(x) * rnorm(100, 0, 0.2)
lp_reg(x = x, y = y, x_eval = 5, h = 0.5, p = 1)

fit_center <- poly_reg(x = x, y = y, x_eval = 5, h = NULL, p = 1)
h_opt <- fit_center$h

x_grid <- seq(min(x), max(x), length.out = 200)
m_hat <- sapply(x_grid, function(x0) {
  fit <- lp_reg(
    x = x,
    y = y,
    x_eval = x0,
    h = h_opt,
    p = 1
  )
  if (is.null(fit)) NA else fit$beta[1]
})
par(mar = c(5, 4, 4, 6))  # chừa chỗ cho legend

plot(
  x, y,
  pch = 16,
  col = rgb(0, 0, 0, 0.5),
  xlab = "x",
  ylab = "y",
  main = "Local Linear Regression (p = 1)"
)

lines(
  x_grid, m_hat,
  col = "red",
  lwd = 2
)

