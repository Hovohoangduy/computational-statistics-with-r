# KNN homework
## kernel regression

# gaussian kernel function
kernel_gauss <- function(u) exp(-u^2/2) / sqrt(2*pi)

# kernel regression
kernel_reg <- function(x, y, x_eval, h, kernel = "gauss") {
  wi <- kernel_gauss((x_eval - x) / h)
  wi_sum <- sum(wi)
  wi <- wi / wi_sum
  y_hat <- sum(y * wi)
  return(y_hat)
}

kernel_reg <- Vectorize(FUN = kernel_reg, vectorize.args = "x_eval")

# cross-validation
cv1_h <- function(x, y, h, kernel = "gauss") {
  n <- length(x)
  err <- numeric(n) 
  for (i in 1:n) {
    x_new = x [-i]
    y_new <- y[-i]
    y_i <- kernel_reg(x = x_new, y = y_new, x_eval = x[i], h = h, kernel = kernel)
    err[i] <- (y[i] - y_i)^2

  }
  cv_hat <- mean(err)
  return(cv_hat)
}

cv1_h <- Vectorize(FUN = cv1_h, vectorize.args = "h")

kernel_reg_all <- function(x, y, h, kernel = "gauss") {
  
  n <- length(x)
  m_hat <- numeric(n)
  
  for (i in 1:n) {
    w <- kernel_gauss((x - x[i]) / h)
    m_hat[i] <- sum(w * y) / sum(w)
  }
  
  return(m_hat)
}

cv2_h <- function(x, y, h, kernel = "gauss") {
  n <- length(x)
  m_hat <- kernel_reg_all(x, y, h)
  s_ii <- numeric(n)
  for (i in 1:n) {
    w <- kernel_gauss((x - x[i]) / h)
    s_ii[i] <- kernel_gauss(0) / sum(w)
  }

  mean(((y - m_hat) / (1 - s_ii))^2)
}

cv2_h <- Vectorize(FUN = cv2_h, vectorize.args = "h")
