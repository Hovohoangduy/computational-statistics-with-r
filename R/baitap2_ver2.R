kernel_gauss <- function(u) {
    kern <- exp(-0.5 * u^2) / sqrt(2 * pi)
    return(kern)
}

kernel_reg <- function(x, y, x_eval, h, kernel = "gauss") {
    wi <- kernel_gauss((x_eval - x) / h)
    wi_sum <- sum(wi)
    w <- wi / wi_sum
    y_hat <- sum(y * w)
    return(y_hat)
}

kernel_reg <- Vectorize(FUN = kernel_reg, vectorize.args = "x_eval")

# cv1 h
cv1_h <- function(x, y, h, kernel = "gauss") {
  n <- length(x)
  err <- numeric(n)
  for(i in 1:n) {
    x_new <- x[-i]
    y_new <- y[-i]
    y_i <- kernel_reg(x = x_new, y = y_new, x_eval = x[i], h = h, kernel = kernel)
    err[i] <- (y[i] - y_i)^2
  }
  cv_hat <- mean(err)
  return(cv_hat)
}

cv1_h <- Vectorize(FUN = cv1_h, vectorize.args = "h")

# cv2 h
cv2_h <- function(x, y, h, kernel = "gauss") {
  n <- length(x)
  wi0 <- numeric(n)
  for (i in 1:n) {
    wi0[i] <- sum(kernel_gauss(x[i] - x) / h)
  }
  w0 <- kernel_gauss(0) / wi0
  m0 <- kernel_reg(x = x, y = y, x_eval = x, h = h, kernel = kernel)
  cv <- mean(((y - m0) / (1 - w0))^2)
  return(cv)
}

cv2_h <- Vectorize(FUN = cv2_h, vectorize.args = "h")

# GCV
gcv_h <- function(x, y, h, kernel = "gauss") {
  n <- length(x)
  wi0 <- numeric(n)
  for (i in 1:n) {
    wi0[i] <- sum(kernel_gauss((x[i] - x) / h))
  }
  w0 <- kernel_gauss(0) / wi0
  m0 <- kernel_reg(x = x, y = y, x_eval = x, h = h, kernel = kernel)
  gcv <- sum((y - m0)^2) / (n - sum(w0))^2
  return(gcv)
}

gcv_h <- Vectorize(FUN = gcv_h, vectorize.args = "h")

# combined CV - GCV
cv_gcv_h <- function(x, y, h, kernel = "gauss", mode = 3) {
  if (mode == 1) {
    return(cv2_h(x = x, y = y,  h = h, kernel = kernel))
  }
  if (mode == 2) {
    return(gcv_h(x = x, y = y, h = h, kernel = kernel))
  }
  return(
    list(
      cv = cv2_h(x = x, y = y,  h = h, kernel = kernel), 
      gcv = gcv_h(x = x, y = y, h = h, kernel = kernel)
      )
    )
}

## ------------
data(mcycle, package = "MASS")
head(mcycle)

accel_kr_1 <- kernel_reg(x = mcycle$times, y = mcycle$accel, x_eval = 30, h = 0.5)

plot(x = mcycle$times, y = mcycle$accel, pch=16)
points(x = 30, y = accel_kr_1, col = "blue", pch = 4)

## obtain the estimate at several x's.
range(mcycle$times)
x_plot <- seq(0, 60, length.out = 201)
y_hat <- kernel_reg(x = mcycle$times, y = mcycle$accel, x_eval = x_plot,
                    h = 0.5)

plot(x = mcycle$times, y = mcycle$accel, pch = 16, cex = 0.7)
lines(x = x_plot, y = y_hat, col = "blue")

## CV ----
cv1_h(x = mcycle$times, y = mcycle$accel, h = 0.5)
cv1_h(x = mcycle$times, y = mcycle$accel, h = 1.5)
cv1_h(x = mcycle$times, y = mcycle$accel, h = 2.5)

h_plot <- seq(0.1, 4, length.out = 41)
cv1_est <- cv1_h(x = mcycle$times, y = mcycle$accel, h = h_plot)

plot(x = h_plot, y = cv1_est, type = "b", pch = 16)

system.time({
  cv1_est <- cv1_h(x = mcycle$times, y = mcycle$accel, h = h_plot)
})
h_plot[which.min(cv1_est)]

system.time({
  cv2_est <- cv2_h(x = mcycle$times, y = mcycle$accel, h = h_plot)
})

h_plot[which.min(cv2_est)]
all.equal(cv1_est[-(1:2)], cv2_est[-(1:2)])

## GCV ----

system.time({
  gcv_est <- gcv_h(x = mcycle$times, y = mcycle$accel, h = h_plot)
})
h_plot[which.min(gcv_est)]

## benchmark in 100 times
res_banchmark <- microbenchmark::microbenchmark(
  cv1 = cv1_h(x = mcycle$times, y = mcycle$accel, h = h_plot),
  cv2 = cv2_h(x = mcycle$times, y = mcycle$accel, h = h_plot),
  gcv = gcv_h(x = mcycle$times, y = mcycle$accel, h = h_plot)
)

ggplot2::autoplot(res_banchmark)

gcv_h(x = mcycle$times, y = mcycle$accel, h = 0.5)

# calc cv and gcv
cv_gcv_h(x = mcycle$times, y = mcycle$accel, h = 0.5)

# calc cv2
cv_gcv_h(x = mcycle$times, y = mcycle$accel, h = 0.5, mode = 1)

# calc gcv
cv_gcv_h(x = mcycle$times, y = mcycle$accel, h = 0.5, mode = 2)
