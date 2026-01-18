source(file = "R/kernel_regression.R")

# mcycle data load
data(mcycle, package = "MASS")
head(mcycle)


u1 <- kernel_reg(x = mcycle$times, y = mcycle$accel, x_eval = 30, h = 0.5)
plot(x = mcycle$times, y = mcycle$accel, pch = 16)
points(x = 30, y = u1, col = "blue")

range(mcycle$time)

x_plot <- seq(0, 69, length.out = 201)
y_hat <- kernel_reg(x = mcycle$times, y = mcycle$accel, x_eval = x_plot, h = 0.5)
plot(x = mcycle$times, y = mcycle$accel, pch = 16, cex = 0.7)
lines(x = x_plot, y = y_hat, col = "blue")

## CV
cv1_h(x = mcycle$times, y = mcycle$accel, h = 0.5)
cv1_h(x = mcycle$times, y = mcycle$accel, h = 1.5)
cv1_h(x = mcycle$times, y = mcycle$accel, h = 2.5)

h_plot <- seq(0.1, 4, length.out = 22)
cv1_est <- cv1_h(x = mcycle$times, y = mcycle$accel, h = h_plot)
plot(x = h_plot, y = cv1_est, type = "o", pch = 16)

system.time({
  cv1_est <- cv1_h(x = mcycle$times, y = mcycle$accel, h = h_plot)
})


## CV2
cv2_h(x = mcycle$times, y = mcycle$accel, h = 0.5)
cv2_h(x = mcycle$times, y = mcycle$accel, h = 1.5)
cv2_h(x = mcycle$times, y = mcycle$accel, h = 2.5)

h_plot <- seq(0.1, 4, length.out = 22)
cv2_est <- cv2_h(x = mcycle$times, y = mcycle$accel, h = h_plot)
plot(x = h_plot, y = cv2_est, type = "o", pch = 16)

system.time({
  cv2_est <- cv1_h(x = mcycle$times, y = mcycle$accel, h = h_plot)
})

