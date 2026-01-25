kernel_gauss <- function(u) {
  exp(-0.5 * u^2) / sqrt(2 * pi)
}

cv_gcv_kernel <- function(x, y, h, method = 3) {
  n <- length(y)

  dist_mat <- outer(x, x, "-")
  K_mat <- kernel_gauss(dist_mat / h)
  row_sums <- rowSums(K_mat)
  S <- K_mat / row_sums
 
  y_hat <- S %*% y

  # Initialize results
  cv_value <- NA
  gcv_value <- NA

  # Diagonal elements of S
  s_ii <- diag(S)

  # calc CV (Leave-One-Out Cross-Validation)
  if (method == 1 || method == 3) {
    # Using the shortcut formula: CV = 1/n * sum( ((y_i - y_hat_i) / (1 - s_ii))^2 )
    # We must handle cases where s_ii is very close to 1
    denom_cv <- 1 - s_ii
    valid_cv <- denom_cv > 1e-10
    if (any(valid_cv)) {
      cv_val_vector <- ((y - y_hat) / denom_cv)^2
      cv_value <- mean(cv_val_vector[valid_cv])
      # If some points are invalid, CV might be misleadingly low or just NA/Inf
      if (any(!valid_cv)) cv_value <- Inf
    } else {
      cv_value <- Inf
    }
  }

  # calc GCV (Generalized Cross-Validation)
  if (method == 2 || method == 3) {
    rss <- sum((y - y_hat)^2)
    trS <- sum(s_ii)
    denom_gcv <- (1 - trS / n)^2
    if (denom_gcv > 1e-10) {
      gcv_value <- (rss / n) / denom_gcv
    } else {
      gcv_value <- Inf
    }
  }

  if (method == 1) {
    return(cv_value)
  }
  if (method == 2) {
    return(gcv_value)
  }
  return(list(cv = cv_value, gcv = gcv_value))
}

set.seed(123)
n <- 30
x <- sort(runif(n))
y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)

# h_grid: Bắt đầu từ h rất nhỏ để chứng minh tính không ổn định của CV
# Khi h nhỏ, ít nhất một điểm x_i sẽ bị cô lập khiến s_ii gần bằng 1.
# CV bị ảnh hưởng bởi từng s_ii riêng lẻ nên dễ bị Inf nhanh hơn.
# GCV dùng trung bình (trS/n) nên ổn định hơn khi h nhỏ.
h_grid <- 10^seq(-3, -0.5, length.out = 100)

cv_vals <- numeric(length(h_grid))
gcv_vals <- numeric(length(h_grid))

for (k in 1:length(h_grid)) {
  res <- cv_gcv_kernel(x, y, h_grid[k], method = 3)
  cv_vals[k] <- res$cv
  gcv_vals[k] <- res$gcv
}

# Kiểm tra dữ liệu
cv_finite <- is.finite(cv_vals)
gcv_finite <- is.finite(gcv_vals)

cat("Số điểm CV hữu hạn:", sum(cv_finite), "/", length(h_grid), "\n")
cat("Số điểm GCV hữu hạn:", sum(gcv_finite), "/", length(h_grid), "\n")

# Lọc các điểm có thể vẽ (GCV hữu hạn và ko quá lớn để giữ Scale đồ thị đẹp)
plot_idx <- gcv_finite & (gcv_vals < 0.5)

if (sum(plot_idx) > 0) {
  h_plot <- h_grid[plot_idx]
  cv_plot <- cv_vals[plot_idx]
  gcv_plot <- gcv_vals[plot_idx]

  # Xác định ylim dựa trên các giá trị hữu hạn
  y_max <- max(c(cv_plot[is.finite(cv_plot)], gcv_plot), na.rm = TRUE)

  plot(h_plot, gcv_plot,
    type = "l",
    log = "x",
    col = "blue",
    lwd = 2,
    ylim = c(0, min(y_max, 0.2)), # Giới hạn lại để nhìn rõ vùng optimal
    xlab = "Bandwidth h (log scale)",
    ylab = "Score",
    main = "CV mất ổn định nhanh hơn GCV khi h nhỏ"
  )

  lines(h_plot, cv_plot,
    type = "l",
    col = "red",
    lwd = 2
  )

  legend("topright",
    legend = c("CV (LOOCV) - Dễ bị Inf/NaN", "GCV - Ổn định hơn"),
    col = c("red", "blue"),
    lty = 1,
    lwd = 2
  )

  # Đánh dấu h tối ưu theo GCV
  opt_h_gcv <- h_plot[which.min(gcv_plot)]
  abline(v = opt_h_gcv, col = "darkgreen", lty = 2)
  cat("Optimal h (GCV):", opt_h_gcv, "\n")
}
