data <- log10(lynx)

n <- length(data)
init_train <- 50
errors_ar1 <- numeric(n - init_train)
errors_ar2 <- numeric(n - init_train)

for (i in init_train:(n - 1)) {
    train_data <- data[1:i]

    actual_val <- data[i + 1]

    fit_ar1 <- arima(train_data, order = c(1, 0, 0), method = "ML")
    fit_ar2 <- arima(train_data, order = c(2, 0, 0), method = "ML")

    pred_ar1 <- predict(fit_ar1, n.ahead = 1)$pred[1]
    pred_ar2 <- predict(fit_ar2, n.ahead = 1)$pred[1]

    errors_ar1[i - init_train + 1] <- actual_val - pred_ar1
    errors_ar2[i - init_train + 1] <- actual_val - pred_ar2
}

rmse_ar1 <- sqrt(mean(errors_ar1^2))
rmse_ar2 <- sqrt(mean(errors_ar2^2))

par(mfrow = c(2, 1), mar = c(4, 4, 2, 2))
plot(errors_ar2, type = "l", main = "Sai số dự báo 1-bước ngẫu nhiên của mô hình AR(2) qua Time Series CV", ylab = "CV Out-of-Sample Error", col = "blue")
abline(h = 0, col = "red", lty = 2)

best_model <- arima(data, order = c(2, 0, 0), method = "ML")
resid <- na.omit(best_model$residuals)

n_ahead <- 15
B <- 1000

sieve_boot_forecast <- function(model, residuals, h) {
    res_boot <- sample(residuals, size = h, replace = TRUE)

    coef_ar <- model$coef[1:2]
    intercept <- model$coef["intercept"]

    sim_path <- numeric(h)

    y_prev <- tail(data, 2)

    for (t in 1:h) {
        y_t <- intercept + coef_ar[1] * (y_prev[2] - intercept) +
            coef_ar[2] * (y_prev[1] - intercept) + res_boot[t]

        sim_path[t] <- y_t

        y_prev[1] <- y_prev[2]
        y_prev[2] <- y_t
    }
    return(sim_path)
}

set.seed(42)
boot_paths <- replicate(B, sieve_boot_forecast(best_model, resid, n_ahead))

lower_bound <- apply(boot_paths, 1, quantile, probs = 0.025)
upper_bound <- apply(boot_paths, 1, quantile, probs = 0.975)
mean_bound <- apply(boot_paths, 1, mean)

par(mfrow = c(1, 1), mar = c(4, 4, 4, 2))

plot_len <- 30
time_axis_past <- (n - plot_len + 1):n
time_axis_fut <- (n + 1):(n + n_ahead)

plot(time_axis_past, data[time_axis_past],
    type = "l", lwd = 2,
    xlim = c(n - plot_len, n + n_ahead), ylim = c(range(data)[1], max(upper_bound) + 0.5),
    ylab = "Số linh miêu (Log10)", xlab = "Thời điểm (Thời gian thực / Chỉ số mẫu)",
    main = "Sieve Bootstrap tạo Dải khoảng dự báo AR(2)"
)

polygon(c(time_axis_fut, rev(time_axis_fut)),
    c(lower_bound, rev(upper_bound)),
    col = rgb(0, 0, 1, 0.2), border = NA
)

lines(time_axis_fut, mean_bound, col = "red", lwd = 2, lty = 1)

legend("topleft",
    legend = c("Dữ liệu gốc (Gần đây)", "Dự báo điểm (Lấy trung bình)", "Khoảng tin cậy thực nghiệm 95% (Boot)"),
    col = c("black", "red", NA), fill = c(NA, NA, rgb(0, 0, 1, 0.2)),
    border = c(NA, NA, NA), lty = c(1, 1, NA), lwd = c(2, 2, NA), pch = NA
)
