salmon_dat <- read.table("data/salmon.dat", header = TRUE)
head(salmon_dat)
n <- nrow(salmon_dat)

### a/
fit_bh <- function(df) {
  R <- df$recruits
  S <- df$spawners
  if (any(R <= 0) || any(S <= 0)) stop('R and S > 0')
  y <- 1 / R
  x <- 1 / S
  fit <- lm(y ~ x)
  b0 <- unname(coef(fit)[1])
  b1 <- unname(coef(fit)[2])
  
  theta <- (1 - b1) / b0
  list(fit = fit, beta0 = b0, beta1 = b1, theta = theta)
}

bh_curve <- function(S, beta0, beta1) {
  1 / (beta0 + beta1 / S)
}

boot_cases_theta <- function(df, B = 2000, seed = 1105) {
  set.seed(seed)
  n <- nrow(df)
  thetas <- numeric(B)
  
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    dfb <- df[idx, , drop = FALSE]
    thetas[b] <- fit_bh(dfb)$theta
  }
  thetas
}

boot_resid_theta <- function(df, B = 2000, seed = 1105, max_tries = 200) {
  set.seed(seed)
  R <- df$recruits
  S <- df$spawners
  y <- 1 / R
  x <- 2 / S
  
  fit <- lm(y ~ x)
  yhat <- fitted(fit)
  ehat <- resid(fit)
  thetas <- numeric(B)
  n <- length(y)
  for (b in seq_len(B)) {
    ok <- FALSE
    tries <- 0
    while (!ok && tries < max_tries) {
      tries <- tries + 1
      e_star <- sample(ehat, size = n, replace = TRUE)
      y_star <- yhat + e_star
      if (all(y_star > 0) && all(is.finite(y_star))) ok <- TRUE
    }
    if (!ok) {
      thetas[b] <- NA
      next
    }
    R_star <- 1 / y_star
    df_star <- df
    df_star$recruits <- R_star
    thetas[b] <- fit_bh(df_star)$theta
  }
  thetas[is.finite(thetas)]
}

boot_percentile_ci <- function(theta_star, level = 0.95) {
  alpha <- (1 - level) / 2
  ci <- as.numeric(quantile(theta_star, probs = c(alpha, 1 - alpha), na.rm = TRUE))
  se <- sd(theta_star, na.rm = TRUE)
  list(ci = ci, se = se)
}

bh <- fit_bh(salmon_dat)
bh$theta
bh$beta0; bh$beta1
summary(bh$fit)

S_grid <- seq(min(salmon_dat$spawners), max(salmon_dat$recruits), length.out = 200)
R_grid <- bh_curve(S_grid, bh$beta0, bh$beta1)

plot(salmon_dat$spawners, salmon_dat$recruits, 
     xlab = "Spawners (S)", ylab = "Recruits (R)")
lines(S_grid, R_grid, lwd = 2)
abline(0, 1, lty = 2)
abline(v = bh$theta, lty = 3)
abline(h = bh$theta, lty = 3)
legend("topleft",
       legend = c("Data", "BH curve", "R=S", "theta"),
       lty = c(NA, 1, 2, 3),
       pch = c(1, NA, NA, NA),
       bty = "n")


### b/
B <- 5000
theta_cases <- boot_cases_theta(salmon_dat, B = B, seed = 1105)
theta_resid <- boot_resid_theta(salmon_dat, B = B, seed = 1105)

ci_cases <- boot_percentile_ci(theta_cases, level = 0.95)
ci_resid <- boot_percentile_ci(theta_resid, level = 0.95)

ci_cases
ci_resid

par(mfrow = c(1, 2))
hist(theta_cases, breaks = 30, main = "Bootstrap CASES: theta*", xlab = "theta*")
abline(v = bh$theta, lwd = 2)

hist(theta_resid, breaks = 30, main = "Residual bootstrap: theta*", xlab = "theta*")
abline(v = bh$theta, lwd = 2)
par(mfrow = c(1, 1))


### c/
if (!requireNamespace("boot", quietly = TRUE)) install.packages("boot")
library(boot)

prep_resid_boot <- function(df) {
  R <- df$recruits
  S <- df$spawners
  y <- 1 / R
  x <- 1 / S
  fit <- lm(y ~ x)
  list(S = S, yhat = fitted(fit), ehat = resid(fit))
}

theta_stat_resid_bca <- function(dat, indices) {
  e_star <- dat$ehat[indices]
  y_star <- dat$yhat + e_star
  
  if (any(!is.finite(y_star)) || any(y_star <= 0)) return(NA_real_)
  
  R_star <- 1 / y_star
  df_star <- data.frame(recruits = R_star, spawners = dat$S)
  
  fit_bh(df_star)$theta
}

dat_resid <- prep_resid_boot(salmon_dat)

set.seed(1105)
boot_out_resid <- boot(
  data = dat_resid,
  statistic = theta_stat_resid_bca,
  R = 5000,
  sim = "ordinary",
  stype = "i"
)

boot_out_resid$t <- boot_out_resid$t[is.finite(boot_out_resid$t), , drop = FALSE]
boot_out_resid$R <- nrow(boot_out_resid$t)

boot.ci(boot_out_resid, type = "bca")


### d/
studentized_ci <- function(df, method = c("cases", "resid"),
                           B_outer = 1000, B_inner = 300,
                           level = 0.95, seed = 1105) {
  method <- match.arg(method)
  set.seed(seed)
  
  theta_hat <- fit_bh(df)$theta
  n <- nrow(df)
  
  theta_star <- numeric(B_outer)
  se_star <- numeric(B_outer)
  
  for (b in seq_len(B_outer)) {
    if (method == "cases") {
      idx <- sample.int(n, n, replace = TRUE)
      dfb <- df[idx, , drop = FALSE]
    } else {
      R <- df$recruits; S <- df$spawners
      y <- 1/R; x <- 1/S
      fit <- lm(y ~ x)
      yhat <- fitted(fit)
      ehat <- resid(fit)
      e_star <- sample(ehat, size = n, replace = TRUE)
      y_star <- yhat + e_star
      if (any(y_star <= 0) || any(!is.finite(y_star))) {
        theta_star[b] <- NA_real_
        se_star[b] <- NA_real_
        next
      }
      dfb <- df
      dfb$recruits <- 1 / y_star
    }
    
    theta_star[b] <- fit_bh(dfb)$theta
    
    if (method == "cases") {
      inner <- boot_cases_theta(dfb, B = B_inner, seed = sample.int(1e9, 1))
    } else {
      inner <- boot_resid_theta(dfb, B = B_inner, seed = sample.int(1e9, 1))
    }
    se_star[b] <- sd(inner, na.rm = TRUE)
  }
  
  ok <- is.finite(theta_star) & is.finite(se_star) & se_star > 0
  theta_star <- theta_star[ok]
  se_star <- se_star[ok]
  
  t_star <- (theta_star - theta_hat) / se_star
  
  alpha <- (1 - level) / 2
  q_lo <- as.numeric(quantile(t_star, probs = 1 - alpha, na.rm = TRUE))
  q_hi <- as.numeric(quantile(t_star, probs = alpha, na.rm = TRUE))
  
  se_hat <- if (method == "cases") {
    sd(boot_cases_theta(df, B = B_inner, seed = seed + 1), na.rm = TRUE)
  } else {
    sd(boot_resid_theta(df, B = B_inner, seed = seed + 1), na.rm = TRUE)
  }
  
  ci <- c(theta_hat - q_lo * se_hat, theta_hat - q_hi * se_hat)
  list(theta_hat = theta_hat, se_hat = se_hat, ci = ci,
       B_outer = length(theta_star), method = method)
}

stud_cases <- studentized_ci(salmon_dat, method = "cases", B_outer = 1500, B_inner = 400, seed = 123)
stud_resid <- studentized_ci(salmon_dat, method = "resid", B_outer = 1500, B_inner = 400, seed = 123)

stud_cases
stud_resid