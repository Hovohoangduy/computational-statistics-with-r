### bootstrap baitap 3

boot_percentile_ci <- function(x, stat_fun, B = 2000, alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(x)
  t0 <- stat_fun(x)
  tstar <- replicate(B, stat_fun(sample(x, n, replace = TRUE)))
  ci <- quantile(tstar, probs = c(alpha/2, 1 - alpha/2), names = FALSE, type = 6)
  list(t0 = t0, ci = ci, tstar = tstar)
}

### a/
sim_cauchy <- function(n = 30, location = 0, scale = 1,
                       B = 2000, R = 500, alpha = 0.05, seed = 1) {
  set.seed(seed)
  
  cover <- logical(R)
  avg_len <- numeric(R)
  example_x <- NULL
  example_tstar <- NULL
  
  for (r in 1:R) {
    x <- rcauchy(n, location = location, scale = scale)
    
    out <- boot_percentile_ci(x, mean, B = B, alpha = alpha)
    ci <- out$ci
    
    cover[r] <- (location >= ci[1] && location <= ci[2])
    avg_len[r] <- (ci[2] - ci[1])
    
    if (r == 1) {
      example_x <- x
      example_tstar <- out$tstar
    }
  }
  
  list(
    coverage = mean(cover),
    mean_len = mean(avg_len),
    cover_vec = cover,
    len_vec = avg_len,
    example_x = example_x,
    example_tstar = example_tstar
  )
}

resA <- sim_cauchy(n = 30, location = 0, scale = 1, B = 4000, R = 500, seed = 42)
cat("=== (3a) Cauchy mean ===\n")
cat("Empirical coverage of 95% percentile CI:", round(resA$coverage, 3), "\n")
cat("Average CI length:", round(resA$mean_len, 3), "\n\n")

par(mfrow = c(1, 3))
hist(resA$example_x, breaks = 30, main = "One Cauchy sample", xlab = "x")
hist(resA$example_tstar, breaks = 30, main = "Bootstrap means (one run)", xlab = "mean*")
hist(resA$len_vec, breaks = 30, main = "CI length over replications", xlab = "CI length")
par(mfrow = c(1, 1))

sim_cauchy_compare_mean_median <- function(n = 30, location = 0, scale = 1,
                                           B = 2000, R = 500, alpha = 0.05, seed = 1) {
  set.seed(seed)
  cov_mean <- cov_med <- logical(R)
  
  for (r in 1:R) {
    x <- rcauchy(n, location = location, scale = scale)
    
    ci_mean <- boot_percentile_ci(x, mean, B = B, alpha = alpha)$ci
    ci_med  <- boot_percentile_ci(x, median, B = B, alpha = alpha)$ci
    
    cov_mean[r] <- (location >= ci_mean[1] && location <= ci_mean[2])
    cov_med[r]  <- (location >= ci_med[1]  && location <= ci_med[2])
  }
  c(coverage_mean = mean(cov_mean), coverage_median = mean(cov_med))
}
cmpA <- sim_cauchy_compare_mean_median(n = 30, B = 3000, R = 400, seed = 99)
cat("Cauchy: coverage(mean) =", round(cmpA["coverage_mean"], 3),
    " | coverage(median) =", round(cmpA["coverage_median"], 3), "\n\n")


### b/
sim_unif_theta <- function(n = 20, theta = 1,
                           B = 4000, R = 500, alpha = 0.05, seed = 1) {
  set.seed(seed)
  
  cover_boot <- logical(R)
  len_boot <- numeric(R)
  
  cover_exact <- logical(R)
  len_exact <- numeric(R)
  
  example_x <- NULL
  example_boot_tstar <- NULL
  
  for (r in 1:R) {
    x <- runif(n, 0, theta)
    
    outB <- boot_percentile_ci(x, max, B = B, alpha = alpha)
    ciB <- outB$ci
    
    cover_boot[r] <- (theta >= ciB[1] && theta <= ciB[2])
    len_boot[r] <- (ciB[2] - ciB[1])

    M <- max(x)
    theta_L <- M / (1 - alpha/2)^(1/n)
    theta_U <- M / (alpha/2)^(1/n)
    
    cover_exact[r] <- (theta >= theta_L && theta <= theta_U)
    len_exact[r] <- (theta_U - theta_L)
    
    if (r == 1) {
      example_x <- x
      example_boot_tstar <- outB$tstar
    }
  }
  
  list(
    coverage_boot = mean(cover_boot),
    mean_len_boot = mean(len_boot),
    coverage_exact = mean(cover_exact),
    mean_len_exact = mean(len_exact),
    len_boot = len_boot,
    len_exact = len_exact,
    example_x = example_x,
    example_boot_tstar = example_boot_tstar
  )
}

resB <- sim_unif_theta(n = 20, theta = 1, B = 6000, R = 500, seed = 202)
cat("=== (3b) Uniform(0, theta) with theta=1, stat=max ===\n")
cat("Empirical coverage of 95% bootstrap percentile CI:", round(resB$coverage_boot, 3), "\n")
cat("Average CI length (bootstrap):", round(resB$mean_len_boot, 3), "\n")
cat("Empirical coverage of 95% exact CI:", round(resB$coverage_exact, 3), "\n")
cat("Average CI length (exact):", round(resB$mean_len_exact, 3), "\n\n")

par(mfrow = c(1, 3))
hist(resB$example_x, breaks = 20, main = "One U(0,θ) sample", xlab = "x")
hist(resB$example_boot_tstar, breaks = 30, main = "Bootstrap max* (one run)", xlab = "max*")
boxplot(resB$len_boot, resB$len_exact, names = c("Boot", "Exact"),
        main = "CI length comparison", ylab = "length")
par(mfrow = c(1, 1))

cat("NOTES:\n")
cat("- (3a) Cauchy: mean has no finite variance; bootstrap CI for mean is unstable, coverage can be poor.\n")
cat("- (3b) U(0,θ): theta is an endpoint; bootstrap from empirical data cannot exceed observed max -> downward bias.\n")
cat("- Exact CI uses the known distribution of max and typically achieves correct coverage.\n")