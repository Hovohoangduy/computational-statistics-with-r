df <- read.table("data/cancersurvival.dat", header = TRUE)
head(df)

df$logtime <- log(df$survivaltime)

boot_mean <- function(x, B = 2000, seed = 1) {
  set.seed(seed)
  n <- length(x)
  bmeans <- replicate(B, mean(sample(x, n, replace = TRUE)))
  list(t0 = mean(x), t = bmeans)
}

### a/
bca_ci <- function(x, B = 2000, alpha = 0.05, seed = 1) {
  set.seed(seed)
  n <- length(x)
  t0 <- mean(x)
  
  t <- replicate(B, mean(sample(x, n, replace = TRUE)))
  
  z0 <- qnorm(mean(t < t0))
  
  tj <- sapply(1:n, function(i) mean(x[-i]))
  tbar <- mean(tj)
  num <- sum((tbar - tj)^3)
  den <- 6 * (sum((tbar - tj)^2))^(3/2)
  a <- ifelse(den == 0, 0, num / den)
  
  zalpha1 <- qnorm(alpha/2)
  zalpha2 <- qnorm(1 - alpha/2)
  
  adj1 <- pnorm(z0 + (z0 + zalpha1) / (1 - a * (z0 + zalpha1)))
  adj2 <- pnorm(z0 + (z0 + zalpha2) / (1 - a * (z0 + zalpha2)))
  
  ci <- quantile(t, probs = c(adj1, adj2), names = FALSE, type = 6)
  list(ci_log = ci, ci_orig = exp(ci), t0_log = t0, t0_orig = exp(t0), boot_t = t)
}

studentized_ci <- function(x, B1 = 2000, B2 = 200, alpha = 0.05, seed = 1) {
  set.seed(seed)
  n <- length(x)
  t0 <- mean(x)
  se0 <- sd(x) / sqrt(n)
  
  tstar <- numeric(B1)
  sest  <- numeric(B1)
  
  for (b in 1:B1) {
    xb <- sample(x, n, replace = TRUE)
    tstar[b] <- mean(xb)
    
    inner_means <- replicate(B2, mean(sample(xb, n, replace = TRUE)))
    sest[b] <- sd(inner_means)
    if (!is.finite(sest[b]) || sest[b] == 0) sest[b] <- sd(xb) / sqrt(n)
  }
  
  tstat <- (tstar - t0) / sest
  qlo <- quantile(tstat, probs = 1 - alpha/2, names = FALSE, type = 6)
  qhi <- quantile(tstat, probs = alpha/2,   names = FALSE, type = 6)
  
  ci_log <- c(t0 - qlo * se0, t0 - qhi * se0)
  list(ci_log = ci_log, ci_orig = exp(ci_log), t0_log = t0, t0_orig = exp(t0),
       tstat = tstat, tstar = tstar)
}


groups <- sort(unique(df$disease))

res <- lapply(groups, function(g) {
  x <- df$logtime[df$disease == g]
  list(
    group = g,
    n = length(x),
    bca  = bca_ci(x, B = 3000, seed = 10 + g),
    stud = studentized_ci(x, B1 = 2000, B2 = 300, seed = 100 + g)
  )
})

# Print summary
for (r in res) {
  cat("\n===== Disease =", r$group, " | n =", r$n, "=====\n")
  cat("Point estimate (mean log):", r$bca$t0_log, "\n")
  cat("Point estimate (exp(mean log)):", r$bca$t0_orig, "\n\n")
  
  cat("BCa 95% CI on log-scale:", paste(round(r$bca$ci_log, 4), collapse = " , "), "\n")
  cat("BCa 95% CI on original scale (exp):", paste(round(r$bca$ci_orig, 2), collapse = " , "), "\n\n")
  
  cat("Studentized 95% CI on log-scale:", paste(round(r$stud$ci_log, 4), collapse = " , "), "\n")
  cat("Studentized 95% CI on original scale (exp):", paste(round(r$stud$ci_orig, 2), collapse = " , "), "\n")
}

par(mfrow = c(length(res), 2))
for (r in res) {
  hist(r$bca$boot_t, breaks = 25,
       main = paste0("BCa bootstrap mean(log) - disease=", r$group),
       xlab = "mean(log(survivaltime))")
  abline(v = r$bca$t0_log, lwd = 2)
  
  hist(exp(r$bca$boot_t), breaks = 25,
       main = paste0("BCa bootstrap exp(mean(log)) - disease=", r$group),
       xlab = "exp(mean(log)) (geometric mean, original scale)")
  abline(v = r$bca$t0_orig, lwd = 2)
}
par(mfrow = c(1, 1))

### b/
perm_test <- function(df, B = 20000, seed = 999, two_sided = TRUE) {
  set.seed(seed)
  
  df <- df[!is.na(df$logtime) & !is.na(df$disease), ]
  df$disease <- as.integer(as.character(df$disease))
  
  gr <- sort(unique(df$disease))
  if (length(gr) != 2) {
    stop("Permutation test can dung CHINH XAC 2 nhom disease (vd 0/1). Hien tai: ",
         paste(gr, collapse = ", "))
  }
  
  g0 <- gr[1]; g1 <- gr[2]
  n0 <- sum(df$disease == g0)
  n1 <- sum(df$disease == g1)
  if (n0 == 0 || n1 == 0) stop("Mot nhom dang rong. Kiem tra du lieu disease.")
  
  obs <- mean(df$logtime[df$disease == g1]) - mean(df$logtime[df$disease == g0])
  
  y <- df$logtime
  perm_stats <- replicate(B, {
    yp <- sample(y, length(y), replace = FALSE)
    mean(yp[1:n1]) - mean(yp[(n1 + 1):(n1 + n0)])
  })
  
  if (two_sided) {
    pval <- mean(abs(perm_stats) >= abs(obs))
  } else {
    pval <- mean(perm_stats >= obs)
  }
  
  list(obs = obs, perm = perm_stats, pval = pval, groups = c(g0, g1), n = c(n0, n1))
}

pt <- perm_test(df, B = 20000, seed = 999)
cat("Groups:", pt$groups[1], "vs", pt$groups[2], " | n0,n1 =", pt$n, "\n")
cat("Permutation p-value =", pt$pval, "\n")

hist(pt$perm, breaks = 40,
     main = "Permutation distribution of diff in means (log-scale)",
     xlab = "mean(log)_group2 - mean(log)_group1")
abline(v = pt$obs, lwd = 2)

### c/
df$survivaltime <- as.numeric(df$survivaltime)
df$disease <- as.integer(as.character(df$disease))
df <- df[is.finite(df$survivaltime) & df$survivaltime > 0 & !is.na(df$disease), ]
df$logtime <- log(df$survivaltime)

breast <- df[df$disease != 1, , drop = FALSE]
if (nrow(breast) < 2) stop("Nhom ung thu vu khong du du lieu (disease != 1).")

boot_percentile_ci <- function(x, stat_fun, B = 20000, alpha = 0.05, seed = 999) {
  set.seed(seed)
  n <- length(x)
  t_star <- replicate(B, stat_fun(sample(x, n, replace = TRUE)))
  ci <- quantile(t_star, probs = c(alpha/2, 1 - alpha/2), names = FALSE, type = 6)
  list(ci = ci, t0 = stat_fun(x), t_star = t_star)
}

B <- 20000

# (c1) Percentile on log scale for mean(logtime), then exp endpoints
res_log <- boot_percentile_ci(
  x = breast$logtime,
  stat_fun = mean,
  B = B, seed = 999
)
ci_log_then_exp <- exp(res_log$ci)
point_log_then_exp <- exp(res_log$t0)

# (c2) Percentile on original scale for mean(survivaltime)
res_orig <- boot_percentile_ci(
  x = breast$survivaltime,
  stat_fun = mean,
  B = B, seed = 999
)
ci_orig <- res_orig$ci
point_orig <- res_orig$t0

cat("=== BAI 2(c) - Ung thu vu (disease != 1) ===\n")
cat("n =", nrow(breast), "\n\n")

cat("[C1] Percentile tren log, sau do exp() 2 can\n")
cat("Point estimate exp(mean(log)) (geometric mean) =", round(point_log_then_exp, 3), "\n")
cat("95% CI =", paste(round(ci_log_then_exp, 3), collapse = " , "), "\n\n")

cat("[C2] Percentile tren thang goc (mean survivaltime)\n")
cat("Point estimate mean(survivaltime) (arithmetic mean) =", round(point_orig, 3), "\n")
cat("95% CI =", paste(round(ci_orig, 3), collapse = " , "), "\n")

par(mfrow = c(1,2))
hist(exp(res_log$t_star), breaks = 30,
     main = "Bootstrap dist: exp(mean(log))",
     xlab = "geometric mean (original scale)")
abline(v = point_log_then_exp, lwd = 2)

hist(res_orig$t_star, breaks = 30,
     main = "Bootstrap dist: mean(original)",
     xlab = "arithmetic mean (days)")
abline(v = point_orig, lwd = 2)
par(mfrow = c(1,1))