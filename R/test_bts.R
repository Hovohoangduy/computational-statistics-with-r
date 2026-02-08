source("R/bootstrap.R")

x <- rnorm(10)
my_boot(data = x, R = 5)

y <- data.frame(x1 = c(1:10), x2 = rnorm(10))
my_boot(data = y, R = 5)

## Gamma distribution

x <- rgamma(n = 20, shape = 3, rate = 0.5)
x_bts <- my_boot(data = x, R = 500)
x_bts
x_mean_bts <- apply(x_bts, MARGIN = 2, FUN = mean)
hist(x_mean_bts)
mean(x)
mean(x_mean_bts)

x_mean_bts_s <- sort(x_mean_bts)
x_mean_bts_s[c(floor(500*0.05/2), floor(500*(1 - 0.05/2)))]

quantile(x = x_mean_bts, probs = c(0.05/2, 1 - 0.05/2), type = 3)


# BCa
x_mean <- mean(x)
b_est <- qnorm(mean(x_mean_bts <= x_mean))
b_est

x_mean_jack <- numeric(length(x))
# for (i in 1:length(x)) {
#  x_mean_jack[i] <- mean[-i]
# }
x_mean_jack <- sapply(1:length(x), function(i) mean(x[-i]))

x_mean_jack2 <- mean(x_mean_jack)
u <- x_mean_jack2 - x_mean_jack
a_est <- ((1/6) * sum(u^3)) * (sum(u^2)^(3/2))
a_est
alpha <- 0.05
Z_alpha <- qnorm(alpha / 2)

beta1_est <- pnorm(b_est + 1/ (1/(b_est + Z_alpha) - alpha))
beta2_est <- pnorm(b_est + 1/ (1/(b_est - Z_alpha) - alpha))

x_mean_bts_s <- sort(x_mean_bts)
x_mean_bts_s[c(floor(500 * beta1_est), floor(500 * beta2_est))]


## claridge
claridge <- read.csv(file = "data/claridge.csv")

cor(claridge, method = "pearson")[1, 2]
cor(x = claridge$dnan, y = claridge$hand, method = "pearson")

plot(x = claridge$dnan, y = claridge$hand, pch = 16)

claridge_bts <- my_boot(data = claridge, R = 500, seed = 28)

my_cor <- function(data, type = "pearson"){
  cor(data, method = type)[1, 2]
}
my_cor(data = claridge)

out_cor_1 <- my_ci_boot(boot_sample = claridge_bts, data = claridge,
                        fun = my_cor)
hist(out_cor_1$est_boot, xlim = c(-0.5, 1))

my_ci_boot(boot_sample = claridge_bts, data = claridge, fun = my_cor,
           type = "kendall")

out_cor_2 <- my_ci_boot(boot_sample = claridge_bts, data = claridge, 
                        fun = my_cor, ci_type = "bca")
