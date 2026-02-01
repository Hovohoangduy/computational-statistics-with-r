source(file = "R/bootstrap.R")

x <- rnorm(10)
my_boot(data = x, R = 5)

replicate(n = 5, expr = sample(1:5, size = 10, replace = TRUE))

y <- data.frame(x1 = c(1:10), x2 = rnorm(10))
my_boot(data = y, R = 5)

x <- rgamma(n = 20, shape = 3, rate = 0.5)
x_bts <- my_boot(data = x, R = 5)
x_bts

x_mean_bts <- apply(x_bts, MARGIN = 2, FUN = mean)
hist(x_mean_bts)
mean(x)
mean(x_mean_bts)

x_mean_bts_s <- sort(x_mean_bts)
x_mean_bts_s[c(floor(500*0.05/2), floor(500*(1 - 0.05/2)))]

quantile(x = x_mean_bts, probs = c(0.05/2, 1 - 0.05/2), type = 3)
