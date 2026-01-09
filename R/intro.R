A <- matrix(rnorm(10000 * 2000), 2000)
B <- matrix(rnorm(10000 * 2000), 10000)

print(system.time(sum(diag(A %*%B))))

system.time(sum(A * t(B)))

logf <- function(x) {
  f <- 0.3*dnorm(x, -1, sqrt(0.1)) + 0.7*dnorm(x, 1, 0.1)
  return(log(f))
}
print(logf(c(--15, 15)))

logf <- function(x) {
  log1 <- log(0.3) + dnorm(x, mean=-1, sd=sqrt(0.1), log=TRUE)
  log2 <- log(0.7) + dnorm(x, mean=1, sd=0.1, log=TRUE)
  m <- pmax(log1, log2)
  m + log(exp(log1 - m) + exp(log2 - m))
  }
logf(c(-15, 15))