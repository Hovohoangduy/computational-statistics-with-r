set.seed(1105)
n <- 800
x <- runif(n, -3, 3)
y <- sin(x) + rnorm(n, sd = 0.25)

X <- cbind(x, y)
plot(X, pch = 16, cex = 0.4, xlab = "x", ylab = "y")


kernel_gauss <- function(u) {
  exp(-0.5 * sum(u^2)) / (2 * pi)
}

kde <- function(x, X, h) {
  n <- nrow(X)
  d <- ncol(X)
  
  s <- 0
  for (i in 1:n) {
    u <- (x - X[i, ]) / h
    s <- s + kernel_gauss(u)
  }
  s / (n * h^d)
}

kde_grad <- function(x, X, h) {
  n <- nrow(X)
  d <- ncol(X)
  
  g <- rep(0, d)
  for (i in 1:n) {
    u <- (x - X[i, ]) / h
    K <- kernel_gauss(u)
    g <- g + K * (X[i, ] - x)
  }
  g / (n * h^(d + 2))
}

kde_hessian <- function(x, X, h) {
  n <- nrow(X)
  d <- ncol(X)
  
  H <- matrix(0, d, d)
  I <- diag(d)
  
  for (i in 1:n) {
    diff <- X[i, ] - x
    u <- diff / h
    K <- kernel_gauss(u)
    H <- H + K * (outer(diff, diff) - h^2 * I)
  }
  H / (n * h^(d + 4))
}

scms <- function(x0, X, h, max_iter = 50, eps = 1e-4) {
  x <- x0
  
  for (iter in 1:max_iter) {
    grad <- kde_grad(x, X, h)
    H <- kde_hessian(x, X, h)
    
    eig <- eigen(H)
    
    # eigenvector ứng với eigenvalue nhỏ nhất
    idx <- which.min(eig$values)
    v2 <- eig$vectors[, idx, drop = FALSE]
    
    # projector lên không gian trực giao ridge
    P <- v2 %*% t(v2)
    
    x_new <- x + as.vector(P %*% grad)
    
    if (sum((x_new - x)^2) < eps) break
    x <- x_new
  }
  x
}

h <- 0.4

init_idx <- sample(1:n, 200)
ridge_pts <- matrix(0, 200, 2)

for (i in 1:200) {
  ridge_pts[i, ] <- scms(X[init_idx[i], ], X, h)
}

plot(X, pch = 16, cex = 0.35, main = "Density Ridge Estimation", xlab = "x", ylab = "y")

points(ridge_pts, col = "red", pch = 16, cex = 0.6)

