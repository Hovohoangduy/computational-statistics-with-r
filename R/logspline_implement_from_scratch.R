library(splines)

# mixture gaussian data
set.seed(123)

x <- c(
  rnorm(300, -2, 0.7),
  rnorm(700, -2, 1.0)
)

n <- length(x)

# domain selection and knots
xmin <- min(x) - 1
xmax <- max(x) + 1

knots <- quantile(x, probs = c(0.2, 0.4, 0.6, 0.8))
print(knots)

# create B-splines
bs_basic <- function(x, knots, xmin, xmax) {
  bs(x, knots = knots, degree = 3, Boundary.knots = c(xmin, xmax), intercept = TRUE)
}

# defined log-density spline
g_fun <- function(x, beta) {
  B <- bs_basic(x, knots, xmin, xmax)
  drop(B %*% beta)
}

# normalize - calc partition function
Z_fun <- function(beta) {
  f <- function(t) {
    # Ensure numerical stability: prevent overflow in exp
    val <- g_fun(t, beta)
    # If values are too large, exp will be Inf
    if (any(val > 700)) {
      # handling large values if needed, or let integrate fail/handle it
    }
    exp(val)
  }

  res <- tryCatch(
    {
      integrate(f, xmin, xmax)$value
    },
    error = function(e) {
      NA
    }
  )

  if (is.na(res) || res <= 0) {
    return(NA)
  }
  return(res)
}

# log-likelihood logspline
loglik <- function(beta) {
  gx <- g_fun(x, beta)
  Z <- Z_fun(beta)

  if (is.na(Z)) {
    return(-1e10)
  } # Return a very small value (large penalty)

  sum(gx) - n * log(Z)
}

# log-likelihood gradient
grad_loglik <- function(beta) {
  Bx <- bs_basic(x, knots, xmin, xmax)
  term1 <- colSums(Bx)

  Z <- Z_fun(beta)
  if (is.na(Z)) {
    return(rep(0, length(beta)))
  } # Should not happen if loglik checks Z

  EB <- function(j) {
    f <- function(t) {
      B <- bs_basic(t, knots, xmin, xmax)
      B[, j] * exp(g_fun(t, beta))
    }
    # Remove Vectorize wrapper for speed, assuming f is vectorized
    integrate(f, xmin, xmax)$value / Z
  }

  term2 <- sapply(1:ncol(Bx), EB)
  term1 - n * term2
}

# MLE optimizer
K <- ncol(bs_basic(x, knots, xmin, xmax))
beta0 <- rep(0, K)

# BFGS usually requires finite values. We minimize negative loglik.
fit <- optim(
  par = beta0, fn = function(b) -loglik(b),
  gr = function(b) -grad_loglik(b), method = "BFGS",
  control = list(maxit = 500)
)

beta_hat <- fit$par

# density logspline
density_logspline <- function(x, beta) {
  Z <- Z_fun(beta)
  exp(g_fun(x, beta)) / Z
}

# results plot
grid <- seq(xmin, xmax, length.out = 400)
dens_hat <- density_logspline(grid, beta_hat)

# Determine ylim to ensure density is visible
max_dens <- max(dens_hat, na.rm = TRUE)
hist_info <- hist(x, plot = FALSE, breaks = 40)
max_hist <- max(hist_info$density)

ylim_val <- c(0, max(max_dens, max_hist) * 1.1)

hist(x,
  prob = TRUE, breaks = 40,
  col = "gray", border = "white",
  main = "Logspline density",
  ylim = ylim_val
)

lines(grid, dens_hat, col = "red", lwd = 2)
box()
