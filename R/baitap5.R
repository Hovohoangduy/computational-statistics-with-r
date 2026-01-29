source(file = "R/local_linear_reg.R")
# bai tap 5
# a/
set.seed(1105)
n <- 100

X <- runif(n, 0, pi)

# ham m(x)
m <- function(x) {
  result <- 1 + sin(x^2)/x^2
  return(result)
}

sd_eps <- abs(m(X)) / 8
eps <- rnorm(n, mean = 0, sd = sd_eps)
Y <- m(X) + eps
Y

x_grid <- seq(0.01, pi, length.out = 300)
plot(X, Y, pch = 16, col = "blue",
     xlab = "X", ylab = "Y")
lines(x_grid, m(x_grid), col = "red", lwd = 2)


# b/
h_grid <- seq(0.05, 1, length.out = 100)

h_CV_gauss <- h_grid[which.min(CV_loclin(X, Y, h_grid, "gauss"))]
h_CV_epan  <- h_grid[which.min(CV_loclin(X, Y, h_grid, "epan"))]

h_GCV_gauss <- h_grid[which.min(GCV_loclin(X, Y, h_grid, "gauss"))]
h_GCV_epan  <- h_grid[which.min(GCV_loclin(X, Y, h_grid, "epan"))]

bandwidth_table <- data.frame(
  Method = c("RT", "CV", "GCV"),
  Gaussian = c(
    h_ll_RT(X, Y, "gauss"),
    h_CV_gauss,
    h_GCV_gauss
  ),
  Epanechnikov = c(
    h_ll_RT(X, Y, "epan"),
    h_CV_epan,
    h_GCV_epan
  )
)

print(bandwidth_table)

# c/
h_RT_gauss  <- h_ll_RT(X, Y, kernel = "gauss")
h_RT_epan  <- h_ll_RT(X, Y, kernel = "epan")

m_RT_gauss  <- loclin_reg(X, Y, x_grid, h_RT_g,  kernel = "gauss")
m_CV_gauss  <- loclin_reg(X, Y, x_grid, h_CV_gauss,  kernel = "gauss")
m_GCV_gauss <- loclin_reg(X, Y, x_grid, h_GCV_gauss, kernel = "gauss")

m_RT_epan  <- loclin_reg(X, Y, x_grid, h_RT_epan,  kernel = "epan")
m_CV_epan  <- loclin_reg(X, Y, x_grid, h_CV_epan,  kernel = "epan")
m_GCV_epan <- loclin_reg(X, Y, x_grid, h_GCV_epan, kernel = "epan")

par(mar = c(5, 4, 4, 8))

plot(X, Y, pch = 16, col = rgb(0, 0, 0, 0.3), xlab = "X", ylab = "Y")

# Hàm m(x) thật
lines(x_grid, m(x_grid), col = "black", lwd = 2)

# RT
lines(x_grid, m_RT_gauss, col = "blue",  lwd = 2, lty = 1)
lines(x_grid, m_RT_epan, col = "blue",  lwd = 2, lty = 2)

# CV
lines(x_grid, m_CV_gauss, col = "red",   lwd = 2, lty = 1)
lines(x_grid, m_CV_epan, col = "red",   lwd = 2, lty = 2)

# GCV
lines(x_grid, m_GCV_gauss, col = "green", lwd = 2, lty = 1)
lines(x_grid, m_GCV_epan, col = "green", lwd = 2, lty = 2)

legend(
  "bottomleft",
  legend = c(
    "True m(x)",
    "RT – Gaussian", "RT – Epanechnikov",
    "CV – Gaussian", "CV – Epanechnikov",
    "GCV – Gaussian", "GCV – Epanechnikov"
  ),
  col = c(
    "black",
    "blue", "blue",
    "red", "red",
    "green", "green"
  ),
  lwd = 2,
  bty = "n",
  cex = 0.7
)


# d/
B <- 100
n <- 100
h_grid <- seq(0.05, 1, length.out = 100)
h_mat <- matrix(NA, nrow = B, ncol = 6)
colnames(h_mat) <- c(
  "RT_Gauss", "CV_Gauss", "GCV_Gauss",
  "RT_Epan",  "CV_Epan",  "GCV_Epan"
)

for (b in 1:B) {
  X <- runif(n, 0, pi)
  m <- function(x) 1 + sin(x^2)/x^2
  
  sd_eps <- abs(m(X)) / 8
  eps <- rnorm(n, 0, sd_eps)
  Y <- m(X) + eps

  h_mat[b, "RT_Gauss"]  <- h_ll_RT(X, Y, "gauss")
  h_mat[b, "CV_Gauss"]  <- h_grid[which.min(CV_loclin(X, Y, h_grid, "gauss"))]
  h_mat[b, "GCV_Gauss"] <- h_grid[which.min(GCV_loclin(X, Y, h_grid, "gauss"))]

  h_mat[b, "RT_Epan"]  <- h_ll_RT(X, Y, "epan")
  h_mat[b, "CV_Epan"]  <- h_grid[which.min(CV_loclin(X, Y, h_grid, "epan"))]
  h_mat[b, "GCV_Epan"] <- h_grid[which.min(GCV_loclin(X, Y, h_grid, "epan"))]
}

par(mar = c(5, 4, 4, 2))

boxplot(
  h_mat,
  col = c("lightblue", "pink", "lightgreen",
          "lightblue", "pink", "lightgreen"),
  ylab = "Optimal bandwidth"
)
