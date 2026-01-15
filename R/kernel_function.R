# include all kernel function
kernel_gaus <- function(x) {
  kern <- exp(-x^2/2)/sqrt(2*pi)
  return(kern)
}

kernel_epan <- function(x) {
  kern <- 0.75 * (1 - x^2) * (abs(x) <= 1)
  return(kern)
}

kernel_tri <- function(x) {
  kern <- (1 - abs(x)) * (abs(x) <= 1)
  return(kern)
}

kernel_consine <- function(x) {
  kern <- (pi / 4) * cos((pi / 2) * x) * (abs(x) <= 1)
  return(kern)
}