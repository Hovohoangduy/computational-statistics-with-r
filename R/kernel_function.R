# include all kernel function
kernel_gaus <- function(x) {
  kern <- exp(-x^2/2)/sqrt(2*pi)
  return(kern)
}