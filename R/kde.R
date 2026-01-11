source("R/kernel_function.R")

kde <- function(data, y, h, kernel="gaus") {
  # n <- length(data)
  if(kernel == "gaus") {
    u <- kernel_gaus((data - y) / h) / h
    kde_out <- mean(u)
  } else {
    kde_out  <- NA
  }
  return(kde_out)
}