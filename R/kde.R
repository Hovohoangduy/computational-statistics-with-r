source("R/kernel_function.R")

kde <- function(data, y, h, kernel = "gaus") {
  # n <- length(data)
  if (kernel == "gaus") {
    u <- kernel_gaus((data - y) / h) / h
    kde_out <- mean(u, na.rm = TRUE)
  } else if (kernel == "epan") {
    u <- kernel_epan((data - y) / h) / h
    kde_out <- mean(u, na.rm = TRUE)
  } else if (kernel == "tri") {
    u <- kernel_tri((data - y) / h) / h
    kde_out <- mean(u, na.rm = TRUE)
  } else if (kernel == "consine") {
    u <- kernel_consine((data - y) / h) / h
    kde_out <- mean(u, na.rm = TRUE)
  }
  else {
    kde_out <- NA
  }
  return(kde_out)
}