integ_trap <- function(x, f, ...) {
  # x: vector points
  # f: class function
  y <- f(x, ...)
  n <- length(x)
  res <- 0

  for(i in 1:(n-1)) {
    res <- res + (x[i + 1] - x[i]) * (y[i] + y[i+1])
  }
  res <- res*0.5
  return(res)
}
