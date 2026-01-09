# load library ...

func_sin2 <- function(x, y) {
  z <- x^2 + y^2
  res <- sin(z)
  return (res)
}

func_sin2(x = 2,  y = 15)

func_sin <- function(x) {
  sin(x)
}

x <- seq(0,1, by=0.1)
y <- func_sin(x)


plot(x, y, type = "l")