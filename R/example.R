x <- 5
x

# vector
x1 <- c(-4, 2, 4.5, 0.1)
x1
x2 <- c("student", "worker", "teacher", "0.01")
x2

# matrix
A <- matrix(data = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 9.6, -0.5, 12),
            nrow = 3, ncol = 4, byrow = TRUE)
A
B <- matrix(data = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 9.6, -0.5, 12),
            nrow = 3, ncol = 4, byrow = FALSE)
B

A <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 9.6, -0.5, 12))
A

## read csv file
data_students <- read.csv("data/students.csv", na = c("", "NA", "N/A"))
data_students

library(janitor)
data_students <- read.csv("data/students.csv", na = c("", "NA", "N/A")) |> clean_names()
data_students

1 + 2.5

x <- c(1, 2, 3, 4, 5)
y <- c(1, 2, 3, 4, 5)
x + y
x/y
x^2

A <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 9.6, 1.2, 1.9))
t(A)
A %*%t(A)

A * A
# lay nghich dao
solve(A %*% t(A))

x <- c(6, 7, 8, 9)
x[3]
x[c(1, 4, 2)]

A <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 9.6, 1.2, 1.9))
A[1, 1]
A[2, ]

# phan tu tren duong cheo
x <- c(1.2, 5.6, 3.2, 0.6, -1)
x
diag(x)
sum(x)

A <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 9.6,-0.5, 12))
sum(diag(A))

x <- c(1.2, 5.6, 3.2, 0.6, -1)
x > 0
x[x > 0]
x[x > 0 & x <= 2]

## function
mean_x <- function(x) {
    n <- length(x)
    res <- sum(x)/n
    return(res)
}

x <- rnorm(30, mean = 1, sd = 2)
mean_x(x)

## plot visualize
g <- function(x) {
    log(x) / (1 + x)
}
x <- seq(from = 1, to = 7, by = 0.1)
g_y <- g(x)
plot(x = x, y = g_y, type='l', lty = 1, xlab = "x", ylab="g(x)")