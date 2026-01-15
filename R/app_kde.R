source("R/kde.R")
# laod data
whale <- read.table(file="data/whalemigration.dat")

head(whale)
range(whale)

kde(data=whale$V1, y=mean(whale$V1), h=10, kernel = "consine")

v <- seq(1, 10, by=0.1)
# seq(1, 10, length.out=101)
vv <- numeric(length(v))

min_whale <- min(whale)
max_whale <- max(whale)

y_val <- seq(min_whale, max_whale, length.out=201)
kde_whale <- numeric(length(v))

for (i in 1:length(y_val)) {
  kde_whale[i] <- kde(data = whale$V1, y = y_val[i], h = 10, kernel = "consine")
}

plot(x = y_val, y = kde_whale, type="l")
rug(x = whale$V1)

hist(x = whale$V1, probability = TRUE, ylim = c(0, 0.015), breaks = 18)
lines(x = y_val, y = kde_whale, col = "blue")

kde_vect <- Vectorize(FUN = kde, vectorize.args = "y")
kde_whale_2 <- kde_vect(data = whale$V1, y = y_val, h = 10)
plot(x = y_val, y = kde_whale, type="l")
lines(x = y_val, y = kde_whale_2, col = "blue")
rug(x = whale$V1)

bw.nrd(x = whale$V1)
bw.SJ(x = whale$V1, method = "ste")
bw.SJ(x = whale$V1, method = "dpi")
bw.ucv(x = whale$V1)
bw.bcv(x = whale$V1)
