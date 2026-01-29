source(file = "R/baitap3.R")
# bai tap 4
# a/
mars_data <- read.table(file = "data/mars.dat", header = TRUE)
head(mars_data)

x <- mars_data$radius
y <- mars_data$temperature

fit_center <- poly_reg(x = x, y = y, x_eval = mean(x), h = NULL, p = 1)
h_use <- fit_center$h

x_grid <- seq(min(x), max(x), length.out = 200)

m_hat <- sapply(x_grid, function(x0) {
  fit <- lp_reg(x, y, x_eval = x0, h = h_use, p = 1)
  if (is.null(fit)) NA else fit$beta[1]
})

plot(x, y, pch = 16, col = "blue", xlab = "Radius (km)", ylab = "Temperature")

lines(x_grid, m_hat, col = "red", lwd = 2)

## b/
marsbig_data <- read.table(file = "data/marsbig.dat", header = TRUE)
head(marsbig_data)
x_big <- marsbig_data$radius
y_big <- marsbig_data$temperature
orbit_big <- marsbig_data$orbit

fit_center_big <- poly_reg(x = x_big, y = y_big, x_eval = mean(x_big), h = NULL, p = 1)
h_big <- fit_center_big$h

x_grid_big <- seq(min(x_big), max(x_big), length.out = 200)

orbits <- sort(unique(orbit_big))
colors <- rainbow(length(orbits))

plot(x_big, y_big, pch = 16, col = rgb(0,0,0,0.2), xlab = "Radius (km)", ylab = "Temperature (K)")

for (k in seq_along(orbits)) {
  ok <- orbit_big == orbits[k]
  
  m_hat_k <- sapply(x_grid_big, function(x0) {
    fit <- lp_reg(
      x = x_big[ok],
      y = y_big[ok],
      x_eval = x0,
      h = h_big,
      p = 1
    )
    if (is.null(fit)) NA else fit$beta[1]
  })
  
  lines(x_grid_big, m_hat_k, col = colors[k], lwd = 2)
}
legend("topright", 
       legend = paste("Orbit", orbits), 
       col = colors, lwd = 2,
       bty = "n", cex = 0.7)

## c/
marsbig_data <- read.table(file = "data/marsbig.dat", header = TRUE)
head(marsbig_data)
x_pre <- marsbig_data$pressure
y_pre <- marsbig_data$temperature
orbit_pre <- marsbig_data$orbit

fit_center_pre <- poly_reg(x = x_pre, y = y_pre, x_eval = mean(x_pre), h = NULL, p = 1)
h_pre <- fit_center_pre$h

x_grid_pre <- seq(min(x_pre), max(x_pre), length.out = 200)

orbits <- sort(unique(orbit_pre))
colors <- rainbow(length(orbits))

plot(x_pre, y_pre, pch = 16, col = rgb(0,0,0,0.2), xlab = "Radius (km)", ylab = "Temperature (K)")

for (k in seq_along(orbits)) {
  ok <- orbit_pre == orbits[k]
  
  m_hat_k <- sapply(x_grid_pre, function(x0) {
    fit <- lp_reg(
      x = x_pre[ok],
      y = y_pre[ok],
      x_eval = x0,
      h = h_pre,
      p = 1
    )
    if (is.null(fit)) NA else fit$beta[1]
  })
  
  lines(x_grid_pre, m_hat_k, col = colors[k], lwd = 2)
}
legend("bottomright",
       legend = paste("Orbit", orbits), 
       col = colors, lwd = 2, bty = "n",
       cex = 0.7)