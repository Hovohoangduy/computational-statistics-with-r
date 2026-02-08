perm_test <- function(x, g, label, fun, R = 1000, ...){
  n <- length(x)
  nA <- sum(g == label)
  # data_perm <- list()
  t_perm <- numeric(R)
  for(i in 1:R){
    x_perm <- sample(x, size = n, replace = FALSE)
    xA_perm <- x_perm[1:nA]
    xB_perm <- x_perm[(nA + 1):n]
    ## test statistics
    xA_mean <- fun(xA_perm, ...)
    xB_mean <- fun(xB_perm, ...)
    t_perm[i] <- xA_mean - xB_mean
    # data_perm[[i]] <- data.frame(x = x_perm, 
    #                              gr = rep(c("A", "B"), times = c(nA, n - nA)))
  }
  t_obs <- mean(x[g == label]) - mean(x[g != label])
  p_value <- mean(abs(t_perm) >= abs(t_obs))
  return(list(t_perm = t_perm, t_obs = t_obs, p_value = p_value))
}
