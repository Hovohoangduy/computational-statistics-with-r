set.seed(1105)
sample(x = 1:10, size = 10, replace = TRUE)

my_boot <- function(data, R = 500, seed = NULL){
  if(is.null(seed)) seed <- 1105
  set.seed(seed)
  if(is.null(dim(data))){
    n <- length(data)
    id_bts <- replicate(n = R, 
                        expr = sample(1:n, size = n, replace = TRUE))
    data_bts <- matrix(0, nrow = n, ncol = R)
    for(i in 1:R){
      data_bts[, i] <- data[id_bts[, i]]
    }
  } else{
    n <- nrow(data) # dim(data)[1] # number's rows 
    id_bts <- replicate(n = R, 
                        expr = sample(1:n, size = n, replace = TRUE))
    data_bts <- list() # list-object
    for(i in 1:R){
      data_bts[[i]] <- data[id_bts[, i], ]
    }
  }
  return(data_bts)
}


my_ci_boot <- function(boot_sample, data, fun, alpha = 0.05, 
                       ci_type = c("quantile", "bca"), ...){
  ci_type <- match.arg(ci_type) ## defaut is quantile
  est_boot <- sapply(boot_sample, function(x) fun(x, ...))
  est_boot_st <- sort(est_boot)
  R <- length(boot_sample)
  alpha2 <- alpha/2
  if (ci_type == "quantile"){
    ci_boot <- est_boot_st[c(floor(R * alpha2), floor(R * (1 - alpha2)))]
  }
  if (ci_type == "bca"){
    est_sample <- fun(data, ...)
    b_est <- qnorm(mean(est_boot <= est_sample))
    if(is.null(dim(data))){
      est_jack <- sapply(1:length(data), function(i){
        fun(data[-i], ...)
      })
    } else{
      est_jack <- sapply(1:nrow(data), function(i){
        fun(data[-i, ], ...)
      })
    }
    est_jack2 <- mean(est_jack)
    u <- est_jack2 - est_jack
    a_est <- ((1/6) * sum(u^3)) * (sum(u^2))^(-3/2)
    z_alp <- qnorm(alpha2)
    beta1_est <- pnorm(b_est + 1/(1/(b_est + z_alp) - a_est))
    beta2_est <- pnorm(b_est + 1/(1/(b_est - z_alp) - a_est))
    ci_boot <- est_boot_st[c(floor(R * beta1_est), floor(R * beta2_est))]
    print(c(a_est, b_est))
  }
  return(list(est_boot = est_boot, ci_type = ci_type, ci_boot = ci_boot))
}
