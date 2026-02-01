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
