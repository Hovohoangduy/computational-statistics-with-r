pkgs <- c("rpart", "pROC", "randomForest")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
library(rpart)
library(pROC)
library(randomForest)

df <- read.csv("data/breast_cancer.csv", header = TRUE, stringsAsFactors = FALSE)
head(df)

if ("X" %in% names(df)) df$X <- NULL
df$diagnosis <- factor(df$diagnosis, levels = c("B", "M"))  # B=0, M=1 (positive)
stopifnot(all(!is.na(df$diagnosis)))
if ("id" %in% names(df)) df$id <- NULL

acc <- function(y_true, prob_M, thr = 0.5) {
  pred <- ifelse(prob_M >= thr, "M", "B")
  mean(pred == as.character(y_true))
}

auc_score <- function(y_true, prob_M) {
  roc_obj <- pROC::roc(response = y_true, predictor = prob_M, levels = c("B","M"), direction = "<", quiet = TRUE)
  as.numeric(pROC::auc(roc_obj))
}

train_test_split <- function(df, train_frac = 0.7, seed = 1) {
  set.seed(seed)
  n <- nrow(df)
  idx_train <- sample.int(n, size = floor(train_frac * n), replace = FALSE)
  list(train = df[idx_train, , drop = FALSE],
       test  = df[-idx_train, , drop = FALSE])
}

### a/
fit_tree <- function(train) {
  rpart(diagnosis ~ ., data = train, method = "class",
        control = rpart.control(cp = 0.01, minsplit = 20, xval = 10))
}

predict_tree_probM <- function(model, newdata) {
  p <- predict(model, newdata = newdata, type = "prob")
  p[, "M"]
}

run_part_a <- function(df, seed_split = 1) {
  spl <- train_test_split(df, train_frac = 0.7, seed = seed_split)
  tr <- spl$train; te <- spl$test
  
  model <- fit_tree(tr)
  probM <- predict_tree_probM(model, te)
  
  list(
    seed = seed_split,
    accuracy = acc(te$diagnosis, probM),
    auc = auc_score(te$diagnosis, probM)
  )
}

seeds_demo <- c(1, 2, 3, 10, 20)
res_a_demo <- lapply(seeds_demo, function(s) run_part_a(df, s))
print(do.call(rbind, lapply(res_a_demo, as.data.frame)))

### b/
bagging_predict_probM <- function(train, test, B = 200, seed = 1) {
  set.seed(seed)
  n <- nrow(train)
  prob_mat <- matrix(NA_real_, nrow = nrow(test), ncol = B)
  
  for (b in 1:B) {
    idx_boot <- sample.int(n, size = n, replace = TRUE)
    boot_train <- train[idx_boot, , drop = FALSE]
    
    m <- fit_tree(boot_train)
    prob_mat[, b] <- predict_tree_probM(m, test)
  }
  
  rowMeans(prob_mat)
}

run_part_b <- function(df, seed_split = 1, B = 200) {
  spl <- train_test_split(df, train_frac = 0.7, seed = seed_split)
  tr <- spl$train; te <- spl$test
  
  probM_bag <- bagging_predict_probM(tr, te, B = B, seed = 1000 + seed_split)
  
  list(
    seed = seed_split,
    accuracy = acc(te$diagnosis, probM_bag),
    auc = auc_score(te$diagnosis, probM_bag)
  )
}

res_b_one <- run_part_b(df, seed_split = 1, B = 200)
print(res_b_one)

### c/
run_part_c <- function(df, R = 50, B = 200, seed0 = 123) {
  set.seed(seed0)
  seeds <- sample.int(1e7, R)
  
  acc_tree <- numeric(R)
  acc_bag  <- numeric(R)
  
  auc_tree <- numeric(R)
  auc_bag  <- numeric(R)
  
  for (i in 1:R) {
    s <- seeds[i]
    spl <- train_test_split(df, train_frac = 0.7, seed = s)
    tr <- spl$train; te <- spl$test
    
    m_tree <- fit_tree(tr)
    p_tree <- predict_tree_probM(m_tree, te)
    acc_tree[i] <- acc(te$diagnosis, p_tree)
    auc_tree[i] <- auc_score(te$diagnosis, p_tree)
    
    p_bag <- bagging_predict_probM(tr, te, B = B, seed = 1000 + s)
    acc_bag[i] <- acc(te$diagnosis, p_bag)
    auc_bag[i] <- auc_score(te$diagnosis, p_bag)
  }
  
  mean_tree <- mean(acc_tree); var_tree <- var(acc_tree)
  mean_bag  <- mean(acc_bag);  var_bag  <- var(acc_bag)
  
  var_reduction_pct <- (var_tree - var_bag) / var_tree * 100
  
  list(
    acc_tree = acc_tree, acc_bag = acc_bag,
    auc_tree = auc_tree, auc_bag = auc_bag,
    summary = data.frame(
      model = c("Tree", "Bagging"),
      mean_accuracy = c(mean_tree, mean_bag),
      var_accuracy  = c(var_tree, var_bag),
      mean_auc      = c(mean(auc_tree), mean(auc_bag)),
      var_auc       = c(var(auc_tree), var(auc_bag))
    ),
    var_reduction_pct = var_reduction_pct
  )
}

res_c <- run_part_c(df, R = 50, B = 200, seed0 = 2026)
print(res_c$summary)
cat("\n% giam phuong sai Accuracy (Bagging vs Tree) =", round(res_c$var_reduction_pct, 2), "%\n")

boxplot(res_c$acc_tree, res_c$acc_bag,
        names = c("Tree", "Bagging"),
        main = "Accuracy over 50 splits", ylab = "Accuracy")

### d/
run_part_d <- function(df, R = 50, ntree = 200, seed0 = 2026) {
  set.seed(seed0)
  seeds <- sample.int(1e7, R)
  
  acc_rf <- numeric(R)
  auc_rf <- numeric(R)
  
  for (i in 1:R) {
    s <- seeds[i]
    spl <- train_test_split(df, train_frac = 0.7, seed = s)
    tr <- spl$train; te <- spl$test
    
    set.seed(5000 + s)
    rf <- randomForest(diagnosis ~ ., data = tr, ntree = ntree)  # mặc định mtry=sqrt(p)
    
    probM <- predict(rf, newdata = te, type = "prob")[, "M"]
    acc_rf[i] <- acc(te$diagnosis, probM)
    auc_rf[i] <- auc_score(te$diagnosis, probM)
  }
  
  data.frame(
    model = "RandomForest",
    mean_accuracy = mean(acc_rf),
    var_accuracy  = var(acc_rf),
    mean_auc      = mean(auc_rf),
    var_auc       = var(auc_rf)
  )
}

res_d <- run_part_d(df, R = 50, ntree = 200, seed0 = 2026)
print(res_d)