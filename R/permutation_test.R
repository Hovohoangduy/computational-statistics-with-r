source(file = "R/permutation.R")
grade <- read.table(file = "data/grades.txt", header = TRUE)
head(grade)

grade$avg <- (grade$e1 + grade$e2)/2
head(grade)

boxplot(avg ~ mf, data = grade)

out_perm <- perm_test(x = grade$avg, g = grade$mf, label = "F", fun = mean)
hist(out_perm$t_perm)

