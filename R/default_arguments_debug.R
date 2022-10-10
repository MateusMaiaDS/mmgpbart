# Running random arguments to debugg GP-BART
library(purrr)
n_train <- n_test <-  100
seed <- 42
set.seed(seed)
sim_data_train <- mlbench::mlbench.friedman1(n = n_train,sd = 0.01) %>% as.data.frame %>% .[,c(1:5,11), drop = FALSE]
sim_data_test <- mlbench::mlbench.friedman1(n = n_train,sd = 0.01) %>% as.data.frame %>% .[,c(1:5,11), drop = FALSE]


# Setting the data
x_train <- sim_data_train[,1:5, drop = FALSE]
x_test <- sim_data_test[,1:5, drop = FALSE]

y_train <- sim_data_train$y
y_test <- sim_data_test$y

# Setting default parameters
numcut <- 100

# Creating the x_cut object
a_min <- NULL
b_max <- NULL

a_min <- min(y_train)
b_max <- max(y_train)

# Cut matrix
xcut <- matrix(NA,ncol = ncol(x_train),nrow = numcut)

# Getting possible x values
for(j in 1:ncol(x_train)){
        xs <- stats::quantile(x_train[ , j], type=7,
                              probs=(0:(numcut+1))/(numcut+1))[-c(1, numcut+2)]

        xcut[,j] <-xs
}# Error of the matrix

# Dummy for a new tree
tree <- new_tree(x_train = x_train,x_test = x_test)


# Setting the parameters

