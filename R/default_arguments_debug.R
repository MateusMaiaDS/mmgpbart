# # Running random arguments to debugg GP-BART
# rm(list=ls())
# # source("R/bart.R")
#
# library(purrr)
# n_train <- n_test <-  100
# seed <- 42
# set.seed(seed)
# sim_data_train <- mlbench::mlbench.friedman1(n = n_train,sd = 0.01) %>% as.data.frame %>% .[,c(1:5,11), drop = FALSE]
# sim_data_test <- mlbench::mlbench.friedman1(n = n_train,sd = 0.01) %>% as.data.frame %>% .[,c(1:5,11), drop = FALSE]
# sim_data_train <- matrix(seq(-pi,pi,length.out = n_train))
# sim_data_train <- cbind(sim_data_train,sin(sim_data_train))
# colnames(sim_data_train) <- c("x","y")
# sim_data_test <- sim_data_train
#
# # Setting the data
# x_train <- sim_data_train[,1, drop = FALSE]
# x_test <- sim_data_test[,1, drop = FALSE]
#
# y_train <- sim_data_train[,2]
# y_test <- sim_data_test[,2]
#
# # Setting default parameters
# numcut <- 100
# alpha <- 0.95
# beta <- 2
# tau <- tau_mu <- 100
# n_tree <- 20
# n_mcmc <- 2000
# n_burn <- 250
# n_min_size <- 5
# df <- 3
# sigquant <- 0.9
# numcut <- 100
# scale_boolean <- TRUE
# K_bart <- 2
#
#
# # Creating the x_cut object
# a_min <- NULL
# b_max <- NULL
#
# a_min <- min(y_train)
# b_max <- max(y_train)
#
# # Cut matrix
# xcut <- matrix(NA,ncol = ncol(x_train),nrow = numcut)
#
# # Getting possible x values
# for(j in 1:ncol(x_train)){
#         xs <- stats::quantile(x_train[ , j], type=7,
#                               probs=(0:(numcut+1))/(numcut+1))[-c(1, numcut+2)]
#
#         xcut[,j] <-xs
# }# Error of the matrix
#
# # Dummy for a new tree
# # tree <- new_tree(x_train = x_train,x_test = x_test)
# #
# # grow_tree <- grow(res_vec = y_train,tree = tree,x_train = x_train,
# #                   x_test = x_test,xcut = xcut,tau = tau,tau_mu = tau_mu,
# #                   alpha = alpha, beta = beta)
# #
# # grow_tree_two <- grow(res_vec = y_train,tree = grow_tree,x_train = x_train,
# #                   x_test = x_test,xcut = xcut,tau = tau,tau_mu = tau_mu,
# #                   alpha = alpha, beta = beta)
# #
# # # Setting the parameters
# # res_vec <- y_train
# # node_min_size <- 5
# #
# # #
# bart_mod <- bart(x_train = x_train,y_train = y_train,x_test = x_train,n_tree = 200,n_mcmc = 2000,n_burn = 250,n_min_size = 1,
#      tau = 1,alpha = 0.95,beta = 2,df = 3,sigquant = 0.9,numcut = 100,scale_boolean = TRUE,K_bart = 2)
# # #
# #
# plot(bart_mod$tau_post^(-1/2), type = "l")
#
# bart_mod$y_hat_post %>% colMeans() %>% unique()
# #
#
# plot(x_train,bart_mod$y_hat_post %>% colMeans())
# points(x_train,y_train,pch=20)
# dbart_mod <- dbarts::bart(x.train = x_train,y.train = y_train,x.test = x_train,nskip = 1,ndpost = 2000,ntree = 200,keeptrees = TRUE)
# # dbart_mod$fit$plotTree(chainNum = 1,treeNum = 1,sampleNum = 1000)
# lines(dbart_mod$sigma, col = "red")
# #
# # oopBART <- oopBART::r_bart(x_train = as.matrix(x_train),y = y_train,
# #                            x_test = as.matrix(x_train),num_cut = 100,n_tree = 1,df = 3,sigquant = 0.9,
# #                            n_mcmc = 2000,n_burn = 250,n_min_size = 1,tau = 1,mu = 1,alpha = 0.95,beta = 2,scale_boolean = TRUE,K_bart = 2,
# #                            tau_linero_bool = FALSE,tau_mu_linero_bool = FALSE)
# # oopBART$bart_obj$tau_post^(-1/2) %>% points(col = "red")
