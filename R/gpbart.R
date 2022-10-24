## GP-Bart
#' @useDynLib mmgpbart
#' @importFrom Rcpp sourceCpp
#'
# Building the GP-BART function
gp_bart <- function(x_train,
                 y_train,
                 x_test,
                 n_tree,
                 n_mcmc,
                 n_burn,
                 n_min_size,
                 tau,
                 alpha, beta,
                 df, sigquant,
                 numcut,
                 scale_boolean = TRUE,
                 K_bart = 2,
                 bart_boolean = FALSE){

        # Saving a_min and b_max
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
        if(is.null(colnames(x_train)) || is.null(colnames(x_test)) ) {
                stop("Insert a valid NAMED matrix")
        }

        if(!is.vector(y_train)) {
                stop("Insert a valid y vector")
        }

        # Scale values
        if(scale_boolean) {
                # Normalizing y
                y_scale <- normalize_bart(y = y_train)

                # Calculating \tau_{\mu} based on the scale of y
                tau_mu <- (4 * n_tree * (K_bart^2))

                # Getting the naive sigma
                nsigma <- naive_sigma(x = x_train,y = y_scale)

                # Getting the shape
                a_tau <- df/2

                # Calculating lambda
                qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
                lambda <- (nsigma*nsigma*qchi)/df
                d_tau <- (lambda*df)/2

                cat("a_tau is " , round(a_tau,digits = 3)," \n")
                cat("d_tau is " , round(d_tau,digits = 3)," \n")

        } else {

                # Not scaling the y
                y_scale <- y_train

                # Calculating \tau_{\mu} based on the scale of y
                # Need to change this value in case of non-scaling
                tau_mu <- (4 * n_tree * (K_bart^2))/((b_max-a_min)^2)
                nsigma <- naive_sigma(x = x_train,y = y_scale)

                # Getting the naive sigma
                nsigma <- naive_sigma(x = x_train,y = y_scale)

                # Getting the shape
                a_tau <- df/2

                # Calculating lambda
                qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
                lambda <- (nsigma*nsigma*qchi)/df
                d_tau <- (lambda*df)/2

                cat("a_tau is " , round(a_tau,digits = 3)," \n")
                cat("d_tau is " , round(d_tau,digits = 3)," \n")

        }

        # a_tau <- 1.5
        # d_tau <- 0.003

        # Defining other quantities
        n_train <- nrow(x_train)
        n_test <- nrow(x_test)
        n_post <- n_mcmc-n_burn

        # Getting the y_hat for train and test
        y_train_hat_post <- matrix(0, ncol = n_train,nrow = n_post)
        y_test_hat_post <- matrix(0, ncol = n_test,nrow = n_post)
        curr <- 0
        y_train_hat_trees <- matrix(0, nrow = n_tree, ncol = n_train)
        y_test_hat_trees <- matrix(0, nrow = n_tree, ncol = n_test)


        # Storing the current partial residuals
        current_partial_residuals_matrix <- matrix(NA,nrow = n_tree, ncol = n_train)
        current_partial_residuals_list <- list()


        # Initialising values for phi_vec, and nu
        phi_vec_matrix <- matrix(1, nrow = n_tree,ncol = ncol(x_train))
        phi_post <- list(n_post)
        nu_vector <- rep(1,n_tree)
        nu_post <- matrix(NA,nrow = n_post, ncol = n_tree)

        tau_post <- numeric(n_post)

        # Getting initial trees
        current_trees <- list()

        for(i in 1:n_tree){
                current_trees[[i]] <- new_tree(x_train = x_train,x_test = x_test)
        }

        # Creating a manual progress bart
        progress_bart_limits <- round(seq(1,n_mcmc,length.out=10))

        #
        for(i in 1:n_mcmc){

                # Small progress bar
                if(i %in% progress_bart_limits){
                        cat(" | ")
                }

                for(t in 1:n_tree){


                        # USING BART BOOLEAN OR NOT
                        if(bart_boolean){
                                partial_residuals <- y_scale - colSums(y_train_hat_trees[-t,,drop = FALSE])

                                # Selecting one verb
                                verb <- sample(c("grow","prune","change"), prob = c(0.3,0.3,0.4),size = 1)

                                if(length(current_trees[[t]])==1){
                                        verb <- "grow"
                                }

                                # Selecting a new tree
                                current_trees[[t]]  <- if(verb=="grow"){
                                        grow(res_vec = partial_residuals,tree = current_trees[[t]],
                                             x_train = x_train,x_test = x_test,xcut = xcut,tau = tau,
                                             tau_mu = tau_mu,alpha = alpha,beta = beta,node_min_size = 1)
                                } else if(verb=="prune"){
                                        prune(tree = current_trees[[t]],res_vec = partial_residuals,
                                              tau = tau,tau_mu = tau_mu,alpha = alpha,beta = beta)
                                } else if(verb=="change"){
                                        change(res_vec = partial_residuals,tree = current_trees[[t]],
                                               x_train = x_train,x_test = x_test,xcut = xcut,
                                               tau = tau,tau_mu = tau_mu,alpha = alpha,beta = beta,node_min_size = 1)
                                }

                                # Updating mu
                                current_trees[[t]] <- update_mu(tree = current_trees[[t]],
                                                                partial_residuals = partial_residuals,
                                                                tau = tau,tau_mu = tau_mu)

                                # Prediction aux
                                pred_obj <- getPrediction(tree = current_trees[[t]])

                                y_train_hat_trees[t,] <- pred_obj$train_pred
                                y_test_hat_trees[t,] <- pred_obj$test_pred

                                # =================================
                        } else {# CALCULATING THE GP-ITERATIONS
                                # =================================

                                # Calculating partial residuals
                                partial_residuals <- y_scale - colSums(y_train_hat_trees[-t,,drop = FALSE])

                                # Storing current partial
                                current_partial_residuals_matrix[t,] <- partial_residuals

                                verb <- sample(x = c("grow","prune","change"),size = 1,prob = c(0.3,0.3,0.4))

                                if(length(current_trees[[t]])==1){
                                        verb <- "grow"
                                }
                                # Selecting one verb movement
                                if(verb == "grow"){
                                        current_trees[[t]] <- grow_gpbart(res_vec = partial_residuals,tree = current_trees[[t]],
                                                                          x_train = x_train,x_test = x_test,xcut = xcut,tau = tau,
                                                                          tau_mu = tau_mu,alpha = alpha,beta = beta,node_min_size = node_min_size,
                                                                          nu = nu_vector[t],phi_vector_p = phi_vec_matrix[t,])
                                } else if( verb == "prune"){
                                        current_trees[[t]] <- prune_gpbart(res_vec = partial_residuals,
                                                                           tree = current_trees[[t]],
                                                                           tau = tau, tau_mu = tau_mu, alpha = alpha, beta = beta,
                                                                           nu = nu_vector[t], phi_vector_p = phi_vec_matrix[t,])

                                } else if( verb == "change"){
                                        current_trees[[t]] <- change_gpbart(res_vec = partial_residuals,tree = current_trees[[t]],
                                                                          x_train = x_train,x_test = x_test,xcut = xcut,tau = tau,
                                                                          tau_mu = tau_mu,alpha = alpha,beta = beta,node_min_size = node_min_size,
                                                                          nu = nu_vector[t],phi_vector_p = phi_vec_matrix[t,])
                                } else {
                                        stop("Error no valid-verb")
                                }



                                # Changing for the current tree
                                # Updating the phi
                                phi_vec_matrix[t,] <- update_phi_gpbart(x_train = x_train,res_vec = partial_residuals,
                                                          phi_vector_p = phi_vec_matrix[t,],nu = nu_vector[t],tau = tau,tau_mu = tau_mu)


                                nu_vector[t] <- update_nu_gpbart(x_train = x_train,res_vec = partial_residuals,
                                                                   phi_vector_p = phi_vec_matrix[t,],nu = nu_vector[t],tau = tau,tau_mu = tau_mu)

                                # Update the mu values
                                current_trees[[t]] <- update_mu_gpbart(tree = tree,x_train = x_train,res_vec = partial_residuals,nu = nu_vector[t],
                                                                  phi_vector_p = phi_vec_matrix[t,],tau = tau,tau_mu = tau_mu)

                                # This one is the most complicated, I need to update the predictions based on the tree structure
                                sample_g_aux <- update_g_gpbart(tree = current_trees[[t]],
                                                                x_train = x_train,
                                                                x_test = x_test,
                                                                res_vec = partial_residuals,
                                                                tau =  tau,
                                                                tau_mu = tau_mu,
                                                                nu = nu_vector[t],
                                                                phi_vector_p = phi_vec_matrix[t,])

                                y_train_hat_trees[t,] <- sample_g_aux$g_sample
                                y_test_hat_trees[t,] <- sample_g_aux$g_sample_test
                        }
                }

                # Storing tau and getting new tau
                tau <- update_tau(y = y_scale,y_hat = colSums(y_train_hat_trees),a_tau = a_tau,d_tau = d_tau)

                tau_post[i] <- tau

                # Storing the posterior elements
                if(i>n_burn){
                        curr = curr + 1
                        current_partial_residuals_list[[curr]] <- current_partial_residuals_matrix
                        y_train_hat_post[curr,] <- colSums(y_train_hat_trees)
                        y_test_hat_post[curr,] <- colSums(y_test_hat_trees)
                        nu_post[curr,] <- nu_vector
                        phi_post[[curr]] <- phi_vec_matrix

                        # Saving the trees
                        if(keeptrees) {
                                post_trees[[curr]] <- current_trees
                        }

                        # Saving the posterior of the hyperparameters
                        phi_post[[curr]] <- phi_vec_matrix

                }


        }


        # Adjusting tau and y_hat for the scale factor
        if(scale_boolean){
                tau_post <- tau_post/((b_max-a_min)^2)
                y_train_hat_post <- unnormalize_bart(z = y_train_hat_post,a = a_min,b = b_max)
                y_test_hat_post <- unnormalize_bart(z = y_test_hat_post,a = a_min,b = b_max)
        }

        # Diagnostic
        plot(y_train,colMeans(y_train_hat_post))

        # Returning the posterior objets
        if(keeptrees){
                post_obj <- list(tau_post = tau_post,
                     y_hat_post = y_train_hat_post,
                     y_test_hat_post = y_test_hat_post,
                     last_trees = current_trees,
                     posterior = list(phi_post = phi_post,
                                      nu_post = nu_post,
                                      partial_residuals = current_partial_residuals_list,
                                      trees = post_trees))
        } else {
                post_obj <- list(tau_post = tau_post,
                     y_hat_post = y_train_hat_post,
                     y_test_hat_post = y_test_hat_post,
                     last_trees = current_trees,
                     posterior = list(phi_post = phi_post,
                                      nu_post = nu_post,
                                      partial_residuals = current_partial_residuals_list))

        }

        return(post_obj)

}
