# Building the GP-BART function
bart <- function(x_train,
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
                 K_bart = 2){

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

        }

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


        tau_post <- numeric()

        # Getting initial trees
        current_trees <- list()

        for(i in 1:n_tree){
                current_trees[[i]] <- new_tree(x_train = x_train,x_test = x_test)
        }


        #
        for(i in 1:n_mcmc){

                for(t in 1:n_tree){

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
                        current_trees[[t]] <- update_mu(tree = current_trees[[t]],partial_residuals = partial_residuals,
                                                        tau = tau,tau_mu = tau_mu)

                        # Prediction aux
                        pred_obj <- getPrediction(tree = current_trees[[t]])

                        y_train_hat_trees[t,] <- pred_obj$train_pred
                        y_test_hat_trees[t,] <- pred_obj$test_pred

                }

                # Storing tau and getting new tau
                tau <- update_tau(y = y_scale,y_hat = colSums(y_train_hat_trees),a_tau = a_tau,d_tau = d_tau)

                tau_post[i] <- tau

                # Storing the posterior elements
                if(i>n_burn){
                        curr = curr + 1
                        y_train_hat_post[curr,] <- colSums(y_train_hat_trees)
                        y_test_hat_post[curr,] <- colSums(y_test_hat_post)
                }


        }


        # Adjusting tau and y_hat for the scale factor
        if(scale_boolean){
                tau_post <- tau_post/((b_max-a_min)^2)
                y_train_hat_post <- unnormalize_bart(z = y_train_hat_post,a = a_min,b = b_max)
                y_test_hat_post <- unnormalize_bart(z = y_test_hat_post,a = a_min,b = b_max)
        }

        # Returning the posterior objets
        return(list(tau_post = tau_post,
                    y_hat_post = y_train_hat_post,
                    y_test_hat_post = y_test_hat_post,
                    last_trees = current_trees))

}
