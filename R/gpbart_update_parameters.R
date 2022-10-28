update_phi_gpbart <- function(tree,
                              x_train,
                              res_vec,
                              nu,
                              phi_vector_p,
                              tau,
                              tau_mu,
                              gp_variables){

        # Getting terminal nodes
        t_nodes <- get_terminals(tree)

        up_crossings <- c(0.001,0.01,0.05,0.1,0.5,1,2,5,10,100)

        for(i in 1:length(phi_vector_p)){
                phi_proposal <- sample(x = 1/(2*pi*up_crossings),size = 1)
                new_phi_vector_p <- phi_vector_p
                new_phi_vector_p[i] <- phi_proposal

                old_log_like <- Reduce("+",lapply(t_nodes, function(node){ node_loglikelihood_gpbart(node = node,res_vec = res_vec,x_train = x_train[,gp_variables],
                                                                                          tau = tau,tau_mu = tau_mu,nu = nu,phi_vector = phi_vector_p)}))
                new_log_like <- Reduce("+",lapply(t_nodes, function(node){ node_loglikelihood_gpbart(node = node,res_vec = res_vec,x_train = x_train[,gp_variables],
                                                                                                     tau = tau,tau_mu = tau_mu,nu = nu,phi_vector = new_phi_vector_p)}))


                # Calculating acceptance
                acceptance <- exp(new_log_like-old_log_like + stats::dgamma(x = phi_proposal,shape = 5,rate = 1,log = TRUE) - stats::dgamma(x = phi_vector_p[i],shape = 5,rate = 1,log = TRUE))

                if(stats::runif(n = 1,min = 0,max = 1)<=acceptance){
                        phi_vector_p <- new_phi_vector_p
                }
        }

        return(phi_vector_p)
}


update_nu_gpbart <- function( tree,
                              x_train,
                              res_vec,
                              nu,
                              phi_vector_p,
                              tau,
                              tau_mu){

        # Getting terminal nodes
        t_nodes <- get_terminals(tree)

        nu_proposal <- stats::runif(n = 1,min = 1,max = 100)

        old_log_like <- Reduce("+",lapply(t_nodes, function(node){ node_loglikelihood_gpbart(node = node,res_vec = res_vec,x_train = x_train,
                                                                                             tau = tau,tau_mu = tau_mu,nu = nu,phi_vector = phi_vector_p)}))
        new_log_like <- Reduce("+",lapply(t_nodes, function(node){ node_loglikelihood_gpbart(node = node,res_vec = res_vec, x_train = x_train,
                                                                                             tau = tau,tau_mu = tau_mu,nu = nu_proposal,phi_vector = phi_vector_p)}))

        # Calculating acceptance
        acceptance <- exp(new_log_like-old_log_like)

        if(stats::runif(n = 1,min = 0,max = 1)<=acceptance){
                return(nu_proposal)
        } else {
                return(nu)
        }

}



# Update mu node
update_mu_node <- function(node,
                           x_train,
                           x_test,
                           res_vec, # The complete res_vector
                           tau,
                           tau_mu,
                           nu,
                           phi_vector_p){

        # Calculating the x_train from node
        x_train_node <- x_train[node$obs_train,,drop = FALSE]
        res_node <- res_vec[node$obs_train]
        # Calculating the v factor from equation 11
        distance_sq_matrix <- symm_distance_matrix(m1 = x_train_node,phi_vector = phi_vector_p)
        omega_plus_tau_diag <- (nu^(-1))*exp(-distance_sq_matrix)+diag(1/tau,nrow = nrow(distance_sq_matrix))
        inv_omega_plus_tau <- chol2inv(chol(omega_plus_tau_diag))
        # Calculating S
        S <- sum(inv_omega_plus_tau)+tau_mu

        # Calculating the mean
        mu_mean <- (S^(-1))*crossprod(res_node,inv_omega_plus_tau)
        mu_var <- S^(-1)

        node$mu <- stats::rnorm(n = 1,mean = mu_mean,sd = sqrt(mu_var))

        return(node)
}


update_mu_gpbart <- function(tree,
                             x_train,
                             res_vec,
                             nu,
                             phi_vector_p,
                             tau,
                             tau_mu){

        # Getting terminal nodes
        t_nodes <- get_terminals(tree)

        new_t_nodes_tree <-lapply(t_nodes, function(node){ update_mu_node(node = node,x_train = x_train,res_vec = res_vec,
                                                                          tau = tau,tau_mu = tau_mu,nu = nu,phi_vector_p = phi_vector_p)})

        # Updating all nodes
        tree[names(t_nodes)] <- new_t_nodes_tree

        return(tree)

}

# Update mu node
update_g_node <- function(node,
                           x_train,
                           x_test,
                           res_vec, # The complete res_vector
                           tau,
                           tau_mu,
                           nu,
                           phi_vector_p,
                           test_only = FALSE){

        # Calculating the x_train from node
        x_train_node <- x_train[node$obs_train,,drop = FALSE]
        x_test_node <- x_test[node$obs_test,,drop = FALSE]

        # Creating the vetors
        g_sample <- numeric(nrow(x_train_node))
        g_sample_test <- numeric(nrow(x_test_node))

        res_node <- res_vec[node$obs_train]
        # Calculating the v factor from equation 11
        distance_sq_matrix <- symm_distance_matrix(m1 = x_train_node,phi_vector = phi_vector_p)
        omega <- (nu^(-1))*exp(-distance_sq_matrix)
        omega_plus_tau_diag <- omega+ diag(1/tau,nrow = nrow(distance_sq_matrix))
        inv_omega_plus_tau <- chol2inv(chol(omega_plus_tau_diag))



        # Calculating S
        S <- sum(inv_omega_plus_tau)+tau_mu

        if(!test_only){
                # Calculating the mean
                g_mean <- node$mu + crossprod(omega,crossprod(inv_omega_plus_tau,(res_node-node$mu)))
                g_var <- omega - crossprod(omega,crossprod(inv_omega_plus_tau,omega))
        }

        # Diagnostics
        # plot(res_node,g_mean)

        # Calculating test quantities
        distance_sq_matrix_test_star <- distance_matrix(m1 = x_train_node,m2 = x_test_node,phi_vector = phi_vector_p)
        distance_sq_matrix_test_star_star <- symm_distance_matrix(m1 = x_test_node,phi_vector = phi_vector_p)

        omega_star <- (nu^(-1))*exp(-distance_sq_matrix_test_star)
        omega_star_star <- (nu^(-1))*exp(-distance_sq_matrix_test_star_star)


        g_test_mean <- node$mu+crossprod(omega_star,crossprod(inv_omega_plus_tau,res_node-node$mu))
        g_test_var <- omega_star_star - crossprod(omega_star,crossprod(inv_omega_plus_tau,omega_star))
        # ====

        # Returning g sample
        if(!test_only){
                g_sample <- mvnfast::rmvn(n = 1,mu = g_mean,sigma = g_var+diag(1e-8,nrow = nrow(g_var)))
        }

        g_sample_test <- mvnfast::rmvn(n = 1,mu = g_test_mean,sigma = g_test_var+diag(1e-6,nrow = nrow(g_test_var)))

        return(list(train_sample = g_sample, test_sample = g_sample_test))
}


update_g_gpbart <- function(tree,
                            x_train,
                            x_test,
                            res_vec,
                            tau,
                            tau_mu,
                            nu,
                            phi_vector_p){

        # new g
        g_sample_train <- numeric(nrow(x_train))
        g_sample_test <- numeric(nrow(x_test))
        # Getting terminal nodes
        t_nodes <- get_terminals(tree)

        # Iterate over terminal nodes
        for(i in 1:length(t_nodes)){
               g_aux <- update_g_node(node = t_nodes[[i]],x_train = x_train,x_test = x_test,
                                      res_vec = res_vec,tau = tau,
                                      tau_mu = tau_mu,nu = nu,phi_vector_p = phi_vector_p)

               g_sample_train[t_nodes[[i]]$obs_train] <- g_aux$train_sample
               g_sample_test[t_nodes[[i]]$obs_test] <- g_aux$test_sample
        }


        if(any(is.na(g_sample_train)) || any(is.na(g_sample_test))){
                stop("error somehwere in the prediction")
        }
        return(list(g_sample = g_sample_train, g_sample_test = g_sample_test))
}

