# Creating the main BART function and it's generalities
new_tree <- function(x_train,
                     x_test){

        # Creating the list of node
        node_0 <- list(index = 0,
                       obs_train = 1:nrow(x_train),
                       obs_test  = 1:nrow(x_test),
                       left = NA,
                       right = NA,
                       parent = NA,
                       terminal = 1,
                       nog = 0,
                       depth = 0,
                       var = NA,
                       var_split_rule = NA,
                       mu = 0)

        # Returning the new tree
        return(list(node_0))
}

# Getting terminal nodes
get_terminals <- function(tree){
       return( tree[unlist(lapply(tree, function(t){is.na(t$left)&is.na(t$right)}))] )
}

get_nonterminals <- function(tree){
        return( tree[unlist(lapply(tree, function(t){!is.na(t$left)&!is.na(t$right)}))] )
}

count_nog <- function(tree){
        return(sum(unlist(lapply(tree,function(t){t$nog==1}))))
}


# Getting all nodes index
get_indexes <- function(tree){
        return(unlist(lapply(tree, function(t)(t$index))))
}

# Calculating the node likelihood
node_loglikelihood <- function(node,
                               res_vec,
                               tau,
                               tau_mu){

        # Slicing the current res_vec
        res_node <- res_vec[node$obs_train]
        n_obs <- length(res_node)

        return(-0.5*tau*crossprod(res_node)[[1]]-0.5*log(tau_mu+(n_obs*tau)) + (0.5*(tau^2)*(sum(res_node)^2))/(tau*n_obs+tau_mu))

}

# Creating a function to grow a tree
grow <- function(res_vec,
                 tree,
                 x_train,
                 x_test,
                 xcut,
                 tau,
                 tau_mu,
                 alpha,
                 beta){

        # Getting the terminal nodes
        terminal_nodes <- get_terminals(tree)

        # Sampling one terminal node
        g_node_position <- sample(1:length(terminal_nodes),size = 1)
        g_node <- terminal_nodes[[g_node_position]]

        # Initializing the sample
        split_var_candidates <- 1:ncol(x_train)
        good_tree_index <- 0


        while(good_tree_index==0){
                # Selecting a valid split
                split_var <- sample(split_var_candidates,size = 1)

                # Getting the min and maximum observed value within the terminal node
                min_node_obs <- sort(x_train[g_node$obs_train,split_var])[node_min_size]
                max_node_obs <- sort(x_train[g_node$obs_train,split_var])[length(x_train[g_node$obs_train,split_var])-node_min_size]

                # Getting the column from xcut
                xcut_valid <- xcut[which(xcut[,split_var]>=min_node_obs & xcut[,split_var]<=max_node_obs),split_var]


                # No valid tree found
                if(length(xcut_valid) == 0 ){

                        split_var_candidates <-  split_var_candidates[-which(split_var==split_var_candidates)]

                        if(length(split_var_candidates)==0){
                                return(tree) # There are no valid candidates for this node
                        }

                } else {
                        good_tree_index <- 1
                }
        }
        # Sampling a x_cut_rule
        split_var_sampled_rule <- sample(xcut_valid,size = 1)

        # Creating the left and the right nodes
        max_index <- max(get_indexes(tree))


        # Creating the vector of new train and test index
        left_train_id <- which(x_train[g_node$obs_train,split_var]<=split_var_sampled_rule)
        right_train_id <- which(x_train[g_node$obs_train,split_var]>split_var_sampled_rule)

        left_test_id <- which(x_test[g_node$obs_test,split_var]<=split_var_sampled_rule)
        right_test_id <- which(x_test[g_node$obs_test,split_var]<=split_var_sampled_rule)



        # Creating the left node
        left_node <- list(index = max_index+1,
                          obs_train = left_train_id,
                          obs_test  = left_test_id,
                          left = NA,
                          right = NA,
                          parent = g_node$index,
                          terminal = 1,
                          nog = 0,
                          depth = g_node$index+1,
                          var = split_var_sampled_rule,
                          var_split_rule = split_var_sampled_rule,
                          mu = 0)

        right_node <- list(index = max_index+1,
                          obs_train = right_train_id,
                          obs_test  = right_test_id,
                          left = NA,
                          right = NA,
                          parent = g_node$index,
                          terminal = 1,
                          nog = 0,
                          depth = g_node$index+1,
                          var = split_var,
                          var_split_rule = split_var_sampled_rule,
                          mu = 0)

        # Modifying the new g_node
        new_g_node <-g_node
        new_g_node$left <- max_index+1
        new_g_node$right <- max_index+2
        new_g_node$terminal <- 0


        if(!is.na(g_node$parent)){new_g_node$nog = 1}


        # Get nog counter ( FOR THE NEW TREE )
        nog_counter <- count_nog(tree = tree[-g_node_position]) + 1

        # Calculating the acceptance for two new nodes
        tree_loglikeli <- node_loglikelihood(res_vec = res_vec,node = left_node,tau = tau,tau_mu = tau_mu) +
                node_loglikelihood(res_vec = res_vec,node = rigth_node, tau = tau,tau_mu = tau_mu) -
                node_loglikelihood(res_vec = res_vec,node = g_node,tau = tau, tau_mu = tau_mu)

        # Calculate the transition
        transition_loglike <- log(0.3/nog_counter)-log(length(terminal_nodes)) # prob of getting from the new tree to the old (PRUNE), minus getting to the old to the new (GROW)

        # Calculate the tree prior contribution
        tree_prior <- 2*log(1-alpha*(1+(g_node$depth+1))^(-beta)) + (-beta)*log(alpha*(1+g_node$depth)) - log(1-alpha*(1+g_node$depth)^(-beta))

        log_acceptance <- tree_loglikeli+transition_loglike+tree_prior

        # Accepting the tree ornot
        if(stats::runif(n = 1)<exp(log_acceptance)){
                tree[[g_node_position]] <- new_g_node
                tree[[max_index+1]] <- left_node
                tree[[max_index+2]] <- right_node
        }

        return(tree)

}

