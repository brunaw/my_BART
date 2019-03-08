#' @name bcart
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title B-CART model. 
#' @description This function runs a B-CART model and ...
#' @param formula The model formula.  
#' @param data The data to be used in the modeling.
#' @param iter The number of iterations for the MCMC. 
#' @return A list containing:
#'  sigma^2_y, the current errors, the final tree, the
#'  ratios used in the MCMC step and the uniform values sampled. 
#' @details 
#' Priors used ----------------------------------------------------------
#' sigma^2 ~ InvGamma(nu/2, lambda*nu/2), with
#' the parameters chosen in a way that there is a high probability mass
#' 'before' the RMSE of a linear model 
#' 
#' mu_i ~ Normal(mu_mu, sd_mu^2), with 
#' mu_mu = 0
#' sd_mu = (max(y))/2
#' 
#' Posteriors (both conjugate) ------------------------------------------
#' sigma^2 ~ InvGamma((nu + n) / 2, 1/2 *(sum_i error^2_i + lambda * nu)))
#' Note: the error is calculated with the current mu_i for each node 
#' 
#' mu_i ~ Normal(mu_post, sigma_mu_post), with 
#' sigma_mu_post = 1/(1/sd_mu + n_in_the_node/sigma^2)
#' mu_post = (0 + (n_in_the_node/sigma^2) * node.avgResponse* sigma_mu_post
#' n_in_the_node = observations in each node
#' node.avgResponse = average of the response in each node

bcart <- function(formula, data, iter = 5000){
  # B-CART Function 
  # --------------------------------------------------------------------
  # Extracting the data and the response from the formula
  # --------------------------------------------------------------------
  # Removing the intercept
  formula <- stats::as.formula(paste(c(formula), "- 1"))
  response_name <- all.vars(formula)[1]
  
  # Extracting the model structure
  m <- stats::model.frame(formula, data = data)
  
  X <- stats::model.matrix(formula, m)
  # Scaling the response variable 
  data[ , response_name] <- scale(data[, response_name]) %>% as.vector()
  
  # Defining the response
  y <- stats::model.response(m)
  # renaming the covariates 
  names(data[ , colnames(X)]) <- paste("X", 1:length(colnames(X)))
  
  # --------------------------------------------------------------------
  # Finding a value for the parameters of the prior to sigma
  #---------------------------------------------------------------------
  my_lm <- stats::lm(formula, data)
  res <- sqrt(stats::deviance(my_lm)/stats::df.residual(my_lm))
  
  nu = 0.5
  lambda = 0.5
  p_inv <- invgamma::pinvgamma(q = res, shape = nu/2, rate = nu*lambda/2)
  
  # Putting high probabilities of the BCART improving a linear model
  while(p_inv < 0.8){
    p_inv <- invgamma::pinvgamma(q = res, shape = nu/2, rate = nu*lambda/2)
    if(p_inv < 0.8){
      nu = abs(nu + stats::rnorm(1))
      lambda = abs(lambda + stats::rnorm(1))
    }
  }
  
  #---------------------------------------------------------------------
  # Defining current distributions parameters
  #---------------------------------------------------------------------
  # Grow and prune probabilities
  p_grow = 0.5
  p_prune = 0.5
  
  # Prior parameters for mu: mean and std. deviation 
  mu <- list()
  mu[[1]] <- mu_mu <- 0
  k <- 2
  sd_mu <- (max(y))/k
  
  # Minimum batch for each node
  keep_node <- 0.005 * nrow(data)
  p = ncol(X)
  #---------------------------------------------------------------------
  # Initializing acessory columns in the data
  #---------------------------------------------------------------------
  data$node <- "root"           # To save the current node
  data$parent <- "root"         # To save the current parent of each node
  data$d = 0                    # To save the current depth of each node
  data$node_index = 1           # To save the current index of each node
  data$criteria = 'left'        # To initializa the root as a 'left' node
  
  #---------------------------------------------------------------------
  # Initializing useful vectors
  #---------------------------------------------------------------------
  # For grow or prune
  action_taken = vector()       # To save the actions taken in the algorithm
  selec_var = vector()          # To save the selected variable when growing
  rule = vector()               # To save the selected splittinh rule when growing
  drawn_node = vector()         # To save the selected node to grow or prune
  
  # For the sampling
  sigma_2 = vector()           # To save the simulated values for sigma
  # Prior to sigma^2
  sigma_2[1] <- invgamma::rinvgamma(n = 1, shape = nu/2, rate = nu*lambda/2)
  
  r <- vector()              # To save the ratios of grow or prune 
  u <- vector()              # To save the sampled uniform values
  
  # For the trees 
  my_trees = list()            # To save each new tree 
  my_trees[[1]] <- data        # Initializing the first tree as the 
  # data (the 'root')
  n <- nrow(data)              # Total number of observations
  
  sigma_mu_post <- vector()    # To save posterior values of 
  mu_mu_post <- vector()       # mu_mu and sigma_mu
  mu_vec <- vector()           # To save the sampled mus
  
  errors <- vector()          # To save the value of the errors used in
  parent_action <- vector()          # the posterior of sigma
  
  results <- list()
  results[[1]] <- data.frame(node = NA, var = NA, rule = NA)
  verif <- list()
  #---------------------------------------------------------------------
  # A simple progress bar 
  pb <- progress::progress_bar$new(
    format = "  Iterations of the BCART model [:bar]  :current/:total (:percent)",
    clear = FALSE, width = 60, total = iter)
  
  
  # Loop to perform the BCART model 
  for(i in 1:iter){
    # Sampling a number from the uniform distribution
    k[i] <- stats::runif(1)
    
    # Checking if the depths of the nodes are bigger than 0 
    #  (if 0, we can't prune the tree)
    depths <- sum(unique(my_trees[[i]]$d) > 0)
    
    # Growing, pruning or staying in the same tree ---------------------
    # Should we grow the tree?
    if(k[i] <  p_grow){
      action_taken[i] = "grow"
      
      # Selecting the node to grow, uniformly
      drawn_node[i] <- sample(unique(my_trees[[i]]$node), size = 1)
      # Selecting the variable and splitting rule, uniformly
      selec_var[i] <- colnames(X)[sample(1:ncol(X), size = 1)]  
      rule[i] <- mybcart::p_rule(selec_var[i], data = my_trees[[i]], 
                                 sel_node = drawn_node[i])
      
      my_trees[[i+1]] <- my_trees[[i]] %>% 
        dplyr::mutate(
          # Increasing the depth of the node giving the grow 
          d =  ifelse(node == drawn_node[i], d + 1, d), 
          # Updating the parent of the splitted node 
          parent = ifelse(node == drawn_node[i], drawn_node[i], parent),
          
          # Changing the node "side" of each observation: left and right
          criteria = ifelse(
            node == drawn_node[i],
            ifelse(!!rlang::sym(selec_var[i]) > rule[i], 
                   "left", "right"), "no split"),
          
          # Updating the node accordingly to the new split
          node = ifelse(node == drawn_node[i], 
                        ifelse(!!rlang::sym(selec_var[i]) > rule[i],
                               paste(node, selec_var[i], "left"), 
                               paste(node, selec_var[i], "right")), node), 
          # Updating the node index
          node_index =  as.numeric(as.factor(node)))
      
      # Checking whether all of the nodes have the minimum
      # of observations required
      temps <-  my_trees[[i + 1]]  %>% 
        dplyr::count(node) %>% 
        dplyr::pull(n)
      
      if(sum(temps <= keep_node) > 0){
        # If all nodes don't have the minimum of observations, 
        # just keep the previous tree 
        my_trees[[i + 1]] <- my_trees[[i]]
        results[[i+1]] <- results[[i]]
      } else{
        
        
        # Saving the parent of the new node to use in the 
        # calculation of the transition ratio
        parent_action[i] <- my_trees[[i + 1]] %>% 
          dplyr::filter(parent == drawn_node[i]) %>% 
          dplyr::distinct(parent) %>% 
          dplyr::pull(parent)
        
        results[[i+1]] <- suppressWarnings(
          dplyr::bind_rows(stats::na.omit(results[[i]]), 
                           data.frame(
                             node = parent_action[i], 
                             var = selec_var[i], 
                             rule = rule[i])) 
        )
        
        
        # Calculating the acceptance ratio for the grow, 
        # which uses the ratios of: 
        # 1. The transition probabilities of the trees
        # 2. The likelihoods of the trees
        # 3. The probability of the tree structures
        
        r[i] <- mybcart::ratio_grow(tree = my_trees[[i + 1]], 
                                    current_node =  parent_action[i],
                                    sigma_2_mu = (sd_mu^2), 
                                    sigma_2_y =  sigma_2[i], 
                                    current_selec_var = selec_var[i],
                                    p = ncol(X))
        
      } 
      # Should we prune the tree?
    } else if(k[i] >  p_grow && depths > 0){
      
      action_taken[i] = "prune"

      # Selecting node to prune, uniformly 
      drawn_node[i] <- sample(unique(my_trees[[i]]$node), size = 1)
      # Detect the two nodes (left and right) to be pruned
      nodes_to_prune <- stringr::str_remove(drawn_node[i], '( right| left)$')
      # Detect the variable that was split in the node that will 
      # be now pruned
      variable_in_question <- stringr::str_extract(drawn_node[i], 
                                                   'X[0-9][^X[0-9]]*$') %>% 
        stringr::str_remove(" left| right")
      # Updating the node that will suffer the prune (returning
      # to the point where it was before the growing)
      new_node <- nodes_to_prune %>% stringr::str_remove('( X[0-9])$')
      
      my_trees[[i + 1]] <- my_trees[[i]] %>% 
        dplyr::mutate(
          # Calculating the reduction in the depth of the pruned 
          # node 
          d_reduction = 
            ifelse(
              stringr::str_detect(node, nodes_to_prune), 
              stringr::str_match(node, 
                                 paste0('(?<=', nodes_to_prune, '\\s).*')) %>% 
                stringr::str_count(pattern = "X[0-9]") + 1, 0),
          
          # Reducing the depth of the nodes
          d = ifelse(
            stringr::str_detect(node, nodes_to_prune), d - d_reduction, d),
          
          # Changing the node for the new_node (with the pruning) 
          temp_node = ifelse(
            stringr::str_detect(node, nodes_to_prune), paste(new_node), node),  
          
          # Removing from the new_node the last occurrence of the name of 
          # a variable + the side of the node, to update the parent
          parent = ifelse(
            stringr::str_detect(node, nodes_to_prune),
            stringr::str_remove(temp_node, 
                                paste0('( X[0-9] right| X[0-9] left)$')), parent),
          
          # Updating the node and node index
          node = temp_node, 
          node_index = as.numeric(as.factor(node))) %>%  
        # Discarding temporary variables
        dplyr::select(-temp_node, -d_reduction)
      
      
      # Saving the node that was pruned to use in the 
      # calculation of the transition ratio
      parent_action[i] <- my_trees[[i]] %>% 
        dplyr::filter(node == drawn_node[i]) %>% 
        dplyr::distinct(parent) %>% 
        dplyr::pull(parent)
      
      parent_prune <- stringr::str_remove(parent_action[i], '( right| left)$')
      
      results[[i+1]] <- results[[i]] %>% 
        dplyr::filter(!stringr::str_detect(node, parent_prune))
      
      # Calculating the acceptance ratio for the prune, 
      # which uses the ratios of: 
      # 1. The transition probabilities of the trees
      # 2. The likelihoods of the trees
      # 3. The probability of the tree structures
      
      r[i] <- mybcart::ratio_prune(old_tree = my_trees[[i]],
                                   tree = my_trees[[i + 1]], 
                                   current_node =  parent_action[i],
                                   var_in_prune = variable_in_question, 
                                   sigma_2_mu = (sd_mu^2), 
                                   sigma_2_y =  sigma_2[i], 
                                   p = ncol(X), 
                                   nodes_to_prune = nodes_to_prune)
      
      # Should we stay in the same tree?   
    } else {
      my_trees[[i + 1]] <- my_trees[[i]]
      results[[i+1]] <- results[[i]]
    }
    
    # Checking if the tree will be accepted or not ---------------------
    if(!identical(my_trees[[i + 1]], my_trees[[i]])){
      
      # Should we accept the new tree?
      u[i] <- stats::runif(1)
      
      # Checking if the tree should be accepted or not, based
      # on the acceptance ratio calculated and a value sampled
      # from a uniform distributioon
      if(u[i] >= r[i]){
        my_trees[[i + 1]] <- my_trees[[i]] 
        results[[i+1]] <- results[[i]]
      }
    }
    
    # Updating posteriors ----------------------------------------------
    
    # Sampling from the posterior distribution of mu_mu
    # Calculating how many nodes we currently have 
    m_mu <- my_trees[[i + 1]] %>% 
      dplyr::distinct(node) %>% 
      nrow()
    
    # Caclulating the average and n of each node
    res_node <- my_trees[[i+1]] %>% 
      dplyr::group_by(node) %>% 
      dplyr::summarise(mu_avg = mean(y), n = dplyr::n())
    
    # Separating the average in each node
    avg_node <- res_node %>%  dplyr::arrange(node) %>% dplyr::pull(mu_avg)
    # Separating the number of observations in each node
    n_node <- res_node %>% dplyr::arrange(node) %>% dplyr::pull(n)
    
    # Updating posteriors
    for(j in 1:m_mu){
      
      sigma_mu_post[j] <- 1/((1/sd_mu^2) + n_node[j]/sigma_2[i])
      mu_mu_post[j] <-  ((n_node[j]/sigma_2[i]) * avg_node[j]) * sigma_mu_post[j]
      
      mu_vec[j] <- stats::rnorm(n = 1, 
                                mean = mu_mu_post[j],
                                sd = sqrt(sigma_mu_post[j]))
    }
    
    # Saving the results for the posterior of mu
    mu[[i + 1]] <- mu_vec
    
    # Resetting the auxiliar vectors 
    sigma_mu_post = mu_mu_post = mu_vec = vector()
    
    # Sampling from the posterior distribution of sigma^2_y
    
    # Calculating the sum of squares, given the current values 
    # sampled for mu in each node
    errors[i] <- my_trees[[i+1]] %>% 
      dplyr::distinct(node) %>% 
      dplyr::arrange(node) %>% 
      dplyr::mutate(mu_samp = mu[[i+1]]) %>% 
      dplyr::right_join(my_trees[[i+1]], by = 'node') %>% 
      dplyr::mutate(err = (y - mu_samp)^2) %>% 
      dplyr::summarise(sum_errors = sum(err)) %>% 
      dplyr::pull(sum_errors)
    
    # Updating posterior for sigma^2_y
    sigma_2[i+1] <- invgamma::rinvgamma(n = 1, shape = (nu+n)/2, 
                              (errors[i] + lambda * nu)/2)
    
    i <- i + 1
    pb$tick()
  }
  
  cat("
# ----------------------------------------------------------  
The tree has finished growing, with a final error 
of:", errors[i-1], "and a final posterior variance of:", sigma_2[i-1], 
"\nThe minimum node depth is", min(my_trees[[i-1]]$d), "and the maximum is ", max(my_trees[[i-1]]$d),   
'\n# ---------------------------------------------------------- ') 
  
  return(list(
    sigma_2 = sigma_2,
    errors = errors, 
    final_tree = my_trees[[i-1]],
    ratios = r,
    samp_unif = u,
    nodes = drawn_node,
    action_taken = action_taken,
    results = results[[i]], 
    mu = mu[[i-1]], 
    model_formula = formula
  ))
}