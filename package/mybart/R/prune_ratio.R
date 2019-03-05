#' @name transition_ratio_prune
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Transition ratio for prune.
#' @description Transition probability of going to the candidate tree, 
#' given that the action step was a prune.
#' @param tree The current tree.
#' @param current_node The current node.
#' @param p The number of predictors still available. 
#' @return The transition ratio 
#' @details When transitioning to a prune, we need the probabilities of:
#' 1. Pruning the tree
#' 2. Selecting  node to prune
#' When transitioning from a prune back to the split, 
#' we need the probabilities of:
#'  1. Growing the tree
#'  2. Growing from the specific pruned node, that has to consider the 
#'  available predictors and available values for that grow.   
#' @example
#' transition_ratio_prune(tree, current_node) 

transition_ratio_prune <- function(tree, current_node, p = ncol(X)){
  # Number of available final nodes to prune -------
  b <-  tree %>% distinct(node_index) %>% nrow()
  # Number of internal node  -----------------------
  w_2 <-  b - 1
  # Probability of pruning -------------------------
  p_t_to_tstar <- p_prune/w_2

  # Probability of splitting a variable ------------
  p_adj <- p
  
  # Available values to split ----------------------
  n_j_adj <-  tree %>% 
    dplyr::filter(parent == current_node) %>% 
    distinct(!!sym(selec_var[i])) %>% nrow()
  
  # Probability of transitioning from the new tree
  # to the old one --------------------------------
  p_tstar_to_t <- p_grow * 1/((b-1) * p_adj * n_j_adj)
  
  # Transition ratio in log scale
  trans_ratio <- log(p_tstar_to_t/p_t_to_tstar)
  return(trans_ratio)
}

#' @name lk_ratio_prune
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Likelihood ratio for prune.
#' @description Likelihood ratio of the candidate tree and the previous 
#' tree, given that the action step was a prune.
#' @param tree The current tree.
#' @param current_node The current pruned node.
#' @param sigma_2_y The current value of sigma^2 for y.
#' @param sigma_2_mu The current valur of sigma^2_mu. 
#' @return The likelihood ratio. 
#' @details For the likelihood ratio of the new pruned tree, we need 
#' to calculate an inversion of the one in the grow version. 
#' The value is based on the likelihood of all regions, given the
#' tree and the parameters, for the new and the previous trees. 
#' @example
#' lk_ratio_prune(tree, current_node, sigma_2_y, sigma_2_mu)


lk_ratio_prune <- function(tree, current_node, sigma_2_y, sigma_2_mu){
  
  filtered_tree <- tree %>% filter(parent == current_node)
  
  # The first node is on the left, the second is on the right,
  # meaning that the left node has the smaller index --------
  
  # Counting how many observations are in each 
  # region (left and right) ---------------------------------
  nl_nr <-  filtered_tree %>% count(node_index) %>% arrange(n) %>% pull(n) 
  
  # Calculating the sums of y in each region ----------------
  sums_nodes <- filtered_tree %>% 
    group_by(node_index) %>% 
    summarise(sum = sum(y)) %>% 
    arrange(node_index) %>% 
    pull(sum)
  
  # Calculating the equation of the lk. ratio ---------------
  first_term <- log(sqrt(
    (((sigma_2_y + sigma_2_mu * nl_nr[1])*
          (sigma_2_y + sigma_2_mu * nl_nr[2])) / 
       (sigma_2_y*(sigma_2_y + sigma_2_mu * sum(nl_nr)))
    )
  ))
  
  # Exponential part -----------------------------------------
  first_term_exp <- sigma_2_mu/(2*sigma_2_y)
  # Left node is in the first position of the objects
  second_term_exp <- (sums_nodes[1]^2)/(sigma_2_y + nl_nr[1] * sigma_2_mu)
  # Right node is in the second position of the objects
  third_term_exp <- (sums_nodes[2]^2)/(sigma_2_y + nl_nr[2] * sigma_2_mu)
  
  fourth_term_exp <- (sum(sums_nodes)^2)/(sigma_2_y + sum(nl_nr) * sigma_2_mu)
  
  
  # The exponential part: it's an inversion of the one 
  # calculated in the grow version --------------------------
  exp_part <- first_term_exp * 
    (fourth_term_exp - third_term_exp - second_term_exp)
  
  lk_ratio <- exp_part + first_term
  return(lk_ratio)
}


#' @name structure_ratio_prune
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Tree structure ratio for prune.
#' @description Tree structure ratio of the candidate tree and 
#' the previous tree, given that the action step was a prune.
#' @param tree The current tree.
#' @param current_node The current pruned node.
#' @return The tree structure ratio. 
#' @details For the tree structure ratio of the new pruned tree, we need 
#' to calculate an inversion of the tree structure ratio for the grow. 
#' We need the probabilities of:
#' 1. Splitting at node n
#' 2. Splitting at the node of the left
#' 3. Splitting at the node of the right
#' 4. Using each rule at node n
#' @example 
#' structure_ratio_prune(tree, current_node)

structure_ratio_prune <- function(tree, current_node){
  
  # Finding the probability of selecting one
  # available predictor -------------------------------------
  p <-  ncol(X)
  p_adj <- 1/p
  
  # Counting the distinct rule options from
  # this available predictor -------------------------------
  n_j_adj <-  tree %>% 
    filter(parent == current_node) %>% 
    distinct(!!sym(selec_var[i])) %>% nrow()
  
  # Calculating the probability of the chosen rule --------
  p_rule <- p_adj * (1/n_j_adj)
  
  # Calculating the probability of split
  terminal_nodes <- tree %>% distinct(node_index) %>% nrow()
  p_split <- 1/terminal_nodes
  
  p_t <- ((1-p_split)^2)*p_split*p_rule
  
  p_t_star <- (1 - p_split)
  
  st_ratio <- log(p_t_star/p_t)
  
  return(st_ratio)
}


#' @name ratio_prune
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Final ratio for a prune step.
#' @description The final ratio is to be used as the acceptance 
#' criteria in the MCMC of the b-cart model.
#' @param tree The current tree.
#' @param current_node The current pruned node.
#' @param sigma_2_y The current value of sigma^2 for y.
#' @param sigma_2_mu The current valur of sigma^2_mu. 
#' @param p The number of available predictors
#' @return The final ratio for the candidate tree. 
#' @example 
#' ratio_prune(tree, current_node, sigma_2_mu, sigma_2)

ratio_prune <- function(tree, current_node, sigma_2_mu, 
                         sigma_2_y, p = ncol(X)){
  # All ratios:
  trans <- transition_ratio_prune(tree, current_node)
  lk <- lk_ratio_prune(tree, current_node, sigma_2_y, sigma_2_mu)
  struct <- structure_ratio_prune(tree, current_node)
  
  r <- min(1, exp(trans+lk+struct))
  return(r)
}


