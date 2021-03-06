#' @name transition_ratio_grow
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Transition ratio for grow. 
#' @description Transition probability of going to the candidate tree, 
#' given that the action step was a grow.
#' @param tree The current tree.
#' @param current_node The current node.
#' @param p The number of predictors still available. 
#' @param current_selec_var The variable selected for the split. 
#' @param results The current results to find the number of second
#' generation internal nodes. 
#' @return The transition ratio 
#' @details When transitioning to a grow, we need the probabilities of:
#'  1. Growing the tree
#'  2. Growing from the specific node, that has to consider the 
#'  available predictors and available values for that grow. 
#' When transitioning from a grow back to the split, 
#' we need the probabilities of:
#' 1. Pruning the tree
#' 2. Selecting  node to prune
#' @example
#' transition_ratio_grow(tree, current_node) 

transition_ratio_grow <- function(tree, current_node, p, 
                                  current_selec_var, results){
  
  p_grow = 0.5
  p_prune = 0.5
  # Number of available final nodes to break -------
  b <-  tree %>% dplyr::distinct(node_index) %>% nrow()
  
  # Probability of splitting a variable ------------
  p_adj = 1/p
  
  # Available values to split given the variable selected for the
  # split and the node selected ----------------------
  n_j_adj <-  tree %>% 
    dplyr::filter(parent == current_node) %>% 
    dplyr::distinct(!!rlang::sym(current_selec_var)) %>% nrow()
  
  # Probability of the transition -------------------
  p_t_to_tstar <- p_grow * (1/b) * (p_adj) * (1/n_j_adj)
  
  
  #  Probability of transitioning from the new tree
  # back to the original -------------------------
  w_2 <-  nrow(results)
  p_tstar_to_t <- p_prune/w_2
  
  trans_ratio <- log(p_t_to_tstar/p_tstar_to_t)
  return(trans_ratio)
}

#' @name lk_ratio_grow
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Likelihood ratio for grow. 
#' @description Likelihood ratio of the candidate tree and the previous 
#' tree, given that the action step was a grow. 
#' @param tree The current tree.
#' @param current_node The current pruned node.
#' @param sigma_2_y The current value of sigma^2 for y.
#' @param sigma_2_mu The current valur of sigma^2_mu. 
#' @return The likelihood ratio. 
#' @details The likelihood ratio for growth needs the 
#' joint distribution of each region (of each node), given the 
#' parameter sigma^2 (here, mu is considered to be integrated out)
#' @example
#' lk_ratio_grow(tree, current_node, sigma_2_y, sigma_2_mu)

lk_ratio_grow <- function(tree, current_node, sigma_2_y, sigma_2_mu){

  filtered_tree <- tree %>% dplyr::filter(parent == current_node)
  
  # The first node is on the left, the second is on the right,
  # meaning thta the left node has the smaller index --------
  
  # Counting how many observations are in each 
  # region (left and right) ---------------------------------
  nl_nr <-  filtered_tree %>% 
      dplyr::count(node_index) %>% 
    dplyr::arrange(n) %>% 
    dplyr::pull(n) 
  
  # Calculating the sums of y in each region ----------------
  sums_nodes <- filtered_tree %>% 
    dplyr::group_by(node_index) %>% 
    dplyr::summarise(sum = sum(y)) %>% 
    dplyr::arrange(node_index) %>% 
    dplyr::pull(sum)
  
  # Calculating the equation of the lk. ratio ---------------
  first_term <- log(sqrt(
    ((sigma_2_y*(sigma_2_y + sigma_2_mu * sum(nl_nr)))/
       ((sigma_2_y + sigma_2_mu * nl_nr[1])*
       (sigma_2_y + sigma_2_mu * nl_nr[2]))
    )
  ))
  
  # Exponential part -----------------------------------------
  first_term_exp <- sigma_2_mu/(2*sigma_2_y)
  # Left node is in the first position of the objects
  second_term_exp <- (sums_nodes[1]^2)/(sigma_2_y + nl_nr[1] * sigma_2_mu)
  # Right node is in the second position of the objects
  third_term_exp <- (sums_nodes[2]^2)/(sigma_2_y + nl_nr[2] * sigma_2_mu)
  
  fourth_term_exp <- (sum(sums_nodes)^2)/(sigma_2_y + sum(nl_nr) * sigma_2_mu)
  
  
  # The exponential part 
  # check if the - sign is there or not
  exp_part <- first_term_exp * 
    (second_term_exp + third_term_exp - fourth_term_exp)
  
  lk_ratio <- exp_part + first_term
  return(lk_ratio)
}

#' @name structure_ratio_grow
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Tree structure ratio for grow
#' @description Tree structure ratio of the candidate tree and 
#' the previous tree, given that the action step was a grow.
#' @param tree The current tree.
#' @param current_node The current grown node.
#' @param current_selec_var The variable selected for the split. 
#' @param p The number of available predictors. 
#' @return The tree structure ratio. 
#' @details For the tree structure ratio of the new tree, we need 
#' the probabilities of:
# 1. Splitting at node n
# 2. Splitting at the node of the left
# 3. Splitting at the node of the right
# 4. Using each rule at node n
#' @example 
#' structure_ratio_grow(tree, current_node)

structure_ratio_grow <- function(tree, current_node, 
                                 current_selec_var, p){

  # Finding the probability of selecting one
  # available predictor -------------------------------------
  p_adj <- 1/p
  
  # Counting the distinct rule options from
  # this available predictor -------------------------------
  n_j_adj <-  tree %>% 
    dplyr::filter(parent == current_node) %>% 
    dplyr::distinct(!!rlang::sym(current_selec_var)) %>% nrow()
  
  # Calculating the probability of the chosen rule --------
  p_rule <- p_adj * (1/n_j_adj)
  
  # Calculating the probability of split
  terminal_nodes <- tree %>% dplyr::distinct(node_index) %>% nrow()
  p_split <- 1/terminal_nodes

  p_t_star <- ((1-p_split)^2)*p_split*p_rule
  
  p_t <- (1 - p_split)
  
  st_ratio <- log(p_t_star/p_t)
    
  return(st_ratio)
}

#' @name ratio_grow
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Final ratio for a growth step.
#' @description The final ratio is to be used as the acceptance 
#' criteria in the MCMC of the b-cart model.
#' @param tree The current tree.
#' @param current_node The current grown node.
#' @param sigma_2_y The current value of sigma^2 for y.
#' @param sigma_2_mu The current valur of sigma^2_mu. 
#' @param p The number of available predictors
#' @param current_selec_var The variable selected for the split. 
#' @param results The current results to find the number of second
#' generation internal nodes. 
#' @return The final ratio for the candidate tree. 
#' @example 
#' ratio_grow(tree, current_node, sigma_2_mu, sigma_2)

ratio_grow <- function(tree, current_node, sigma_2_mu, 
                       sigma_2_y, p, current_selec_var, 
                       results){
  # All ratios:
  trans <- transition_ratio_grow(tree, current_node, 
                                 current_selec_var = current_selec_var,
                                 p = p, results = results)
  lk <- lk_ratio_grow(tree, current_node, sigma_2_y, sigma_2_mu)
  struct <- structure_ratio_grow(tree, current_node, 
                                 current_selec_var = current_selec_var,  p = p)
  
  r <- min(1, exp(trans + lk + struct))
  return(r)
}
