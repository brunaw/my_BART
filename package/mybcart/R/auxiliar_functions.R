#' @name p_rule
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Rule selection. 
#' @description Selects a value to split the tree in a grow step.
#' @param variable_index The variable to create the split. 
#' @param data The current tree.
#' @param sel_node The node to break from. 
#' @return The selected splitting value.  

p_rule <- function(variable_index, data, sel_node){
  selected_rule <- data %>% 
    dplyr::filter(node == sel_node) %>% 
    dplyr::mutate(var = !!rlang::sym(variable_index)) %>% 
    dplyr::distinct(var) %>% 
    dplyr::filter(var > stats::quantile(var, 0.15), 
                  var < stats::quantile(var, 0.85)) %>% 
    dplyr::pull(var) %>% 
    # selecting the cut point
    base::sample(size = 1)
  
  return(selected_rule)
}

#' Pipe operator
#'
#' See \code{\link[magrittr]{\%>\%}} for more details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @examples{
#'
#'   iris %>% as.matrix()
#'}
NULL