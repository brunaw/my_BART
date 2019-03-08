#' @name predict_bcart
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Predictions for the B-CART model. 
#' @description This function predicts for the final tree of a 
#' Bayesian CART model. 
#' @param model The model object.  
#' @param newdata The new data to predict. 
#' @return A dataframe with the "prediction" column. 
predict_bcart <- function(model, newdata){
  
  formula <- model$model_formula
  response_name <- all.vars(formula)[1]
  # Scaling the response variable 
  newdata[ , response_name] <- scale(newdata[, response_name]) %>% as.vector()
  
  m <- stats::model.frame(formula, data = newdata)
  X <- stats::model.matrix(formula, m)
  # renaming the covariates 
  names(newdata[ , colnames(X)]) <- paste("X", 1:length(colnames(X)))
  
  res <- model$results %>% dplyr::mutate_if(is.factor, as.character)
  newdata$node <- "root"
  pred <- newdata

  # Creating the nodes in the newdata
  for(i in 1:nrow(res)){
    
    pred <- pred %>% 
      dplyr::mutate(node = 
                      ifelse(node == res$node[i], 
                             ifelse(!!rlang::sym(res$var[i]) > res$rule[i],
                                    paste(node, res$var[i], "left"), 
                                    paste(node, res$var[i], "right")),
                             node))
  }
  # Inserting the mu's and returning the data.frame
  nodes_and_mus <- data.frame(node = model$final_tree %>% 
                                dplyr::distinct(node) %>% 
                                dplyr::arrange(node) ,
                              prediction = model$mu)
  
  pred <- pred %>% 
    dplyr::left_join(nodes_and_mus, by = "node") 
  return(pred)
}
