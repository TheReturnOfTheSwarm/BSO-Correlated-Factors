################## BSO subfunctions #############################
#' @title Bee Swarm Optimization (BSO) Helper Functions
#' @description This file contains a collection of sub-functions which are utilized during BSO algorithm execution.
#' These functions handle model initialization, scout bee operations, onlooker bee operations, and various helper utilities.

#' @title Initialize Bee Solutions
#' @description Generates the random start models in the first iteration of the BSO algorithm.
#' 
#' @param item_names Vector of item names.
#' @param data Data frame containing the data set.
#' @param max_iter Maximum number of iterations.
#' @param n_bees Number of bees.
#' @param n_items Number of items (inferred from name vector by default).
#' @param verbose Logical; if TRUE, print detailed information.
#' @param debug_fit_mode Logical; if TRUE, return fit object for debugging.
#' @param fit_crit Vector of fit criteria to extract.
#' @param logistic_weights Weights for logistic transformation.
#' @param nu_weights Weights for the aggregation function.
#' @param nu_min Threshold value for fit criteria.
#' @param n_start_items How many items should be in the start model (allows starting with a subset of items).
#' @param balance_n_fac Logical; if TRUE, start models are sampled such that the number of factors is uniformly distributed, if FALSE, random item assignments are used resulting in more candidate models with fewer factors 
#' @param ... Additional arguments passed to lavaan.
#'
#' @return A vector containing item assignments, seed, iteration number, factor count, age, and fit statistics.
#'
#' @details
#' Item assignment is coded as follows:
#' - -1: Item is not in the model
#' -  1: max_nest_fac: Item loads on this correlated factor (nest is the name for consistency with the original function)
#'
#' The function ensures that each factor has at least 3 items.
initialize.bees <- function(item_names, # a vector of item names (length must equal the length of item assignments)
                            data, # a data frame containing the data set
                            max_iter, # maximum number of iterations
                            n_bees, # number of bees
                            n_items = length(item_names), # number of items (inferred from name vector)
                            verbose = FALSE, # verbose mode (passed to fit function)
                            debug_fit_mode = FALSE, # debug mode (passed to fit function)
                            fit_crit = c("cfi", "rmsea", "min_omega2", "min_loading"), # see fit-function()
                            logistic_weights, # see fit-function()
                            nu_weights, # see fit-function()
                            nu_min, # see fit-function()
                            n_start_items = length(item_names), # How many items should be in the start model? (allows start with a subselection of items)
                            balance_n_fac = FALSE, # if TRUE, the start models are sampled such that the number of factors is uniformly distributed
                            ...){ # further arguments to be passed to lavaan
  

#####
#  Function start for the initialization of the bees
#####
  # item_assignment is a vector with numbers with the following coding scheme:
  # -1                ... Item is not in the model
  
  # 1:max_nest_fac    ... Item loads on this correlated factor (nest is the name for consistency with the original function)
  
  
  if (!balance_n_fac){ 
    # If not balancing factor counts, create random assignment
    # if we start with a subset of the items, compute how many are excluded
    n_items_excluded <- n_items - n_start_items
    
    # return error for impossible argument combinations
    if(n_items_excluded < 0 | (n_items - n_items_excluded) < 0) stop("Invalid value for n_start_items")
    
    # assign random item allocations
    item_assignment <- rep(-1, n_items_excluded)
    item_assignment <- c(item_assignment, resamp(1:max_nest_fac, n_items - n_items_excluded, replace = TRUE)) # generate models
    
    # Note that across all start models this procedure will result in a higher number of models
    # with fewer factors (because these are more frequent among all possible random assignments)

  } else {
    
    # This approach tries to balance the number of factors across all start models. 
    
    # draw a random number of factors first
    n_fac_tmp <- sample(min_nest_fac:max_nest_fac, size = 1)
    item_assignment <- c(rep(1:n_fac_tmp, each = 3)) # place at least 3 items for each factor
    
    # the other items can be assigned to any of the factors 
    # items_in_model <- n_items - length()
    
      n_items_excluded <- n_items - n_start_items
      n_items_remaining <- n_items - length(item_assignment) - n_items_excluded
      # return error for impossible argument combinations
      if(n_items_excluded < 0 | n_items_remaining < 0) stop("Invalid value for n_start_items")
      item_assignment <- c(item_assignment, rep(-1, n_items_excluded), sample(1:n_fac_tmp, size = n_items_remaining, replace = TRUE))
      
    # randomize order
    item_assignment <- sample(item_assignment, size = length(item_assignment))
  }
  
 
  # remove correlated factors with less than 3 items
  item_counts <- table(item_assignment[item_assignment %in% 1:max_nest_fac])
  fac_2b_removed <- as.numeric(names(item_counts)[item_counts < 3])
  # remove the factor by removing these items 
  item_assignment[item_assignment %in% fac_2b_removed] <- -1 #Items will be removed
  
  # example: it is possible that a factor 5 occurs but factor 2 was removed
  # in that case, the factors should be relabeled 1:4
  item_assignment <- renew_factor_numbers(item_assignment = item_assignment)
  
  # Reorder Factors here
  item_assignment <-reorder_factors(item_assignment = item_assignment)

  
  # Name assignment so that the names are exported
  names(item_assignment) <- item_names
  
  # n_nest_fac is the max number of items (normally less due to min requirement of 3 items) 
  n_nest_fac <- max(item_assignment)
  
  # fit the the model
  fit <- fit.function(item_assignment = item_assignment, 
                      item_names = item_names,
                      dat = data,
                      verbose = verbose,
                      fit_crit = fit_crit,
                      logistic_weights = logistic_weights,
                      nu_weights = nu_weights,
                      nu_min = nu_min,
                      debug_fit_mode = debug_fit_mode,
                      ...) #evaluate model
  
  # save items - factor allocation, seed, iteration, number of factors, "age" and quality in solution object
  new_sol <- c(item_assignment, 
               seed = seed, 
               iteration = 0, 
               n_nest_fac = n_nest_fac, 
               age = 0, 
               fit)
  new_sol 
}


#' @title Scout Bee Function
#' @description Performs major model modifications as a scout bee in the BSO algorithm.
#' 
#' @param item_assignment Current item assignment vector.
#' @param max_nest_fac Maximum number of correlated factors allowed.
#'
#' @return A modified item assignment vector after applying a major operation.
#'
#' @details Scout bees perform major operations on models:
#' 1. Add a new factor
#' 2. Split an existing factor
#' 3. Remove a factor
#' 4. Merge two factors
#' 
#' The function randomly selects one of these operations based on what's possible 
#' with the current model configuration.
scout.bee <- function(item_assignment, max_nest_fac){
  n_items_per_fac <- table(factor(item_assignment, levels = 1:max_nest_fac))
  n_nest_fac   <- max(item_assignment)
  
  # Check which major changes are allowed and choose one randomly 
  allowed_operations <- check.operations.scouts(item_assignment = item_assignment,
                                                max_nest_fac = max_nest_fac)

  chosen_operation <- resamp(allowed_operations, size = 1) 

  switch (chosen_operation,
          { # add new factor
            free_items <- which(item_assignment <= 0) # Which items are available?
            n_new_items <- resamp(3:length(free_items), size = 1) # How many items should the new factor have?
            new_fac_items <- resamp(free_items, size = n_new_items) # Choose items randomly
            item_assignment[new_fac_items] <- n_nest_fac + 1}, # Assign these items to the new factor
          { #split factors
            big_fac <- which(n_items_per_fac >= 6) # Which factors could be split?
            split_fac <- resamp(big_fac, size = 1) # Decide which factor should be split
            split_fac_items <- which(item_assignment == split_fac) # Which items belong to the factors?
            n_new_fac_items <- resamp(3:(length(split_fac_items) - 3), size = 1) # Split the factor randomly such that each factor has at least 3 items
            new_fac_items <- resamp(split_fac_items, size = n_new_fac_items) # Sample items for new factor
            
            item_assignment[new_fac_items]<- n_nest_fac + 1},# Assign these items to the new factor
          { #remove factor
            old_fac <- resamp(1:n_nest_fac, size = 1)
            item_assignment[item_assignment == old_fac] <- -1 }, # Remove the factor including Items
          { #merge factors
            old_fac <- resamp(1:n_nest_fac, size = 2)
            item_assignment[item_assignment == old_fac[1]] <- old_fac[2] }
  )
  
  # renew factor numbers again
  item_assignment <- renew_factor_numbers(item_assignment = item_assignment)   
  item_assignment
  

}

#' @title Onlooker Bee Function
#' @description Performs minor model modifications as an onlooker bee in the BSO algorithm.
#' 
#' @param item_assignment Current item assignment vector.
#' @param max_nest_fac Maximum number of correlated factors allowed.
#'
#' @return A modified item assignment vector after applying a minor operation.
#'
#' @details Onlooker bees perform minor operations on models:
#' 1. Add an item to a factor
#' 2. Remove an item from a factor
#' 3. Swap items between factors
#' 4. Delete an item from the item pool (redundant for correlated factors with operation 2)
#' 
#' The function randomly selects one of these operations based on what's possible 
#' with the current model configuration.
onlooker.bee <- function(item_assignment, max_nest_fac){
  n_items_per_fac <- table(factor(item_assignment,levels=1:max_nest_fac))
  n_nest_fac   <- max(item_assignment)

  # Check which minor changes are allowed and choose one randomly
  allowed_operations <- check.operations.onlookers(item_assignment = item_assignment,
                                                max_nest_fac = max_nest_fac)

  chosen_operation <- resamp(allowed_operations, size = 1) 

  switch (chosen_operation,
          { # add item to correlated factor
            which_item    <- resamp(which(item_assignment <= 0), size = 1) # look for items not assigned to a correlated factor
            to_which_factor   <- resamp(1:n_nest_fac, 1) # sample a new factor for this item
            item_assignment[which_item] <- to_which_factor}, # change assignment
          { # remove item from correlated factor
            which_factor   <- resamp(which(n_items_per_fac > 3), size = 1) # look for a factor with more than the minimal item count
            which_item   <- resamp(which(item_assignment == which_factor), size = 1) # sample any item from this factor
            item_assignment[which_item] <- -1}, # remove item from factor
          { # swap item between correlated factors
            item_1 <- resamp(which(item_assignment > 0), size = 1) # sample any item assigned to a correlated factor
            item_2 <- resamp(which((item_assignment > 0) & (item_assignment != item_assignment[[item_1]])), size = 1) # sample any item from another correlated factor
            item_assignment[c(item_1, item_2)] <- item_assignment[c(item_2, item_1)]},
          { # delete item from item pool #redundant with 2
            candidate_factors   <- c(which(n_items_per_fac > 3)) # look for a factor with more than the minimal item count 
            which_item   <- resamp(which(item_assignment %in% candidate_factors), size = 1) # sample any item from this factor 
            item_assignment[which_item] <- -1} # remove item from item pool completely
  )
  item_assignment # return
}

#' @title Scout Bee Operation Selection
#' @description Decides which scout operations are possible based on an item assignment and the number of factors.
#' 
#' @param item_assignment Current item assignment vector.
#' @param max_nest_fac Maximum number of correlated factors allowed.
#'
#' @return A vector of indices corresponding to allowed operations.
#'
#' @details Checks which of these major operations can be performed:
#' 1. Add a new factor
#' 2. Split a factor
#' 3. Remove a factor
#' 4. Merge factors
check.operations.scouts <- function(item_assignment, max_nest_fac){
  
  n_items_per_fac <- table(factor(item_assignment,levels=1:max_nest_fac))
  n_nest_fac   <- max(item_assignment)
  
  operations <- c(
    add_fac_allowed = (n_nest_fac < max_nest_fac) & (sum(item_assignment <= 0) >= 3),     #can scouts add factors?
    spl_fac_allowed = n_nest_fac < max_nest_fac & any(n_items_per_fac >= 6),  #can scouts split factors?
    rmv_fac_allowed = n_nest_fac > min_nest_fac,                      #can scouts remove factors?
    mer_fac_allowed = n_nest_fac > min_nest_fac                      #can scouts merge factors?
  )
   which(operations)
}

#' @title Onlooker Bee Operation Selection
#' @description Decides which onlooker operations are possible based on an item assignment and the number of factors.
#' 
#' @param item_assignment Current item assignment vector.
#' @param max_nest_fac Maximum number of correlated factors allowed.
#'
#' @return A vector of indices corresponding to allowed operations.
#'
#' @details Checks which of these minor operations can be performed:
#' 1. Add an item to a factor
#' 2. Remove an item from a factor
#' 3. Swap items between factors
check.operations.onlookers <- function(item_assignment, max_nest_fac){
  
  n_items_per_fac <- table(factor(item_assignment,levels=1:max_nest_fac))
  n_nest_fac   <- max(item_assignment)
  
  operations <- c(
    add_item_allowed = any(item_assignment <= 0),         # can onlookers add item-factor allocations?
    rmv_item_allowed = any(n_items_per_fac > 3),   # can onlookers remove item-correlated-factor allocations?
    swap_item_allowed = n_nest_fac > 1#,         # can onlookers swap item-correlated-factor allocations?
  # not considered because its redundant with 2
  #del_item_allowed = any(item_assignment != -1) # can onlookers delete items from the model?
  )
  
   which(operations)
}


#### Minor helper functions ####

#' @title Renew Factor Numbers
#' @description Ensures correct sequential factor numbering.
#' 
#' @param item_assignment Current item assignment vector.
#'
#' @return A vector with corrected sequential factor numbering.
#'
#' @details
#' Example: if a factor 5 occurs but factor 2 was removed,
#' the factors should be relabeled as 1:4 for consistency.
renew_factor_numbers <- function(item_assignment){
  # example: it is possible that a factor 5 occurs but factor 2 was removed
  # in that case, the factors should be relabeled 1:4
  # which items load on which correlated factors?
  items_nested <- item_assignment > 0
  item_assignment[items_nested] <- match(item_assignment[items_nested], sort(unique(item_assignment[items_nested])))
  item_assignment  
}

#' @title Robust Sample Function
#' @description Sample function that will also work correctly with vectors of length 1.
#' 
#' @param y Vector to sample from.
#' @param ... Additional arguments passed to sample().
#'
#' @return Sampled value(s) from the input vector.
resamp <- function(y,...){if(length(y)==1) y else sample(y,...)} 

#' @title Debug Helper
#' @description Debugging function to copy all temporary objects into the global environment.
#' 
#' @param tmp_env The environment to copy from (default is the current environment).
#'
#' @return Nothing, but copies all objects from the environment to the global environment.
allglobal <- function(tmp_env = environment()) {
  lss <- ls(envir = tmp_env)
  for (i in lss) {
    assign(i, get(i, envir = tmp_env), envir = globalenv())
  }
}

#' @title Reorder Factors
#' @description Reorders factor numbers based on their average positions in the item vector.
#' 
#' @param item_assignment Current item assignment vector.
#'
#' @return A vector with reordered factor numbering.
#'
#' @details
#' This improved function reorders factor numbers so that factors with items 
#' appearing earlier in the item vector tend to get lower factor numbers.
#' This creates more consistent labeling across multiple runs.
reorder_factors <- function(item_assignment) {
  # Remove -1 and other special values from consideration
  factors <- unique(item_assignment[item_assignment != -1])
  
  # Calculate average indices for each factor
  index_means <- sapply(factors, function(f) {
    mean(which(item_assignment == f))
  }, simplify = TRUE, USE.NAMES = TRUE)
  
  # Order factors by their average indices
  ordered_factors <- factors[order(index_means)]
  
  # Create a mapping from old factor numbers to new ordered factor numbers
  factor_mapping <- setNames(seq_along(ordered_factors), ordered_factors)
  
  # Apply the mapping to item_assignment to reorder factor labels
  item_assignment <- vapply(item_assignment, function(x) {
    if (x %in% names(factor_mapping)) {
      as.integer(factor_mapping[[as.character(x)]])
    } else {
      as.integer(x)  # Convert special values to integer explicitly
    }
  }, FUN.VALUE = integer(1))
  
  return(item_assignment)
}
