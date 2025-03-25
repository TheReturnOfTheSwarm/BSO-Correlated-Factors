# BSOs Main-Function

# Install packages if required 
install_and_load <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Load all required packages
install_and_load("psych")        # For psychometric functions and factor analysis
install_and_load("lavaan")       # For structural equation modeling 
install_and_load("parallel")     # For parallel computing capabilities
install_and_load("GPArotation")  
install_and_load("OpenMx")       # For extended structural equation modeling (and the Holzinger-Dataset for Testing Purposes)
install_and_load("semTools")     # For additional SEM tools and utilities
install_and_load("ggplot2")      # For visualization
install_and_load("scales")       # For scale transformations and formatting
install_and_load("tidyverse")    # For data manipulation and visualization

# Load necessary helper functions from external files
source("tool_functions.R")  # Contains subfunctions for BSO operations
source("fit_function.R")    # Contains the fitness evaluation function (may need individual changes for optimization functions which aren´t included yet)

########################## Main BSO Function ################################################
#' @title Bee Swarm Optimization (BSO) for Correlated Factor Analysis
#' @description 
#' Implements the Bee Swarm Optimization algorithm for discovering optimal factor structures
#' in psychometric data. The algorithm mimics the behavior of bee colonies, where scout bees
#' explore new solutions (major model changes) and onlooker bees refine existing solutions
#' (minor model changes). The BSO approach avoids local optima through its exploratory nature.
#'
#' @param item_names Vector of item names to be used in the factor model
#' @param data Data frame containing the variables
#' @param max_iter Maximum number of iterations without improvement before stopping
#' @param n_bees Total number of bees (Number of Solutions to explore in one Iteration) to use in the optimization
#' @param n_start_bees Number of random starting models to generate (defaults to n_bees)
#' @param n_start_items Number of items to include in the initial models (defaults to all items)
#' @param percent_scouts Proportion of bees assigned as scouts (major model changes)
#' @param top_best Number of best solutions to focus modification efforts on
#' @param min_nest_fac Minimum number of correlated factors allowed
#' @param max_nest_fac Maximum number of correlated factors allowed
#' @param depletion Maximum "age" of a solution before it's considered depleted (exhausted)
#' @param summaryfile File path to save detailed results of all explored models
#' @param summaryfile_fin File path to save only the best final model
#' @param write_solutions Logical; if TRUE, write solutions to files (default is TRUE)
#' @param scouts Vector specifying number of scout bees per top solution (calculated from percent_scouts)
#' @param onlookers Vector specifying number of onlooker bees per top solution (calculated from percent_scouts)
#' @param fit_crit Vector of fit criteria to extract and optimize
#' @param logistic_weights List of parameters for logistic transformation of fit criteria
#' @param nu_weights Weights for the aggregation function that combines fit criteria
#' @param nu_min Threshold value for minimum acceptability of individual fit criteria
#' @param n_items_min Minimum number of items required in a factor
#' @param verbose Logical; if TRUE, print detailed progress information
#' @param debug_fit_mode Logical; if TRUE, return fit object for debugging
#' @param parallel Logical; if TRUE, use parallel processing
#' @param nCores Number of CPU cores to use when parallel=TRUE
#' @param seed Random seed for reproducibility (in case of parallel processing the number of cores must be the same as before for reproducibility)
#' @param fun Logical; if TRUE, display ASCII art (Easter egg feature)
#' @param plot_nectar Logical; if TRUE, create plots of optimization progress
#' @param plot_list List of aesthetic parameters for the nectar plot
#' @param cluster_mode Logical; if TRUE, disable plotting (for use on clusters)
#' @param balance_n_fac Logical; if TRUE, start with uniform factor distribution
#' @param ... Additional arguments passed to lavaan
#'
#' @return A ggplot object showing optimization progress (if plot_nectar=TRUE)
#' @details
#' The BSO algorithm works as follows:
#' 1. Generate initial random factor models as starting points
#' 2. Evaluate each model using fit criteria
#' 3. Select top solutions for further exploration
#' 4. Apply scout bees to make major changes (add/remove/merge/split factors)
#' 5. Apply onlooker bees to make minor changes (add/remove/swap items)
#' 6. Re-evaluate new models and keep the best solutions
#' 7. Continue as long as a better model is found within max_iter iterations
#'
#' The "correlated factors" terminology refers to factor models where factors are allowed to
#' correlate with each other; nest is used for correlated factors for consistency with the original function.
BSO = function(item_names,      # Vector of item names to be analyzed
               data,            # Data frame containing the data set
               max_iter,        # Maximum iterations without improvement before stopping
               n_bees,          # Total number of bees (computational agents)
               n_start_bees = n_bees,  # Number of random starting models (defaults to n_bees)
               n_start_items = length(item_names),  # Number of items in starting models
               percent_scouts,  # Proportion of scout bees (major model changes)
               top_best,        # Number of best solutions to focus on
               min_nest_fac,    # Minimum number of correlated factors allowed
               max_nest_fac,    # Maximum number of correlated factors allowed
               depletion,       # Maximum "age" of a solution before it's considered exhausted and won´t be longer used as top solution
               summaryfile,     # File path for detailed results output
               summaryfile_fin, # File path for final best model output
               write_solutions = TRUE,  # Whether to write solutions to files
               # Calculate scout and onlooker distribution based on solution quality (weighted by rank)
               scouts = round(prop.table(top_best:1) * n_bees * percent_scouts),     # More scouts for better solutions
               onlookers = round(prop.table(top_best:1) * n_bees * (1 - percent_scouts)),  # More onlookers for better solutions
               fit_crit = c("cfi", "rmsea", "min_omega2", "min_loading"),  # Fit criteria to optimize (for other implemented criteria see fit_function)
               # Parameters for logistic transformation of fit indices (shape and scale parameters)
               logistic_weights = list(c(d = 0.9, a = 70),    # CFI transformation (higher is better)
                                       c(d = 0.06, a = -70),  # RMSEA transformation (lower is better)
                                       c(d = 0.4, a = 70),    # min_omega2 transformation (higher is better)
                                       c(d = 0.33, a = 70)),  # min_loading transformation (higher is better)
               nu_weights = c(1,1,1,1),  # Equal weights for combining fit criteria by default
               nu_min = 0,               # Minimum threshold for any single fit criterion
               n_items_min = 3,          # Minimum items per factor
               verbose = FALSE,          # Whether to print detailed information (not compatible with parallel mode)
               debug_fit_mode = FALSE,   # Whether to return full fit object for debugging
               parallel = TRUE,          # Whether to use parallel processing
               nCores = parallel::detectCores() - 2,  # Number of cores to use (leaves 2 for system)
               seed = 1,                 # Random seed for reproducibility
               fun = TRUE,               # Easter egg: display ASCII art
               plot_nectar = TRUE,       # Whether to create convergence plot
               # Aesthetic parameters for the nectar (fitness) plot
               plot_list = list(xlim = c(0,max_iter*2), 
                                ylim = c(0,sum(nu_weights)),
                                ylab = "Overall Nectar Value",
                                xlab = "Iteration",
                                jitter_width = 0.5,
                                alpha = 0.2,
                                size = 1,
                                width = 10,
                                height = 8,
                                dpi = 450),
               cluster_mode = FALSE,     # For use on computing clusters (disables plots)
               balance_n_fac = FALSE,    # Whether to balance factor counts in initial models
               ...){                     # Additional arguments passed to lavaan
  
  # Display ASCII art if fun=TRUE (Easter egg feature)
  if(fun){
    cat(readLines("fun1.txt", warn = F), sep = "\n")
  }
  
  # Initialize the nectar (fitness) plot if enabled
  if (plot_nectar){
    # Create base plot with proper axis limits and labels
    conv_plot <- ggplot() + 
      xlim(plot_list$xlim[1], plot_list$xlim[2]) +
      ylim(plot_list$ylim[1], plot_list$ylim[2]) +
      xlab(plot_list$xlab) +
      ylab(plot_list$ylab) 
    # Display empty plot at start (unless in cluster mode)
    if(!cluster_mode) plot(conv_plot)
  }
  
  # Validate bee distribution: total scouts + onlookers must equal n_bees
  if (sum(scouts + onlookers) != n_bees) {
    # Adjust n_bees if parameters don't match and warn user
    n_bees <- sum(scouts + onlookers)
    warning(paste0("The current combination of n_bees, percent_scouts and top_best does not work out, I changed n_bees to "), n_bees)
  }
  
  # Set up parallel processing environment if enabled
  if (parallel) {
    # Create PSOCK cluster for Windows compatibility
    myCL <- makePSOCKcluster(names = nCores, outfile = ifelse(verbose, paste0("parallel_verbose",Sys.time(),".txt"), "")) # the same number of cores is necessary to reproduce the same results 
    
    # Set up random number generation for reproducibility in parallel mode
    RNGkind("L'Ecuyer-CMRG") # Use a parallel-safe RNG
    clusterSetRNGStream(cl = myCL, iseed = seed) 
    
    # Export all variables and functions to the cluster nodes
    clusterExport(cl = myCL, varlist = ls(envir = environment()), envir = environment())
    clusterExport(cl = myCL, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
    
    # Load necessary packages on each cluster node
    clusterEvalQ(cl = myCL, expr = {
      library(psych)
      library(lavaan)
      library(semTools)
      library(ggplot2)
      library(scales)
    })
  } else {
    # For single-core mode, just set the seed
    set.seed(seed, kind = "L'Ecuyer-CMRG")  # Use same RNG as parallel for consistency
  }
  
  # Store total number of items
  n_items <- length(item_names)
  
  ########################### Phase 1: Generate Initial Solutions ##############################
  # Generate n_start_bees random starting models using initialize.bees function
  
  if (!parallel){
    # Single-core implementation
    solutions <- sapply(1:n_start_bees, function(i_scout){ 
      if (verbose) print(i_scout)
      
      initialize.bees(item_names = item_names,
                      data = data,
                      max_iter = max_iter,
                      n_bees = n_bees,
                      n_start_items = n_start_items,
                      fit_crit = fit_crit,
                      logistic_weights = logistic_weights,
                      nu_weights = nu_weights,
                      nu_min = nu_min,
                      balance_n_fac = balance_n_fac,
                      debug_fit_mode = debug_fit_mode,
                      verbose = verbose,
                      ...)
    }, simplify = "matrix")
  } else if (parallel){
    # Parallel implementation
    solutions <- parSapply(myCL, 1:n_start_bees, function(i_scout){ 
      if (verbose) print(i_scout)
      initialize.bees(item_names = item_names,
                      data = data,
                      max_iter = max_iter,
                      n_bees = n_bees,
                      n_start_items = n_start_items,
                      fit_crit = fit_crit,
                      logistic_weights = logistic_weights,
                      nu_weights = nu_weights,
                      nu_min = nu_min,
                      balance_n_fac = balance_n_fac,
                      debug_fit_mode = debug_fit_mode,
                      verbose = verbose,
                      ...
      )
    }, simplify = "matrix")
  }
  
  # Transform solutions to a more convenient data frame format
  solutions <- data.frame(t(solutions))
  
  # Sort solutions by overall fitness (nectar value) in descending order
  solutions <- solutions[order(solutions$fit_overall, decreasing = TRUE),]
  
  ################################################ Phase 2: Main Optimization Loop ##############################################
  # In this phase, scout bees perform major modifications and onlooker bees perform minor modifications
  # to the top solutions in each iteration until convergence
  
  iter <- 1        # Overall iteration counter
  counter <- 1     # Counter for iterations without improvement (resets when improvement found)
  best_fit <- solutions$fit_overall[1]  # Track best fitness seen so far
  best_solution <- solutions[1,]        # Store best solution seen so far
  
  # Main optimization loop
  while (counter <= max_iter){ # Continue until max_iter iterations without improvement
    start <- Sys.time() # Measure runtime
    
    # Select top solutions to explore in this iteration
    tmp_solutions <- solutions[1:top_best,]
    
    # Check which solutions are still "fresh" (not depleted)
    fresh <- solutions$age <= depletion
    
    # Increase age of best/top solutions that are still fresh
    solutions$age[fresh][1:top_best] <- solutions$age[fresh][1:top_best] + 1
    
    # Prepare a data frame that defines which bees will modify which solutions
    allBees <- data.frame(i_bee_ind = 1:n_bees)
    # Assign bees to solutions based on scouts and onlookers vectors
    allBees$solution_focus = c(rep(c(1:top_best), times = scouts), rep(c(1:top_best), times = onlookers)) 
    # Assign bee roles - scouts first, then onlookers
    allBees$role = c(rep("scout", sum(scouts)), rep("onlooker", sum(onlookers)))
    # Add the current model parameters for each solution
    allBees <- data.frame(allBees, tmp_solutions[allBees$solution_focus,item_names]) 
    
    # Process each bee's operation (either in serial or parallel)
    if (!parallel){
      # Single-core implementation
      new_solutions <- sapply(1:nrow(allBees), function(i_bee){
        
        if (verbose) print(allBees[i_bee,])
        
        # Get the current item assignment for this bee to modify
        item_assignment_orig <- allBees[i_bee, item_names]
        
        # Apply appropriate bee operation based on the bee's role
        if (allBees$role[i_bee] == "scout"){
          # Scout bees make major model changes (add/remove/merge/split factors)
          item_assignment <- scout.bee(item_assignment = item_assignment_orig,
                                       max_nest_fac = max_nest_fac)
          # Ensure consistent factor numbering
          item_assignment <- reorder_factors(item_assignment = item_assignment)
          
        } else if (allBees$role[i_bee] == "onlooker") {
          # Onlooker bees make minor model changes (add/remove/swap items)
          item_assignment <- onlooker.bee(item_assignment = item_assignment_orig, 
                                          max_nest_fac = max_nest_fac)
          
          # Ensure consistent factor numbering
          item_assignment <- reorder_factors(item_assignment = item_assignment)
        } else {
          stop("Fatal error: Bee swarm out of control.") # If an unknown role occurs
        }
        
        # Check if this model configuration has been tested before
        compare_item_assignment <- apply(solutions[, item_names], 1, function(i_solution) 
          identical(item_assignment, i_solution))
        
        if (any(compare_item_assignment)){ 
          # If this model was seen before, retrieve its existing evaluation
          if (verbose) print ("Empty Flower")
          new_sol <- solutions[max(which(compare_item_assignment)),]
          # Increase age of retrieved solution
          new_sol$age <- max(solutions$age[which(compare_item_assignment)]) + 1 
        } else { 
          # If this is a new model, evaluate it
          
          # Calculate fit metrics for this model
          fit <- fit.function(item_assignment = item_assignment, 
                              item_names = item_names,
                              dat = data,
                              verbose = verbose,
                              fit_crit = fit_crit,
                              logistic_weights = logistic_weights,
                              nu_weights = nu_weights,
                              nu_min = nu_min,
                              debug_fit_mode = debug_fit_mode,
                              ...) # evaluate model
          
          
          # Save model results: item-factor allocation, seed, iteration, 
          # number of factors, age, and fit statistics
          n_nest_fac <- max(item_assignment)
          
          new_sol <- c(unlist(item_assignment), 
                       seed = seed, 
                       iteration = iter, 
                       n_nest_fac = n_nest_fac, 
                       age = 0,  # New solutions start with age 0
                       fit)
          new_sol
        }
        # Return the solution (either new or existing)
        new_sol
      }, simplify = "matrix")
      
    } else if (parallel){
      # Parallel implementation - similar to single-core but distributed across nodes
      
      # Export updated environment variables to all cluster nodes
      clusterExport(cl = myCL, varlist = ls(envir = environment()), envir = environment())
      clusterExport(cl = myCL, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
      
      # Process each bee's operation in parallel
      new_solutions <- parSapply(myCL, 1:nrow(allBees), function(i_bee){
        
        if (verbose) print(allBees[i_bee,])
        
        # Get the current item assignment for this bee to modify
        item_assignment_orig <- allBees[i_bee, item_names]
        
        # Apply appropriate bee operation based on the bee's role
        if (allBees$role[i_bee] == "scout"){
          # Scout bees make major model changes
          item_assignment <- scout.bee(item_assignment = item_assignment_orig,
                                       max_nest_fac = max_nest_fac)
          
          # Ensure consistent factor numbering
          item_assignment <- reorder_factors(item_assignment = item_assignment)
          
        } else if (allBees$role[i_bee] == "onlooker") {
          # Onlooker bees make minor model changes
          item_assignment <- onlooker.bee(item_assignment = item_assignment_orig, 
                                          max_nest_fac = max_nest_fac)
          
          # Ensure consistent factor numbering
          item_assignment <- reorder_factors(item_assignment = item_assignment)
        } else {
          stop("Fatal error: Bee swarm out of control.") # If an unknown role occurs
        }
        
        # Check if this model configuration has been tested before
        compare_item_assignment <- apply(solutions[, item_names], 1, function(i_solution) 
          identical(item_assignment, i_solution))
        
        if (any(compare_item_assignment)){ 
          # If this model was seen before, retrieve its existing evaluation
          if (verbose) print ("Empty Flower")
          new_sol <- solutions[max(which(compare_item_assignment)),]
          new_sol$age <- new_sol$age + 1 
        } else { 
          # If this is a new model, evaluate it
          
          # Calculate fit metrics for this model
          fit <- fit.function(item_assignment = item_assignment, 
                              item_names = item_names,
                              dat = data,
                              verbose = verbose,
                              fit_crit = fit_crit,
                              logistic_weights = logistic_weights,
                              nu_weights = nu_weights,
                              nu_min = nu_min,
                              debug_fit_mode = debug_fit_mode,
                              ...
          ) # evaluate model
          
          
          # Store model information and results
          n_nest_fac <- max(item_assignment)
          
          
          new_sol <- c(unlist(item_assignment), 
                       seed = seed, 
                       iteration = iter, 
                       n_nest_fac = n_nest_fac, 
                       age = 0, 
                       fit)
          
          new_sol
        }
      }, simplify = "matrix")
    }
    
    # Transform results to a data frame for easier processing
    new_solutions <- data.frame(t(new_solutions))
    
    
    # Update visualization if plotting is enabled
    if (plot_nectar) {
      # Create a copy of solutions for plotting
      new_solutions2 <- new_solutions
      # Add information needed for the plot
      new_solutions2$iter <- iter 
      new_solutions2$role <- as.factor(allBees$role)
      new_solutions2$n_nest_fac <- as.factor(new_solutions2$n_nest_fac)
      # Update plot and suppress NA-warnings (e.g., when fit is not converged)
      suppressWarnings(expr = {
        conv_plot <- conv_plot + 
          geom_jitter(data = new_solutions2, 
                      aes(x = iter,
                          y = fit_overall,
                          color = n_nest_fac,
                          shape = role),
                      width = plot_list$jitter_width,
                      alpha = plot_list$alpha,
                      size = plot_list$size) 
        if(!cluster_mode) plot(conv_plot)
      }) 
      
    }
    
    
    # Verify that at least some models have valid fit values
    if (all(is.na(new_solutions$fit_overall))) {
      stop("All initial fit values are NA. Please check fit-function. Try setting verbose = TRUE or debug_fit_mode = TRUE.")
    }
    
    # Sort new solutions by fit quality
    new_solutions <- new_solutions[order(new_solutions$fit_overall, decreasing = TRUE),]
    
    # Check if we've found a new best solution
    if (max(new_solutions$fit_overall, na.rm = TRUE) > best_fit) { 
      # Update best solution information
      best_fit <- max(new_solutions$fit_overall, na.rm = TRUE)
      best_solution <- new_solutions[1,]
      best_solution$nCores <- nCores # Record number of cores used (Only necessary if cluster with varying core numbers are used)
      counter <- 0 # Reset counter because we found improvement
    } 
    
    
    # Combine old and new solutions, then sort by quality
    solutions <- rbind(solutions, new_solutions)
    solutions <- solutions[order(solutions$fit_overall, decreasing = T),]
    
    # Update iteration counters
    counter <- counter + 1 # Count iterations without improvement
    iter <- iter + 1       # Count total iterations
    end <- Sys.time()      # Record end time for this iteration
    
  } # End of main optimization loop
  
  # Clean up: stop cluster if parallel processing was used
  if (parallel) stopCluster(myCL)
  
  # Save results to files if requested
  if (write_solutions) {
    # Save all explored solutions
    write.table(solutions, summaryfile, sep=";", row.names = F, col.names = T, append = F)
  }
  
  # Capture session information for reproducibility
  session_info <- capture.output(sessionInfo())
  
  # Write session information and best solution to final summary file
  writeLines(c("# Session Information:", session_info, "# End of Session Information", ""), con = summaryfile_fin)
  write.table(best_solution, summaryfile_fin, sep=";", row.names = F, col.names = T, append = F)
  
  # Return convergence plot if created
  if(plot_nectar) {
    conv_plot + xlim(0, iter)
  } else {
    NULL
  }
  
}


############################################# End BSO #################################################
