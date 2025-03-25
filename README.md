# Bee Swarm Optimization (BSO) for Model Specification Search

Bee Swarm Optimization (BSO) is a metaheuristic algorithm developed for model specification search in psychometric test development. Inspired by the foraging behavior of bees, BSO simultaneously explores the optimal factor structure and performs item selection. In the algorithm, the labor is splitted to two types of bees:
 
- **Scout Bees** perform a global search determening the overall structure;  they perform major modifications to the factor structure (adding, splitting, merging or removing factors).
- **Onlooker Bees** perform a more localized and targeted search; they perform minor more localized adjustments to the item-factor assignment (adding, reassigning or removing items from a factor).

This dual search strategy helps prevent premature convergence on local optima and allows the simultaneous optimization of multiple fit criteria (e.g., global and local model fit, factor loading, maximum factor correlations, reliability). For a detailed explanation of the algorithm and its simulation performance, please refer to the associated publication (see the corresponding reference).

## Algorithm Overview

The BSO algorithm is controlled by several hyperparameters (1-3 are the key hyperparameters):

1. **Number of Bees (`n_bees`)**: Total candidate solutions evaluated per iteration.
2. **Percentage of Scouts (`percent_scouts`)**: Proportion of scout bees that perform major modifications in proportion to the number of onlooker bees.
3. **Top Solutions (`top_best`)**: Number of best-performing solutions selected for further refinement (those will be modified by scouts & onlookers in the next iteration).
4. **Factor Range Constraints (`min_nest_fac` and `max_nest_fac`)**: Minimum and maximum number of correlated factors allowed in the model (search space restriction).
5. **Maximum Iterations (`max_iter`)**: Maximum number of iterations without improvement before termination.
6. **Depletion (`depletion`)**: Number of iterations after which a solution is considered exhausted and is updated (top solutions will only be explored `depletion` iterations).
7. **Number of Start Bees (`n_start_bees`)**: The Number of solutions evaluated in the first iteration. It´s often a good idea to set it as a multiple of n_bees especially with convergence issues to improve the inital exploration.   

The algorithm evaluates candidate models using several fit criteria (e.g., CFI, RMSEA, maximum factor correlation, item retention, and local fit critera based on the resiudals). Each criterion is transformed to a common scale via logistic functions; these transformed values are then aggregated into an overall fit score (the “nectar” value).

## Special Features

- **Simultaneous Optimization:** BSO integrates the search for the optimal factor structure while simultaneously performing item selection in one process.
- **Global and Local Search:** The distinct roles of scout and onlooker bees facilitate both wide-ranging exploration and targeted refinement.
- **Multiple Fit Criteria:** The algorithm balances various fit indices (e.g. global fit, local fit, reliability, and parsimony) to yield a robust and interpretable model.


## Usage

1. **Prepare Your Data:** Load your dataset and modify it if necessary (e.g. standardization needed?).
2. **Define the Item Pool:** Create a vector of item names, only items contained in this vector are considered for the Model Specification Search.
3. **Set Hyperparameters:** Adjust the hyperparameters as needed (we would recommend using at least n_bees = 200, with percent_scouts = .25 and top_solutions = 25% of n_bees; for a decent compromise between efficiency and performance).
4. **Run BSO:** Call the `BSO` function with the specified parameters.
5. **Review Outputs:** Results are saved to CSV files and, if enabled, a convergence plot is generated.

## Example

The following example uses the Holzinger & Swineford dataset from the `OpenMx` package to run the BSO algorithm:

```r
# Clear the workspace
rm(list = ls())

# Load the Bee Swarm Optimization main function
source("BSO_main_function.R")

# Load the Holzinger & Swineford dataset from OpenMx
data("HS.ability.data", package = "OpenMx")
HS_data <- HS.ability.data

# Extract item names from columns 7 to 30
item_names <- colnames(HS_data[, 7:30])

# Standardize the selected items (may not be necessary depending on your data)
HS_data <- as.data.frame(apply(HS_data[, item_names], 2, scale))

# Set hyperparameters for the BSO algorithm
max_iter <- 10             # Maximum iterations without improvement
n_bees <- 250              # Total number of bees (candidate solutions) per iteration
percent_scouts <- 0.25     # Proportion of bees designated as scouts (global modifications)
top_best <- 10             # Number of top solutions selected for further refinement (this is an absolute number here!)
min_nest_fac <- 1          # Minimum number of correlated factors allowed
max_nest_fac <- 5          # Maximum number of correlated factors allowed
depletion <- 5             # Iterations after which a solution is considered exhausted
cluster_mode <- FALSE      # Set to TRUE when running on a computing cluster

# Define file paths for saving output results
res.name <- paste0("BSO_Holz_s1_b", n_bees, "_s", percent_scouts * n_bees, "_t", top_best, Sys.Date())
summaryfile <- paste0(res.name, ".csv") # File to save all solutions to (the entire convergence process)
summaryfile_fin <- paste0(res.name, "_final.csv")	# File to save the best solution to

# Run the BSO algorithm
conv_plot <- BSO(
  item_names = item_names, # Candidate items for model specification
  data = HS_data, # Your dataset
  max_iter = max_iter, # Stopping Criterion: Maximum iterations without improvement before termination 
  n_bees = n_bees, # Number of bees (candidate solutions) per iteration
  
  # Hyperparameters:
  percent_scouts = percent_scouts, # Proportion of bees designated as scouts (global modifications)
  top_best = top_best, 	# Number of top solutions selected for further refinement
  min_nest_fac = min_nest_fac, 	 # Minimum number of correlated factors allowed
  max_nest_fac = max_nest_fac, 	 # Maximum number of correlated factors allowed
  depletion = depletion, 	# Iterations after which a solution is considered exhausted
  
  # Output settings:
  summaryfile = summaryfile,	# File to save all solutions to (the entire convergence process)
  summaryfile_fin = summaryfile_fin,# File to save the best solution to
  
  seed = 1,                  # Random seed for reproducibility
  verbose = FALSE,           # Set TRUE for detailed output
  parallel = TRUE,           # Enable parallel processing
  nCores = 16,               # Number of CPU cores for parallel processing
  plot_nectar = TRUE,        # Enable convergence plot
  
  # Fit criteria and their logistic transformation parameters
  fit_crit = c("cfi", "rmsea", "max_fac_cor1", "n_items", "residual"),
  logistic_weights = list(
    c(d = 0.90, a = 30),     # CFI: Higher values are better
    c(d = 0.06, a = -50),    # RMSEA: Lower values are better
    c(d = 0.35, a = 15),     # Maximum Factor Correlation: Penalizes high factor correlations
    c(d = 0.95, a = 50),     # Number of Items: Balances item retention 
    c(d = NA, a = NA)        # Residual: Local fit criterion that penalizes the highest 10% of values in the standardized residual covariance matrix -> No logistic transformation is applied, as this is already done in the fit_function
  ),
  nu_weights = c(0.5, 0.5, 1, 1, 1),
  nu_min = 1e-10,
  
  # Additional settings:
  balance_n_fac = FALSE, # Balances the number of factors considered at the start
  n_start_bees = n_bees * 3, # More starting solutions to improve initial exploration
  n_start_items = 24,	# Number of Items as starting point (all items vs. a selcted number of items)
  bounds = "pos.var",	# Bounds to prevent convergence & run time issues
  cluster_mode = cluster_mode, # Set to TRUE when running on a computing cluster (prevents plotting after each iteration)
  plot_list = list(	# ggplot arguments for the convergence plot 
    xlim = c(0, max_iter * 10),
    ylim = c(0, sum(c(1, 1, 0.5, 1, 1, 1))),
    ylab = "Overall Nectar Value",
    xlab = "Iteration",
    jitter_width = 0.5,
    alpha = 0.2,
    size = 1
  )
)

# Save the convergence plot as a PDF (if created)
if (!is.null(conv_plot)) {
  suppressWarnings({
    ggsave(plot = conv_plot, filename = paste0(res.name, ".pdf"), device = "pdf")
  })
}
