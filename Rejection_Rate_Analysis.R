rm(list = ls())
setwd("~/ETH/Convexity_test/Strict_convexity_test")
source("rejection_rate_parallel.R")


###########################################
## 1. Regression functions
f_list <- list(
  f0 = function(x){4*(x-0.5)^2},                       # super strictly convex
  f1 = function(x){(exp(x-1)-exp(-1))/(1-exp(-1))},    # slightly strictly convex
  f2 = function(x){x},                                 # convex (not strictly)
  f3 = function(x){-4*(x-0.5)^2},                      # super non convex
  f4 = function(x){1-(exp(x-1)-exp(-1))/(1-exp(-1))},  # slightly non convex
  f5 = function(x){((x^3-x^2)+(4/27))*27/4}            # small convexity violation
)

## 2. Error standard-deviations
sigma_list <- c(0.05, 0.1, 0.4, 1)

## 3. Sample sizes
n_list <- c(25, 50, 100)

## 4. Monte-Carlo settings
N      <- 501 
alpha  <- 0.05
set.seed(987)

## 5. Storage
results_list <- list()
counter <- 1

## 6. Main loop
for (f_name in names(f_list)) {
  for (sigma in sigma_list) {
    for (n in n_list) {
      
      ## equi‐spaced design on [0,1]
      x <- seq(0, 1, length.out = n)
      f <- f_list[[f_name]]
      m <- f(x)
      
      message(sprintf("Simulation %d:  f = %s, sigma = %g, n = %d",
                      counter, f_name, sigma, n))
      
      # returns a named numeric vector: c(decision_leq=…, decision_le=…, …)
      result_values <- rejection_rate_parallel(x, m, sigma, n, N, alpha, counter)
      
      # bind your 3 “metadata” columns with the 4 decision‐rates in one row
      results_list[[counter]] <- cbind(
        data.frame(
          regression_function    = f_name,
          sigma                  = sigma,
          n                      = n,
          stringsAsFactors       = FALSE
        ),
        # turn the named vector into a one‐row data.frame
        as.data.frame(
          as.list(result_values),
          stringsAsFactors = FALSE
        )
      )
      
      counter <- counter + 1
    }
  }
}

# Combine all rows into one data frame
results_df <- do.call(rbind, results_list)
results_df

#write.csv(results_df, file = "results_second_order_bootstrap_no_continuity_correction_1.csv", row.names = FALSE)

