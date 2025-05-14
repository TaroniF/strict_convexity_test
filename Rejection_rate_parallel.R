rejection_rate_parallel <- function(x, m, sigma, 
                                    n = length(x), 
                                    N = 50, 
                                    alpha = 0.05, 
                                    seed = 0) {
  library(foreach)
  library(parallel)
  library(doParallel)
  
  # ---- 1) spin up the cluster ----
  ncores <- detectCores() - 1
  cl <- makeCluster(ncores, outfile = "")
  registerDoParallel(cl)
  # make sure each worker sees your working directory (adjust path as needed)
  clusterEvalQ(cl, setwd("C:/Users/franc/OneDrive - Politecnico di Milano/Documenti/ETH/Convexity_test/Strict_convexity_test"))
  
  # ---- 2) parallel loop, returning one row per i ----
  iter_results <- foreach(i = seq_len(N), 
                          .combine = rbind) %dopar% {
                            set.seed(167 * (i + 1) + seed)
                            
                            # placeholder for error message
                            err_msg <- NA_character_
                            
                            # try to compute the four decisions; on error, record it and return NAs
                            decisions <- tryCatch({
                              source("main.R")
                              eps <- rnorm(n, 0, sigma)
                              y   <- m + eps
                              
                              test <- convexity_test(x = x, y = y, alpha = alpha)
                              c(test$decision_leq,
                                test$decision_le,
                                test$decision_leq_corrected,
                                test$decision_le_corrected)
                            }, error = function(e) {
                              err_msg <<- conditionMessage(e)
                              rep(NA_real_, 4)
                            })
                            
                            # build a one‐row data.frame
                            data.frame(
                              iteration               = i,
                              decision_leq            = decisions[1],
                              decision_le             = decisions[2],
                              decision_leq_corrected  = decisions[3],
                              decision_le_corrected   = decisions[4],
                              error                   = err_msg,
                              stringsAsFactors        = FALSE
                            )
                          }
  
  # ---- 3) tear down the cluster ----
  stopCluster(cl)
  message("Cluster stopped.")
  
  # ---- 4) report any errors ----
  errs <- subset(iter_results, !is.na(error))
  if (nrow(errs) > 0) {
    message("Errors occurred in these iterations:")
    print(errs)
  } else {
    message("No errors in any iteration.")
  }
  
  # ---- 5) compute the four rejection‐rates ----
  rejection_rate <- with(iter_results, c(
    decision_leq           = mean(decision_leq,           na.rm = TRUE),
    decision_le            = mean(decision_le,            na.rm = TRUE),
    decision_leq_corrected = mean(decision_leq_corrected, na.rm = TRUE),
    decision_le_corrected  = mean(decision_le_corrected,  na.rm = TRUE)
  ))
  class(rejection_rate) <- "rejection_rate"
  
  return(rejection_rate)
}
