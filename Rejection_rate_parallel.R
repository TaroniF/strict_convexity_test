rejection_rate_parallel <- function(x, m, sigma, n = length(x), N = 50, alpha = 0.05, seed = 0) {
  library(foreach)
  library(parallel)
  library(doParallel)
  
  ncores <- parallel::detectCores() - 1
  cl <- NULL
  
  result <- tryCatch({
    cl <- makeCluster(ncores, outfile = "")
    registerDoParallel(cl)
    
    # Set working directory on each worker
    clusterEvalQ(cl, setwd("C:/Users/franc/OneDrive - Politecnico di Milano/Documenti/ETH/Convexity_test"))
    
    # Parallel loop: return a data.frame with both decision and error message
    iter_results <- foreach(i = 1:N, 
                            .combine = rbind,
                            .errorhandling = "pass") %dopar% {
                              set.seed(167 * (i + 1) + seed)
                              
                              decision <- NA
                              err_msg <- ""
                              
                              # Wrap iteration in tryCatch to capture errors and log them
                              tryCatch({
                                source("main-v5.R")
                                eps <- rnorm(n, 0, sigma)
                                y <- m + eps
                                
                                # Assume convexity_test returns a list with an element named 'decision'
                                decision <- convexity_test(x = x, y = y, alpha = alpha)$decision
                              }, error = function(e) {
                                err_msg <<- conditionMessage(e)
                              })
                              
                              data.frame(iteration = i, decision = decision, error = err_msg, stringsAsFactors = FALSE)
                            }
    
    # Print error log: show only those iterations where an error occurred.
    error_log <- iter_results[iter_results$error != "", ]
    if(nrow(error_log) > 0){
      message("Errors occurred in the following iterations:")
      print(error_log)
    } else {
      message("No errors in any iteration.")
    }
    
    # Compute and return the mean decision, ignoring NA values.
    mean(iter_results$decision, na.rm = TRUE)
    
  }, error = function(e) {
    message("ERROR in rejection_rate_parallel: ", conditionMessage(e))
    NA_real_
    
  }, finally = {
    if(!is.null(cl)){
      stopCluster(cl)
      message("Cluster stopped correctly")
    }
  })
  
  return(result)
}