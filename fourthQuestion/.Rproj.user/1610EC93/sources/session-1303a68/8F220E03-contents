#######Now we will implement the EM algorithm to estimate the parameters of a Gaussian Mixture Model (GMM)#######
dextreme_value <- function(y, mu, beta) {
  z <- (y - mu) / beta
  (1 / beta) * exp(-z) * exp(-exp(-z))
}

# Assumes you have a vector Y of observed counts

# Step 1: specify convergence criteria
EM_algorithm <- function(Y, pis, thetas, normal = TRUE) {
  obs <- length(Y)
  
  iter <- 1
  maxiter <- 1000
  lik0 <- 10000
  normdiff <- 10000
  tolerance <- 1e-6
  
  while (normdiff > tolerance && iter <= maxiter) {
    
    # === E-step ===
    EZ <- matrix(0, nrow = obs, ncol = length(pis))
    
    if (normal) {
      temp_pdf <- sapply(thetas, function(th) {
        dnorm(Y, mean = th[1], sd = th[2])
      })
    } else {
      temp_pdf <- sapply(thetas, function(th) {
        dextreme_value(Y, mu = th[1], beta = th[2])  # <- you must define this
      })
    }
    
    for (j in 1:length(pis)) {
      EZ[, j] <- pis[j] * temp_pdf[, j]
    }
    
    row_sums <- rowSums(EZ)
    EZ <- EZ / row_sums
    
    # === M-step ===
    theta_fun <- function(theta_vec) {
      # Convert flat vector back to list of duples
      thetas <- split(theta_vec, ceiling(seq_along(theta_vec) / 2))
      
      pdf_mat <- if (normal) {
        sapply(thetas, function(th) dnorm(Y, mean = th[1], sd = th[2]))
      } else {
        sapply(thetas, function(th) dextreme_value(Y, mu = th[1], beta = th[2]))
      }
      
      -sum(log(pdf_mat) * EZ)
    }
    
    theta_init <- unlist(thetas)
    
    opt_result <- optim(par = theta_init, fn = theta_fun)
    theta_vec <- opt_result$par
    thetas <- split(theta_vec, ceiling(seq_along(theta_vec) / 2))
    
    # Update pi
    pis <- colMeans(EZ)
    
    # === Compute likelihood for convergence ===
    pdf_mat <- if (normal) {
      sapply(thetas, function(th) dnorm(Y, mean = th[1], sd = th[2]))
    } else {
      sapply(thetas, function(th) dextreme_value(Y, mu = th[1], beta = th[2]))
    }
    
    lik1 <- sum(log(pdf_mat) * EZ) + sum(log(pis) * EZ)
    normdiff <- (lik0 - lik1)^2
    
    lik0 <- lik1
    iter <- iter + 1
  }
  
  return(list(theta = thetas, pi = pis, iter = iter))
}


plot_mixture_fit <- function(Y, pis, thetas, normal=TRUE, x_range = NULL, bins = 30) {
  # Y: vector of observed data
  # pis: vector of mixture weights
  # thetas: list of parameter vectors (e.g., list(c(4,2.5), c(6,2.5)))
  # x_range: optional x-axis values
  # bins: number of bins for histogram (if continuous)
  
  print(paste("Plotting mixture fit for", length(thetas), "distributions."))
  
  if (is.null(x_range)) {
      x_range <- seq(min(Y), max(Y), length.out = 300)
  }
  
  # Define PDF generator
  get_pdf <- function(x, theta, normal) {
    if (normal == TRUE) {
      dnorm(x, mean = theta[1], sd = theta[2])
    } else {
      z <- (x - theta[1]) / theta[2]
      (1 / theta[2]) * exp(-(z + exp(-z)))
    }
  }
  
  # Compute fitted mixture PDF
  #The sample codeonly did the sum but in this case we dont know if we need to sum 2 or 3 distributions
  Y_fit <- rowSums(sapply(1:length(pis), function(j) {
    pis[j] * get_pdf(x_range, thetas[[j]], normal)
  }))
  
  # Plot mixture
  plot(x_range, Y_fit, type = "l", col = "blue", lwd = 2,
       ylab = "Probability", xlab = "Y",
       main = paste("Fitted Mixture vs Empirical", sep=""))
  
  # Add empirical histogram
  hist(Y, breaks = bins, freq = FALSE,
         col = rgb(0.4, 0.4, 0.4, 0.5), add = TRUE)
  
  legend("topright", legend = c("Estimated mixture", "Empirical distribution"),
         col = c("blue", rgb(0.4, 0.4, 0.4, 0.5)), lwd = 2, pch = c(NA, 15))
}

