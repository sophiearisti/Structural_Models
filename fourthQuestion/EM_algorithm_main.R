##upload the file

global_path <- "/Users/sophiaaristizabal/Documents/GitHub/structural models/fourthQuestion"

doc_name <- "/Q4_data.csv"
#open csv

data <- read.csv(file.path(global_path, doc_name), header = TRUE)

#######first we will graph the data's distribution just to know how it looks #################
names(data)

hist(as.numeric(data[[1]]), prob = TRUE, 
     main = "empirical distribution", 
     xlab = "X-Axis Label", 
     ylab = "Probability",
     col = "lightblue")

#EM algorithm is implemented in the file EM_algorithm_functions.R

#just as the sample code
pisX3<-c(0.5,0.5,0.5);
pisX2<-c(0.5,0.5);

###################################First Assumption:two Normal distributions###################################
#initial guesses 
#miu and sigma repectively 
theta_0 <- c(4, 2.5) #first distribution parameters
theta_1 <- c(6, 2.5) #second distribution parameters

result <- EM_algorithm(Y = as.numeric(data[[1]]), 
             pis = pisX2, 
             thetas = list(theta_0, theta_1), 
             normal = TRUE)

# Print results
cat("The first assumption is two Normal distributions.\n\n")

cat("Final parameter estimates (mean, sd):\n")
for (i in seq_along(result$theta)) {
  cat(sprintf("  Component %d: mean = %.4f, sd = %.4f\n",
              i, result$theta[[i]][1], result$theta[[i]][2]))
}

cat("\nFinal mixture probabilities:\n")
for (i in seq_along(result$pi)) {
  cat(sprintf("  Pi[%d] = %.4f\n", i, result$pi[i]))
}

cat("\nNumber of iterations:\n")
cat(result$iter, "\n")

png("mixture_fit_plot_2Normal.png", width = 800, height = 600)

#create a plot of the fitted mixture
plot_mixture_fit(Y = as.numeric(data[[1]]), 
                 pis = result$pi, 
                 thetas = result$theta, 
                 normal = TRUE)

dev.off()

###################################Second Assumption:two Extreme Value distributions###################################

theta_0 <- c(8, 2.5) #first distribution parameters
theta_1 <- c(6, 2.5) #second distribution parameters

result <- EM_algorithm(Y = as.numeric(data[[1]]), 
             pis = pisX2, 
             thetas = list(theta_0, theta_1), 
             normal = FALSE)

# Print results
cat("The second assumption is two Extreme Value distributions.\n")
cat("Final parameter estimates:\n")
print(result$theta)
cat("Final mixture probabilities:\n")
print(result$pi)
cat("Number of iterations:\n")
print(result$iter)

png("mixture_fit_plot_2EVD.png", width = 800, height = 600)

# Create plot for the second assumption
plot_mixture_fit(Y = as.numeric(data[[1]]), 
                 pis = result$pi, 
                 thetas = result$theta, 
                 normal = FALSE)
dev.off()

###################################Third Assumption:three Normal distributions###################################
theta_0 <- c(4, 2.5)  #first distribution parameters
theta_2 <- c(8, 2.5)  #second distribution parameters
theta_1 <- c(6, 2.5)  #Third distribution parameters

result <- EM_algorithm(Y = as.numeric(data[[1]]), 
             pis = pisX3, 
             thetas = list(theta_0, theta_1, theta_2), 
             normal = TRUE)

# Print results
cat("The third assumption is three Normal distributions.\n")
cat("Final parameter estimates:\n")
print(result$theta)
cat("Final mixture probabilities:\n")
print(result$pi)
cat("Number of iterations:\n")
print(result$iter)


png("mixture_fit_plot_3Normal.png", width = 800, height = 600)

# Create plot for the second assumption
plot_mixture_fit(Y = as.numeric(data[[1]]), 
                 pis = result$pi, 
                 thetas = result$theta, 
                 normal = FALSE)
dev.off()

###################################Fourth Assumption:three Extreme Value distributions###################################
theta_0 <- c(4, 2.5)  #first distribution parameters
theta_2 <- c(8, 2.5)  #second distribution parameters
theta_1 <- c(6, 2.5)  #Third distribution parameters

result <- EM_algorithm(Y = as.numeric(data[[1]]), 
             pis = pisX3, 
             thetas = list(theta_0, theta_1, theta_2), 
             normal = FALSE)

# Print results
cat("The third assumption is two Normal distributions.\n")
cat("Final parameter estimates:\n")
print(result$theta)
cat("Final mixture probabilities:\n")
print(result$pi)
cat("Number of iterations:\n")
print(result$iter)

# Create plot for the fourth assumption
png("mixture_fit_plot_3EVD.png", width = 800, height = 600)

# Create plot for the second assumption
plot_mixture_fit(Y = as.numeric(data[[1]]), 
                 pis = result$pi, 
                 thetas = result$theta, 
                 normal = FALSE)
dev.off()