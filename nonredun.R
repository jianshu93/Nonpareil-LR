library(ggplot2) # For plotting
library(stats) # For numerical integration and optimization

# Function to calculate kappa given alpha, beta, and LR
calculate_kappa <- function(alpha, beta, LR) {
  # Handling the edge case when LR is exactly 0
  if (LR == 0) {
    return(0)
  } else {
  numerator <- integrate(function(t) exp(-t) * t^(alpha - 1), lower = 0, upper = beta * log(LR + 1))$value
  denominator <- gamma(alpha)
  kappa <- numerator / denominator
  return(kappa)
  }
}

# Objective function for optimization
# params is a vector that contains the values of alpha and beta to be optimized
objective_function <- function(params, LR, target_kappa) {
  alpha <- params[1]
  beta <- params[2]
  calculated_kappa <- calculate_kappa(alpha, beta, LR)
  # The objective is to minimize the squared difference between calculated kappa and target kappa
  return((calculated_kappa - target_kappa)^2)
}

# Target kappa value for LR = 10^8
target_kappa_init <- 0.43
LR_init <- 1.6e6
target_kappa <- 0.43
LR <- 1.6e6

# Initial guesses for alpha and beta
initial_guess <- c(alpha = 3.5, beta = 0.2)

# Run the optimization to find the best alpha and beta
optim_result <- optim(initial_guess, objective_function, LR = LR, target_kappa = target_kappa, method = "L-BFGS-B", lower = c(0.001,0.001))

# Extract the optimized alpha and beta
optimized_alpha <- optim_result$par[1]
optimized_beta <- optim_result$par[2]

cat("Optimized alpha:", optimized_alpha, "\n")
cat("Optimized beta:", optimized_beta, "\n")

N_d = (optimized_alpha-1)/optimized_beta

cat("Nonpareil diversity:", N_d, "\n")
# Generating a sequence of logarithmically spaced numbers from 1e-5 to 1e8
LR_range <- exp(seq(log(1e-12), log(1e24), length.out = 100))

# Calculating kappa for the range of LR values using placeholder values for alpha and beta
alpha_opt <- optimized_alpha # Placeholder for the optimized alpha value
beta_opt <-  optimized_beta# Placeholder for the optimized beta value
kappa_values <- sapply(LR_range, function(LR) calculate_kappa(alpha_opt, beta_opt, LR))

# Plotting kappa vs. LR
data <- data.frame(LR = LR_range, Kappa = kappa_values)
ggplot(data) +
  geom_line(aes(x = LR, y = Kappa*100), size=1.5) +
  scale_x_log10() +
  labs(x = "Sequencing efforts (bp)", y = "Estimated average coverage(%)") +
  theme_classic() +
  theme(panel.grid.minor=element_blank())+
  theme(axis.text.x=element_text(size=13))+
  theme(axis.text=element_text(colour="black",size=14))+
  theme(legend.text=element_text(color="black",size=13),
        legend.title=element_text(color="black",size=13))+
  theme(axis.title = element_text(size=14))+guides(fill=guide_legend(ncol=1))+
  
  geom_point(aes(x=LR_init,y=target_kappa_init*100),colour="darkred",size=3, shape=17)+
  geom_hline(yintercept=100, linetype="dashed", color = "red", size=1.5) + 
  geom_hline(yintercept=95, linetype="dashed", color = "pink", size=1.5) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=1.5)

