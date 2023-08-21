# Preamble
library(MCMCpack)
library(tibble)
library(ggplot2)
library(pilot)
library(progress)
set.seed(221)
y <- read.csv("data/obs.csv")[, 'x']
n <- length(y)
y <- c(0, y)
iterations <- 2000
B <- 200
y_matrix <- t(matrix(rep(y, iterations), nrow = length(y)))
d <- 3
I_0 <- 1
z_0 <- 0

# Initialization
mult_helper <- function(p) {
  sample <- rmultinom(1, 1, p)
  return(which(sample == max(sample)))
}
gamma <- c(0)
sigma_squared <- matrix(c(1, 6, 15, 2), nrow = 1)
tau <- matrix(rep(1 / d, d * d), nrow = d)

# Progress bar
pb <- progress::progress_bar$new(
  format = "[:bar] :percent :elapsedfull", 
  total = B, 
  clear = FALSE, 
  width = 60
)

# Run Monte Carlo EM
for (l in 2:B) {
  gamma_EM <- gamma[l - 1]
  tau_EM <- tau[((l - 1) * 3 - 2):((l - 1) * 3), ]
  sigma_squared_EM <- sigma_squared[l - 1, ]
  
  # Draw I and z
  I <- matrix(rep(0, iterations * (n + 1)), nrow = iterations)
  I[, 1] <- I_0
  z <- matrix(rep(0, iterations * (n + 1)), nrow = iterations)
  z[, 1] <- z_0
  for (t in 2:(n + 1)) {
    I[1, t] <- mult_helper(tau_EM[I[1, t - 1], ])
    z[1, t] <- rnorm(1, gamma_EM * z[1, t - 1], sqrt(sigma_squared_EM[I[1, t]]))
  }
  for (k in 2:iterations) {
    I_EM <- I[k - 1, ]
    z_EM <- z[k - 1, ]
    
    # Draw new I and z
    for (t in 2:(n + 1)) { 
      p <- rep(0, d) 
      for (i in 1:d) {
        p[i] <- tau_EM[I_EM[t - 1], i] * 
          ifelse(t + 1 > (n + 1), 1, tau_EM[i, I_EM[t + 1]]) / sqrt(sigma_squared_EM[i]) * 
          dnorm((z_EM[t] - gamma_EM * z_EM[t  -1]) / sqrt(sigma_squared_EM[i]))
      }
      I_EM[t] <- mult_helper(p / sum(p))
      if (t + 1 > (n + 1)) {
        s <- (1 / sigma_squared_EM[I_EM[t]] + 1 / sigma_squared_EM[d + 1])^{-1}
        u <- (gamma_EM * z_EM[t - 1] / sigma_squared_EM[I_EM[t]] + y[t] / sigma_squared_EM[d + 1]) * s
      } else {
        s <- (1 / sigma_squared_EM[I_EM[t]] + 
                (gamma_EM)^2 / sigma_squared_EM[I_EM[t + 1]] + 1 / 
                sigma_squared_EM[d + 1])^{-1}
        u <- (gamma_EM * z_EM[t - 1] / sigma_squared_EM[I_EM[t]] + gamma_EM * 
                z_EM[t+1] / sigma_squared_EM[I_EM[t + 1]] + y[t] / sigma_squared_EM[d + 1]) * s
      }
      z_EM[t] <- rnorm(1, u, sqrt(s))
    }
    I[k, ] <- I_EM
    z[k, ] <- z_EM
  }
  
  # Update tau
  tau_counter <- matrix(rep(0, d * d), nrow = d)
  for (i in 1:d) {
    for (j in 1:d) {
      tau_counter[i, j] <- sum((I[, 1:n] == i) * (I[, 2:(n+1)] == j))
    }
  }
  tau_EM <- tau_counter / rowSums(tau_counter)
  tau <- rbind(tau, tau_EM)
  
  # Iteratively update sigma_squared and gamma
  for (m in 1:B) {
    sigma_vector <- matrix(rep(0, iterations * (n + 1)), nrow = iterations)
    for (k in 1:iterations) {
      for (t in 2:(n + 1)) {
        sigma_vector[k, t] <- sigma_squared_EM[I[k, t]]
      }
    }
    
    # Update sigma
    for (j in 1:d) {
      sigma_squared_EM[j] <- sum((I[, 2:(n + 1)] == j) * 
                                   (z[, 2:(n + 1)] - gamma_EM * z[, 1:n])^2) / 
        sum(I[, 2:(n + 1)] == j)
    }
    sigma_squared_EM <- c(sort(sigma_squared_EM[1:d]), sigma_squared_EM[d + 1])
    sigma_squared_EM[d + 1] <- sum((y_matrix[, 2:(n + 1)] - z[, 2:(n + 1)])^2) / (iterations * n)
    
    # Update gamma
    gamma_EM <- sum(z[, 2:(n + 1)] * z[, 1:n] / sigma_vector[, 2:(n + 1)]) / 
      sum(z[, 1:n]^2 / sigma_vector[, 2:(n + 1)])
  }
  
  # Save values
  gamma <- c(gamma, gamma_EM)
  sigma_squared <- rbind(sigma_squared, sigma_squared_EM)
  rownames(sigma_squared) <- NULL
  
  # Update progress bar
  pb$tick()
}

# Output MLE
tau[(dim(tau)[1] - 2):dim(tau)[1], ]
gamma[length(gamma)]
sqrt(sigma_squared[dim(sigma_squared)[1], ])

# Trace plots
gamma_frame <- tibble(x = seq_along(gamma), y = gamma)
ggplot(gamma_frame, aes(x = x, y = y)) +
  geom_line(color = pilot_color("navy")) +
  xlab("Samples") +
  ylab(expression(gamma)) +
  ggtitle(expression(paste("Trace Plot for ", gamma))) +
  theme_pilot()
sigma_frame <- data.frame(sqrt(sigma_squared))
sigma_frame_1 <- tibble(x = seq_along(sigma_frame$X1), y = sigma_frame$X1)
ggplot(sigma_frame_1, aes(x = x, y = y)) +
  geom_line(color = pilot_color("brown")) +
  xlab("Samples") +
  ylab(expression(sigma[1])) +
  ggtitle(expression(paste("Trace Plot for ", sigma[1]))) +
  theme_pilot()
sigma_frame_2 <- tibble(x = seq_along(sigma_frame$X2), y = sigma_frame$X2)
ggplot(sigma_frame_2, aes(x = x, y = y)) +
  geom_line(color = pilot_color("orange")) +
  xlab("Samples") +
  ylab(expression(sigma[2])) +
  ggtitle(expression(paste("Trace Plot for ", sigma[2]))) +
  theme_pilot()
sigma_frame_3 <- tibble(x = seq_along(sigma_frame$X3), y = sigma_frame$X3)
ggplot(sigma_frame_3, aes(x = x, y = y)) +
  geom_line(color = pilot_color("yellow")) +
  xlab("Samples") +
  ylab(expression(sigma[3])) +
  ggtitle(expression(paste("Trace Plot for ", sigma[3]))) +
  theme_pilot()
sigma_frame_4 <- tibble(x = seq_along(sigma_frame$X4), y = sigma_frame$X4)
ggplot(sigma_frame_4, aes(x = x, y = y)) +
  geom_line(color = pilot_color("purple")) +
  xlab("Samples") +
  ylab(expression(sigma[4])) +
  ggtitle(expression(paste("Trace Plot for ", sigma[4]))) +
  theme_pilot()

# ACF plots
acf_gamma_values <- acf(gamma, plot = FALSE, lag.max = 300)
acf_gamma_df <- data.frame(
  Lag = seq_along(acf_gamma_values$acf) - 1, 
  ACF = acf_gamma_values$acf
)
ggplot(acf_gamma_df, aes(x = Lag, y = ACF)) +
  geom_line(color = pilot_color("navy")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Lag", 
    y = expression(gamma), 
    title = expression(paste("ACF for ", gamma))
  ) +
  theme_pilot()
acf_sigma_1_values <- acf(sigma_frame$X1, plot = FALSE, lag.max = 300)
acf_sigma_1_df <- data.frame(
  Lag = seq_along(acf_sigma_1_values$acf) - 1,
  ACF = acf_sigma_1_values$acf
)
ggplot(acf_sigma_1_df, aes(x = Lag, y = ACF)) +
  geom_line(color = pilot_color("brown")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Lag",
    y = expression(sigma[1]),
    title = expression(paste("ACF for ", sigma[1]))
  ) +
  theme_pilot()
acf_sigma_2_values <- acf(sigma_frame$X2, plot = FALSE, lag.max = 300)
acf_sigma_2_df <- data.frame(
  Lag = seq_along(acf_sigma_2_values$acf) - 1,
  ACF = acf_sigma_2_values$acf
)
ggplot(acf_sigma_2_df, aes(x = Lag, y = ACF)) +
  geom_line(color = pilot_color("orange")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Lag",
    y = expression(sigma[2]),
    title = expression(paste("ACF for ", sigma[2]))
  ) +
  theme_pilot()
acf_sigma_3_values <- acf(sigma_frame$X3, plot = FALSE, lag.max = 300)
acf_sigma_3_df <- data.frame(
  Lag = seq_along(acf_sigma_3_values$acf) - 1,
  ACF = acf_sigma_3_values$acf
)
ggplot(acf_sigma_3_df, aes(x = Lag, y = ACF)) +
  geom_line(color = pilot_color("yellow")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Lag",
    y = expression(sigma[3]),
    title = expression(paste("ACF for ", sigma[3]))
  ) +
  theme_pilot()
acf_sigma_4_values <- acf(sigma_frame$X4, plot = FALSE, lag.max = 300)
acf_sigma_4_df <- data.frame(
  Lag = seq_along(acf_sigma_4_values$acf) - 1,
  ACF = acf_sigma_4_values$acf
)
ggplot(acf_sigma_4_df, aes(x = Lag, y = ACF)) +
  geom_line(color = pilot_color("purple")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Lag",
    y = expression(sigma[4]),
    title = expression(paste("ACF for ", sigma[4]))
  ) +
  theme_pilot()