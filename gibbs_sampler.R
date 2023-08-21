# Preamble
library(MCMCpack)
library(truncnorm)
library(tibble)
library(ggplot2)
library(pilot)
library(reshape2)
library(plyr)
library(progress)
set.seed(221)
y <- read.csv("data/obs.csv")[, 'x']
n <- length(y)
y <- c(0, y)
d <- 3
I_0 <- 1
z_0 <- 0
nu_0 <- 4
tau_0 <- 2

# Helper multinomial function
mult_helper <- function(p) {
  sample <- rmultinom(1, 1, p)
  return(which(sample == max(sample)))
}

# Initialization
iterations <- 10000
I <- matrix(rep(0, iterations * (n + 1)), nrow = iterations)
I[, 1] <- I_0
z <- matrix(rep(0, iterations * (n + 1)), nrow = iterations)
z[, 1] <- z_0
tau <- matrix(rep(0, iterations * d * d), nrow = iterations * d)
tau[1:3, ] <- 1 / d
gamma <- rep(0, iterations)
gamma[1] <- 0
sigma_squared <- matrix(rep(0, iterations * (d + 1)), nrow = iterations)
sigma_squared[1, ] <- 8
for (t in 2:(n + 1)) {
  I[1, t] <- mult_helper(tau[I[1, t - 1], ])
  z[1, t] <- rnorm(1, gamma[1] * z[1, t - 1], sqrt(sigma_squared[1, I[1, t]]))
}

# Progress bar
pb <- progress::progress_bar$new(
  format = "[:bar] :percent :elapsedfull", 
  total = iterations, 
  clear = FALSE, 
  width = 60
)

# Implement Gibbs sampler
for (k in 2:iterations) {
  I_gibbs <- I[k - 1, ]
  z_gibbs <- z[k - 1, ]
  tau_gibbs <- tau[((k - 1) * 3 - 2):((k - 1) * 3), ]
  gamma_gibbs <- gamma[k - 1]
  sigma_squared_gibbs <- sigma_squared[k - 1, ]
  
  # Draw I and z
  for (t in 2:(n + 1)) {
    p <- rep(0, d)
    for (i in 1:d) {
      p[i] <- tau_gibbs[I_gibbs[t - 1], i] * 
        ifelse(t + 1 > (n + 1), 1, tau_gibbs[i, I_gibbs[t + 1]]) / sqrt(sigma_squared_gibbs[i]) * 
        dnorm((z_gibbs[t] - gamma_gibbs * z_gibbs[t - 1]) / sqrt(sigma_squared_gibbs[i]))
    }
    I_gibbs[t] <- mult_helper(p / sum(p))
    if (t + 1 > (n + 1)) {
      s <- (1 / sigma_squared_gibbs[I_gibbs[t]] + 1 / sigma_squared_gibbs[d + 1])^{-1}
      u <- (gamma_gibbs * z_gibbs[t - 1] / sigma_squared_gibbs[I_gibbs[t]] + 
              y[t] / sigma_squared_gibbs[d + 1]) * s
    } else {
      s <- (1 / sigma_squared_gibbs[I_gibbs[t]] + (gamma_gibbs)^2 / 
              sigma_squared_gibbs[I_gibbs[t + 1]] + 1 / sigma_squared_gibbs[d + 1])^{-1}
      u <- (gamma_gibbs * z_gibbs[t - 1] / sigma_squared_gibbs[I_gibbs[t]] + 
              gamma_gibbs * z_gibbs[t + 1] / sigma_squared_gibbs[I_gibbs[t + 1]] + y[t] / 
              sigma_squared_gibbs[d + 1]) * s
    }
    z_gibbs[t] <- rnorm(1, u, sqrt(s))
  }
  
  # Draw gamma
  placeholder <- rep(0, n + 1)
  for (t in 2:(n + 1)) {
    placeholder[t] <- sigma_squared_gibbs[I_gibbs[t]]
  }
  o <- (sum(z_gibbs[1:n] ** 2 / placeholder[2:(n + 1)]))^{-1}
  l <- sum(z_gibbs[1:n] * z_gibbs[2:(n + 1)] / placeholder[2:(n+1)]) * o
  gamma_gibbs <- rtruncnorm(1, -1, 1, l, sqrt(o))
  
  # Draw tau
  for (i in 1:d) {
    c <- rep(0, d)
    for (j in 1:d) {
      c[j] <- sum((I_gibbs[1:n] == i) * (I_gibbs[2:(n + 1)] == j))
    }
    tau_gibbs[i, ] <- rdirichlet(1, c + 1)
  }
  
  # Draw sigma_squared
  while (TRUE) {
    for (j in 1:d) {
      shape <- (nu_0 + sum(I_gibbs[2:(n + 1)] == j)) / 2
      scale <- (nu_0 * (tau_0)^2 + 
                  sum((I_gibbs[2:(n + 1)] == j) * 
                        (z_gibbs[2:(n + 1)] - gamma_gibbs * z_gibbs[1:n]) ** 2)) / 2
      sigma_squared_gibbs[j] <- rinvgamma(1, shape, scale)
    }
    if ((sigma_squared_gibbs[1] <= sigma_squared_gibbs[2]) && 
        (sigma_squared_gibbs[2] <= sigma_squared_gibbs[3])) {
      break
    }
  }
  shape <- (nu_0 + n) / 2
  scale <- (nu_0 * (tau_0)^2 + sum((y[2:(n + 1)] - z_gibbs[2:(n + 1)]) ** 2)) / 2
  sigma_squared_gibbs[d + 1] <- rinvgamma(1, shape, scale) 
  
  # Update values
  I[k, ] <- I_gibbs
  z[k, ] <- z_gibbs
  tau[(k * 3 - 2):(k * 3), ] <- tau_gibbs
  gamma[k] <- gamma_gibbs
  sigma_squared[k, ] <- sigma_squared_gibbs
  
  # Update progress bar
  pb$tick()
}

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

# Posterior means and standard deviations
colMeans(tau[seq(1, nrow(tau), 3), ])
colMeans(tau[seq(2, nrow(tau), 3), ])
colMeans(tau[seq(3, nrow(tau), 3), ])
sd(tau[seq(1, nrow(tau), 3), ])
sd(tau[seq(2, nrow(tau), 3), ])
sd(tau[seq(3, nrow(tau), 3), ])
mean(gamma)
sd(gamma)
sapply(sigma_frame, mean)
sapply(sigma_frame, sd)

# Marginal posterior distributions
ecdf_gamma_df <- data.frame(
  x = sort(unique(gamma)),
  y = ecdf(gamma)(sort(unique(gamma)))
)
ggplot(ecdf_gamma_df, aes(x = x, y = y)) +
  geom_line(color = pilot_color("navy")) +
  labs(
    x = expression(gamma),
    y = "ECDF",
    title = expression(paste("Marginal Posterior Distribution for ", gamma))
  ) +
  theme_pilot()
ecdf_sigma_1_df <- data.frame(
  x = sort(unique(sigma_frame$X1)),
  y = ecdf(sigma_frame$X1)(sort(unique(sigma_frame$X1)))
)
ggplot(ecdf_sigma_1_df, aes(x = x, y = y)) +
  geom_line(color = pilot_color("brown")) +
  labs(
    x = expression(sigma[1]),
    y = "ECDF",
    title = expression(paste("Marginal Posterior Distribution for ", sigma[1]))
  ) +
  theme_pilot()
ecdf_sigma_2_df <- data.frame(
  x = sort(unique(sigma_frame$X2)),
  y = ecdf(sigma_frame$X2)(sort(unique(sigma_frame$X2)))
)
ggplot(ecdf_sigma_2_df, aes(x = x, y = y)) +
  geom_line(color = pilot_color("orange")) +
  labs(
    x = expression(sigma[2]),
    y = "ECDF",
    title = expression(paste("Marginal Posterior Distribution for ", sigma[2]))
  ) +
  theme_pilot()
ecdf_sigma_3_df <- data.frame(
  x = sort(unique(sigma_frame$X3)),
  y = ecdf(sigma_frame$X3)(sort(unique(sigma_frame$X3)))
)
ggplot(ecdf_sigma_3_df, aes(x = x, y = y)) +
  geom_line(color = pilot_color("yellow")) +
  labs(
    x = expression(sigma[3]),
    y = "ECDF",
    title = expression(paste("Marginal Posterior Distribution for ", sigma[3]))
  ) +
  theme_pilot()
ecdf_sigma_4_df <- data.frame(
  x = sort(unique(sigma_frame$X4)),
  y = ecdf(sigma_frame$X4)(sort(unique(sigma_frame$X4)))
)
ggplot(ecdf_sigma_4_df, aes(x = x, y = y)) +
  geom_line(color = pilot_color("purple")) +
  labs(
    x = expression(sigma[1]),
    y = "ECDF",
    title = expression(paste("Marginal Posterior Distribution for ", sigma[4]))
  ) +
  theme_pilot()

# Duration at each state 
first_state <- apply(I, 2, function(x) sum(x == 1) / iterations)
second_state <- apply(I, 2, function(x) sum(x == 2) / iterations)
third_state <- apply(I, 2, function(x) sum(x == 3) / iterations)
time_series_df <- data.frame(
  t = 1:length(first_state), 
  first_state = first_state, 
  second_state = second_state, 
  third_state = third_state
)
melted_df <- melt(
  time_series_df, 
  id.vars = "t", 
  variable.name = "state", 
  value.name = "duration"
)
melted_df$state <- revalue(
  melted_df$state, 
  c(
    "first_state" = "First State", 
    "second_state" = "Second State", 
    "third_state" = "Third State"
  )
)
ggplot(melted_df, aes(x = t, y = duration, color = state)) +
  geom_line() +
  facet_wrap(~state, scales = "free", ncol = 1) +
  labs(
    x = "Time",
    y = "Duration",
    title = "Duration for Each State Over Time",
    color = "State"
  ) +
  theme_pilot() +
  scale_color_pilot()