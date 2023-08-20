# Sampling and Point Estimation for Three-Level Hierarchical Model
- We study a three-level state-space hierarchical model with observation $y$, hidden variables $(z,I)$, and parameter $\theta = (T,\gamma,\sigma^2)$, where $T$ is the transition matrix for the Markov chain corresponding to $I$, $\gamma$ is a scalar in $(-1,1)$, and $\sigma^2$ is a vector of variances.
- We find the complete-data loglikelihood and reasonable priors for $(T,\gamma,\sigma^2)$, and compute the full logposterior for $(z,I,\theta)$.
- Next, we design and implement a Gibbs sampler to sample from the posterior for $(z,I,\theta)$, then use this to compute marginal posterior distributions as well as posterior means and variances for parameters in $\theta$; diagnostic plots show MCMC chain convergence.
- Finally, we design and implement a Monte Carlo EM algorithm to compute the maximum likelihood estimate for parameters in $\theta$; here as well, diagnostic (trace and autocorrelation) plots show that the proposed Monte Carlo EM algorithm runs to convergence.

# Code
The code for all the computations and plots in this paper can be found at `main.R`.

# Data
The observed data can be accessed under the data folder [here](data/obs.csv).

# Paper
The paper can be accessed under the docs folder [here](docs/hierarchical-mcmc.pdf).

# Figures
The figures for the paper generated by `main.R` can be found in the `figures` folder.
