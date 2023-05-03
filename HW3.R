# Setting up the number of observations and starting value
Totalnumber <- 100
Thetatemp <- 0

# Import data (here it is generated, can be switched to any custom function)
Vt <- rep.int(1, Totalnumber)
Wt <- rep.int(1, Totalnumber)

# Generate playground
df <- data.frame(Vt, Wt)
df['w'] <- 0
df['v'] <- 0
df['theta'] <- 0 # state
df['Yt'] <- 0 # observation

for (x in 1:Totalnumber) {
  df[x,'w'] <- rnorm(n = 1, mean = 0, sd = df[x,'Wt'])
  df[x,'v'] <- rnorm(n = 1, mean = 0, sd = df[x,'Vt'])
  df[x,'theta'] <- Thetatemp + df[x,'w']
  Thetatemp <- df[x,'theta']
  df[x,'Yt'] <- df[x,'theta'] + df[x,'v']
}

# Kalman filter implementation
P <- 1
theta_kalman <- numeric(Totalnumber)
theta_kalman[1] <- 0
for (t in 2:Totalnumber) {
  # Prediction
  theta_pred <- theta_kalman[t - 1]
  P_pred <- P + Wt[t]
  
  # Update
  K <- P_pred / (P_pred + Vt[t]) # Kalman gain
  theta_kalman[t] <- theta_pred + K * (df$Yt[t] - theta_pred)
  P <- (1 - K) * P_pred
}

df['Theta_est'] <- theta_kalman

# SIS (Sequential Importance Sampling) using Normal importance distribution
N_particles <- 1000
particles <- matrix(0, nrow = N_particles, ncol = Totalnumber)
weights <- rep(1 / N_particles, N_particles)
theta_SIS <- numeric(Totalnumber)
theta_SIS[1] <- 0

for (t in 2:Totalnumber) {
  # Importance sampling step
  particles[, t] <- rnorm(N_particles, mean = particles[, t - 1], sd = sqrt(Wt[t]))
  
  # Weight update step
  weights <- weights * dnorm(df$Yt[t], mean = particles[, t], sd = sqrt(Vt[t]))
  weights <- weights / sum(weights)
  
  # Resampling step
  index <- sample(1:N_particles, size = N_particles, replace = TRUE, prob = weights)
  particles <- particles[index, ]
  weights <- rep(1 / N_particles, N_particles)
  
  # SIS estimate
  theta_SIS[t] <- mean(particles[, t])
}

df['Theta_SIS'] <- theta_SIS

# Plot results
library(ggplot2)

ggplot(df, aes(x = 1:Totalnumber)) +
  geom_line(aes(y = Yt, color = "Observation")) +
  geom_line(aes(y = Theta_est, color = "Kalman Filter Estimate")) +
  geom_line(aes(y = Theta_SIS, color = "SIS Estimate")) +
  labs(title = "Local Level Model: Observation, Kalman Filter, and SIS Estimates",
       x = "Time", y = "Value", color = "Legend") +
  theme_minimal()

