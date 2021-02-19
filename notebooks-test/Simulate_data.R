rm(list = ls())
library(simstudy)
library(glmnet)
library(stats)
library(rdd)
library(ggplot2)

set.seed(1)
# Number of Observations
N <- 1e3
total.covar <- 50 + 1e3
# Number of covariates (excluding W and unobservable)
p <- total.covar - 2
# Simulate Data
mu.vector <- rep(0, total.covar)
variance.vector <- abs(rnorm(total.covar, mean = 1, sd = .5))
simulated.data <- as.data.frame.matrix(genCorGen(n = N, nvars = total.covar, params1 = mu.vector, params2 = variance.vector, dist = 'normal',  rho = .5,
                                                 corstr = 'ar1', wide='True'))[2:(total.covar+1)]
colnames(simulated.data)[1] <- 'W' # Running variable for RDD
colnames(simulated.data)[total.covar] <- 'C' # Unobservable variable
# Random assignment for A/B test
T <- rep(0, N)
T[0:(N/2)] <- 1
T <- sample(T)
X <- simulated.data[, 2:(total.covar-1)]
covariate.names <- colnames(X)
# Independent error terms
error <- rnorm(n = N)
# Make W a function of the X's and unobservable (for RDD)
simulated.data$W <- simulated.data$W + .5 * simulated.data$C + 3 * simulated.data$V80 - 6 * simulated.data$V81
# Assign treatment, based on a threshold along W
treated <- (simulated.data$W > 0) * 1.0
# True coefficients on controls
beta.true.linear <- rnorm(p, mean = 5, sd = 5)
beta.true.linear[30:p] <- 0
# Functional form of Y for A/B test:
Y.ab <- 2.0 * T + data.matrix(X) %*% beta.true.linear + .6 * simulated.data$C + error
# Functional form of Y for RDD (function of treatment, W, X's, and unobservable C)
Y.rdd <- 1.2 * treated - 4 * simulated.data$W  + data.matrix(X) %*% beta.true.linear + .6 * simulated.data$C + error
df <- cbind(Y.ab, Y.rdd, T, simulated.data)
colnames(df)[1:3] <- c('Y.ab', 'Y.rdd', 'T')
X.colnames <- colnames(X)