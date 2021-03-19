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
## A/B test
# Use Lasso of Y on X to select H

lasso.fit.outcome <- cv.glmnet(data.matrix(X), df$Y.ab, alpha=1)
coef <- predict(lasso.fit.outcome, type = "nonzero")
H <- X.colnames[unlist(coef)]
# Variables selected by LASSO:
H

# Use Lasso of T on X to select K
lasso.fit.propensity <- cv.glmnet(data.matrix(X), df$T, alpha=1)
coef <- predict(lasso.fit.propensity, type = "nonzero")
K <- X.colnames[unlist(coef)]
# Variables selected by LASSO:
K

# Perform OLS of Y on T, controlling for H union K
# Union of selected variables:
H_union_K.names <- unique(c(H, K))
H_union_K.names
sum.H_union_K <- paste(H_union_K.names, collapse = " + ")
eq.H_union_K <- paste("Y.ab ~ T + ", sum.H_union_K)

# OLS regression, using all covariates selected by double selection
fit.double <- lm(eq.H_union_K, data = df)
T.double <- fit.double$coefficients[2]
ci.double <- confint(fit.double, 'T', level = 0.95)

# Results:
T.double
ci.double
summary(fit.double)

## comparison plots
# Alternate methods:
#    OLS of Y on T
fit.simple <- lm('Y.ab ~ T', data = df)
#summary(fit.simple)
T.simple <- fit.simple$coefficients[2]
ci.simple <- confint(fit.simple, 'T', level = 0.95)

#    OLS of Y on T,controlling for all of X
sum.X <- paste(X.colnames, collapse = " + ")
eq.control.all <- paste("Y.ab ~ T + ", sum.X)
fit.allX <- lm(eq.control.all, data = df)
#summary(fit.allX)
T.allX <- fit.allX$coefficients[2]
ci.allX <- confint(fit.allX, 'T', level = 0.95)

#    OLS of Y on T,controlling for (almost) all of X
sum.Xmost <- paste(X.colnames[0:(N-10)], collapse = " + ")
eq.control.almost <- paste("Y.ab ~ T + ", sum.Xmost)
fit.mostX <- lm(eq.control.almost, data = df)
#summary(fit.allX)
T.mostX <- fit.mostX$coefficients[2]
ci.mostX <- confint(fit.mostX, 'T', level = 0.95)

#    OLS of Y on T, with a subset of X
sum.X.subset <- paste(c(X.colnames[5:15], X.colnames[80:90]), collapse = " + ")
eq.control.subset <- paste("Y.ab ~ T + ", sum.X.subset)
fit.subsetX <- lm(eq.control.subset, data = df)
#summary(fit.subsetX)
T.subsetX <- fit.subsetX$coefficients[2]
ci.subsetX <- confint(fit.subsetX, 'T', level = 0.95)

# Vector of T coefficients and confidence intervals
coefs.ab <- c(2.0, T.simple, T.allX, T.mostX, T.subsetX, T.double)
ci.low <- c(NaN, ci.simple[1], ci.allX[1], ci.mostX[1], ci.subsetX[1], ci.double[1])
ci.high <- c(NaN, ci.simple[2], ci.allX[2], ci.mostX[2], ci.subsetX[2], ci.double[2])
dat <- cbind(coefs.ab, ci.low, ci.high, c('True Effect', 'No Controls', 'All X', 'Largest Subset of X', 'Limited Controls', 'Double Selection'))
colnames(dat)[4] <- 'Model'
dat <- data.frame(dat)
dat$coefs.ab <- as.double(levels(dat$coefs.ab))[dat$coefs.ab]
dat$ci.low <- as.double(levels(dat$ci.low))[dat$ci.low]
dat$ci.high <- as.double(levels(dat$ci.high))[dat$ci.high]
dat$Model <- factor(dat$Model, levels = dat$Model)

# Create bar graph
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = dat$ci.high,
              ymin = dat$ci.low)
p <- ggplot(data = dat, aes(x = Model, y = coefs.ab, fill = Model))
p + geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  ylab("Coefficient on Treatment") +
  ggtitle("Estimated Causal Effect of T on Y, for various models") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


## RDD
# Lasso of Y on X to select H
lasso.fit.outcome <- cv.glmnet(data.matrix(X), df$Y.rdd, alpha=1)
coef <- predict(lasso.fit.outcome, type = "nonzero")
H <- X.colnames[unlist(coef)]
# Variables selected by LASSO:
H

# Lasso of W on X to select K
lasso.fit.propensity <- cv.glmnet(data.matrix(X), df$W, alpha=1)
coef <- predict(lasso.fit.propensity, type = "nonzero")
K <- X.colnames[unlist(coef)]
# Variables selected by LASSO:
K

# ### Perform RDD of Y on W, Controlling for H union K

# Union of selected variables:
H_union_K.names <- unique(c(H, K))
H_union_K.names
sum.H_union_K <- paste(H_union_K.names, collapse = " + ")
eq.H_union_K <- paste("Y.rdd ~ W | ", sum.H_union_K)

# RDD, using all covariates selected by double selection
fit.rdd <- RDestimate(eq.H_union_K, data = df)
summary(fit.rdd)
