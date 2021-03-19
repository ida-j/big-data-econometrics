library(glmnet)
library(stats)
library(ggplot2)
library(AER)

#Simulate Data
set.seed(1)
N <- 1e3 # Number of Observations
# Simulate Data for Example in Intro
weak <- rnorm(n = N)
good <- rnorm(n = N)
correl <- rnorm(n = N) 
C <- 3 * correl + rnorm(n = N)
X <- -.000001 * round(weak) + 2.5 * correl + 2.5 * good + .5 * C + rnorm(n = N)
Y <- 2.0 * X + 1.1 * C + rnorm(n = N)
df <- cbind(Y, X, correl, good, weak)
colnames(df) <- c('Y', 'X', 'correl', 'good', 'weak')
df <- data.frame(df)

# No Instrument, naive regression
fit.none <- lm(Y ~ X, data = df)
summary(fit.none)
X.none <- fit.none$coefficients[2]
ci.none <- confint(fit.none, 'X', level = 0.95)

# Fit using a good instrument (strong first stage and satisfies exclusion restriction)
fit.good <- ivreg(Y ~ X | good, data = df)
summary(fit.good, vcov = sandwich, df = Inf, diagnostics = TRUE)
X.good <- fit.good$coefficients[2]
ci.good <- confint(fit.good, 'X', level = 0.95)


# Fit using a weak instrument (satisfies exclusion restriction)
fit.weak <- ivreg(Y ~ X | weak, data = df)
summary(fit.weak, vcov = sandwich, df = Inf, diagnostics = TRUE)
X.weak <- fit.weak$coefficients[2]
ci.weak <- confint(fit.weak, 'X', level = 0.95)

# Fit using an instrument that doesn't satisfy the exclusion restriction (but with a strong first stage)
fit.correl <- ivreg(Y ~ X | correl, data = df)
summary(fit.correl, vcov = sandwich, df = Inf, diagnostics = TRUE)
X.correl <- fit.correl$coefficients[2]
ci.correl <- confint(fit.correl, 'X', level = 0.95)


# Vector of T coefficients and confidence intervals
coefs.ab <- c(2.0, X.good, X.none, X.weak, X.correl)
ci.low <- c(NaN, ci.good[1], ci.none[1], ci.weak[1], ci.correl[1])
ci.high <- c(NaN, ci.good[2], ci.none[2], ci.weak[2], ci.correl[2])
dat <- cbind(coefs.ab, ci.low, ci.high, c('True Effect', 'Good Instrument', 'No Instrument', 'Weak Instrument', 'Failed Exclusion Restriction'))
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
  ylab("Coefficient on X") +
  ggtitle("Estimated Causal Effect of X on Y") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5))


# Clear workspace
rm(list=ls())
invisible(gc()) 

#Simulate Data
set.seed(1)
N <- 1e3 # Number of Observations
p <- 500 # Number of instruments (excluding unobservable)

# Generate matrix M of p randomly assigned AB tests
M <- matrix(0, N, p)
M[0:(N/2),] <- 1
M <- apply(t(M), 1, function(d) sample(d, N))
# Create an unobservable C
C <- rnorm(n = N)

# Create X as a function of weak instruments, good instruments, instruments that don't satisfy the exclusion restriction, and C
beta.X <- rnorm(p, mean = 0, sd = 5) * sign(rnorm(n = p))
beta.X[1:20] <- rnorm(20, mean = 1, sd = 4)
beta.X[21:25] <- 100
beta.X[201:p] <- 0
X <- M %*% beta.X + 120 * C + rnorm(n = N)

# Create Y as a function of X, C, and other variables that affect X
beta.Y <- rnorm(20, mean = 5, sd = 3) * sign(rnorm(n = 20))
Y <- 2.0 * X + 200 * C + rnorm(n = N) + M[, 1:20] %*% beta.Y

# Create data frames with and without variables that don't satisfy the exclusion restriction
colnames(M) <- colnames(data.frame('M' = M))
M.excl <- M[, 21:p]
df.allM <- data.frame('Y' = Y, 'X' = X, M)
df.excl <- data.frame('Y' = Y, 'X' = X, M.excl)

lasso.fit <- cv.glmnet(M.excl, df.excl$X, alpha=1)
coef <- predict(lasso.fit, type = "nonzero")
M.excl.names <- colnames(M.excl)
Z <- M.excl.names[unlist(coef)]
# Instruments selected by LASSO:
length(Z)
Z

# Run the IV
Z.list <- paste("~ ", paste(Z, collapse = " + "))
fit.lasso.excl <- ivreg(Y ~ X, Z.list, data = df.excl)
summary(fit.lasso.excl, vcov = sandwich, df = Inf, diagnostics = TRUE)
X.lasso.excl <- fit.lasso.excl$coefficients[2]
ci.lasso.excl <- confint(fit.lasso.excl, 'X', level = 0.95, vcov = sandwich)

# 1st Stage
firststage <- paste("X ~ ", paste(Z, collapse = " + "))
fit.firststage <- lm(firststage, data = df.excl)
summary(fit.firststage)

# No instruments
fit.noinst <- lm(Y ~ X, data = df.allM)
X.noinst <- fit.noinst$coefficients[2]
ci.noinst <- confint(fit.noinst, 'X', level = 0.95)

# Largest subset of all possible IV
M.names <- colnames(M)
M.largest.list <- paste("~ ", paste(M.names[1:p], collapse = " + "))
fit.largest <- ivreg(Y ~ X, M.largest.list, data = df.allM)
X.largest <- fit.largest$coefficients[2]
ci.largest <- confint(fit.largest, 'X', level = 0.95, vcov = sandwich)


# Vector of T coefficients and confidence intervals
coefs.ab <- c(2.0, X.lasso.excl, X.noinst, X.largest)
ci.low <- c(NaN, ci.lasso.excl[1], ci.noinst[1], ci.largest[1])
ci.high <- c(NaN, ci.lasso.excl[2], ci.noinst[2], ci.largest[2])
dat <- cbind(coefs.ab, ci.low, ci.high, c('True Effect', 'LASSO','No Instruments',
                                          'All Instruments'))
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
  ylab("Coefficient on X") +
  ggtitle("Estimated Causal Effect of X on Y") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5))
