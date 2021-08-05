library(rstan)
SLR <- "
data {
  int<lower=0> N;   // number of observations
  int<lower=0> P; // number of covariates
  matrix[N, P] X;   // predictor matrix
  vector[N] y;      // outcome vector
}
parameters {
  real alpha;           // intercept
  vector[P] beta;       // coefficients for predictors
  real<lower=0> sigma;  // error scale
}
model {
  y ~ normal(alpha + X * beta, sigma);  // target density
}
"

weakSLR <- "
data {
  int<lower=0> N;   // number of observations
  int<lower=0> P; // number of covariates
  matrix[N, P] X;   // predictor matrix
  vector[N] y;      // outcome vector
}
parameters {
  real alpha;           // intercept
  vector[P] beta;       // coefficients for predictors
  real<lower=0> sigma;  // error scale
}
model {
alpha ~ normal(0, 5);
beta ~ normal(0, 1);
sigma ~ cauchy(0,1);
y ~ normal(alpha + X * beta, sigma);  // target density
}
"

N <- 60
P <- 50
alpha <- 0
set.seed(357203)
beta <- rbinom(n = P, size = 1, prob = .1)

x <- rnorm(N*P, 0, 1) %>% matrix(nrow = N)
y <- as.vector(alpha + x%*%beta + rnorm(N, 0, 1))
data <- list(X = x, y = y, N = length(y), P = ncol(x))

SLR_model <- stan(model_code = SLR, data = data)


weakSLR_model <- stan(model_code = weakSLR, data = data)




beta_samps <- extract(SLR_model, par = "beta")$beta

weakbeta_samps <- extract(weakSLR_model, par = "beta")$beta

ggplot() +
  geom_histogram(aes(x = beta_samps[,as.logical(beta)]), alpha = I(.5)) +
  geom_histogram(aes(x = weakbeta_samps[,as.logical(beta)]), alpha = I(.5), fill = "red")

ggplot() +
  geom_histogram(aes(x = beta_samps[,!as.logical(beta)]), alpha = I(.5)) +
  geom_histogram(aes(x = weakbeta_samps[,!as.logical(beta)]), alpha = I(.5), fill = "red")

sd(beta_samps)

sd(weakbeta_samps)


plot(apply(beta_samps, 2, mean), beta)
plot(apply(weakbeta_samps, 2, mean), beta)

sqrt(mean((apply(beta_samps, 2, mean)- beta)^2))

sqrt(mean((apply(weakbeta_samps, 2, mean)- beta)^2))
