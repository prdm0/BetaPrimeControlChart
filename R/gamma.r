# Distribuição Gamma
densidade_gamma_mean <- function(x, k, theta) {
  # k é shape
  # theta é scale
  # mu = k * theta
  # var = k * theta^2/n

  n <- length(x)
  x <- mean(x)

  k <- n * k
  theta <- theta / n
  x^(k - 1) * exp(-x / theta) / (theta^k * gamma(k))
}

densidade_gamma_mean <- Vectorize(FUN = densidade_gamma_mean, vectorize.args = c("x"))

# Dados de uma v.a. com distribuição gamma
random_gamma_mean <- function(n_sample, n, k, theta) {
  observacoes_mean <- function() {
    rgamma(n = n, shape = k, scale = theta)
  }
  lapply(X = 1L:n_sample, FUN = \(i) observacoes_mean())
}

n <- 10

dados <- random_gamma_mean(n_sample = 1000, n = n, k = 0.1, theta = 0.9)

log_like_gamma <- function(par, x) {
  fdp <- function(par, x) {
    k <- par[1L]
    theta <- par[2L]
    densidade_gamma_mean(x = x, k = k, theta = theta)
  }
  -sum(log(fdp(par, x)))
}

optim(par = c(1.5, 1.5), fn = log_like_gamma, x = dados, method = "Nelder-Mead")

# integrate(f = \(x) densidade_gamma_mean(x = x, k = 1, theta = 1), lower = 0, upper = 100)
# integrate(f = \(x) x * densidade_gamma_mean(x = x, k = 1, theta = 3), lower = 0, upper = 100)

# E <-
# integrate(
#   f = \(x) x * densidade_gamma_mean(x = x, k = 11, theta = 0.7),
#   lower = 0,
#   upper = Inf)$value

# E2 <-
# integrate(
#   f = \(x) x^2 * densidade_gamma_mean(x = x, k = 11, theta = 0.7),
#   lower = 0,
#   upper = Inf)$value
# E2 - E^2
