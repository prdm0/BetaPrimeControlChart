# Distribuição Gamma
densidade_gamma_mean <- function(x, k, theta, n) {
  # k é shape
  # theta é scale
  # mu = k * theta
  # var = k * theta^2/n
  k <- n * k
  theta <- theta / n
  x^(k - 1) * exp(-x / theta) / (theta^k * gamma(k))
}

# Dados de uma v.a. com distribuição gamma
random_gamma_mean <- function(n_sample, n, k, theta) {
  observacoes_mean <- function() {
    mean(rgamma(n = n, shape = k, scale = theta))
  }
  sapply(X = 1L:n_sample, FUN = \(i) observacoes_mean())
}

n <- 10

dados <- random_gamma_mean(n_sample = 10000, n = n, k = 3, theta = 5.5)

log_like_gamma <- function(par, x, n) {
  fdp <- function(par, x) {
    k <- par[1L]
    theta <- par[2L]
    densidade_gamma_mean(x = x, k = k, theta = theta, n = n)
  }
  -sum(log(fdp(par, x)))
}

optim(par = c(1.5, 1.5), fn = log_like_gamma, n = n, x = dados, method = "Nelder-Mead")

integrate(f = \(x) densidade_gamma_mean(x = x, k = 1, theta = 1, n = 1), lower = 0, upper = 100)


integrate(f = \(x) x * densidade_gamma_mean(x = x, alpha = 1, beta = 3, n = 2), lower = 0, upper = 100)
