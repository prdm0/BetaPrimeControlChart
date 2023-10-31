densidade_gamma_mean <- function(x, alpha, beta, n) {
  alpha <- alpha / n
  beta <- n * beta
  beta^(n * alpha) / gamma(n * alpha) * x^(n * alpha - 1) * exp(-beta * x)
}

# Dados de uma v.a. com distribuição gamma
random_gamma_mean <- function(n_sample, n, alpha, beta) {
  observacoes_mean <- function() {
    mean(rgamma(n = n, shape = alpha / n, rate = beta))
  }
  sapply(X = 1L:n_sample, FUN = \(i) observacoes_mean())
}

n <- 20

dados <- random_gamma_mean(n_sample = 10000, n = n, alpha = 20, beta = 5)

log_like_gamma <- function(par, x, n) {
  fdp <- function(par, x) {
    alpha <- par[1L]
    beta <- par[2L]
    densidade_gamma_mean(x = x, alpha = alpha, beta = beta, n = n)
  }
  -sum(log(fdp(par, x)))
}

optim(par = c(1.5, 1.5), fn = log_like_gamma, n = n, x = dados, method = "Nelder-Mead")
