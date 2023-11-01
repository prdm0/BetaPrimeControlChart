pak::pkg_install("statmod")

# Distribuição Gamma
densidade_inversa_gaussiana_mean <- function(x, mu, lambda) {
  # mu é escala
  # lambda é forma
  # E(\olverline{X}) = \mu
  # Var(\olverline{X}) = \mu^3 / \lambda
  n <- length(x)
  x <- mean(x)

  lambda <- n * lambda

  sqrt(lambda / (2 * pi * x^3)) * exp(-(lambda * (x - mu)^2) / (2 * mu^2 * x))
}

densidade_inversa_gaussiana_mean <- Vectorize(FUN = densidade_inversa_gaussiana_mean, vectorize.args = c("x"))

# Dados de uma v.a. com distribuição gamma
random_inversa_gaussiana <- function(n_sample, n, mu, lambda) {
  observacoes <- function() {
    statmod::rinvgauss(n = n, mean = mu, shape = lambda)
  }
  lapply(X = 1L:n_sample, FUN = \(i) observacoes())
}

n <- 3

dados <- random_inversa_gaussiana(n_sample = 1000, n = n, mu = 10, lambda = 0.9)

log_like_inversa_gaussiana <- function(par, x) {
  fdp <- function(par, x) {
    mu <- par[1L]
    lambda <- par[2L]
    densidade_inversa_gaussiana_mean(x = x, mu = mu, lambda = lambda)
  }
  -sum(log(fdp(par, x)))
}

optim(par = c(1.5, 1.5), fn = log_like_inversa_gaussiana, x = dados, method = "Nelder-Mead")

# integrate(f = \(x) densidade_gamma_mean(x = x, k = 1, theta = 1), lower = 0, upper = 100)

# E <- integrate(f = \(x) x * densidade_inversa_gaussiana_mean(x = x, mu = 10, lambda = 1), lower = 0, upper = Inf)$value
# E2 <- integrate(f = \(x) x^2 * densidade_inversa_gaussiana_mean(x = x, mu = 10, lambda = 1), lower = 0, upper = Inf)$value
# E2 - E^2
