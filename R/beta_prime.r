rbp <- function(n, a, b, ...) {
  x <- qbeta(
    runif(n = n, min = 0, max = 1),
    shape1 = a,
    shape2 = b,
    ...
  )
  x / (1 - x)
}

# Distribuição Beta Prime
densidade_beta_prime_mean <- function(x, mu, phi) {
  n <- length(x)
  x <- sum(x)

  a <- mu * (phi + 1)
  b <- phi + 2

  phi <- (n * a * (a + b^2 - 2 * b + n * a * b - 2 * n * a + 1)) / ((b - 1) * (a + b - 1))
  delta <- (2 * a + b^2 - b + n * a * b - 2 * n * a) / (a + b - 1)

  n * (1 / beta(phi, delta) * x^(phi - 1) / (1 + x)^(phi + delta))
}

densidade_beta_prime_mean <- Vectorize(FUN = densidade_beta_prime_mean, vectorize.args = c("x"))


random_beta_prime_mean <- function(n_sample, n, mu, phi) {
  a <- mu * (phi + 1)
  b <- phi + 2

  phi <- (n * a * (a + b^2 - 2 * b + n * a * b - 2 * n * a + 1)) / ((b - 1) * (a + b - 1))
  delta <- (2 * a + b^2 - b + n * a * b - 2 * n * a) / (a + b - 1)

  observacoes_mean <- function() {
    rbp(n = n, a = phi, b = delta)
  }

  lapply(X = 1L:n_sample, FUN = \(i) observacoes_mean())
}


dados <- rbp(n = 1000L, a = 0.5, b = 0.3)

log_like <- function(par, x) {
  a <- par[1]
  b <- par[2]
  -sum(log(densidade_beta_prime(x = x, a = a, b = b)))
}

optim(
  par = c(1, 1),
  fn = log_like,
  x = dados,
  method = "Nelder-Mead"
)
