rbp <- function(n, mu, phi, ...) {
  a <- mu * (phi + 1)
  b <- phi + 2

  x <- qbeta(
    runif(n = n, min = 0, max = 1),
    shape1 = a,
    shape2 = b,
    ...
  )
  x / (1 - x)
}

densidade_beta_prime <- function(x, mu, phi) {
  1 / beta(mu * (phi + 1), phi + 2) * x^(mu * (phi + 1) - 1) * (1 + x)^(-(mu + 1) * (phi + 1) - 1)
}

densidade_beta_prime_media <- function(x, mu, phi) {
  x <- sum(x)
  n <- length(x)

  lambda <- (n * mu * (mu + phi^2 - 2 * phi + n * mu * phi - 2 * n * mu + 1)) / ((phi - 1) * (mu + phi - 1))
  delta <- (2 * mu + phi^2 - phi + n * mu * phi - 2 * n * mu) / (mu + phi - 1)

  n * densidade_beta_prime(x = x, mu = lambda, phi = delta)
}


densidade_beta_prime_media <- Vectorize(FUN = densidade_beta_prime_media, vectorize.args = c("x"))

random_beta_prime <- function(n_sample, n, mu, phi) {
  observacoes <- function() {
    rbp(n = n, mu = mu, phi = phi)
  }
  lapply(X = 1L:n_sample, FUN = \(i) observacoes())
}

n <- 100
dados <- random_beta_prime(n_sample = 200, n = n, mu = 0.07, phi = 10.5)

log_like <- function(par, x) {
  mu <- par[1]
  phi <- par[2]
  -sum(log(densidade_beta_prime_media(x = unlist(x), mu = mu, phi = phi)))
}

optim(
  par = c(1.5, 1.5),
  fn = log_like,
  x = dados,
  method = "Nelder-Mead"
)
