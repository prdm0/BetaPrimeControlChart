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

d_beta_prime <- function(x, mu, phi) {
  1 / beta(mu * (phi + 1), phi + 2) * x^(mu * (phi + 1) - 1) * (1 + x)^(-(mu + 1) * (phi + 1) - 1)
}

#' @export
d_xbarra_beta_prime <- function(x, mu, phi) {
  x <- sum(x)
  n <- length(x)

  lambda <- (n * mu * (mu + phi^2 - 2 * phi + n * mu * phi - 2 * n * mu + 1)) / ((phi - 1) * (mu + phi - 1))
  delta <- (2 * mu + phi^2 - phi + n * mu * phi - 2 * n * mu) / (mu + phi - 1)

  n * d_beta_prime(x = x, mu = lambda, phi = delta)
}

#' @export
r_beta_prime <- function(n_lotes, n, mu, phi, lotes = TRUE) {
  observacoes <- function() {
    rbp(n = n, mu = mu, phi = phi)
  }

  obs <- lapply(X = 1L:n_lotes, FUN = \(i) observacoes())

  if (lotes) {
    return(obs)
  } else {
    return(sapply(X = obs, FUN = mean))
  }
}

#' @export
log_like_beta_prime <- function(par, x) {
  mu <- par[1]
  phi <- par[2]
  -sum(log(d_xbarra_beta_prime(x = unlist(x), mu = mu, phi = phi)))
}
