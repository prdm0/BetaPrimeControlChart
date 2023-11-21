# Distribuição Gamma A densidade da variavel aleatoria x barra quando x segue a
# distribuicao gamma com k e theta, em que k eh parametro de forma e theta eh
# escala.
#' @export
d_xbarra_gamma <- function(x, mu, k) {
  n <- length(x)
  x <- mean(x)
  dgamma(x = x, shape = n * k, scale = mu / (n * k))
}

# As observacoes sao geradas de uma distribuicao gamma reparametrizada com theta
# = mu/k.
#' @export
r_gamma <- function(n_lotes, n, mu, k, lotes = TRUE, ...) {
  observacoes <- function() {
    rgamma(n = n, shape = k, scale = mu / k, ...)
  }

  obs <- lapply(X = 1L:n_lotes, FUN = \(i) observacoes())

  if (lotes) {
    return(obs)
  } else {
    return(sapply(X = obs, FUN = mean))
  }
}

#' @export
log_like_gamma <- function(par, x) {
  fdp <- function(par, x) {
    mu <- par[1L]
    k <- par[2L]
    d_xbarra_gamma(x = unlist(x), mu = mu, k = k)
  }
  -sum(log(fdp(par, x)))
}
