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
r_gamma <- function(n_lotes, n, mu, k, ...) {
  observacoes <- function() {
    rgamma(n = n, shape = k, scale = mu / k, ...)
  }
  lapply(X = 1L:n_lotes, FUN = \(i) observacoes())
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

#' @export
limites_gamma <- function(dados = NULL, alpha = 0.0027, mu, k) {
  if (!is.list(dados)) {
    stop("O objeto dados deve ser uma lista.")
  }
  n <- length(dados[[1]])
  r <- qgamma(p = c(alpha / 2, 1 - alpha / 2), shape = n * k, scale = mu / (n * k))
  return(list(li = r[1L], ls = r[2L]))
}

arl0_gamma <- function(dados, alpha = 0.0027, mu = 1, k = 0.7) {
  vec_mean <- sapply(X = dados, FUN = mean)
  limites <- limites_gamma(dados = dados, alpha = alpha, mu = mu, k = k)
  fora <- sum(vec_mean < limites$li) + sum(vec_mean > limites$ls)
  alpha_hat <- fora / length(vec_mean)
  list(alpha_hat = alpha_hat, arl0 = 1 / alpha_hat)
}

dados <- r_gamma(n_lotes = 100e3L, n = 50L, mu = 1, k = 1.7)
arl0_gamma(dados = dados, alpha = 0.0027, mu = 1, k = 1.7)
