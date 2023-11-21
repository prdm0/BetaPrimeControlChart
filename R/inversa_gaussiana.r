box::use(
  statmod
)

# Densidade da variável aleatória x barra quando x segue a distribuicao inversa
# gaussiana com mu e lambda, em que mu eh parametro de localizacao e lambda eh
# parametro de escala.
#' @export
d_xbarra_inversa_gaussiana <- function(x, mu, lambda) {
  n <- length(x)
  x <- mean(x)

  lambda <- n * lambda

  sqrt(lambda / (2 * pi * x^3)) * exp(-(lambda * (x - mu)^2) / (2 * mu^2 * x))
}

# Gera observacoes de uma variavel aleatoria com distribuicao inversa gaussiana
#' @export
r_inversa_gaussiana <- function(n_lotes, n, mu, lambda, lotes = TRUE) {
  observacoes <- function() {
    statmod$rinvgauss(n = n, mean = mu, shape = lambda)
  }

  obs <- lapply(X = 1L:n_lotes, FUN = \(i) observacoes())

  if (lotes) {
    return(obs)
  } else {
    return(sapply(X = obs, FUN = mean))
  }
}

#' @export
log_like_inversa_gaussiana <- function(par, x) {
  fdp <- function(par, x) {
    mu <- par[1L]
    lambda <- par[2L]
    d_xbarra_inversa_gaussiana(x = unlist(x), mu = mu, lambda = lambda)
  }
  -sum(log(fdp(par, x)))
}
