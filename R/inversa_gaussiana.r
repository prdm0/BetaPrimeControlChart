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
  statmod$dinvgauss(x = x, mean = mu, shape = n * lambda)
}

# Gera observacoes de uma variavel aleatoria com distribuicao inversa gaussiana
#' @export
r_inversa_gaussiana <- function(n_lotes, n, mu, lambda, lotes = TRUE) {
  observacoes <- function() {
    statmod$rinvgauss(n = n, mean = mu, shape = lambda)
  }
  lapply(X = 1L:n_lotes, FUN = \(i) observacoes())
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

#' @export
limites_inversa_gaussiana <- function(dados = NULL, alpha = 0.0027, mu, lambda) {
  if (!is.list(dados)) {
    stop("O objeto dados deve ser uma lista.")
  }
  n <- length(dados[[1]])
  r <- statmod$qinvgauss(p = c(alpha / 2, 1 - alpha / 2), mean = mu, shape = n * lambda)
  return(list(li = r[1L], ls = r[2L]))
}

arl0_inversa_gaussiana <- function(dados, alpha = 0.0027, mu = 1, lambda = 0.7) {
  vec_mean <- sapply(X = dados, FUN = mean)
  limites <-
    limites_inversa_gaussiana(
      dados = dados,
      alpha = alpha,
      mu = mu,
      lambda = lambda
    )
  fora <- sum(vec_mean < limites$li) + sum(vec_mean > limites$ls)
  alpha_hat <- fora / length(vec_mean)
  list(alpha_hat = alpha_hat, arl0 = 1 / alpha_hat)
}

dados <- r_inversa_gaussiana(n_lotes = 100e3L, n = 50L, mu = 1, lambda = 1.7)
arl0_inversa_gaussiana(dados = dados, alpha = 0.0027, mu = 1, lambda = 1.7)
