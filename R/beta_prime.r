box::use(
  extraDistr
)

#' @export
d_xbarra_beta_prime <- function(x, mu, phi) {
  n <- length(x)
  x <- sum(x)

  alpha <- mu * (phi + 1)
  beta <- phi + 2

  lambda <- n * alpha * (alpha + beta^2 - 2 * beta + n * alpha * beta - 2 * n * alpha + 1) / ((beta - 1) * (alpha + beta - 1))
  delta <- (2 * alpha + beta^2 - beta + n * alpha * beta - 2 * n * alpha) / (alpha + beta - 1)

  n * extraDistr$dbetapr(x = x, shape1 = lambda, shape2 = delta)
}

#' @export
r_beta_prime <- function(n_lotes, n, mu, phi, ...) {
  observacoes <- function() {
    alpha <- mu * (phi + 1)
    beta <- phi + 2
    extraDistr$rbetapr(n = n, shape1 = alpha, shape2 = beta, ...)
  }
  lapply(X = 1L:n_lotes, FUN = \(i) observacoes())
}

#' @export
log_like_beta_prime <- function(par, x) {
  mu <- par[1]
  phi <- par[2]
  -sum(log(d_xbarra_beta_prime(x = unlist(x), mu = mu, phi = phi)))
}

qbp <- function(p, n, mu, phi, epsilon = 1e-4) {
  fdp <- function(x, n, mu, phi) {
    alpha1 <- mu * (phi + 1)
    beta <- phi + 2

    lambda <- n * alpha1 * (alpha1 + beta^2 - 2 * beta + n * alpha1 * beta - 2 * n * alpha1 + 1) / ((beta - 1) * (alpha1 + beta - 1))
    delta <- (2 * alpha1 + beta^2 - beta + n * alpha1 * beta - 2 * n * alpha1) / (alpha1 + beta - 1)

    n * extraDistr$dbetapr(x = x, shape1 = lambda, shape2 = delta)
  }
  repeat{
    r <-
      integrate(
        f = function(x) fdp(x = n * x, n = n, mu = mu, phi = phi), lower = 0, upper = q
      )$value
    if (abs(r - p) < epsilon) {
      return(q)
    }
    q <- q + epsilon
  }
}

qbp <- Vectorize(qbp, vectorize.args = "p")

limites_beta_prime <- function(data = NULL, alpha = 0.0027, mu, phi) {
  if (!is.list(data)) {
    stop("O objeto dados deve ser uma lista.")
  }

  n <- length(data[[1L]])

  r <- qbp(p = c(alpha / 2, 1 - alpha / 2), n = n, mu = mu, phi = phi)
  return(list(li = r[1L], ls = r[2L]))
}

peformance_beta_prime <- function(data, alpha = 0.0027, mu = 1, phi = 0.7) {
  n <- data[[1L]] |> length()
  vec_sum <- sapply(X = data, FUN = \(i) mean(i))
  limites <- limites_beta_prime(data = data, alpha = alpha, mu = mu, phi = phi)
  fora <- sum(vec_sum < limites$li) + sum(vec_sum > limites$ls)
  alpha_hat <- fora / length(vec_sum)
  list(alpha_hat = alpha_hat, arl0 = 1 / alpha_hat)
}

set.seed(0)
dados <- r_beta_prime(n_lotes = 100e3L, n = 1000, mu = 1, phi = 1.7)
peformance_beta_prime(data = dados, alpha = 0.0027, mu = 1, phi = 1.7)

# #' @export
# limites_beta_prime <- function(dados = NULL, n = NULL, alpha = 0.0027, mu, phi) {
#   if (!is.null(dados) && is.list(dados)) {
#     n <- length(dados[[1L]])
#   } else {
#     if (is.null(n)) {
#       stop("Dados ou n devem ser fornecidos")
#     }
#   }

#   qbp <- function(p, n, mu, phi) {
#     fdp <- function(x, n, mu, phi) {
#       alpha1 <- mu * (phi + 1)
#       beta <- phi + 2

#       lambda <- n * alpha1 * (alpha1 + beta^2 - 2 * beta + n * alpha1 * beta - 2 * n * alpha1 + 1) / ((beta - 1) * (alpha1 + beta - 1))
#       delta <- (2 * alpha1 + beta^2 - beta + n * alpha1 * beta - 2 * n * alpha1) / (alpha1 + beta - 1)

#       n * extraDistr$dbetapr(x = x, shape1 = lambda, shape2 = delta)
#     }
#     q <- 1e-4
#     repeat{
#       r <-
#         integrate(
#           f = function(x) fdp(x = x, n = n, mu = mu, phi = phi), lower = 0, upper = q
#         )$value
#       if (abs(r - p) < 1e-4) {
#         return(q)
#       }
#       q <- q + 1e-4
#     }
#   }

#   qbp <- Vectorize(qbp, vectorize.args = "p")

#   r <- qbp(p = c(alpha / 2, 1 - alpha / 2), n = n, mu = mu, phi = phi)
#   return(list(li = r[1L], ls = r[2L]))
# }

# dados <- r_beta_prime(n_lotes = 1000, n = 100, mu = 0.6, phi = 0.7)


# fora <- function(dados, alpha = 0.0027, mu = 0.6, phi = 0.7) {
#   vec_mean <- sapply(X = dados, FUN = mean)

#   n <- length(dados[[1L]])

#   limites <- limites_beta_prime(dados = dados, alpha = alpha, mu = mu, phi = phi)
#   return(sum(vec_mean < limites[[1L]] | vec_mean > limites[[2L]]))
# }

# fora(dados = dados, alpha = 0.0027, mu = 0.6, phi = 0.7)
# limites_beta_prime(dados = dados, alpha = 0.0027, mu = 1, phi = 0.7)
