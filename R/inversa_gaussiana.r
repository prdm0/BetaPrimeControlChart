# Densidade da variável aleatória x barra quando x segue a distribuicao inversa
# gaussiana com mu e lambda, em que mu eh parametro de localizacao e lambda eh
# parametro de escala.
#' @export
d_xbarra_inverse_gaussian <- function(x, n, mu, lambda) {
  statmod::dinvgauss(x = x, mean = mu, shape = n * lambda)
}

# Gera observacoes de uma variavel aleatoria com distribuicao inversa gaussiana
#' @export
r_inverse_gaussian <- function(n_lotes, n, mu, lambda, lotes = TRUE) {
  observacoes <- function() {
    statmod::rinvgauss(n = n * n_lotes, mean = mu, shape = lambda)
  }
  m <- matrix(data = observacoes(), nrow = n_lotes, ncol = n, byrow = TRUE)
  colnames(m) <- paste("n", 1L:ncol(m), sep = "_")
  rownames(m) <- paste("sample", 1L:nrow(m), sep = "_")
  return(m)
}

#' @export
mle_inverse_gaussian <-
  function(data,
           start = c(1, 1),
           epsilon = 1e-6,
           linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
           ...) {
    x <- data
    n <- ncol(x)
    x <- apply(X = x, FUN = \(i) mean(i), MARGIN = 1L)

    log_like_inverse_gaussian <- function(x, n, mu, lambda) {
      -sum(log(d_xbarra_inverse_gaussian(x = x, n = n, mu = mu, lambda = lambda)))
    }

    grad_inverse_gaussian <- function(x, n, mu, lambda) {
      numDeriv::grad(
        func = \(par) log_like_inverse_gaussian(x = x, n = n, mu = par[1L], lambda = par[2L]),
        x = c(mu, lambda)
      )
    }

    lbfgs::lbfgs(
      call_eval = \(par) log_like_inverse_gaussian(x = x, n = n, mu = par[1L], lambda = par[2L]),
      call_grad = \(par) grad_inverse_gaussian(x = x, n = n, mu = par[1L], lambda = par[2L]),
      vars = start,
      invisible = 1,
      linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
      epsilon = 1e-6,
      ...
    )
  }

limits_inverse_gaussian <-
  function(data = NULL,
           alpha = 0.0027,
           mu = NULL,
           lambda = NULL,
           ...) {
    if (!is.matrix(data)) {
      stop("The object data must be a matrix!")
    }

    if (is.null(mu) || is.null(lambda)) {
      mle <- mle_inverse_gaussian(data = data, ...)
      mu <- mle$par[1L]
      lambda <- mle$par[2L]
    }
    n <- ncol(data)
    r <- statmod::qinvgauss(p = c(alpha / 2, 1 - alpha / 2), mean = mu, shape = n * lambda)
    return(list(li = r[1L], ls = r[2L], mu_hat = mu, lambda_hat = lambda))
  }

stats_inverse_gaussian <- function(data, alpha = 0.0027, mu = NULL, lambda = NULL, ...) {
  if (!is.matrix(data)) {
    stop("The object data must be a matrix!")
  }
  n <- ncol(data)
  vec_sum <- apply(X = data, FUN = \(i) mean(i), MARGIN = 1L)
  limites <- limits_inverse_gaussian(data = data, alpha = alpha, mu = mu, lambda = lambda, ...)
  fora <- sum(vec_sum < limites$li) + sum(vec_sum > limites$ls)
  alpha_hat <- fora / length(vec_sum)
  list(
    mu_hat = limites$mu_hat,
    lambda_hat = limites$lambda_hat,
    alpha_hat = alpha_hat,
    ARL = 1 / alpha_hat,
    MRL = log(0.5) / log(1 - alpha_hat),
    SDRL = sqrt((1 - alpha_hat) / (alpha_hat^2)),
    li = limites$li,
    ls = limites$ls
  )
}

set.seed(0)
dados <- r_inverse_gaussian(n_lotes = 100000, n = 200, mu = 1, lambda = 1.7)

# Estimando os parâmetros por maxima verossimilhança
stats_inverse_gaussian(data = dados, alpha = 0.0027)

# Utilizando os parametros verdadeiros
stats_inverse_gaussian(data = dados, alpha = 0.0027, mu = 1, lambda = 1.7)
