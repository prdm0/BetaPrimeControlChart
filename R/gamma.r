d_xbarra_gamma <- function(x, n, mu, k) {
  dgamma(x = x, shape = n * k, scale = mu / (n * k))
}

# As observacoes sao geradas de uma distribuicao gamma reparametrizada com theta
# = mu/k.
#' @export
r_gamma <- function(n_lotes, n, mu, k, ...) {
  observacoes <- function() {
    rgamma(n = n * n_lotes, shape = k, scale = mu / k, ...)
  }
  m <- matrix(data = observacoes(), nrow = n_lotes, ncol = n, byrow = TRUE)
  colnames(m) <- paste("n", 1L:ncol(m), sep = "_")
  rownames(m) <- paste("sample", 1L:nrow(m), sep = "_")
  return(m)
}

#' @export
mle_gamma <-
  function(data,
           start = c(1, 1),
           epsilon = 1e-6,
           linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
           ...) {
    x <- data
    n <- ncol(x)
    x <- apply(X = x, FUN = \(i) mean(i), MARGIN = 1L)

    log_like_gamma <- function(x, n, mu, k) {
      -sum(log(d_xbarra_gamma(x = x, n = n, mu = mu, k = k)))
    }

    grad_gamma <- function(x, n, mu, k) {
      numDeriv::grad(
        func = \(par) log_like_gamma(x = x, n = n, mu = par[1L], k = par[2L]),
        x = c(mu, k)
      )
    }

    lbfgs::lbfgs(
      call_eval = \(par) log_like_gamma(x = x, n = n, mu = par[1L], k = par[2L]),
      call_grad = \(par) grad_gamma(x = x, n = n, mu = par[1L], k = par[2L]),
      vars = start,
      invisible = 1,
      linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
      epsilon = 1e-6,
      ...
    )
  }

limits_gamma <-
  function(data = NULL,
           alpha = 0.0027,
           mu = NULL,
           k = NULL) {
    if (!is.matrix(data)) {
      stop("The object data must be a matrix!")
    }

    if (is.null(mu) || is.null(k)) {
      mle <- mle_gamma(data = data)
      mu <- mle$par[1L]
      k <- mle$par[2L]
    }

    n <- ncol(data)

    r <- qgamma(p = c(alpha / 2, 1 - alpha / 2), shape = n * k, scale = mu / (n * k))
    return(list(li = r[1L], ls = r[2L], mu_hat = mu, k_hat = k))
  }

stats_gamma <- function(data, alpha = 0.0027, mu = NULL, k = NULL) {
  if (!is.matrix(data)) {
    stop("The object data must be a matrix!")
  }
  n <- ncol(data)
  vec_sum <- apply(X = data, FUN = \(i) mean(i), MARGIN = 1L)
  limites <- limits_gamma(data = data, alpha = alpha, mu = mu, k = k)
  fora <- sum(vec_sum < limites$li) + sum(vec_sum > limites$ls)
  alpha_hat <- fora / length(vec_sum)
  list(
    mu_hat = limites$mu_hat,
    k_hat = limites$k_hat,
    alpha_hat = alpha_hat,
    ARL = 1 / alpha_hat,
    MRL = log(0.5) / log(1 - alpha_hat),
    SDRL = sqrt((1 - alpha_hat) / (alpha_hat^2)),
    li = limites$li,
    ls = limites$ls
  )
}

set.seed(0)
dados <- r_gamma(n_lotes = 10000, n = 100, mu = 1, k = 1.7)

# Estimando os parâmetros por maxima verossimilhança
stats_gamma(data = dados, alpha = 0.0027)

# Utilizando os parametros verdadeiros
stats_gamma(data = dados, alpha = 0.0027, mu = 1, k = 1.7)
