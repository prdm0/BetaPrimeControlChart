d_xbarra_beta_prime <- function(x, n, mu, phi) {
  alpha <- mu * (phi + 1)
  beta <- phi + 2

  lambda <- n * alpha * (alpha + beta^2 - 2 * beta + n * alpha * beta - 2 * n * alpha + 1) / ((beta - 1) * (alpha + beta - 1))
  delta <- (2 * alpha + beta^2 - beta + n * alpha * beta - 2 * n * alpha) / (alpha + beta - 1)

  n * extraDistr::dbetapr(x = n * x, shape1 = lambda, shape2 = delta)
}

r_beta_prime <- function(n_lotes, n, mu, phi, ...) {
  observacoes <- function() {
    alpha <- mu * (phi + 1)
    beta <- phi + 2
    extraDistr::rbetapr(n = n * n_lotes, shape1 = alpha, shape2 = beta, ...)
  }
  m <- matrix(data = observacoes(), nrow = n_lotes, ncol = n, byrow = TRUE)
  colnames(m) <- paste("n", 1L:ncol(m), sep = "_")
  rownames(m) <- paste("sample", 1L:nrow(m), sep = "_")
  return(m)
}

#' @export
mle_beta_prime <-
  function(data,
           start = c(1, 1),
           epsilon = 1e-6,
           linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
           ...) {
    x <- data
    n <- ncol(x)
    x <- apply(X = x, FUN = \(i) mean(i), MARGIN = 1L)

    log_like_beta_prime <- function(x, n, mu, phi) {
      -sum(log(d_xbarra_beta_prime(x = x, n = n, mu = mu, phi = phi)))
    }

    grad_beta_prime <- function(x, n, mu, phi) {
      numDeriv::grad(
        func = \(par) log_like_beta_prime(x = x, n = n, mu = par[1L], phi = par[2L]),
        x = c(mu, phi)
      )
    }

    lbfgs::lbfgs(
      call_eval = \(par) log_like_beta_prime(x = x, n = n, mu = par[1L], phi = par[2L]),
      call_grad = \(par) grad_beta_prime(x = x, n = n, mu = par[1L], phi = par[2L]),
      vars = start,
      invisible = 1,
      linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
      epsilon = 1e-6,
      ...
    )
  }

qbp <- function(p, n, mu, phi) {
  alpha1 <- mu * (phi + 1)
  beta <- phi + 2
  lambda <- n * alpha1 * (alpha1 + beta^2 - 2 * beta + n * alpha1 * beta - 2 * n * alpha1 + 1) / ((beta - 1) * (alpha1 + beta - 1))
  delta <- (2 * alpha1 + beta^2 - beta + n * alpha1 * beta - 2 * n * alpha1) / (alpha1 + beta - 1)
  extraDistr::qbetapr(p = p, shape1 = lambda, shape2 = delta) / n
}

qbp <- Vectorize(qbp, vectorize.args = "p")

limits_beta_prime <-
  function(data = NULL,
           alpha = 0.0027,
           mu = NULL,
           phi = NULL,
           ...) {
    if (!is.matrix(data)) {
      stop("The object data must be a matrix!")
    }

    if (is.null(mu) || is.null(phi)) {
      mle <- mle_beta_prime(data = data, ...)
      mu <- mle$par[1L]
      phi <- mle$par[2L]
    }

    n <- ncol(data)

    r <- qbp(p = c(alpha / 2, 1 - alpha / 2), n = n, mu = mu, phi = phi)
    return(list(li = r[1L], ls = r[2L], mu_hat = mu, phi_hat = phi))
  }

stats_beta_prime <- function(data, alpha = 0.0027, mu = NULL, phi = NULL, ...) {
  if (!is.matrix(data)) {
    stop("The object data must be a matrix!")
  }
  n <- ncol(data)
  vec_sum <- apply(X = data, FUN = \(i) mean(i), MARGIN = 1L)
  limites <- limits_beta_prime(data = data, alpha = alpha, mu = mu, phi = phi, ...)
  fora <- sum(vec_sum < limites$li) + sum(vec_sum > limites$ls)
  alpha_hat <- fora / length(vec_sum)
  list(
    mu_hat = limites$mu_hat,
    phi_hat = limites$phi_hat,
    alpha_hat = alpha_hat,
    ARL = 1 / alpha_hat,
    MRL = log(0.5) / log(1 - alpha_hat),
    SDRL = sqrt((1 - alpha_hat) / (alpha_hat^2)),
    li = limites$li,
    ls = limites$ls
  )
}

set.seed(0)
dados <- r_beta_prime(n_lotes = 10e3L, n = 2050, mu = 1, phi = 1.7)

# Estimando os parâmetros por maxima verossimilhança
stats_beta_prime(data = dados, alpha = 0.0027)

# Utilizando os parametros verdadeiros
stats_beta_prime(data = dados, alpha = 0.0027, mu = 1, phi = 1.7)
