#' @ export
d_xbarra_normal <- function(x, mu, var) {
  n <- length(x)
  x <- mean(x)

  dnorm(x, mean = mu, sd = var / sqrt(n))
}

#' @export
r_normal <- function(n_lotes, n, mu, var, lotes = TRUE) {
  observacoes <- function() {
    rnorm(n = n, mean = mu, sd = sqrt(var))
  }

  obs <- lapply(X = 1L:n_lotes, FUN = \(i) observacoes())

  if (lotes) {
    return(obs)
  } else {
    return(sapply(X = obs, FUN = mean))
  }
}

#' @export
log_like_normal <- function(par, x) {
  fdp <- function(par, x) {
    mu <- par[1L]
    var <- par[2L]
    d_xbarra_normal(x = unlist(x), mu = mu, var = var)
  }
  -sum(log(fdp(par, x)))
}
