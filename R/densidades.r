# Função densidade de probabilidade da média amostral de uma sequência de
# variáveis aleatórias iids com distribuição \Gamma(\alpha, \beta). Temos   que
# \overline{x} \sim \Gamma(n \alpha, \beta / n). Temos que \alpha > 0 e \beta >
# 0, em que \alpha é parâmetro de forma e \beta é parâmetro de escala.
densidade_gamma_xbarra <- function(x, mu, beta, ...) {
  # Densidade de \overline{X} sem reparametrização(beta / n)^(n * alpha) * x^(n
  # * alpha - 1) * exp(-beta * x / n) / gamma(n * alpha)
  n <- length(x)
  1 / gamma(n * mu) * (n * mu / beta)^(n * mu) * x^(n * mu - 1) * exp(-n * mu * x / beta)
}

# Função densidade de probabilidade da média amostral de uma sequência de
# variáveis aleatórias com distribuição Gaussiana Inversa (IG ~ (\mu, \lambda)).
# x > 0, \mu > 0 e \lambda > 0, em que \mu é parâmetro de locação e \lambda é
# parâmetro de forma. Temos que \overline{x} ~ IG(\mu, n\lambda).
densidade_gaussiana_inversa <- function(x, mu, lambda, ...) {
  n <- length(x)
  lambda <- n * lambda
  sqrt(lambda / (2 * pi * x^3)) * exp(-lambda * (x - mu)^2 / (2 * mu^2 * x))
}

# integrate(densidade_gamma_xbarra, lower = 0, upper = 10, mu = 0.4, beta = 0.6)
# integrate(densidade_gaussiana_inversa, lower = 0, upper = 100, mu = 1, lambda = 1)
