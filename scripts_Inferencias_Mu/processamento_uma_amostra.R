# ==============================================================================
# INFERÊNCIA SOBRE UM VETOR DE MÉDIAS MULTIVARIADO
# Módulo de Funções Procedurais em R
# ==============================================================================

# ------------------------------------------------------------------------------
# TESTE DE NORMALIDADE MULTIVARIADA 
# ------------------------------------------------------------------------------

#' Avalia a Normalidade Multivariada via Shapiro-Wilk Multivariado
#' 
#' @param X Matriz ou data.frame de dados (n x p)
#' @return Objeto htest com o resultado do teste de Shapiro-Wilk Multivariado
teste_normalidade_multivariada <- function(X) {
  if (!requireNamespace("mvnormtest", quietly = TRUE)) {
    stop("O pacote 'mvnormtest' é necessário. Instale-o executando: install.packages('mvnormtest')")
  }
  
  X <- as.matrix(X)
  
  # mshapiro.test exige variáveis nas linhas e observações nas colunas (p x n)
  resultado <- mvnormtest::mshapiro.test(t(X))
  
  return(resultado)
}
# ------------------------------------------------------------------------------
# QUANTIL DE T² DE HOTELLING
# ------------------------------------------------------------------------------

#' Retorna o quantil superior da distribuição T² de Hotelling usando a distribuição F
#' c = [p * (n - 1) / (n - p)] * F_{p, n-p}(1 - alpha)
#' @param alpha Nível de significância (ex: 0.05)
#' @param p Número de variáveis
#' @param n Tamanho da amostra
quantil_t2_hotelling <- function(alpha, p, n) {
  if (n <= p) {
    stop("O tamanho da amostra (n) deve ser maior que o número de variáveis (p).")
  }
  
  f_critico <- qf(1 - alpha, df1 = p, df2 = n - p)
  c <- (p * (n - 1) / (n - p)) * f_critico
  return(c)
}

# ------------------------------------------------------------------------------
# TESTE T² DE HOTELLING PARA UMA AMOSTRA (H0: mu = mu0)
# ------------------------------------------------------------------------------

#' Executa o Teste T² de Hotelling de 1 amostra
#' @param X Matriz ou data.frame de dados (n x p)
#' @param mu0 Vetor de médias sob a hipótese nula H0 (dimensão p x 1)
#' @param alpha Nível de significância (padrão 0.05)
teste_hotelling_1amostra <- function(X, mu0, alpha = 0.05) {
  X <- as.matrix(X)
  mu0 <- matrix(mu0, ncol = 1)
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (length(mu0) != p) {
    stop("A dimensão do vetor mu0 deve ser igual ao número de colunas de X.")
  }
  
  xbar <- colMeans(X)
  S <- cov(X)
  
  # Estatística de teste T² observada
  diff_mu <- xbar - mu0
  T2_obs <- as.numeric(n * t(diff_mu) %*% solve(S) %*% diff_mu)
  
  # Quantil superior T²
  T2_crit <- quantil_t2_hotelling(alpha, p, n)
  
  # Estatística F equivalente e P-Valor
  F_obs <- T2_obs * (n - p) / (p * (n - 1))
  p_valor <- 1 - pf(F_obs, df1 = p, df2 = n - p)
  
  rejeita_H0 <- T2_obs > T2_crit
  
  resultado <- list(
    xbar = xbar,
    S = S,
    n = n,
    p = p,
    T2_obs = T2_obs,
    T2_crit = T2_crit,
    F_obs = F_obs,
    p_valor = p_valor,
    alpha = alpha,
    decisao = ifelse(rejeita_H0, "Rejeita-se H0", "Não se rejeita H0")
  )
  
  class(resultado) <- "teste_hotelling"
  return(resultado)
}

#' Método de exibição limpa para o resultado do teste T²
print.teste_hotelling <- function(x, ...) {
  cat("\n=======================================================\n")
  cat("     TESTE T² DE HOTELLING PARA UMA AMOSTRA            \n")
  cat("=======================================================\n")
  cat(sprintf(" Tamanho Amostral (n): %d  |  Variáveis (p): %d\n", x$n, x$p))
  cat(sprintf(" Nível de Significância (alpha): %.3f\n", x$alpha))
  cat("-------------------------------------------------------\n")
  cat(sprintf(" Estatística T² Observada: %.4f\n", x$T2_obs))
  cat(sprintf(" Quantil Crítico T² (1-alpha): %.4f\n", x$T2_crit))
  cat(sprintf(" Estatística F Equivalente: %.4f (gl1 = %d, gl2 = %d)\n", 
              x$F_obs, x$p, x$n - x$p))
  cat(sprintf(" p-valor: %.5e\n", x$p_valor))
  cat("-------------------------------------------------------\n")
  cat(sprintf(" Decisão: %s\n", x$decisao))
  cat("=======================================================\n\n")
}


# ------------------------------------------------------------------------------
# TESTE PARA AMOSTRAS GRANDES (Teorema do Limite Central - Qui-Quadrado)
# ------------------------------------------------------------------------------

#' Executa o Teste Assintótico (TLC) para o Vetor de Médias de 1 Amostra (H0: mu = mu0)
#' Estatística de Teste: X² = n * (xbar - mu0)' * S^(-1) * (xbar - mu0) ~ Chi-Quadrado(p)
#' 
#' @param X Matriz ou data.frame de dados (n x p)
#' @param mu0 Vetor de médias sob a hipótese nula H0 (dimensão p x 1)
#' @param alpha Nível de significância (padrão 0.05)
teste_chi_quadrado <- function(X, mu0, alpha = 0.05) {
  X <- as.matrix(X)
  mu0 <- matrix(mu0, ncol = 1)
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (length(mu0) != p) {
    stop("A dimensão do vetor mu0 deve ser igual ao número de colunas de X.")
  }
  
  xbar <- colMeans(X)
  S <- cov(X)
  
  # Estatística Qui-Quadrado Observada (mesmo valor numérico de T²_obs)
  diff_mu <- xbar - mu0
  chisq_obs <- as.numeric(n * t(diff_mu) %*% solve(S) %*% diff_mu)
  
  # Quantil crítico da distribuição Qui-Quadrado com p graus de liberdade
  chisq_crit <- qchisq(1 - alpha, df = p)
  
  # P-valor via distribuição Qui-Quadrado
  p_valor <- 1 - pchisq(chisq_obs, df = p)
  
  rejeita_H0 <- chisq_obs > chisq_crit
  
  resultado <- list(
    xbar = xbar,
    S = S,
    n = n,
    p = p,
    chisq_obs = chisq_obs,
    chisq_crit = chisq_crit,
    p_valor = p_valor,
    alpha = alpha,
    decisao = ifelse(rejeita_H0, "Rejeita-se H0", "Não se rejeita H0")
  )
  
  class(resultado) <- "teste_chi_quadrado"
  return(resultado)
}

#' Método de exibição para o Teste Assintótico (Qui-Quadrado)
print.teste_chi_quadrado <- function(x, ...) {
  cat("\n=======================================================\n")
  cat("   TESTE QUI-QUADRADO PARA O VETOR DE MÉDIAS (TLC) \n")
  cat("=======================================================\n")
  cat(sprintf(" Tamanho Amostral (n): %d  |  Variáveis (p): %d\n", x$n, x$p))
  cat(sprintf(" Nível de Significância (alpha): %.3f\n", x$alpha))
  cat("-------------------------------------------------------\n")
  cat(sprintf(" Estatística Qui-Quadrado Observada (X²): %.4f\n", x$chisq_obs))
  cat(sprintf(" Quantil Crítico Qui-Quadrado (gl = %d): %.4f\n", x$p, x$chisq_crit))
  cat(sprintf(" p-valor: %.5e\n", x$p_valor))
  cat("-------------------------------------------------------\n")
  cat(sprintf(" Decisão: %s\n", x$decisao))
  cat("=======================================================\n\n")
}

# ------------------------------------------------------------------------------
# VISUALIZAÇÃO DE ELIPSES DE CONFIANÇA MULTIVARIADAS (p = 2)
# Suporta aproximação Exata (F / Hotelling T²) e Assintótica (Qui-Quadrado)
# ------------------------------------------------------------------------------

#' Plot da Elipse de Confiança (p = 2) usando ggplot2
#' 
#' @param X Matriz ou data.frame (n x 2)
#' @param mu0 Vetor numérico de tamanho 2 para a hipótese nula (opcional)
#' @param alpha Nível de significância (default: 0.05)
#' @param distribuicao Tipo de distribuição para a região de confiança: "F" (Hotelling T²) ou "chisq" (Amostra grande / TLC)
#' @param col_elipse Cor da linha da elipse
#' @param col_centro Cor do ponto da média amostral
#' @param main Título do gráfico
#' @return Objeto ggplot editável
plot_elipse_confianca_2d <- function(X, mu0 = NULL, alpha = 0.05, 
                                     distribuicao = c("F", "chisq"),
                                     col_elipse = "steelblue", col_centro = "forestgreen",
                                     main = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ellipse", quietly = TRUE)) {
    stop("Os pacotes 'ggplot2' e 'ellipse' são necessários.")
  }
  
  distribuicao <- match.arg(distribuicao)
  X <- as.matrix(X)
  if (ncol(X) != 2) {
    stop("A função de plotagem de elipse 2D suporta apenas p = 2 variáveis.")
  }
  
  var_names <- colnames(X)
  if (is.null(var_names)) var_names <- c("X1", "X2")
  
  n <- nrow(X)
  p <- ncol(X)
  
  xbar <- colMeans(X)
  S <- cov(X)
  
  # Fator do raio t segundo a distribuição escolhida
  if (distribuicao == "F") {
    c_valor <- quantil_t2_hotelling(alpha, p, n)
    sub_titulo <- sprintf("Nível de Confiança %.0f%% (Hotelling T² / F)", (1 - alpha) * 100)
    if (is.null(main)) main <- "Elipse de Confiança T² de Hotelling"
  } else {
    c_valor <- qchisq(1 - alpha, df = p)
    sub_titulo <- sprintf("Nível de Confiança de %.0f%% (Qui-Quadrado)", (1 - alpha) * 100)
    if (is.null(main)) main <- "Elipse de Confiança Assintótica"
  }
  
  raio_t <- sqrt(c_valor / n)
  
  # Coordenadas da elipse
  elipse_mat <- ellipse::ellipse(S, centre = as.numeric(xbar), t = raio_t)
  
  df_obs <- data.frame(x = X[, 1], y = X[, 2])
  df_elipse <- data.frame(x = elipse_mat[, 1], y = elipse_mat[, 2])
  df_centro <- data.frame(x = xbar[1], y = xbar[2])
  
  cores <- c(
    "Observações" = "gray60",
    "Vetor Médio (xbar)" = col_centro,
    "Elipse de Confiança" = col_elipse
  )
  
  g <- ggplot2::ggplot() +
    ggplot2::geom_point(data = df_obs, ggplot2::aes(x = x, y = y, color = "Observações"), 
                        alpha = 0.6, size = 2) +
    ggplot2::geom_path(data = df_elipse, ggplot2::aes(x = x, y = y, color = "Elipse de Confiança"), 
                       linewidth = 1) +
    ggplot2::geom_point(data = df_centro, ggplot2::aes(x = x, y = y, color = "Vetor Médio (xbar)"), 
                        size = 3.5, shape = 19)
  
  if (!is.null(mu0)) {
    df_mu0 <- data.frame(x = mu0[1], y = mu0[2])
    cores["H0 (mu0)"] <- "firebrick"
    
    g <- g + ggplot2::geom_point(data = df_mu0, ggplot2::aes(x = x, y = y, color = "H0 (mu0)"), 
                                 size = 4, shape = 17)
  }
  
  g <- g +
    ggplot2::scale_color_manual(name = "Elementos", values = cores) +
    ggplot2::labs(
      title = main,
      subtitle = sub_titulo,
      x = var_names[1],
      y = var_names[2]
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray30"),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  return(g)
}