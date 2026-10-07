# ==============================================================================
# FICHEIRO: processamento_duas_amostras.R
# Módulo de Funções para Inferência Multivariada sobre Duas Amostras
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. AMOSTRAS PAREADAS (Teste T² de Hotelling Pareado)
# ------------------------------------------------------------------------------

#' Teste T² de Hotelling para Amostras Pareadas
#'
#' @param X1 Matriz/data.frame do grupo 1 (n x p) ou vetor de dados
#' @param X2 Matriz/data.frame do grupo 2 (n x p) ou vetor de dados
#' @param delta0 Vetor hipotético de diferenças H0 (default: rep(0, p))
#' @param alpha Nível de significância (default: 0.05)
#' @param amostra_grande TRUE para aproximação assintótica Qui-Quadrado (n - p grande)
teste_pareadoT2 <- function(X1, X2, delta0 = NULL, alpha = 0.05, amostra_grande = FALSE) {
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  if (!all(dim(X1) == dim(X2))) {
    stop("Para testes pareados, X1 e X2 devem ter as mesmas dimensões (n x p).")
  }
  
  n <- nrow(X1)
  p <- ncol(X1)
  
  var_names <- colnames(X1)
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  
  if (is.null(delta0)) delta0 <- rep(0, p)
  delta0 <- matrix(delta0, ncol = 1)
  
  # Matriz de diferenças D = X1 - X2
  D <- X1 - X2
  dbar <- colMeans(D)
  names(dbar) <- var_names
  Sd <- cov(D)
  
  diff_vec <- matrix(dbar, ncol = 1) - delta0
  T2_obs <- as.numeric(n * t(diff_vec) %*% solve(Sd) %*% diff_vec)
  
  if (!amostra_grande) {
    # Caso 1: Pareado, n - p pequeno/moderado (Exato via F)
    tipo_str <- "Pareado (Exato - F)"
    c_crit <- (p * (n - 1) / (n - p)) * qf(1 - alpha, df1 = p, df2 = n - p)
    F_obs <- T2_obs * (n - p) / (p * (n - 1))
    p_valor <- 1 - pf(F_obs, df1 = p, df2 = n - p)
    stat_test <- T2_obs
  } else {
    # Caso 2: Pareado, n - p grande (Assintótico via Qui-Quadrado)
    tipo_str <- "Pareado (Assintótico - X²)"
    c_crit <- qchisq(1 - alpha, df = p)
    F_obs <- NA
    p_valor <- 1 - pchisq(T2_obs, df = p)
    stat_test <- T2_obs
  }
  
  resultado <- list(
    tipo = tipo_str,
    dbar = dbar,
    S = Sd,
    n1 = n, n2 = n, p = p,
    T2_obs = stat_test,
    chisq_obs = if (amostra_grande) stat_test else NULL,
    c_crit = c_crit,
    F_obs = F_obs,
    p_valor = p_valor,
    alpha = alpha,
    amostra_grande = amostra_grande,
    raio_elipse = sqrt(c_crit / n),
    decisao = ifelse(stat_test > c_crit, "Rejeita-se H0", "Não se rejeita H0")
  )
  class(resultado) <- c("teste_pareado", "teste_duas_amostras")
  return(resultado)
}

# Alias para facilitar o uso sem a sigla T2 no final
teste_pareado <- teste_pareadoT2


# ------------------------------------------------------------------------------
# 2. AMOSTRAS INDEPENDENTES - COVARIÂNCIAS IGUAIS
# ------------------------------------------------------------------------------

#' Teste T² para Amostras Independentes com Covariâncias Iguais (Homocedástico)
teste_independente_homogeneo <- function(X1 = NULL, X2 = NULL, 
                                         xbar1 = NULL, S1 = NULL, n1 = NULL,
                                         xbar2 = NULL, S2 = NULL, n2 = NULL,
                                         delta0 = NULL, alpha = 0.05,
                                         amostra_grande = FALSE) {
  if (!is.null(X1) && !is.null(X2)) {
    X1 <- as.matrix(X1); X2 <- as.matrix(X2)
    n1 <- nrow(X1); n2 <- nrow(X2); p <- ncol(X1)
    xbar1 <- colMeans(X1); xbar2 <- colMeans(X2)
    S1 <- cov(X1); S2 <- cov(X2)
    var_names <- colnames(X1)
  } else {
    p <- length(xbar1)
    xbar1 <- as.numeric(xbar1); xbar2 <- as.numeric(xbar2)
    S1 <- as.matrix(S1); S2 <- as.matrix(S2)
    var_names <- names(xbar1)
  }
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  names(xbar1) <- var_names; names(xbar2) <- var_names
  
  if (is.null(delta0)) delta0 <- rep(0, p)
  delta0 <- matrix(delta0, ncol = 1)
  
  v <- n1 + n2 - 2
  Sp <- ((n1 - 1) * S1 + (n2 - 1) * S2) / v
  dbar <- xbar1 - xbar2
  names(dbar) <- var_names
  
  diff_vec <- matrix(dbar, ncol = 1) - delta0
  n12 <- (n1 * n2) / (n1 + n2)
  
  T2_obs <- as.numeric(n12 * t(diff_vec) %*% solve(Sp) %*% diff_vec)
  
  if (!amostra_grande) {
    # Caso 3: Independentes Sigma1 = Sigma2, n - p pequeno/moderado (Exato via F)
    tipo_str <- "Independentes Homocedástico (Exato - F)"
    c_crit <- (p * v / (v - p + 1)) * qf(1 - alpha, df1 = p, df2 = v - p + 1)
    F_obs <- T2_obs * (v - p + 1) / (p * v)
    p_valor <- 1 - pf(F_obs, df1 = p, df2 = v - p + 1)
    stat_test <- T2_obs
  } else {
    # Caso 4: Independentes Sigma1 = Sigma2, n - p grande (Assintótico via X²)
    tipo_str <- "Independentes Homocedástico (Assintótico - X²)"
    c_crit <- qchisq(1 - alpha, df = p)
    F_obs <- NA
    p_valor <- 1 - pchisq(T2_obs, df = p)
    stat_test <- T2_obs
  }
  
  resultado <- list(
    tipo = tipo_str,
    dbar = dbar,
    S = Sp,
    n1 = n1, n2 = n2, p = p, n12 = n12,
    T2_obs = stat_test,
    chisq_obs = if (amostra_grande) stat_test else NULL,
    c_crit = c_crit,
    F_obs = F_obs,
    p_valor = p_valor,
    alpha = alpha,
    amostra_grande = amostra_grande,
    raio_elipse = sqrt(c_crit / n12),
    decisao = ifelse(stat_test > c_crit, "Rejeita-se H0", "Não se rejeita H0")
  )
  class(resultado) <- c("teste_independente_homogeneo", "teste_duas_amostras")
  return(resultado)
}


# ------------------------------------------------------------------------------
# 3. AMOSTRAS INDEPENDENTES - COVARIÂNCIAS DIFERENTES
# 
# ------------------------------------------------------------------------------

#' Teste para Amostras Independentes Heterocedásticas (Sigma1 != Sigma2)
teste_independente_heterogeneo <- function(X1 = NULL, X2 = NULL,
                                           xbar1 = NULL, S1 = NULL, n1 = NULL,
                                           xbar2 = NULL, S2 = NULL, n2 = NULL,
                                           delta0 = NULL, alpha = 0.05,
                                           amostra_grande = FALSE) {
  if (!is.null(X1) && !is.null(X2)) {
    X1 <- as.matrix(X1); X2 <- as.matrix(X2)
    n1 <- nrow(X1); n2 <- nrow(X2); p <- ncol(X1)
    xbar1 <- colMeans(X1); xbar2 <- colMeans(X2)
    S1 <- cov(X1); S2 <- cov(X2)
    var_names <- colnames(X1)
  } else {
    p <- length(xbar1)
    xbar1 <- as.numeric(xbar1); xbar2 <- as.numeric(xbar2)
    S1 <- as.matrix(S1); S2 <- as.matrix(S2)
    var_names <- names(xbar1)
  }
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  names(xbar1) <- var_names; names(xbar2) <- var_names
  
  if (is.null(delta0)) delta0 <- rep(0, p)
  delta0 <- matrix(delta0, ncol = 1)
  
  # Matriz de covariância da diferença de médias: S_dbar = S1/n1 + S2/n2
  S_dbar <- (S1 / n1) + (S2 / n2)
  S_dbar_inv <- solve(S_dbar)
  
  dbar <- xbar1 - xbar2
  names(dbar) <- var_names
  diff_vec <- matrix(dbar, ncol = 1) - delta0
  
  T2_obs <- as.numeric(t(diff_vec) %*% S_dbar_inv %*% diff_vec)
  
  if (!amostra_grande) {
    # Caso 5: Heterocedástico, Normalidade / Nel & Van der Merwe (Johnson & Wichern)
    tipo_str <- "Independentes Heterocedástico (Nel & Van der Merwe / F)"
    
    # Matrizes M1 e M2 para o cálculo dos traços
    M1 <- (1 / n1) * S1 %*% S_dbar_inv
    M2 <- (1 / n2) * S2 %*% S_dbar_inv
    
    tr_M1_sq <- sum(diag(M1 %*% M1))
    tr_M1    <- sum(diag(M1))
    tr_M2_sq <- sum(diag(M2 %*% M2))
    tr_M2    <- sum(diag(M2))
    
    # Denominador do cálculo dos graus de liberdade nu
    denom <- (1 / (n1 - 1)) * (tr_M1_sq + tr_M1^2) + (1 / (n2 - 1)) * (tr_M2_sq + tr_M2^2)
    nu <- (p + p^2) / denom
    
    # Valor crítico F e estatística F observada
    c_crit <- (nu * p / (nu - p + 1)) * qf(1 - alpha, df1 = p, df2 = nu - p + 1)
    F_obs  <- T2_obs * (nu - p + 1) / (nu * p)
    p_valor <- 1 - pf(F_obs, df1 = p, df2 = nu - p + 1)
    gl_info <- nu
  } else {
    # Caso 6: Heterocedástico, n1 e n2 grandes (Assintótico TLC / X²)
    tipo_str <- "Independentes Heterocedástico (Assintótico TLC - X²)"
    c_crit <- qchisq(1 - alpha, df = p)
    F_obs  <- NA
    p_valor <- 1 - pchisq(T2_obs, df = p)
    gl_info <- NULL
  }
  
  resultado <- list(
    tipo = tipo_str,
    dbar = dbar,
    S = S_dbar,
    n1 = n1, n2 = n2, p = p,
    chisq_obs = T2_obs,
    T2_obs = T2_obs,
    nu = gl_info,
    c_crit = c_crit,
    F_obs = F_obs,
    p_valor = p_valor,
    alpha = alpha,
    amostra_grande = amostra_grande,
    raio_elipse = sqrt(c_crit),
    decisao = ifelse(T2_obs > c_crit, "Rejeita-se H0", "Não se rejeita H0")
  )
  class(resultado) <- c("teste_independente_heterogeneo", "teste_duas_amostras")
  return(resultado)
}


# ------------------------------------------------------------------------------
# 4. IMPRESSÕES FORMATADAS (MÉTODOS S3)
# ------------------------------------------------------------------------------

print.teste_pareado <- function(x, ...) {
  cat("\n=======================================================\n")
  cat("     TESTE T² DE HOTELLING PARA AMOSTRAS PAREADAS      \n")
  cat("=======================================================\n")
  cat(sprintf(" Tipo: %s\n", x$tipo))
  cat(sprintf(" Tamanho da Amostra (n): %d  |  Variáveis (p): %d\n", x$n1, x$p))
  cat(sprintf(" Nível de Significância (alpha): %.3f\n", x$alpha))
  cat("-------------------------------------------------------\n")
  cat(" Vetor de Diferenças Médias (dbar):\n")
  print(round(x$dbar, 4))
  cat("\n Matriz de Covariância das Diferenças (Sd):\n")
  print(round(x$S, 4))
  cat("-------------------------------------------------------\n")
  cat(sprintf(" Estatística Observada: %.4f\n", x$T2_obs))
  cat(sprintf(" Valor Crítico: %.4f\n", x$c_crit))
  if (!is.na(x$F_obs)) {
    cat(sprintf(" Estatística F Equivalente: %.4f (gl1 = %d, gl2 = %d)\n", 
                x$F_obs, x$p, x$n1 - x$p))
  }
  cat(sprintf(" p-valor: %.5e\n", x$p_valor))
  cat("-------------------------------------------------------\n")
  cat(sprintf(" Decisão: %s\n", x$decisao))
  cat("=======================================================\n\n")
}

print.teste_independente_homogeneo <- function(x, ...) {
  v <- x$n1 + x$n2 - 2
  gl2 <- v - x$p + 1
  
  cat("\n=======================================================\n")
  cat("   TESTE T² PARA AMOSTRAS INDEPENDENTES (HOMOCEDÁSTICO)\n")
  cat("=======================================================\n")
  cat(sprintf(" Tipo: %s\n", x$tipo))
  cat(sprintf(" n1: %d  |  n2: %d  |  n_efetivo (n12): %.2f  |  p: %d\n", 
              x$n1, x$n2, x$n12, x$p))
  cat(sprintf(" Nível de Significância (alpha): %.3f\n", x$alpha))
  cat("-------------------------------------------------------\n")
  cat(" Vetor de Diferença de Médias (xbar1 - xbar2):\n")
  print(round(x$dbar, 4))
  cat("\n Matriz de Covariância Agrupada (Sp):\n")
  print(round(x$S, 4))
  cat("-------------------------------------------------------\n")
  cat(sprintf(" Estatística Observada: %.4f\n", x$T2_obs))
  cat(sprintf(" Valor Crítico: %.4f\n", x$c_crit))
  if (!is.na(x$F_obs)) {
    cat(sprintf(" Estatística F Equivalente: %.4f (gl1 = %d, gl2 = %d)\n", 
                x$F_obs, x$p, gl2))
  }
  cat(sprintf(" p-valor: %.5e\n", x$p_valor))
  cat("-------------------------------------------------------\n")
  cat(sprintf(" Decisão: %s\n", x$decisao))
  cat("=======================================================\n\n")
}

print.teste_independente_heterogeneo <- function(x, ...) {
  cat("\n=======================================================\n")
  cat("   TESTE PARA DUAS AMOSTRAS INDEPENDENTES (HETEROCEDÁSTICO)\n")
  cat("=======================================================\n")
  cat(sprintf(" Tipo: %s\n", x$tipo))
  cat(sprintf(" n1: %d  |  n2: %d  |  Variáveis (p): %d\n", x$n1, x$n2, x$p))
  if (!is.null(x$nu)) {
    cat(sprintf(" Graus de Liberdade Aproximados (nu): %.4f\n", x$nu))
  }
  cat(sprintf(" Nível de Significância (alpha): %.3f\n", x$alpha))
  cat("-------------------------------------------------------\n")
  cat(" Vetor de Diferença de Médias (xbar1 - xbar2):\n")
  print(round(x$dbar, 4))
  cat("\n Matriz de Covariância da Diferença (S1/n1 + S2/n2):\n")
  print(round(x$S, 4))
  cat("-------------------------------------------------------\n")
  cat(sprintf(" Estatística Observada: %.4f\n", x$T2_obs))
  cat(sprintf(" Valor Crítico: %.4f\n", x$c_crit))
  if (!is.na(x$F_obs)) {
    cat(sprintf(" Estatística F Equivalente: %.4f (gl1 = %d, gl2 = %.2f)\n", 
                x$F_obs, x$p, x$nu - x$p + 1))
  }
  cat(sprintf(" p-valor: %.5e\n", x$p_valor))
  cat("-------------------------------------------------------\n")
  cat(sprintf(" Decisão: %s\n", x$decisao))
  cat("=======================================================\n\n")
}

print.teste_duas_amostras <- function(x, ...) {
  cat("\n=======================================================\n")
  cat(sprintf("   TESTE MULTIVARIADO PARA DUAS AMOSTRAS\n   Tipo: %s\n", x$tipo))
  cat("=======================================================\n")
  cat(sprintf(" n1: %d | n2: %d | Variáveis (p): %d\n", x$n1, x$n2, x$p))
  cat(sprintf(" Nível de Significância (alpha): %.3f\n", x$alpha))
  cat("-------------------------------------------------------\n")
  
  est_val <- if (!is.null(x$T2_obs)) x$T2_obs else x$chisq_obs
  cat(sprintf(" Estatística Observada: %.4f\n", est_val))
  cat(sprintf(" Valor Crítico: %.4f\n", x$c_crit))
  cat(sprintf(" p-valor: %.5e\n", x$p_valor))
  cat("-------------------------------------------------------\n")
  cat(sprintf(" Decisão: %s\n", x$decisao))
  cat("=======================================================\n\n")
}


# ------------------------------------------------------------------------------
# 5. TABELA RESUMO COMPARATIVA
# ------------------------------------------------------------------------------

#' Gera uma Tabela Resumo com os Resultados dos Testes Executados
resumo_testes <- function(...) {
  testes <- list(...)
  
  df_resumo <- data.frame(
    Teste = character(),
    Tipo = character(),
    n1 = integer(),
    n2 = integer(),
    p = integer(),
    Estatistica = numeric(),
    Valor_Critico = numeric(),
    p_valor = character(),
    Decisao = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(testes)) {
    t <- testes[[i]]
    
    is_chisq <- !is.null(t$chisq_obs) && isTRUE(t$amostra_grande)
    est_val  <- if (is_chisq) t$chisq_obs else t$T2_obs
    
    nome_teste <- if (!is.null(names(testes)[i]) && names(testes)[i] != "") {
      names(testes)[i]
    } else {
      paste("Teste", i)
    }
    
    df_resumo <- rbind(df_resumo, data.frame(
      Teste = nome_teste,
      Tipo = t$tipo,
      n1 = t$n1,
      n2 = ifelse(is.null(t$n2), t$n1, t$n2),
      p = t$p,
      Estatistica = round(est_val, 4),
      Valor_Critico = round(t$c_crit, 4),
      p_valor = formatC(t$p_valor, format = "e", digits = 4),
      Decisao = t$decisao,
      stringsAsFactors = FALSE
    ))
  }
  
  return(df_resumo)
}


# ------------------------------------------------------------------------------
# 6. PLOTAGEM DE ELIPSES DE CONFIANÇA PARA A DIFERENÇA (p = 2)
# ------------------------------------------------------------------------------

#' Plot da Elipse de Confiança para o Vetor de Diferenças de Médias
plot_elipse_diferenca_2d <- function(res_teste, delta0 = c(0, 0), L = NULL, main = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ellipse", quietly = TRUE)) {
    stop("Os pacotes 'ggplot2' e 'ellipse' são necessários.")
  }
  
  if (res_teste$p != 2) {
    stop("A função de plotagem de elipse suporta apenas p = 2 variáveis.")
  }
  
  dbar <- res_teste$dbar
  S <- res_teste$S
  raio <- res_teste$raio_elipse
  
  elipse_mat <- ellipse::ellipse(S, centre = as.numeric(dbar), t = raio)
  df_elipse <- data.frame(x = elipse_mat[, 1], y = elipse_mat[, 2])
  df_centro <- data.frame(x = dbar[1], y = dbar[2])
  df_h0     <- data.frame(x = delta0[1], y = delta0[2])
  
  if (is.null(main)) {
    main <- sprintf("Elipse de Região Crítica (%.0f%%) para a Diferença de Médias", (1 - res_teste$alpha) * 100)
  }
  
  cores <- c(
    "Diferença Amostral (dbar)" = "forestgreen",
    "Elipse de Região Crítica"   = "steelblue",
    "Hipótese Nula (delta0)"    = "firebrick"
  )
  
  g <- ggplot2::ggplot() +
    ggplot2::geom_path(data = df_elipse, ggplot2::aes(x = x, y = y, color = "Elipse de Região Crítica"), linewidth = 1) +
    ggplot2::geom_point(data = df_centro, ggplot2::aes(x = x, y = y, color = "Diferença Amostral (dbar)"), size = 3.5, shape = 19) +
    ggplot2::geom_point(data = df_h0, ggplot2::aes(x = x, y = y, color = "Hipótese Nula (delta0)"), size = 4, shape = 17) +
    ggplot2::scale_color_manual(name = "Elementos", values = cores) +
    ggplot2::labs(
      title = main,
      subtitle = sprintf("Tipo: %s", res_teste$tipo),
      x = ifelse(is.null(names(dbar)[1]), "d1", names(dbar)[1]),
      y = ifelse(is.null(names(dbar)[2]), "d2", names(dbar)[2])
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray30"),
      legend.position = "right"
    )
  
  if (!is.null(L)) {
    g <- g + ggplot2::xlim(dbar[1] - L, dbar[1] + L) + ggplot2::ylim(dbar[2] - L, dbar[2] + L)
  }
  
  return(g)
}