# ==============================================================================
# FICHEIRO: processamento_manova.R
# Processamento de MANOVA Unifatorial via Estatística Lambda de Wilks (Lambda)
# Estrutura: Pipeline Procedural Modularizada
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. MÓDULO DE ENTRADA E PREPARAÇÃO DOS DADOS
# ------------------------------------------------------------------------------
manova_preparar_dados <- function(X, grupo = NULL) {
  if (!is.list(X) || is.data.frame(X)) {
    if (is.null(grupo)) {
      stop("Forneça 'grupo' ou passe 'X' como uma lista de matrizes/grupos.")
    }
    X_list <- lapply(split(as.data.frame(X), as.factor(grupo)), as.matrix)
  } else {
    X_list <- lapply(X, as.matrix)
  }
  
  g <- length(X_list)
  if (g < 2) stop("O teste requer pelo menos 2 grupos (g >= 2).")
  
  p_vec <- sapply(X_list, ncol)
  if (length(unique(p_vec)) > 1) {
    stop("Todos os grupos devem conter o mesmo número de variáveis (colunas).")
  }
  
  p   <- p_vec[1]
  n_k <- sapply(X_list, nrow)
  n   <- sum(n_k)
  
  if (n <= g + p) {
    stop("Tamanho amostral total insuficiente (n <= g + p).")
  }
  
  return(list(
    X_list = X_list,
    n_k    = n_k,
    n      = n,
    p      = p,
    g      = g
  ))
}


# ------------------------------------------------------------------------------
# 2. MÓDULO DE CÁLCULO DAS MATRIZES E ESTATÍSTICA LAMBDA
# ------------------------------------------------------------------------------
manova_calcular_matrizes <- function(X_list, n_k, n, p, g) {
  medias_k    <- do.call(rbind, lapply(X_list, colMeans))
  xbar_global <- colSums(medias_k * n_k) / n
  
  W <- matrix(0, nrow = p, ncol = p)
  B <- matrix(0, nrow = p, ncol = p)
  
  # Loop 1: Matriz Dentro (W) via covariância amostral por grupo
  for (k in 1:g) {
    W <- W + (n_k[k] - 1) * cov(X_list[[k]])
  }
  
  # Loop 2: Matriz Entre (B) via desvio da média do grupo
  for (k in 1:g) {
    diff <- medias_k[k, ] - xbar_global
    B <- B + n_k[k] * (diff %*% t(diff))
  }
  
  T_mat <- B + W
  det_W <- det(W)
  det_T <- det(T_mat)
  
  if (det_T == 0) stop("A matriz total T é singular (determinante nulo).")
  
  lambda <- det_W / det_T
  
  return(list(
    medias_k    = medias_k,
    xbar_global = xbar_global,
    W           = W,
    B           = B,
    T_mat       = T_mat,
    lambda      = lambda
  ))
}


# ------------------------------------------------------------------------------
# 3. MÓDULO DE TESTE DE HIPÓTESE E APROXIMAÇÕES (Exato F, Rao F e Bartlett Chi2)
# ------------------------------------------------------------------------------
manova_testar_hipotese <- function(lambda, n, p, g, alpha, metodo, n_corte_bartlett = 30) {
  q <- g - 1
  
  # Função interna para calcular a aproximação de Bartlett
  calc_bartlett <- function() {
    stat <- -(n - 1 - (p + g) / 2) * log(lambda)
    df1  <- p * q
    df2  <- NA
    p_v  <- pchisq(stat, df = df1, lower.tail = FALSE)
    crit <- qchisq(1 - alpha, df = df1)
    
    list(
      stat       = stat,
      tipo_teste = "Qui-Quadrado (Bartlett)",
      df1        = df1,
      df2        = df2,
      p_valor    = p_v,
      critico    = crit,
      caso       = "Aproximação Assintótica de Bartlett (n grande, Chi-Quadrado)"
    )
  }
  
  # Caso 1: Forçar teste de Bartlett
  if (metodo == "bartlett") {
    res_teste <- calc_bartlett()
    
    # Caso 2: Seleção Automática ou Aproximação de Rao
  } else {
    if (metodo == "auto" && p == 1) {
      caso <- "Exato: p = 1 (ANOVA Univariada)"
      df1  <- g - 1; df2 <- n - g
      stat <- ((1 - lambda) / lambda) * (df2 / df1)
      tipo <- "F"
      
    } else if (metodo == "auto" && p == 2) {
      caso <- "Exato: p = 2"
      df1  <- 2 * (g - 1); df2 <- 2 * (n - g - 1)
      stat <- ((1 - sqrt(lambda)) / sqrt(lambda)) * (df2 / df1)
      tipo <- "F"
      
    } else if (metodo == "auto" && g == 2) {
      caso <- "Exato: g = 2 (Hotelling T^2)"
      df1  <- p; df2 <- n - p - 1
      stat <- ((1 - lambda) / lambda) * (df2 / df1)
      tipo <- "F"
      
    } else if (metodo == "auto" && g == 3) {
      caso <- "Exato: g = 3"
      df1  <- 2 * p; df2 <- 2 * (n - p - 2)
      stat <- ((1 - sqrt(lambda)) / sqrt(lambda)) * (df2 / df1)
      tipo <- "F"
      
    } else if (metodo == "auto" && n >= n_corte_bartlett) {
      # Para amostras grandes (n >= 30) onde nenhum caso exato se aplica
      res_teste <- calc_bartlett()
      res_teste$rejeita_H0 <- res_teste$p_valor < alpha
      return(res_teste)
      
    } else {
      # Aproximação F de Rao (p > 2 e g > 3 para n moderado)
      caso <- "Aproximação F de Rao"
      m    <- n - 1 - (p + g) / 2
      denom_t <- p^2 + q^2 - 5
      t_val   <- if (denom_t > 0) sqrt((p^2 * q^2 - 4) / denom_t) else 1
      
      df1  <- p * q
      df2  <- m * t_val - (p * q - 2) / 2
      stat <- ((1 - lambda^(1 / t_val)) / lambda^(1 / t_val)) * (df2 / df1)
      tipo <- "F"
    }
    
    critico <- qf(1 - alpha, df1 = df1, df2 = df2)
    p_val   <- pf(stat, df1 = df1, df2 = df2, lower.tail = FALSE)
    
    res_teste <- list(
      stat       = stat,
      tipo_teste = tipo,
      df1        = df1,
      df2        = df2,
      p_valor    = p_val,
      critico    = critico,
      caso       = caso
    )
  }
  
  res_teste$rejeita_H0 <- res_teste$p_valor < alpha
  return(res_teste)
}

# ------------------------------------------------------------------------------
# 4. FUNÇÃO PRINCIPAL / ORQUESTRADORA DA PIPELINE
# ------------------------------------------------------------------------------
#' Executa MANOVA Unifatorial usando a estatística Lambda de Wilks
#'
#' @param X Matriz, data.frame ou lista de matrizes por grupo.
#' @param grupo Vetor de fatores/grupos (caso X seja matriz/df).
#' @param alpha Nível de significância (padrão 0.05).
#' @param metodo Método de teste: "auto" "bartlett" ou "Rao".
#' @return Objeto da classe 'manova_wilks'.
manova_wilks <- function(X, grupo = NULL, alpha = 0.05, 
                         metodo = c("auto", "bartlett", "rao")) {
  metodo <- match.arg(metodo)
  
  # Passo 1: Preparação e validação
  dados <- manova_preparar_dados(X, grupo)
  
  # Passo 2: Álgebra matricial
  matrizes <- manova_calcular_matrizes(
    X_list = dados$X_list,
    n_k    = dados$n_k,
    n      = dados$n,
    p      = dados$p,
    g      = dados$g
  )
  
  # Passo 3: Inferência estatística
  teste <- manova_testar_hipotese(
    lambda = matrizes$lambda,
    n      = dados$n,
    p      = dados$p,
    g      = dados$g,
    alpha  = alpha,
    metodo = metodo
  )
  
  # Passo 4: Consolidação do objeto de retorno
  res <- list(
    lambda      = matrizes$lambda,
    stat        = teste$stat,
    tipo_teste  = teste$tipo_teste,
    df1         = teste$df1,
    df2         = teste$df2,
    p_valor     = teste$p_valor,
    critico     = teste$critico,
    rejeita_H0  = teste$rejeita_H0,
    caso        = teste$caso,
    W           = matrizes$W,
    B           = matrizes$B,
    T_mat       = matrizes$T_mat,
    medias_k    = matrizes$medias_k,
    xbar_global = matrizes$xbar_global,
    n           = dados$n,
    g           = dados$g,
    p           = dados$p,
    alpha       = alpha
  )
  
  class(res) <- "manova_wilks"
  return(res)
}


# ------------------------------------------------------------------------------
# 5. MÉTODOS DE IMPRESSÃO 
# ------------------------------------------------------------------------------
print.manova_wilks <- function(x, ...) {
  cat("\n=== TESTE MANOVA: LAMBDA DE WILKS ===\n")
  cat("Caso / Método:", x$caso, "\n")
  cat("Amostra total (n):", x$n, "| Grupos (g):", 
      x$g, "| Variáveis (p):", x$p, "\n\n")
  cat(sprintf("Lambda de Wilks: %.6f\n", x$lambda))
  cat(sprintf("Estatística %s: %.4f\n", x$tipo_teste, x$stat))
  
  if (is.na(x$df2)) {
    cat(sprintf("Graus de Liberdade (Qui-Quadrado): %d\n", x$df1))
  } else {
    cat(sprintf("Graus de Liberdade (F): df1 = %.2f, df2 = %.2f\n", x$df1, x$df2))
  }
  
  cat(sprintf("Valor crítico (alpha = %.2f): %.4f\n", x$alpha, x$critico))
  cat(sprintf("p-valor: %.6e\n", x$p_valor))
  cat("Decisão (H0: Médias iguais):", 
      if (x$rejeita_H0) "REJEITA H0" else "NÃO REJEITA H0", "\n\n")
}
