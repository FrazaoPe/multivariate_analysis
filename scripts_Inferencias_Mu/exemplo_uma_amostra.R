# ==============================================================================
# FICHEIRO: exemplo.R
# Aplicação Prática dos Dados de Radiação da Tabela 4.1 (Johnson & Wichern)
# ==============================================================================

options(scipen = 999)

# 1. Carregar as funções do módulo procedural
source("processamento.R")

# 2. Importação e Estruturação dos Dados da Tabela 4.1 (n = 42)
t4_1 <- c(
  0.15, 0.09, 0.18, 0.10, 0.05, 0.12, 0.08, 0.05, 0.08, 0.10, 0.07, 0.02, 0.01, 0.10, 0.10,
  0.10, 0.02, 0.10, 0.01, 0.40, 0.10, 0.05, 0.03, 0.05, 0.15, 0.10, 0.15, 0.09, 0.08, 0.18,
  0.10, 0.20, 0.11, 0.30, 0.02, 0.20, 0.20, 0.30, 0.30, 0.40, 0.30, 0.05
)

# Data frame com os dados brutos
dados_radiacao <- data.frame(
  Forno = 1:42,
  Radiacao_Fechada = t4_1
)

# 3. Transformação Quarta Raiz (x^(1/4))
# Conforme o Exemplo 5.3 de Johnson & Wichern (estabiliza variância e aproxima a normalidade)
x1_transf <- dados_radiacao$Radiacao_Fechada^(1/4)

set.seed(42)
x2_transf <- x1_transf + rnorm(42, mean = 0.025, sd = 0.04)

X <- cbind(
  Radiacao_Fechada = x1_transf,
  Radiacao_Aberta  = x2_transf
)

# 4. Estatísticas Descritivas
cat("\n--- Vetor de Médias Amostrais (xbar) ---\n")
xbar <- colMeans(X)
print(round(xbar, 4))

cat("\n--- Matriz de Covariância Amostral (S) ---\n")
S <- cov(X)
print(round(S, 6))

# 5. Avaliação da Normalidade Multivariada (Pacote leve: mvnormtest)
res_norm <- teste_normalidade_multivariada(X)
cat("\n--- Teste de Normalidade Multivariada (Shapiro-Wilk) ---\n")
print(res_norm)

# ==============================================================================
# 6. INFERÊNCIA MULTIVARIADA (H0: mu = mu0)
# ==============================================================================
# Vetor de médias hipotético (norma ou padrão de referência)
mu0_hipotese <- c(0.562, 0.589)
alpha_nivel  <- 0.05

# 6.1 Teste Exato: T² de Hotelling (Distribuição F)
cat("\n--- 1. TESTE EXATO: T² DE HOTELLING ---\n")
res_hotelling <- teste_hotelling_1amostra(X = X, mu0 = mu0_hipotese, alpha = alpha_nivel)
print(res_hotelling)

# 6.2 Teste Assintótico: Amostra Grande (Distribuição Qui-Quadrado / TLC)
cat("\n--- 2. TESTE ASSINTÓTICO: AMOSTRA GRANDE (TLC) ---\n")
res_tlc <- teste_chi_quadrado(X = X, mu0 = mu0_hipotese, alpha = alpha_nivel)
print(res_tlc)

# ==============================================================================
# 7. VISUALIZAÇÃO DAS ELIPSES DE CONFIANÇA
# ==============================================================================

# 7.1 Elipse de Confiança Exata (Distribuição F / Hotelling T²)
grafico_elipse_f <- plot_elipse_confianca_2d(
  X = X, 
  mu0 = mu0_hipotese, 
  alpha = alpha_nivel,
  distribuicao = "F"
)

# 7.2 Elipse de Confiança Assintótica (Distribuição Qui-Quadrado / TLC)
grafico_elipse_chisq <- plot_elipse_confianca_2d(
  X = X, 
  mu0 = mu0_hipotese, 
  alpha = alpha_nivel,
  distribuicao = "chisq"
)

# Exibir os gráficos
print(grafico_elipse_f)
print(grafico_elipse_chisq)
