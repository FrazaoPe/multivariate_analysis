# ==============================================================================
# FICHEIRO: exemplos_johnson_wichern.R
# Aplicação das funções aos Exemplos 6.1 e 6.4 (Johnson & Wichern)
# ==============================================================================

source("processamento_duas_amostras.R")

# ==============================================================================
# EXEMPLO 6.1: Teste Pareado (Effluent Data - Laboratório Comercial vs Estado)
# ==============================================================================
cat("\n#################################################################\n")
cat(" EXEMPLO 6.1: Diferença de Médias com Observações Pareadas (n = 11)\n")
cat("#################################################################\n")

# 1. Entrada dos Dados da Tabela 6.1
# Grupo 1: Commercial Lab (BOD, SS)
X1_comercial <- matrix(c(
  6,  27,
  6,  23,
  18,  64,
  8,  44,
  11,  30,
  34,  75,
  28,  26,
  71, 124,
  43,  54,
  33,  30,
  20,  14
), ncol = 2, byrow = TRUE)
colnames(X1_comercial) <- c("BOD", "SS")

# Grupo 2: State Lab of Hygiene (BOD, SS)
X2_estatal <- matrix(c(
  25,  15,
  28,  13,
  36,  22,
  35,  29,
  15,  31,
  44,  64,
  42,  30,
  54,  64,
  34,  56,
  29,  20,
  39,  21
), ncol = 2, byrow = TRUE)
colnames(X2_estatal) <- c("BOD", "SS")

# 2. Execução do Teste Pareado
res_ex61 <- teste_pareado(
  X1 = X1_comercial, 
  X2 = X2_estatal, 
  delta0 = c(0, 0), 
  alpha = 0.05
)

# Imprimir Resultados do Exemplo 6.1
print(res_ex61)

# 3. Gráfico da Elipse de Confiança para a Diferença das Médias (Exemplo 6.1)
g_ex61 <- plot_elipse_diferenca_2d(
  res_teste = res_ex61, 
  delta0 = c(0, 0), 
  main = "Exemplo 6.1: Elipse de Confiança (95%) para (Commercial - State)"
)
print(g_ex61)

# ==============================================================================
# EXEMPLO 6.4: Duas Amostras Independentes (Consumo de Energia Com/Sem Ar)
# ==============================================================================
cat("\n#################################################################\n")
cat(" EXEMPLO 6.4: Duas Amostras Independentes (n1 = 45, n2 = 55)\n")
cat("#################################################################\n")

# 1. Estatísticas Descritivas Fornecidas
# Grupo 1: Com Ar Condicionado (n1 = 45)
xbar1 <- c(204.4, 556.6)
S1    <- matrix(c(13825.3, 23823.4, 
                  23823.4, 73107.4), nrow = 2, byrow = TRUE)
n1    <- 45

# Grupo 2: Sem Ar Condicionado (n2 = 55)
xbar2 <- c(130.0, 355.0)
S2    <- matrix(c(8632.0, 19616.7, 
                  19616.7, 55964.5), nrow = 2, byrow = TRUE)
n2    <- 55

names(xbar1) <- c("On_Peak", "Off_Peak")
names(xbar2) <- c("On_Peak", "Off_Peak")

# 2. Execução - Suposição de Covariâncias Iguais (Homocedástico)
res_ex64_homo <- teste_independente_homogeneo(
  xbar1 = xbar1, S1 = S1, n1 = n1,
  xbar2 = xbar2, S2 = S2, n2 = n2,
  alpha = 0.05
)

# 3. Execução - Suposição de Covariâncias Diferentes (Heterocedástico)
res_ex64_hetero <- teste_independente_heterogeneo(
  xbar1 = xbar1, S1 = S1, n1 = n1,
  xbar2 = xbar2, S2 = S2, n2 = n2,
  alpha = 0.05
)

# Imprimir Resultados Numéricos
print(res_ex64_homo)
print(res_ex64_hetero)

# ------------------------------------------------------------------------------
# 4. PLOTAGEM DAS ELIPSES DE CONFIANÇA (95%) PARA A DIFERENÇA DE MÉDIAS
# ------------------------------------------------------------------------------

# Elipse 1: Suposição Homocedástica (Sigma1 = Sigma2)
g_ex64_homo <- plot_elipse_diferenca_2d(
  res_teste = res_ex64_homo, 
  delta0 = c(0, 0), 
  main = "Exemplo 6.4: Elipse de Confiança (95%) - Homocedástico"
)
print(g_ex64_homo)

# Elipse 2: Suposição Heterocedástica (Sigma1 != Sigma2 / Nel & Van der Merwe)
g_ex64_hetero <- plot_elipse_diferenca_2d(
  res_teste = res_ex64_hetero, 
  delta0 = c(0, 0), 
  main = "Exemplo 6.4: Elipse de Confiança (95%) - Heterocedástico"
)
print(g_ex64_hetero)

# Opcional: Plot Sobreposto das duas elipses para comparação direta
df_elipse_homo <- data.frame(
  ellipse::ellipse(res_ex64_homo$S, centre = as.numeric(res_ex64_homo$dbar), t = res_ex64_homo$raio_elipse),
  Modelo = "Homocedástico (Sp)"
)
df_elipse_hetero <- data.frame(
  ellipse::ellipse(res_ex64_hetero$S, centre = as.numeric(res_ex64_hetero$dbar), t = res_ex64_hetero$raio_elipse),
  Modelo = "Heterocedástico (S1/n1 + S2/n2)"
)
colnames(df_elipse_homo)[1:2] <- c("x", "y")
colnames(df_elipse_hetero)[1:2] <- c("x", "y")

df_elipses_ambas <- rbind(df_elipse_homo, df_elipse_hetero)

g_ex64_comparativo <- ggplot2::ggplot() +
  ggplot2::geom_path(data = df_elipses_ambas, ggplot2::aes(x = x, y = y, color = Modelo, linetype = Modelo), linewidth = 1) +
  ggplot2::geom_point(aes(x = res_ex64_homo$dbar[1], y = res_ex64_homo$dbar[2]), color = "forestgreen", size = 3.5) +
  ggplot2::geom_point(aes(x = 0, y = 0), color = "firebrick", size = 4, shape = 17) +
  ggplot2::labs(
    title = "Exemplo 6.4: Comparação das Elipses de Confiança (95%)",
    subtitle = "Diferença do Consumo de Energia (Com Ar - Sem Ar)",
    x = "Diferença On-Peak (kWh)",
    y = "Diferença Off-Peak (kWh)"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

print(g_ex64_comparativo)
