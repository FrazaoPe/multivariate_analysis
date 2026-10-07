# ==============================================================================
# FICHEIRO: exemplo_johnson.R
# Teste da MANOVA via Lambda de Wilks - Exemplo 6.9 (Johnson & Wichern)
# ==============================================================================

# Carrega o módulo de processamento MANOVA
source("processamento_MANOVA.R")

# 1. Construção das matrizes de observações por grupo
# Grupo 1 (n1 = 3)
grupo1 <- matrix(c(9, 3,
                   6, 2,
                   9, 7), ncol = 2, byrow = TRUE)

# Grupo 2 (n2 = 2)
grupo2 <- matrix(c(0, 4,
                   2, 0), ncol = 2, byrow = TRUE)

# Grupo 3 (n3 = 3)
grupo3 <- matrix(c(3, 8,
                   1, 9,
                   2, 7), ncol = 2, byrow = TRUE)

# ------------------------------------------------------------------------------
# Execução 1: Entrada como Lista de Grupos Desbalanceados
# ------------------------------------------------------------------------------
dados_lista <- list(Grupo1 = grupo1, Grupo2 = grupo2, Grupo3 = grupo3)

res_lista <- manova_wilks(X = dados_lista, alpha = 0.05, metodo = "auto")
print(res_lista)

# ------------------------------------------------------------------------------
# Execução 2: Entrada Tradicional (Matriz X + Vetor de Grupos)
# ------------------------------------------------------------------------------
X_mat <- rbind(grupo1, grupo2, grupo3)
fator_grupos <- factor(rep(c("G1", "G2", "G3"), times = c(3, 2, 3)))

res_matriz <- manova_wilks(X = X_mat, grupo = fator_grupos, alpha = 0.05, metodo = "auto")

# ------------------------------------------------------------------------------
# Exibição das Matrizes Intermediárias do Livro
# ------------------------------------------------------------------------------
cat("\n--- DETALHAMENTO MATRICIAL (Johnson & Wichern) ---\n")

cat("\nMédias dos Grupos:\n")
print(res_lista$medias_k)

cat("\nMédia Global:\n")
print(res_lista$xbar_global)

cat("\nMatriz Dentro dos Grupos (W):\n")
print(res_lista$W)

cat("\nMatriz Entre Grupos (B):\n")
print(res_lista$B)

cat("\nMatriz Total (T = B + W):\n")
print(res_lista$T_mat)
