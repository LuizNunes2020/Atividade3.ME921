# Carregar bibliotecas
library(pheatmap)
library(igraph)
library(mixer)
library(ergm)
library(network)
library(GGally)
library(gridExtra)
library(sna)
library(rgl)
library(lda)
library(compositions)
library(sbm)
library(blockmodels)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(mclust)
library(nimble)
library(coda)

# Carregar os dados
dados <- read.csv(file.choose(), row.names = 1)
# Converter os dados em uma matriz de adjacência
adj_matrix <- as.matrix(dados)
data_matrix <- as.matrix(dados)


# Calcular a distância entre as linhas e colunas da matriz de adjacência
dist_rows <- dist(adj_matrix)
dist_cols <- dist(t(adj_matrix))

# Realizar o agrupamento hierárquico
hclust_rows <- hclust(dist_rows)
hclust_cols <- hclust(dist_cols)

# Plotar o dendrograma para as linhas
plot(hclust_rows, main = "Dendrograma das Linhas", xlab = "", sub = "", ylab = "Altura")

# Plotar o dendrograma para as colunas
plot(hclust_cols, main = "Dendrograma das Colunas", xlab = "", sub = "", ylab = "Altura")

#Heatmap
pheatmap(adj_matrix, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         display_numbers = FALSE, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

# Ajustar o modelo SBM usando o pacote mixer
set.seed(12345)
fit <- mixer(adj_matrix, qmin = 1, qmax = 10, method = "variational")
# Obter o melhor modelo
mod <- getModel(fit)



# Converter a matriz de adjacência em um objeto de rede
net <- network(adj_matrix, directed=TRUE)

# Gerar o layout da rede usando Fruchterman-Reingold
set.seed(12345)
layout <- network.layout.fruchtermanreingold(net, layout.par=NULL)

z <- t(mod$Taus) #Função map
cluster_membership <- map(z)

color <- brewer.pal(5,"Accent")
ccolor <- color[cluster_membership] #etiquetas de acordo com o cluster

plot(net, displaylabels=TRUE, mode="fruchtermanreingold", coord = layout, vertex.col = ccolor)

# Adicionar a legenda com os tamanhos dos pontos aumentados
legend("topright", col = color, pch = 20, legend = 1:5, title = "Clusters", inset = c(-0.05, 0.2), pt.cex = 2)

# Identificar os jogadores de cada cluster
nnames <- c("ZywOo", "NiKo", "ropz", "m0NESY", "Spinx", "SunPayus", "s1mple", "sh1ro", "stavn", "broky", "device", "frozen", "huNter", "NertZ", "jabbi", "blameF", "Magisk", "cadiaN", "KSCERATO", "Twistzz")
cluster_membership_names <- split(nnames, cluster_membership)

# Imprimir os jogadores de cada cluster
for (i in 1:length(cluster_membership_names)) {
  cat(paste("Cluster", i, ":\n"))
  cat(paste(cluster_membership_names[[i]], collapse=", "), "\n\n")
}




# Definir o modelo em nimble
code <- nimbleCode({
  for (i in 1:N) {
    for (j in 1:N) {
      Y[i, j] ~ dbern(p[i, j])
      logit(p[i, j]) <- theta[Z[i], Z[j]]
    }
  }
  
  for (k in 1:N) {
    Z[k] ~ dcat(tau[1:K])
  }
  
  tau[1:K] ~ ddirch(delta[1:K])
  
  for (g in 1:K) {
    for (h in 1:K) {
      theta[g, h] ~ dbeta(alpha, beta)
    }
  }
})

# Configurar os dados e os parâmetros iniciais
N <- nrow(adj_matrix)
K <- 5 # número de clusters
constants <- list(N = N, K = K, delta = rep(1, K), alpha = 1, beta = 1, Z = sample(1:K, N, replace = TRUE))
data <- list(Y = adj_matrix)
inits <- list(tau = rep(1/K, K), theta = matrix(runif(K * K), nrow = K, ncol = K))

# Criar o modelo NIMBLE
Rmodel <- nimbleModel(code, constants = constants, data = data, inits = inits)

# Configurar o amostrador de Gibbs
mcmcConf <- configureMCMC(Rmodel, monitors = c("Z", "tau", "theta"))
Rmcmc <- buildMCMC(mcmcConf)

# Compilar o modelo e o MCMC
Cmodel <- compileNimble(Rmodel, showCompilerOutput = TRUE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Executar o MCMC
set.seed(123)
samples <- runMCMC(Cmcmc, niter = 10000, nburnin = 2000, thin = 5)

  # Converter os resultados do MCMC para um dataframe
posterior_df <- as.data.frame(as.matrix(samples))

# Selecionar parâmetros específicos para visualização
selected_parameters <- posterior_df[, c("theta[1, 1]", "theta[2, 2]", "theta[3, 3]", "tau[1]", "tau[2]", "tau[3]")]

# Plotar as amostras da distribuição posterior conjunta
pairs(selected_parameters, main = "Distribuição Posterior Conjunta")

# Visualizar as distribuições marginais
melted_df <- melt(selected_parameters)
ggplot(melted_df, aes(x = value)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  facet_wrap(~variable, scales = "free") +
  labs(title = "Distribuições Marginais das Amostras da Distribuição Posterior",
       x = "Valor", y = "Frequência")





g <- graph_from_adjacency_matrix(adj_matrix, mode = "directed", diag = FALSE)
# Converter para um grafo não direcionado
g_undirected <- as.undirected(g, mode = "collapse")

# Detectar comunidades usando o algoritmo Louvain
louvain_clusters <- cluster_louvain(g_undirected)

# Visualizar a rede com clusters do Louvain
V(g_undirected)$cluster_louvain <- factor(membership(louvain_clusters))
plot(g_undirected, vertex.color = V(g_undirected)$cluster_louvain, main = "Clusters Identificados pelo Louvain")
# Adicionar legenda
legend("topright", legend = levels(V(g_undirected)$cluster_louvain), 
       col = rainbow(length(levels(V(g_undirected)$cluster_louvain))), 
       pch = 19, title = "Clusters")
# Comparar os clusters
comparison <- table(sbm_clusters, membership(louvain_clusters))
print(comparison)

# Visualização comparativa
layout(matrix(1:2, nrow = 1))
plot(g, vertex.color = V(g)$cluster_sbm, main = "SBM Clusters")
plot(g_undirected, vertex.color = V(g_undirected)$cluster_louvain, main = "Louvain Clusters")


# Calcular medidas de similaridade
library(clue)
sbm_vs_louvain <- cl_agreement(as.cl_hard_partition(sbm_clusters), as.cl_hard_partition(membership(louvain_clusters)))
cat("Agreement between SBM and Louvain clusters: ", sbm_vs_louvain, "\n")

comparison_df <- as.data.frame(as.table(comparison))
ggplot(comparison_df, aes(x = sbm_clusters, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Valor de Concordância entre Louvain e SBM = 0.4522",x = "Clusters SBM", y = "Clusters Louvain", fill = "Frequency") +
  theme_minimal()



# Identificar os jogadores de cada cluster
player_names <- c("ZywOo", "NiKo", "ropz", "m0NESY", "Spinx", "SunPayus", "s1mple", "sh1ro", "stavn", "broky", 
                  "device", "frozen", "huNter", "NertZ", "jabbi", "blameF", "Magisk", "cadiaN", "KSCERATO", "Twistzz")
vertex_names <- V(g_undirected)$name

# Mapear jogadores para os clusters
player_clusters <- split(player_names, membership(louvain_clusters))

# Imprimir os jogadores de cada cluster
for (i in seq_along(player_clusters)) {
  cat(paste("Cluster", i, ":\n"))
  cat(paste(player_clusters[[i]], collapse = ", "), "\n\n")
}



# Assumindo que 'mod$Pis' contém a matriz \Theta do modelo ajustado
theta_matrix <- mod$Pis

# Mostrar a matriz \Theta
print(theta_matrix)

# Verificar se \Theta é diagonalmente dominante
is_diagonally_dominant <- all(diag(theta_matrix) > rowSums(theta_matrix) - diag(theta_matrix))

# Imprimir o resultado
if (is_diagonally_dominant) {
  cat("A matriz \\Theta é diagonalmente dominante, indicando mistura assortativa.\n")
} else {
  cat("A matriz \\Theta não é diagonalmente dominante, indicando mistura dissortativa.\n")
}