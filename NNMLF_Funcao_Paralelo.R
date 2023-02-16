# Nearest Neighbor Method for Linear Features (NNMLF) - En
# Método do Vizinho Mais Próximo para Feições Lineares (VMPFL) - Br

options(digits = 20)

# Carregar bibliotecas necessárias

library(raster)
library(sf)
library(magrittr)
library(dplyr)
library(maptools)
library(rgdal)
library(rgeos)
library(foreach)
library(doParallel)

NNMLF <- function(linhas_entrada, A){  # Entrar com nível de confiança? Default = 95%
  # Criando campo ID_Unico
  linhas_entrada$Id_Unico <- as.integer(seq(1, dim(linhas_entrada@data)[1], by = 1))
  
  # Extrair os vértices das linhas de entrada
  ## Linhas e pontos possuem o mesmo "Id" com esse comando. OBS: Não é o campo único
  vertices_linhas_entrada <- as(linhas_entrada, "SpatialPointsDataFrame")
  rownames(vertices_linhas_entrada@data) <- seq(from = 1, to = length(vertices_linhas_entrada))
  
  # Preparando um data frame para receber as distâncias
  tabela_dist <- data.frame(matrix(NA_real_, nrow = dim(linhas_entrada@data)[1], ncol = dim(linhas_entrada@data)[1]))
  names(tabela_dist)  <- seq.int(1, dim(linhas_entrada)[1], 1)
  row.names(tabela_dist)  <- seq.int(1, dim(linhas_entrada)[1], 1)
  
  # Checa quantos núcleos existem
  ncl<-detectCores()
  ncl
  
  # Registra os clusters a serem utilizados
  cl <- makeCluster(ncl)
  registerDoParallel(cl)

  # CALCULANDO A DISTANCIA - UTILIZANDO PROGRAMACAO EM PARALELO
  tabela_dist <- foreach (i = 1:dim(linhas_entrada)[1])  %:% # Pontos
    foreach (j = 1:dim(linhas_entrada)[1], .combine = 'c') %dopar%{ # Linhas
      if (i != j){ # Otimizar, para não calcular distâncias desnecessárias
        mean(rgeos::gDistance(vertices_linhas_entrada[vertices_linhas_entrada@data$Id_Unico == i,], linhas_entrada[linhas_entrada@data$Id_Unico == j,], byid = TRUE))
      } else {
        0 # Distância dos pontos para a mesma linha de origem
      }
    }
  
  tabela_dist <- as.data.frame(tabela_dist)
  colnames(tabela_dist) <- seq(1, dim(tabela_dist)[2], by = 1)
  
  # Calcular a Distância de Hausdorff
  Dist_Hausdorff <- data.frame(matrix(NA_real_, nrow = dim(linhas_entrada@data)[1], ncol = dim(linhas_entrada@data)[1]))
  names(Dist_Hausdorff)  <- seq.int(1,dim(linhas_entrada)[1],1)
  row.names(Dist_Hausdorff)  <- seq.int(1,dim(linhas_entrada)[1],1)
  
  for (i in 1:dim(linhas_entrada)[1]) { # Pontos
    for (j in 1:dim(linhas_entrada)[1]) { # Linhas
      if (i != j){
        if (tabela_dist[i,j] < tabela_dist[j,i]){
          Dist_Hausdorff[i,j] <- tabela_dist[j,i]
        } else {
          Dist_Hausdorff[i,j] <- tabela_dist[i,j]
        }
      } else {
        Dist_Hausdorff[i,j] <- 0
      }
    }
  }
  
  NNMLF <- foreach (i = 1:dim(linhas_entrada)[1], .combine = 'rbind') %dopar%{
    Dist_Hausdorff[i, order(Dist_Hausdorff[i,], decreasing = FALSE)[2]]
  }
  NNMLF_2 <- foreach (i = 1:dim(linhas_entrada)[1], .combine = 'rbind') %dopar%{
    Dist_Hausdorff[i, order(Dist_Hausdorff[i,], decreasing = FALSE)[3]]
  }
  NNMLF_3 <- foreach (i = 1:dim(linhas_entrada)[1], .combine = 'rbind') %dopar%{
    Dist_Hausdorff[i, order(Dist_Hausdorff[i,], decreasing = FALSE)[4]]
  }
  
  # R Observado   ######## NÃO PRECISA ALTERAR ######## 
  R_OBS <- sum(NNMLF)/dim(linhas_entrada@data)[1]
  R_OBS_2 <- sum(NNMLF_2)/dim(linhas_entrada@data)[1]
  R_OBS_3 <- sum(NNMLF_3)/dim(linhas_entrada@data)[1]
  
  # R Esperado
  Gamma_1 <- c(0.5000, 0.7500, 0.9375, 1.0937, 1.2305, 1.3535)
  Gamma_2 <- c(0.2613, 0.2722, 0.2757, 0.2775, 0.2784, 0.2789)
  Gamma <- data.frame(Gamma_1, Gamma_2) 
  
  R_ESP <- Gamma[1,1]*sqrt(A/dim(linhas_entrada@data)[1])
  R_ESP_2 <- Gamma[2,1]*sqrt(A/dim(linhas_entrada@data)[1])
  R_ESP_3 <- Gamma[3,1]*sqrt(A/dim(linhas_entrada@data)[1])
  
  # Calcular a estatística R
  R <- R_OBS/R_ESP
  R_2 <- R_OBS_2/R_ESP_2
  R_3 <- R_OBS_3/R_ESP_3
  
  # Aplicar o teste Z
  SE <- Gamma[1,2]*sqrt(A/(dim(linhas_entrada@data)[1])^2)
  SE_2 <- Gamma[2,2]*sqrt(A/(dim(linhas_entrada@data)[1])^2)
  SE_3 <- Gamma[3,2]*sqrt(A/(dim(linhas_entrada@data)[1])^2)
  
  Z <- (R_OBS - R_ESP)/SE
  Z_2 <- (R_OBS_2 - R_ESP_2)/SE_2
  Z_3 <- (R_OBS_3 - R_ESP_3)/SE_3
  
  # Considerando um nível de confiança de 95% (nível de significância de 5%), temos:
  Z_tab <- 1.96
  
  # 1 - AGRUPADO
  # 2 - ALEATORIO
  # 3 - DISPERSO
  # NA - NÚMERO DE LINHAS INSUFICIENTE PARA A ORDEM EM QUESTÃO
  
  # Inferir sobre o padrão de distribuição espacial das linhas
  # Primeira Ordem
  if (dim(linhas_entrada@data)[1] > 1){
    if (abs(Z) < Z_tab){ # Hipotese Nula não rejeitada
      indice_NNMLF <- 2
    } else {
      if (R < 1){
        indice_NNMLF <- 1
      } else {
        indice_NNMLF <- 3
      }
    }
  } else {
    indice_NNMLF <- NA
  }
  
  # Segunda Ordem
  if (dim(linhas_entrada@data)[1] > 2){
    if (abs(Z_2) < Z_tab){ # Hipotese Nula não rejeitada
      indice_NNMLF_2 <- 2
    } else {
      if (R_2 < 1){
        indice_NNMLF_2 <- 1
      } else {
        indice_NNMLF_2 <- 3
      }
    }
  } else {
    indice_NNMLF_2 <- NA
  }
  
  # Terceira Ordem
  if (dim(linhas_entrada@data)[1] > 3){
    if (abs(Z_3) < Z_tab){ # Hipotese Nula não rejeitada
      indice_NNMLF_3 <- 2
    } else {
      if (R_3 < 1){
        indice_NNMLF_3 <- 1
      } else {
        indice_NNMLF_3 <- 3
      }
    }
  } else {
    indice_NNMLF_3 <- NA
  }
  
  # Stop clusters
  stopCluster(cl)
  
  saida_NNMLF <- cbind(indice_NNMLF, indice_NNMLF_2, indice_NNMLF_3)
  
  return(saida_NNMLF)
}


