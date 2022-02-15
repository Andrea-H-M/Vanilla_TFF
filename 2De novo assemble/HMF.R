#######
# ISSUE: Heatmap
# Author: Olga Andrea Hernandez Miranda, Miranda H 
# Date: 17/08/2021
# Note: Heatmap with colorful cluster
# 
#https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html
#https://www.cienciadedatos.net/documentos/37_clustering_y_heatmaps
#######

# Heatmap
library("gplots")
library("RColorBrewer")
library("viridis")
library(lattice)
library(dendextend)

#Crear matriz y normalizar valores
alpha <- 0.01
directorio <- "C:/Users/andii/OneDrive/Documents/CursoED/Denovo/HMv2"
setwd(directorio)
data <- read.table("HMa.csv", sep = ",", header = T)
row.names(data) <- data[,1]
mat_data <- data.matrix(data[,1:ncol(data)]) 
mat_data2 <- mat_data[,-1]
mat_data2 
countTable.kept <- log2(mat_data2) #log forma de reducir diferencias tecnicas entre los datos
dim(countTable.kept)

datos <- countTable.kept

datos <- scale(datos)

colores <- viridis(256)

colores2 <- magma(256)


#Version 1
dend_r <- datos %>% dist(method = "euclidean") %>%
        hclust(method = "average") %>% 
        as.dendrogram %>% ladderize %>%
        color_branches(k=3)

dend_c <- t(datos) %>% dist(method = "euclidean") %>%
        hclust(method = "average") %>% 
        as.dendrogram %>% ladderize%>%
        color_branches(k=3)


# some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
# some_col_func <- colorspace::diverge_hcl
# some_col_func <- colorspace::sequential_hcl
some_col_func <- function(n) (colorspace::diverge_hcl(n, h = c(246, 40), c = 96, l = c(65, 90)))

# par(mar = c(3,3,3,3))
# library(gplots)
heatmap.2(as.matrix (datos),
          cexCol = 0.5,
          cexRow= 1,
          Rowv = dend_r,
          Colv = dend_c,
          trace="none", hline = NA,         
          margins =c(7,11),
          denscol = "grey",
          density.info = "density",
          key.title=NA,
          col = colores
)

#Version 2
dend_r <- datos %>% dist(method = "euclidean") %>%
        hclust(method = "average") %>% 
        as.dendrogram %>% ladderize %>%
        color_branches(k=3)

dend_c <- t(datos) %>% dist(method = "euclidean") %>%
        hclust(method = "average") %>% 
        as.dendrogram %>% ladderize%>%
        color_branches(k=3)

# some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
# some_col_func <- colorspace::diverge_hcl
# some_col_func <- colorspace::sequential_hcl
some_col_func <- function(n) (colorspace::diverge_hcl(n, h = c(246, 40), c = 96, l = c(65, 90)))

# par(mar = c(3,3,3,3))
# library(gplots)
heatmap.2(as.matrix (datos),
          cexCol = 1,
          cexRow= 0.5,
          Rowv = dend_r,
          Colv = dend_c,
          trace = "none",
          lhei=c(2,8),
          lwid=c(2,3),
          keysize=0.1,
          key.par = list(cex=1),
          margins =c(11,11),
          key.title=NA,
          col = colores
)

