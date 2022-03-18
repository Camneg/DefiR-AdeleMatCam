# Importation des données
expData <- read.table("Mito_Genes.txt", row.names = 1, sep = "\t", header = T)

# Fonction Gaëlle représentation expression gene
plotGenes <- function(expData, title = "", yMin = 0, yMax = NULL, meanProfile = TRUE){
  
  # Check function parameters
  if(is.null(yMax)){
    
    print("You must specify a maximal value for Y axis")
    
  }else{
    
    # Representation of the first expression profile
    plot(1:ncol(expData), expData[1,], col = "grey", type = "l",
         ylim = c(floor(yMin), ceiling(yMax)),
         xlab = "Time point", ylab = "Gene expression level",
         main = title)
    
    # Add expression profile for other genes
    for(i in 2:nrow(expData)){
      
      lines(1:ncol(expData), expData[i,], col = "grey")
      
      # end of for()  
    }
    
    # Average expression profile
    if(meanProfile == TRUE){
      expMean = apply(expData, 2, mean)
      lines(1:ncol(expData), expMean, col = "red", 
            lwd = 1.5, lty = "dashed")
    }
    
    # end of else()   
  }
  
  # end of function plotGenes()  
}

################### méthode k-means/Distance euclidienne #######################

# Clustering avec la méthode k-means/Distance euclidienne = A

N <- 6

#Reset paramettres graphiques
#dev.off()

#kmeans du tableau (par défaut distance euclidiennes)
resKmeansEucl <- kmeans(expData, centers = N)

#Condition de partition feuille graph (à optimiser)
if(N<=9){
  par(mfrow = c(ceiling(sqrt(N)),ceiling(sqrt(N))))
}else{
  par(mfrow = (c(3,3)))
}

#Générer un graph profile expression clusters
for(i in 1:N){
  cluster <- expData[which(resKmeansEucl$cluster == i),]
  plotGenes(cluster, yMax = 100)
}

################### méthode  k-maens/correlation #######################
# k-maens/correlation = B
#dev.off()
# Creer une matrice de distance à partir d'une matrice de corrélation
matDistCorr <- as.dist(1 - cor(t(expData)))

resKmeansCorr <- kmeans(matDistCorr, centers = N)

if(N<=9){
  par(mfrow = c(ceiling(sqrt(N)),ceiling(sqrt(N))))
}else{
  par(mfrow = (c(3,3)))
}

for(i in 1:N){
  cluster <- expData[which(resKmeansCorr$cluster == i),]
  plotGenes(cluster, yMax = 100)
}

################### méthode  Hierarchie / euclidienne #######################
# Hierarchie / euclidienne = C
#dev.off()
matDistEucl <- dist(expData, method = "euclidean")

# clustering en coupant l'arbre et ward.D2 = methode de clustering la plus robuste (Warning)
resHierEucl <- hclust(matDistEucl, method = "ward.D2")

for(i in 1:N){
clusterC <- expData[which(cutree(resHierEucl, k = N) == i),]
plotGenes(clusterC, yMax = 100)
}

################### méthode  Hierarchie / Correlation #######################

resHierCorr <- hclust(matDistCorr, method = "ward.D2")

for(i in 1:N){
  clusterD <- expData[which(cutree(resHierCorr, k = N) == i),]
  plotGenes(clusterD, yMax = 100)
}

heatmap(as.matrix(clusterD))


###################### MatCode ggplot2 ############
#Appel de la laibrairie dplyr pour remodeller des data en tidy
library(dplyr)
N=2
i=2
#Ajouter un nom de colonne aux "noms de lignes"
clusterB<-rownames_to_column(cluster2)
#Transformer les colonnes des temps d'experience en lignes (une ligne = nom de gène + nom timepoint + valeur)
ClustC<-clusterB %>% pivot_longer(c(2:ncol(clusterB)))
#Transformer les noms de gènes en variables pour pouvoir utiliser certans types de graphiques (geom_smooth)
ClustD<-ClustC %>%
  group_by(name) %>%
  mutate(group_id = cur_group_id())
#Graph des profils d'expressions avec moyennes
ggplot(ClustD,aes(group_id,value))+geom_line(aes(group=rowname),color="grey")+geom_smooth(color="red")+
  #graphique ligne +moyenne (ne pas prendre en compte le message)
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Time point")+ ylab("Gene expression level")+  ggtitle(paste("Expression Profile Cluster ",i))+
  scale_x_continuous(expand = c(0, 0), limits = c(NA,NA))

###################### MatCode pour fun ############

Graph_clust<-function(data,N=1,ggsave=FALSE,width=15,height=10){
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  data %>% rownames_to_column() %>% pivot_longer(c(2:ncol(.)))%>%
    group_by(name) %>%
    mutate(group_id = cur_group_id())%>%
    ggplot(aes(group_id,value))+geom_line(aes(group=rowname),color="grey")+
    stat_summary(aes(y = value,group=1), fun=mean, colour="red", geom="line",group=1,linetype="dashed")+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    xlab("Time point")+ ylab("Gene expression level")+  ggtitle(paste("Expression Profile Cluster ",N))+
    scale_x_continuous(expand = c(0, 0), limits = c(NA,NA))
  if(ggsave==TRUE){
    ggsave(paste("Cluster",N,".png"),width=width,height=height)
  }
}

