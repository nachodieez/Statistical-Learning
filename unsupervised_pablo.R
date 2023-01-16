require(tidyverse)
require(mice)
require(dbscan)
require(MASS)
require(cluster)
require(mclust)
require(corpcor)
require(RSpectra)
require(factoextra)

data <- read.csv("diabetes.csv")
data$Insulin[data$Insulin == 0] <- NA
dataIm <- mice(data, m = 1, method = "pmm")
data <- complete(dataIm)

c1 <- "deeppink4"
c2 <- "gold2"
c3 <- "darkorchid4"
c4 <- "seagreen2"
c5 <- "orange2"
c6 <- "firebrick2"

#### DBSCAN ####

X <- data[,-9]
y <- data[,9]
X2 <- scale(X)
pca_fit <- princomp(X2)
X3 <- pca_fit$scores[,1:4]
minPts <- 5
kNNdistplot(X3, k = minPts - 1)
abline(h = 1.2, lwd = 2, lty = 2, col = c1)

db_fit <- dbscan(X3, eps = 1.2, minPts = 5)
colors_db <- c(c1, c2, c3)[db_fit$cluster + 1]
plot(X3[,1:2], pch = 19, col = colors_db, xlab = "First PC",
     ylab = "Second PC")

#### Hierarchical Clustering #### 

euc_dist_X <- daisy(X3, metric = "euclidean")

single_X <- hclust(euc_dist_X, method = "single")
ward_X <- hclust(euc_dist_X, method = "ward.D")
complete_X <- hclust(euc_dist_X, method = "complete")
average_X <- hclust(euc_dist_X, method = "average")

plot(single_X, main = "Single linkage", cex = 0.8)
rect.hclust(single_X, k = 5, border = c1)
#absolutamente terrible

cl_single_X <- cutree(single_X, 2)
cl_single_X
table(cl_single_X)
#absolutamente terrible

sil_single_X <- silhouette(cl_single_X, euc_dist_X)
plot(sil_single_X, col = c1)
#que cojones

plot(ward_X, main = "Ward linkage", cex = 0.8)
rect.hclust(ward_X, k = 2, border = c1)
#bastante mejor

cl_ward_X <- cutree(ward_X, 2)
cl_ward_X
table(cl_ward_X)
table(cl_ward_X, -y+1)
#bastante mejor 

sil_ward_X <- silhouette(cl_ward_X, euc_dist_X)
plot(sil_ward_X, col = c1)
#que cojones

plot(complete_X, main = "Complete linkage", cex = 0.8)
rect.hclust(complete_X, k = 2, border = c1)
#bien

cl_complete_X <- cutree(complete_X, 2)
cl_complete_X
table(cl_complete_X)
table(cl_complete_X, y)
#bastante razonable

sil_complete_X <- silhouette(cl_complete_X, euc_dist_X)
plot(sil_complete_X, col = c1)
#que cojones

plot(average_X, main = "Average linkage", cex = 0.8)
rect.hclust(average_X, k = 2, border = c1)
#mal

cl_average_X <- cutree(average_X, 2)
cl_average_X
table(cl_average_X)
table(cl_average_X, y)
#fatal

sil_average_X <- silhouette(cl_average_X, euc_dist_X)
plot(sil_average_X, col = c1)
#que cojones

#### Partitional Clustering ####

fviz_nbclust(X3, kmeans, method = "silhouette", k.max = 10)
#indicates 2

gap_stat <- clusGap(X3, FUN = kmeans, K.max = 10, B = 100)
fviz_gap_stat(gap_stat, linecolor = c1, maxSE = 
                list(method = "firstSEmax", SE.factor = 1))
#this indicates 2!!!! "oleee"

kmeans_X <- kmeans(X3, centers = 2, iter.max = 1000, nstart = 100)
kmeans_X$cluster

kmeans_col <- c(c1, c2)[kmeans_X$cluster]
plot(X3[,1:2], col = kmeans_col, main = "First two PCs", 
     xlab = "First PC", ylab = "Second PC")
#de super puta madre
sil_kmeans_X <- silhouette(kmeans_X$cluster, euc_dist_X)
plot(sil_kmeans_X, col = c1)

pam_X <- pam(X3, k = 2, metric = "euclidean", stand = FALSE)
pam_col <- c(c1, c2)[pam_X$cluster]

plot(X3[,1:2], col = pam_col, main = "First two PCs", 
     xlab = "First PC", ylab = "Second PC")

#igual que el kmeans
table(pam_X$cluster, kmeans_X$cluster)

sil_pam_X <- silhouette(pam_X$cluster, euc_dist_X)
plot(sil_pam_X, col = c1)
#que cojones

clara_X <- clara(X3, k = 2, metric = "euclidean", stand=FALSE)
clara_col <- c(c1, c2)[clara_X$cluster]

plot(X3[,1:2], col = clara_col, main = "First two PCs", 
     xlab = "First PC", ylab = "Second PC")
#igual que el kmeans y pam
table(clara_X$cluster, kmeans_X$cluster)
table(clara_X$cluster, pam_X$cluster)

#### Model-based Clustering ####

BIC_X <- mclustBIC(Z, G = 1:5)
print(BIC_X)
plot(BIC_X)
summary(BIC_X)

#### K-Medoids ####

X_K <- matrix(NA, nrow = 1, ncol = 19)
for (i in 1:19){
  pam_X_euc_mat <- pam(euc_dist_X, k = i+1, diss = TRUE)
  X_K[i] <- pam_X_euc_mat$silinfo$avg.width
}

plot(2:20, X_K, pch = 19, col = "deepskyblue2", xlab = "Number of clusters",ylab="Average silhouette")
which.max(X_K)+1

# dos clusters!!! oleee

pam_X_euc_mat <- pam(euc_dist_X, k = 2, diss = TRUE)

# Medoids

X[pam_X_euc_mat$medoids,]

# Have a look at the silhouette

sil_pam_X_euc_mat <- silhouette(pam_X_euc_mat$cluster, euc_dist_X)
plot(sil_pam_X_euc_mat, col = c1)
summary(sil_pam_X_euc_mat)
