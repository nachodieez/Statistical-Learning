## Read data and set NAs

data <- read.csv("diabetes.csv")
data$Outcome <- factor(data$Outcome, c(0,1))

setNAs <- function(data, fields){
  percentage <- list()
  for (field in fields){
    data[[field]][data[[field]] == 0] <- NA
    percentage[[field]] <- 100*sum(is.na(data[[field]]))/nrow(data)
  }
  return (list(data = data, percentage = percentage))
}

NAfields    <- c("Glucose", "BloodPressure", "SkinThickness", "Insulin", "BMI")
dataNA      <- setNAs(data, NAfields)
data        <- dataNA$data
percentages <- dataNA$percentage

## Outliers

findOutliers <- function(data, fields){
  outliers <- list()
  for (field in fields){
    qs  <- quantile(data[[field]], c(0.25, 0.75), na.rm = TRUE)
    iqr <- qs[2] - qs[1]
    lq  <- qs[1] - 1.5*iqr
    hq  <- qs[2] + 1.5*iqr
    outliers[[field]] <- which((data[[field]] < lq) & (data[[field]] > hq))
  }
  return (outliers)
}

outliers <- findOutliers(data, names(data)[names(data) != "Outcome"])

## Fill NAs (Predictive Mean Matching)

# Drop individuals with missing insulin (48%)

sum(is.na((data[is.na(data$SkinThickness),]$Insulin))) == nrow(data[is.na(data$SkinThickness),])
dim(data[is.na(data$Insulin),])

data <- data[!is.na(data$Insulin),]
require(mice)
dataIm <- mice(data, m = 1, method = "pmm")
data2 <- complete(dataIm)

## Univariate descriptive analysis

require(tidyverse)
devtools::source_gist("2a1bb0133ff568cbe28d", filename = "geom_flat_violin.R")
require(plyr) # necesario tenerlo instalado
univariate_description <- function(x){
  #summary
  {
    s <- summary(x)
    sd_x <- c(St.dev = sd(x))
    s2 <- append(s, sd_x)
    df <- data.frame(round(as.numeric(s2),2)) %>% t() %>% as.data.frame()
    colnames(df) <- names(s2)
    rownames(df) <- NULL
    df2 <- df
  }
  #plot
  {
    pos <- position_jitter(width = 0.15, seed = 1)
    plot <- ggplot(data = data.frame(y=x), aes(x = "", y = y)) +
      geom_flat_violin(position = position_nudge(x = .3, y = 0), alpha = .8,
                       fill = "pink1") +
      theme_classic() + 
      geom_point(position = pos, size = 1.5, alpha = 0.8, col = "pink") +
      geom_boxplot(width = .1, show.legend = FALSE, outlier.shape = NA, 
                   alpha = 0.5) + xlab("") + ylab(names(x))
  }
  return(list(summary = df2, plot = plot))
}
# la elección del pink es discutible, se ve demasiado clarito

a <- univariate_description(data2$Glucose)
a$plot
a$summary

## Bivariate descriptive analysis

# Scatterplot Matrix

require(car)

colnames(data2)[1:8] <- c("Pregnancies", "Glucose",
                         "Blood\nPressure", "Skin\nThickness", 
                         "Insulin", "BMI", 
                         "Diabetes\nPedigree\nFunction",
                         "Age")
scatterplotMatrix(data2[,-9], regLine = F, diagonal = F, 
                  smooth = F, col = "blue", cex.labels = 1.25)

require(scatterPlotMatrix)

scatterPlotMatrix(data2[,-9], corrPlotType = "Text") #????? el profe quiere plots interactivos ??????

# Correlation plot

require(corrplot)

correlation <- cor(data2[,-9])
colnames(correlation) <- c("Pregnancies", "Glucose",
                         "Blood\nPressure", "Skin\nThickness", 
                         "Insulin", "BMI", 
                         "Diabetes\nPedigree\nFunction",
                         "Age")

corrplot.mixed(correlation, lower = "number", upper = "color",
         order = "AOE", diag = "n", tl.col = "black", tl.cex = 0.75) 

## Multidimensional plots

# Andrew's plot

require(pracma)

andrewsplot(as.matrix(data2[,-9]), data2[,9],
            style = "cart")

# PCP

require(MASS)

colors <- c("pink2", "darkblue")
col1 <- colors[1]
col2 <- colors[2]
vec_col <- as.character(data2$Outcome)
vec_col[vec_col==0] <- col1 # esta línea y la siguiente no van
vec_col[vec_col==1] <- col2

par(las=2)
parcoord(data2[,-9], col = vec_col, var.label = T)

legend("topright", legend = c("No diabetes", "Diabetes"),
       col = colors, lty = 1, lwd = 2)













## vamos a hacer cosas del PCA (nacho)

options(digits=4)

color_1 <- "hotpink2"
color_2 <- "forestgreen"

head(data2)


##################################################################################################################
# Define an auxiliary variable for the categorical variable Private 
# and redefine the data matrix without this variable

Y <- data2[,9]
head(Y)
X <- data2[,-9]
head(X)

##################################################################################################################
# Define sample size and the dimension of the data matrix

n <- nrow(X)
p <- ncol(X)
c(n,p)

##################################################################################################################
# Plot the variables individually

##################################################################################################################
# Transform some of the variables and take logs

X$BMI <- log(X$BMI)
colnames(X)[6] <- "log_BMI"
X$Insulin <- log(X$Insulin)
colnames(X)[5] <- "log_Insulin"
X$`Diabetes\nPedigree\nFunction` <- log(X$`Diabetes\nPedigree\nFunction`)
colnames(X)[7] <- "log_Diabetes\nPedigree\nFunction"


###### preguntar a pablo cual log transformó para que haya consistencia

##################################################################################################################
# Two possibilities: Obtain PCs on the whole data set or on the two groups
# See what happens if we obtain PCs on the whole data set
# Use the sample correlation matrix

X_pcs <- prcomp(X, scale = TRUE)

##################################################################################################################
# PC scores

dim(X_pcs$x)
head(X_pcs$x)

##################################################################################################################
# Make a plot of the first two PCs

colors_X <- c(color_1,color_2)[1*(Y=="Yes")+1]
par(mfrow=c(1,1))
plot(X_pcs$x[,1:2],pch=19,col=colors_X)

# The first two PCs show the presence of the two different groups

##################################################################################################################
# Interpretation of the first PC: Loadings for the first PC

plot(1:p,X_pcs$rotation[,1],pch=19,col=color_1,main="Loadings for the first PC",
     xlab="Variables",ylab="Score", ylim = c(0,0.5))
abline(h=0)
text(1:p, X_pcs$rotation[,1], labels = colnames(X),
     pos = 1, col = color_2, cex = 0.75)

##################################################################################################################
# Interpretation of the second PC: Loadings for the second PC

plot(1:p,X_pcs$rotation[,2], pch = 19, col = color_1,
     main = "Loadings for the second PC",
     xlab = "Variables", ylab = "Score")
abline(h = 0)
text(1:p, X_pcs$rotation[,2], labels = colnames(X), pos = 1, 
     col = color_2, cex = 0.75)

##################################################################################################################
# Have a look at the important variables in the first two PCs
# Note the different groups in the data
# The radius is arbitrary

plot(X_pcs$rotation[,1:2], pch = 19, col = color_1, 
     main = "Weights for the first two PCs")
abline(h = 0, v = 0)
text(X_pcs$rotation[,1:2],labels = colnames(X), pos = 1, col = color_2,cex=0.75)
library(plotrix)
#draw.circle(0,0,0.3,border=color_2,lwd=3)

##################################################################################################################
# The biplot is an alternative way to plot points and the first two PCs together 
# However, it is only useful when the data set is not too large

biplot(X_pcs, col = c(color_1,color_2), cex = c(0.5,0.8))

##################################################################################################################
# The eigenvalues of the sample correlation matrix of X_trans, i.e., the variances of the PCs are the square of sdev

X_pcs$sdev^2

##################################################################################################################
# Have a look at these eigenvalues

# Screeplot

library(factoextra)
fviz_eig(X_pcs,ncp=17,addlabels=T,barfill=color_1,barcolor=color_2)

##################################################################################################################
# How many PCs are important?

# Have a look at the proportion of explained variance and the cumulative proportion of explained variance

get_eigenvalue(X_pcs)

# We need 4 PCs to have the 73% of the total variability of the data set
# The sample mean of the eigenvalues is 1 because we are using the sample correlation matrix and the fourth 
# eigenvalue is 0.9942
# Then, the dimension of the data set is reduced from 17 to 4 that represents either the 23.52% of the number of 
# variables in the original data set keeping around the 73% of the information inside

##################################################################################################################
# Plot the scores (the four PCs)

pairs(X_pcs$x[,1:4],col=
        ifelse(Y==1,color_1, color_2)
      ,pch=19,main="The first four PCs")

plot(X_pcs$x[,1], X_pcs$x[,2], col = ifelse(Y==1,color_1, color_2))

# The PCs show that the second PC is the key to show the two groups

# The first one is a PC that creates a ranking of colleges. See the first ten colleges

sort(X_pcs$x[,1], decreasing = TRUE)[1:10]

##################################################################################################################
# Plot the correlations between the original data set and the scores
# This is useful to understand also which are the important variables in the data set in terms of variability
# First, plot all the PCs
# Second, plot only the first four PCs

library(corrplot)
corrplot(cor(X,X_pcs$x),is.corr=T)
corrplot(cor(X,X_pcs$x[,1:4]),is.corr=T)


# Independent component analysis ------- NACHONACHONACHO


##################################################################################################################
# Two possibilities: Obtain ICs on the whole data set or on the two groups
# See what happens if we obtain ICs on the whole data set

# Load the ica package and obtain the ICs with the College data set after transformations

library("ica")
X_trans_ica <- icafast(X,nc=p,alg="par")

##################################################################################################################
# The ICA scores, i.e., the matrix Z, can be found in the object S 

Z <- X_trans_ica$S
dim(Z)
head(Z)
colnames(Z) <- sprintf("IC-%d",seq(1,8))

# However, we need to re-scale to have S_z=I (icafast considers the sample covariance matrix
# by dividing by n instead of n-1)

cov(Z)
Z <- Z * sqrt((n-1)/n)
cov(Z)

##################################################################################################################
# Histograms of the ICA scores obtained

par(mfrow=c(3,3))
sapply(colnames(Z),function(cname){hist(as.data.frame(Z)[[cname]],
                                        main=cname,col=color_1,xlab="")})

##################################################################################################################
# Compute the neg-entropy of the columns in Z and sort them in decreasing order of neg-entropy

neg_entropy <- function(z){1/12 * mean(z^3)^2 + 1/48 * mean(z^4)^2}
Z_neg_entropy <- apply(Z,2,neg_entropy)
ic_sort <- sort(Z_neg_entropy,decreasing=TRUE,index.return=TRUE)$ix
ic_sort

##################################################################################################################
# Plot the sorted neg-entropy and define the ICs sorted by neg-entropy

par(mfrow=c(1,1))
plot(Z_neg_entropy[ic_sort],type="b",col=color_1,pch=19,
     ylab="Neg-entropy",main="Neg-entropies",lwd=3)
Z_ic_imp <- Z[,ic_sort]

##################################################################################################################
# Plot the two ICs with largest neg-entropy

colors_X <- c(color_1,color_2)[1*(Y=="Yes")+1]
par(mfrow=c(1,1))
plot(Z_ic_imp[,1:2],pch=19,col=colors_X,
     xlab="First IC scores",ylab="Second IC scores",main="First two IC scores")

# The first two ICs (in terms of neg-entropy) show the presence of some outliers

# Let's identify the points with large negative values of the first IC

which(Z_ic_imp[,1]>2)

# Let's identify the points with large positive values of the second IC

which(Z_ic_imp[,2]<(-4))

# It will be a good idea to identify why these points outlie with minimum covariance determinant

##################################################################################################################
# Plot the ICs with the highest and lowest entropy

### change index

pairs(Z_ic_imp[,c(1:5,16:17)],col=colors_X,pch=19,
      main="ICs with the highest and lowest entropy")

# Note that the last ICs are able to distinguish the two groups
# Let's see the two ICs with smallest neg-entropy

plot(Z_ic_imp[,16:17],pch=19,col=colors_X,
     xlab="8-th IC scores",ylab="1-th IC scores",main="Last two IC scores")

##################################################################################################################
# Plot the correlations between the data matrix and the ICs sorted by neg-entropy

library(corrplot)
corrplot(cor(X,Z_ic_imp),is.corr=T)

# Note that the ICs are not closely related with a large set of variables, except the last
# two ICs in terms of neg-entropy

##################################################################################################################
# Plot the correlations between the PCs and the ICs sorted by neg-entropy

corrplot(cor(prcomp(X)$x,Z),is.corr=T)

# In general, the linear relationship between PCs and ICs is not strong


### FACTOR ANALYSIS NACHO NACHO NACHO NACHO NACHO NACHO NACHO


n <- nrow(X)
p <- ncol(X)
c(n,p)

##################################################################################################################
# Brief description of the variables
##################################################################################################################

# Obtain a summary of the data set and a barplot of the variables. Note that some of the variables are skewed
# The kurtosis is measured with respect to 3, i.e., the kurtosis coefficient of the Gaussian
# The differences with respect to 3 are not very large

library(psych)
describe(X)
par(mfrow=c(6,6),mar=c(2,2,2,2))
par(mfrow=c(3,3))
for (i in 1:p){barplot(table(X[,i]),col=color_1,main=colnames(X)[i])}
# la parálisis de sueño tal noseque


##################################################################################################################
# These are large matrices, then have a look the the correlation plot

library(corrplot)
par(mfrow=c(1,1))
corrplot(cor(X),order="hclust")

# There are groups of correlated variables that may suggest a factor structure

##################################################################################################################
# Estimation of the factor model
##################################################################################################################

##################################################################################################################
# Standarize variables and obtain PCs

Y <- scale(X)
Y_pcs <- prcomp(Y)

# Screeplot with all the eigenvalues

library(factoextra)
fviz_eig(Y_pcs,ncp=p,addlabels=T,barfill=color_1,barcolor=color_2)
get_eigenvalue(Y_pcs)

# Let focus on the first five PCs

r <- 5

##################################################################################################################
# Initial estimation of M and Sigma_nu

M_0 <- Y_pcs$rotation[,1:r] %*% diag(Y_pcs$sdev[1:r])
M_0
S_y <- cov(Y)
Sigma_nu_0 <- diag(diag(S_y - M_0 %*% t(M_0)))
Sigma_nu_0

##################################################################################################################
# Estimation of M without varimax rotation

MM <- S_y - Sigma_nu_0
MM_eig <- eigen(MM)
MM_values <- MM_eig$values
MM_vectors <- MM_eig$vectors
M_1 <- MM_eig$vectors[,1:r] %*% diag(MM_eig$values[1:r])^(1/2)
M_1

##################################################################################################################
# Final estimation of M and Sigma_nu after varimax rotation for interpretability

M <- varimax(M_1)
M <- loadings(M)[1:p,1:r]
M
Sigma_nu <- diag(diag(S_y - M %*% t(M)))
Sigma_nu

##################################################################################################################
# Understanding the relationship between observable variables and factors

# The first factor appears to be an index of extraversion 

plot(1:p,M[,1],pch=19,col=color_1,xlab="",ylab="Loadings",main="Loadings for the first factor")
abline(h=0)
text(1:p,M[,1],labels=colnames(X),pos=1,col=color_2,cex=0.75)

# The second factor appears to be an index of professional conscientiousness 

plot(1:p,M[,2],pch=19,col=color_1,xlab="",ylab="Loadings",main="Loadings for the second factor")
abline(h=0)
text(1:p,M[,2],labels=colnames(X),pos=1,col=color_2,cex=0.75)

# The third factor appears to be an index of roughness 

plot(1:p,M[,3],pch=19,col=color_1,xlab="",ylab="Loadings",main="Loadings for the third factor")
abline(h=0)
text(1:p,M[,3],labels=colnames(X),pos=1,col=color_2,cex=0.75)

# The fourth factor appears to be an index of restlessness 

plot(1:p,M[,4],pch=19,col=color_1,xlab="",ylab="Loadings",main="Loadings for the fourth factor")
abline(h=0)
text(1:p,M[,4],labels=colnames(X),pos=1,col=color_2,cex=0.75)

# The fifth factor appears to be an index of friendliness

plot(1:p,M[,5],pch=19,col=color_1,xlab="",ylab="Loadings",main="Loadings for the fifth factor")
abline(h=0)
text(1:p,M[,5],labels=colnames(X),pos=1,col=color_2,cex=0.75)

##################################################################################################################
# Communalities and uniquenesses

comM <- diag(M %*% t(M))
comM
names(comM) <- colnames(X)
plot(1:p,sort(comM,decreasing=TRUE),pch=20,col=color_1,xlim=c(0,34),
     xlab="Variables",ylab="Communalities",main="Communalities")
text(1:p,sort(comM,decreasing=TRUE),labels=names(sort(comM,decreasing=TRUE)),
     pos=4,col=color_2,cex=0.75)

uniqM <- 1 - comM
uniqM
plot(1:p,sort(uniqM,decreasing=TRUE),pch=20,col=color_1,xlim=c(0,34),
     xlab="Variables",ylab="Uniquenesses",main="Uniquenesses")
text(1:p,sort(uniqM,decreasing=TRUE),labels=names(sort(uniqM,decreasing=TRUE)),
     pos=4,col=color_2,cex=0.75)

# The variables better explained by the factors are Outgoing, Quiet and Tense
# The variables worst explained by the factors are Lax, Approving and Perseverant

##################################################################################################################
# Estimate the factor scores

Fact <- Y %*% solve(Sigma_nu) %*% M %*% solve(t(M) %*% solve(Sigma_nu) %*% M)
colnames(Fact) <- c("Factor 1","Factor 2","Factor 3","Factor 4","Factor 5")

# See that the factors are uncorrelated

pairs(Fact,pch=19,col=color_1)
corrplot(cor(Fact),order="hclust")

##################################################################################################################
# Estimate the residuals

Nu <- Y - Fact %*% t(M)
corrplot(cor(Nu),order="hclust")

# The residuals show some very minor correlations that the factor model is not able to explain
# but this is something expected



