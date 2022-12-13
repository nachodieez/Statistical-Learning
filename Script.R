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
