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
library(mice)
dataIn <- mice(data, m = 5, method = "pmm")
