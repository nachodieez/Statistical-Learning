library(class)
library(MASS)
library(moments)
library(nnet)
library(naivebayes)
library(pROC)

set.seed(100487766)

load_data <- function(){
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
  
  #data <- data[!is.na(data$Insulin),]
  require(mice)
  dataIm <- mice(data, m = 1, method = "pmm")
  df <- complete(dataIm)
  return(df)
}

### Colors to be used in the script

color_1 <- "deepskyblue2"
color_2 <- "darkorchid4"
color_3 <- "seagreen2"
color_4 <- "indianred2"

### Data things
{
df = load_data()
X <- df[,1:8]
Y <- df[,9]
X <- scale(X, center=F) #center = T, tal vez cambiar a F

n <- nrow(X)
p <- ncol(X)
c(n,p)


n_no <- sum(Y==0)
n_yes <- sum(Y==1)
c(n_no,n_yes)

pr_no <- n_no/n
pr_yes <- n_yes/n
c(pr_no,pr_yes)

n_train <- floor(.7*n) ## 70/30 partition for train and test
n_test <- n - n_train
c(n_train,n_test)
  
i_train <- sort(sample(1:n,n_train))

X_train <- X[i_train,]
X_test <- X[-i_train,]
Y_train <- Y[i_train]
Y_test <- Y[-i_train]

sum(Y_train==0)/n_train
sum(Y_train==1)/n_train
sum(Y_test==0)/n_test
sum(Y_test==1)/n_test
}

### Visual analysis
{
colors_Y <- c(color_1,color_2)[1*(Y==1)+1]
parcoord(X,col=colors_Y,var.label=F,main="PCP for the diabetes data set")

X_skewness <- apply(X,2,skewness)
sort(X_skewness,decreasing=TRUE)
plot(1:p,apply(X,2,skewness),type="h",pch=19,col=color_1,xlab="Variables",
     ylab="Skewness coefficients",main="Skewness coefficients")

# We have some highly skewed variables
#  The methods that we are going to consider for supervised classification works better for symmetric variables
# Then, we take logarithms of all the variables plus 1, because most of the variables contains 0s 

X_log <- log(X+1)
parcoord(X_log,col=colors_Y,var.label=FALSE,main="PCP for the logs of the diabetes data set")
X_log_skewness <- apply(X_log,2,skewness)
sort(X_log_skewness,decreasing=TRUE)
plot(1:p,apply(X_log,2,skewness),type="h",pch=19,col=color_1,
     xlab="Variables",ylab="Skewness coefficients",main="Skewness coefficients")

# PCA
X_log_pcs <- prcomp(X_log,scale=TRUE)
summary(X_log_pcs)

# Make a plot of the first two PCs

plot(X_log_pcs$x[,1:2],pch=20,col=colors_Y,xlim=c(-5,5),ylim=c(-5,5),
     main="First two PCs for the logs of diabetes data set")
#there is a lot of overlapping between the groups
}

### Logs? (no parecen afectar mucho al resultado)
{
X <- log(X+1)
X_train <- X[i_train,]
X_test <- X[-i_train,]
}

### KNN
{
LER <- rep(NA,40)
for (i in 3 : 40){
  knn_output <- knn.cv(X_train,Y_train,k=i)
  LER[i] <- 1 - mean(knn_output==Y_train)
}
LER
plot(1:40,LER,pch=20,col=color_1,type="b",
     xlab="k",ylab="LER",main="LER for logs of spam data set")
k <- which.min(LER)
k

knn_Y_test <- knn(X_train,X_test,Y_train,k=k,prob=T)
head(knn_Y_test)
summary(knn_Y_test)
table(Y_test,knn_Y_test)
knn_TER <- mean(Y_test!=knn_Y_test)
knn_TER

prob_knn_Y_test <- attributes(knn_Y_test)$prob

{
head(prob_knn_Y_test)
prob <- 2*ifelse(knn_Y_test == "-1", prob_knn_Y_test, 1-prob_knn_Y_test) - 1

pred_knn <- prediction(prob, Y_test)
pred_knn <- performance(pred_knn, "tpr", "fpr")
plot(pred_knn, avg= "threshold", colorize=T, lwd=3, main="Voilà, a ROC curve!", asp=1)
}

# Make a plot of the probabilities of the winner group
# In green, good classifications, in red, wrong classifications

colors_errors <- c(color_4,color_3)[1*(Y_test==knn_Y_test)+1]
plot(1:n_test,prob_knn_Y_test,col=colors_errors,pch=20,type="p",
     xlab="Test sample",ylab="Diabetes probabilities",main="Diabetes probabilities")
}

### Logistic regression
{
lr_train <- multinom(Y_train~.,data=as.data.frame(X_train))

# Have a look at the estimated coefficients and their standard errors
summary(lr_train)

# To see which are the most significant coefficients, we use the t-test that is computed next

t_test_lr_train <- summary(lr_train)$coefficients/summary(lr_train)$standard.errors
t_test_lr_train

# Sort the absolute value of the t-test in decreasing order

sort(abs(t_test_lr_train),decreasing=TRUE)

# We can see that the most relevant variables appear to be Glucose and BMI

##################################################################################################################
# Classify the responses in the test data set and estimate the test error rate

lr_test <- predict(lr_train,newdata=as.data.frame(X_test))
head(lr_test)

# Number of instances classified in each group

summary(lr_test)

# Table with good and bad classifications

table(Y_test,lr_test)

# Obtain the Test Error Rate (TER)

lr_TER <- mean(Y_test!=lr_test)
lr_TER

# Obtain the probabilities of non-diabetes

prob_lr_test <- 1 - predict(lr_train,newdata=X_test,type ="probs")
head(prob_lr_test)

# Make a plot of the probabilities of being nonspam
# In green, good classifications, in red, wrong classifications

colors_errors <- c(color_4,color_3)[1*(Y_test==lr_test)+1]
plot(1:n_test,prob_lr_test,col=colors_errors,pch=20,type="p",
     xlab="Test sample",ylab="Probabilities of Diabetes",
     main="Probabilities of Diabetes")
abline(h=0.5)

##################################################################################################################
# Try to improve the test error rate by deleting predictors without discriminatory power
# NOTE: THIS IS TIME CONSUMING

step_lr_train <- step(lr_train,direction="backward",trace=0)
step_lr_train

### BIC
final_df <- as.data.frame(cbind(X_train, Y_train))
final_df$Y_train <- final_df$Y_train - 1
mod_zero  <- glm(Y_train ~ 1, family = binomial, data = final_df)
mod_all   <- glm(Y_train ~ ., family = binomial, data = final_df)
model_glm <- MASS::stepAIC(mod_zero, scope = list(lower = mod_zero, upper = mod_all),
                      direction = "both", trace = 0, k = log(nrow(final_df)))
X_test_2 <- as.data.frame(X_test[, names(model_glm$coefficients)[-1]])
pred <- predict(model_glm, X_test_2) > 0.1
cm <- table(Y_test, pred)

BAC <- ((cm[1,1]/(cm[1,1] + cm[1,2])) + (cm[2,2]/(cm[2,1] + cm[2,2])))/2
BAC


##### voy a hacer una ROC
get_logistic_pred = function(mod, data, res = "y", pos = 1, neg = 0, cut = 0.5) {
  probs = predict(mod, newdata = data, type = "response")
  ifelse(probs > cut, pos, neg)
}

X_test_2 <- as.data.frame(X_test[, names(model_glm$coefficients)[-1]])

test_pred_10 = get_logistic_pred(model_glm, data = X_test_2, res = "default", 
                                 pos = 1, neg = 0, cut = 0.1)
test_pred_50 = get_logistic_pred(model_glm, data = X_test_2, res = "default", 
                                 pos = 1, neg = 0, cut = 0.5)
test_pred_90 = get_logistic_pred(model_glm, data = X_test_2, res = "default", 
                                 pos = 1, neg = 0, cut = 0.9)

test_tab_10 = table(predicted = test_pred_10, actual = Y_test)
test_tab_50 = table(predicted = test_pred_50, actual = Y_test)
test_tab_90 = table(predicted = test_pred_90, actual = Y_test)

test_con_mat_10 = confusionMatrix(test_tab_10, positive = "1")
test_con_mat_50 = confusionMatrix(test_tab_50, positive = "1")
test_con_mat_90 = confusionMatrix(test_tab_90, positive = "1")

metrics = rbind(
  
  c(test_con_mat_10$overall["Accuracy"], 
    test_con_mat_10$byClass["Sensitivity"], 
    test_con_mat_10$byClass["Specificity"]),
  
  c(test_con_mat_50$overall["Accuracy"], 
    test_con_mat_50$byClass["Sensitivity"], 
    test_con_mat_50$byClass["Specificity"]),
  
  c(test_con_mat_90$overall["Accuracy"], 
    test_con_mat_90$byClass["Sensitivity"], 
    test_con_mat_90$byClass["Specificity"])
  
)



rownames(metrics) = c("c = 0.10", "c = 0.50", "c = 0.90")
colnames(metrics) = c("Accuracy", "Sensitivity", "Specificity")

metrics

test_prob = predict(model_glm, newdata = X_test_2, type = "response")
test_roc = roc(Y_test ~ test_prob, plot = TRUE, print.auc = TRUE)

ROC = plot.roc(Y_test, test_prob,
         main="Confidence interval of a threshold", percent=TRUE,
         ci=TRUE, of="thresholds", # compute AUC (of threshold)
         thresholds="best", # select the (best) threshold
         print.thres="best",
         print.auc=TRUE) # also highlight this threshold on the plot
ROC

pred <- predict(model_glm, X_test_2, type = "response") > as.numeric(rownames(ROC$ci$specificity))
pred <- pred * 1
cm <- table(Y_test, pred)

BAC <- ((cm[1,1]/(cm[1,1] + cm[1,2])) + (cm[2,2]/(cm[2,1] + cm[2,2])))/2
BAC

###

# Have a look at the important variables that have been retained in the model

# Classify the responses in the test data set with this new model

step_lr_test <- predict(step_lr_train,newdata=X_test)
head(step_lr_test)

# Number of emails classified in each group

summary(step_lr_test)


# Table with good and bad classifications

table(Y_test,step_lr_test)

# Obtain the Test Error Rate (TER)

step_lr_TER <- mean(Y_test!=step_lr_test)
step_lr_TER # same performance

# The performance of similar to LR with all predictors, but now we know the important predictors

# Obtain the probabilities of nondiabetes

prob_step_lr_test <- 1 - predict(step_lr_train,newdata=X_test,type ="probs")
head(prob_step_lr_test)

# Make a plot of the errors made and the probabilities to see if they are complicated cases
# In green, good classifications, in red, wrong classifications

colors_errors <- c(color_4,color_3)[1*(Y_test==step_lr_test)+1]
plot(1:n_test,prob_step_lr_test,col=colors_errors,pch=20,type="p",
     xlab="Test sample",ylab="Probabilities of non diabetes",
     main="Probabilities of non diabetes")
abline(h=0.5)
}

### LDA
{
lda_train <- lda(Y_train~.,data=as.data.frame(X_train))

# Estimated prior probabilities for the two groups

lda_train$prior

# Estimated sample mean vectors

t(lda_train$means)
# Classify the observations in the test sample

lda_test <- predict(lda_train,newdata=as.data.frame(X_test))

# The vector of classifications made can be found here

lda_Y_test <- lda_test$class
lda_Y_test

# Number of emails classified in each group

table(lda_Y_test)

# Contingency table with good and bad classifications

table(Y_test,lda_Y_test)

# Test Error Rate (TER)

lda_TER <- mean(Y_test!=lda_Y_test)
lda_TER

# Conditional probabilities of the classifications made with the test sample

prob_lda_Y_test <- lda_test$posterior
head(prob_lda_Y_test)

# Make a plot of the probabilities of being nonspam
# In green, good classifications, in red, wrong classifications

colors_errors <- c(color_4,color_3)[1*(Y_test==lda_Y_test)+1]
plot(1:n_test,prob_lda_Y_test[,1],col=colors_errors,pch=20,type="p",
     xlab="Test sample",ylab="Probabilities of non diabetes",
     main="Probabilities of non diabetes")
abline(h=0.5)

test_prob = predict(lda_train, newdata = as.data.frame(X_test), type = "response")
test_prob <- test_prob$posterior[,1]
test_roc = roc(Y_test ~ test_prob, plot = TRUE, print.auc = TRUE)

ROC = plot.roc(Y_test, test_prob,
               main="Confidence interval of a threshold", percent=TRUE,
               ci=TRUE, of="thresholds", # compute AUC (of threshold)
               thresholds="best", # select the (best) threshold
               print.thres="best",
               print.auc=TRUE) # also highlight this threshold on the plot
ROC


}

### QDA
{
qda_train <- qda(Y_train~.,data=as.data.frame(X_train))

qda_train$prior
t(qda_train$means)
qda_test <- predict(qda_train,newdata=as.data.frame(X_test))
qda_Y_test <- qda_test$class
qda_Y_test

table(qda_Y_test)

# Contingency table with good and bad classifications

table(Y_test,qda_Y_test)

# Test Error Rate (TER)

qda_TER <- mean(Y_test!=qda_Y_test)
qda_TER
prob_qda_Y_test <- qda_test$posterior
head(prob_qda_Y_test)

# Make a plot of the probabilities of being non diabetes
# In green, good classifications, in red, wrong classifications

colors_errors <- c(color_4,color_3)[1*(Y_test==qda_Y_test)+1]
plot(1:n_test,prob_qda_Y_test[,1],col=colors_errors,pch=20,type="p",
     xlab="Test sample",ylab="Probabilities of nonspam",
     main="Probabilities of nonspam")
abline(h=0.5)

test_prob = predict(qda_train, newdata = as.data.frame(X_test), type = "response")
test_prob <- test_prob$posterior[,1]
test_roc = roc(Y_test ~ test_prob, plot = TRUE, print.auc = TRUE)

ROC = plot.roc(Y_test, test_prob,
               main="Confidence interval of a threshold", percent=TRUE,
               ci=TRUE, of="thresholds", # compute AUC (of threshold)
               thresholds="best", # select the (best) threshold
               print.thres="best",
               print.auc=TRUE) # also highlight this threshold on the plot
ROC


}

### NB
{
nb_train <- gaussian_naive_bayes(X_train,Y_train)

# Estimated unconditional probabilities

nb_train$prior

# Estimated sample mean vectors

t(nb_train$params$mu)

# The function does not return the estimated covariance matrices

t(nb_train$params$sd)

##################################################################################################################
# Classify the observations in the test sample

nb_test <- predict(nb_train,newdata=X_test,type="prob")

# The vector of classifications made can be found here

nb_Y_test <- c(0,1)[apply(nb_test,1,which.max)]
nb_Y_test

# Number of instances classified in each group

table(nb_Y_test)

# Contingency table with good and bad classifications

table(Y_test,nb_Y_test)

# Test Error Rate (TER)

nb_TER <- mean(Y_test!=nb_Y_test)
nb_TER

# Conditional probabilities of the classifications made with the test sample

prob_nb_Y_test <- nb_test
head(prob_nb_Y_test)

# Make a plot of the probabilities of being nonspam
# In green, good classifications, in red, wrong classifications

colors_errors <- c(color_4,color_3)[1*(Y_test==nb_Y_test)+1]
plot(1:n_test,prob_nb_Y_test[,1],col=colors_errors,pch=20,type="p",
     xlab="Test sample",ylab="Probabilities of non diabetes",
     main="Probabilities of non diabetes")
abline(h=0.5)

test_prob = predict(nb_train, newdata = X_test, type = "prob")
test_prob <- test_prob[,2]
test_roc = roc(Y_test ~ test_prob, plot = TRUE, print.auc = TRUE)

ROC = plot.roc(Y_test, test_prob,
               main="Confidence interval of a threshold", percent=TRUE,
               ci=TRUE, of="thresholds", # compute AUC (of threshold)
               thresholds="best", # select the (best) threshold
               print.thres="best",
               print.auc=TRUE) # also highlight this threshold on the plot
ROC
}

