---
  title: 'Feature Selection + Classification'
author: 'Saarika'
date: '12 April 2020'
---

library(ggplot2)
library(caret)
library(mlbench)
library(Boruta)
library(randomForest)
library(class)
library(pROC)
library(dplyr)
library(psych)
root_dir = getwd()

#Read files
clinical <- read.table(file = file.path(root_dir, 'Desktop', 'ML', 'train_clinical.tsv'), header=TRUE)
dmy <- dummyVars(" ~ Subgroup", data = clinical)
Y <- data.frame(predict(dmy, newdata =  clinical[c('Subgroup')]))
df1 <- within (clinical,"g" <- match(clinical$Subgroup,unique(clinical$Subgroup)))


z <- read.delim(file.path(root_dir, 'Desktop', 'ML', 'z.tsv'), sep = "\t")


#Feature Selection with Boruta
#---------------------------------------------------#
set.seed(111)
boruta <- Boruta(X ~ ., data = z, doTrace = 2, maxRuns = 1000)

#Tentative feature fix
borutaTfix <- TentativeRoughFix(borat)

#Train and Test split up
set.seed(222)
ind <- sample(2, nrow(z), replace = T, prob = c(0.6, 0.4))

train <- z[ind == 1,]
test <- z[ind == 2,]

train$X <- as.factor(train$X)
test$X <- as.factor(test$X)

#Classification Algorithms
#---------------------------------------------------#

#1. KNN
#On Hyperparameter tuning, the iseal k value was seen to be 6

set.seed(1234)
trControl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
grid <- expand.grid(k=(1:10))
set.seed(222)

#Accuracy : 0.3488 - All Classes
fit <- train(X~., data = train, method = 'knn', trControl = trControl, preProc = c("center", "scale"), tuneGrid=grid)
pknn <- predict(fit, test)
confusionMatrix(pknn, test$X)

#g <- knn(cl = train$X,train = train, test = test, k=6) <- Use this to skip grid search application for each of the cases

set.seed(222)
#Accuracy : 0.7209 - Confirmed + Tentative Classes
fit62 <- train(X ~ X463 + X464 + X465 + X472 + X474 + X475 + X479 + X484 + X485 + 
                 X486 + X487 + X488 + X489 + X490 + X696 + X773 + X840 + X842 + 
                 X850 + X852 + X853 + X854 + X855 + X1000 + X1001 + X1002 + 
                 X1650 + X1652 + X1656 + X1657 + X1661 + X1669 + X1671 + X1678 + 
                 X1679 + X1680 + X1684 + X1687 + X1688 + X2057 + X2105 + X2106 + 
                 X2109 + X2113 + X2114 + X2124 + X2125 + X2127 + X2155 + X2184 + 
                 X2185 + X2186 + X2207 + X2211 + X2212 + X2214 + X2215 + X2220 + 
                 X2221 + X2222 + X2224 + X2225, data = train, method = 'knn', trControl = trControl, preProc = c("center", "scale"), tuneGrid=grid)
pknn62 <- predict(fit62, test)
confusionMatrix(pknn62, test$X)

set.seed(222)
#Accuracy : 0.7907 - Confirmed Classes
fit51 <- train(X ~ X463 + X464 + X465 + X475 + X479 + X485 + X486 + X487 + X488 + 
                 X489 + X490 + X696 + X773 + X840 + X850 + X852 + X853 + X854 + 
                 X855 + X1001 + X1002 + X1650 + X1652 + X1656 + X1657 + X1661 + 
                 X1669 + X1671 + X1678 + X1679 + X1680 + X1684 + X1687 + X1688 + 
                 X2057 + X2105 + X2113 + X2114 + X2125 + X2184 + X2185 + X2186 + 
                 X2207 + X2211 + X2212 + X2214 + X2215 + X2220 + X2221 + X2222 + 
                 X2224, data = train, method = 'knn', trControl = trControl, preProc = c("center", "scale"), tuneGrid=grid)
pknn51 <- predict(fit51, test)
confusionMatrix(pknn51, test$X)

#Neural Network
set.seed(444)
nnetGrid <- expand.grid(decay = c(0.01),  size = c(2))

maxSize <- max(nnetGrid$size)

set.seed(1000)
ctrl <- trainControl(method = "cv",  # corss-validation
                     number = 5,  # 10 folds
                     classProbs = FALSE , # report class probability
                     
)

#Accuracy : 0.6279 - All Classes ; 0.01   2     0.6514286  0.46468838
xx <- data.frame(subset(train, select = -c(X)))
numWts <- 1*(maxSize * (length(xx) + 1) + maxSize + 1)
n <- train(xx, train$X,
           method = "nnet", # train neural network using `nnet` package 
           tuneGrid = nnetGrid, # tuning grid
           trControl = ctrl, # process customization set before
           preProc = c("center", "scale"), # standardize data
           trace = FALSE,  # hide the training trace
           MaxNWts = numWts,  # maximum number of weight
           maxit = 100 # maximum iteration
)


output <- compute(n,xx)
pnn <- output$net.result
predpnn <- apply(pnn,1,function(x) which(x==max(x)))
predpnn <- as.factor(predpnn)
confusionMatrix(predpnn, test$X)

set.seed(444)
#Accuracy : 0.7907 - Confirmed + Tentative Classes
n62 <- neuralnet(X ~ X463 + X464 + X465 + X472 + X474 + X475 + X479 + X484 + X485 + 
                   X486 + X487 + X488 + X489 + X490 + X696 + X773 + X840 + X842 + 
                   X850 + X852 + X853 + X854 + X855 + X1000 + X1001 + X1002 + 
                   X1650 + X1652 + X1656 + X1657 + X1661 + X1669 + X1671 + X1678 + 
                   X1679 + X1680 + X1684 + X1687 + X1688 + X2057 + X2105 + X2106 + 
                   X2109 + X2113 + X2114 + X2124 + X2125 + X2127 + X2155 + X2184 + 
                   X2185 + X2186 + X2207 + X2211 + X2212 + X2214 + X2215 + X2220 + 
                   X2221 + X2222 + X2224 + X2225, data = train, hidden = 1000, err.fct = "ce", linear.output = FALSE)
a <- as.factor(subset(train, select = -c(X)))
output62 <- compute(n62,a)
pnn62 <- output62$net.result
predpnn62 <- apply(pnn62,1,function(x) which(x==max(x)))
predpnn62 <- as.factor(predpnn62)
confusionMatrix(predpnn62, test$X)

set.seed(444)
#Accuracy : 0.7907 - Confirmed only
n51 <- neuralnet(X ~ X463 + X464 + X465 + X475 + X479 + X485 + X486 + X487 + X488 + 
                   X489 + X490 + X696 + X773 + X840 + X850 + X852 + X853 + X854 + 
                   X855 + X1001 + X1002 + X1650 + X1652 + X1656 + X1657 + X1661 + 
                   X1669 + X1671 + X1678 + X1679 + X1680 + X1684 + X1687 + X1688 + 
                   X2057 + X2105 + X2113 + X2114 + X2125 + X2184 + X2185 + X2186 + 
                   X2207 + X2211 + X2212 + X2214 + X2215 + X2220 + X2221 + X2222 + 
                   X2224, data = train, hidden = 1000, err.fct = "ce", linear.output = FALSE)
a <- as.factor(subset(train, select = -c(X)))

output51 <- compute(n62,a)
pnn51 <- output51$net.result
predpnn51 <- apply(pnn51,1,function(x) which(x==max(x)))
predpnn51 <- as.factor(predpnn51)
confusionMatrix(predpnn51, test$X)



