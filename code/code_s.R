---
  title: 'Exploratory PCA analysis'
author: 'Saarika'
date: '12 April 2020'
output: html_document
---
  
library(ggplot2)
library(caret)
library(mlbench)
library(Boruta)
library(randomForest)
library(class)
library(pROC)
library(naivebayes)
library(dplyr)
library(psych)
library(neuralnet)
library(keras)
library(tensorflow)
library(tfruns)
root_dir = getwd()


clinical <- read.table(file = "../saarika/Desktop/ML/Train_clinical.tsv", header=TRUE)
dmy <- dummyVars(" ~ Subgroup", data = clinical)
Y <- data.frame(predict(dmy, newdata =  clinical[c('Subgroup')]))
df1 <- within (clinical,"gay" <- match(clinical$Subgroup,unique(clinical$Subgroup)))

z <- read.delim(file = "../saarika/Desktop/ML/z.tsv", sep = "\t")
#z <- cbind(data, df1$gay)

set.seed(111)
borat <- Boruta(X ~ ., data = z, doTrace = 2, maxRuns = 1000)

naruto <- TentativeRoughFix(borat)

set.seed(222)
ind <- sample(2, nrow(z), replace = T, prob = c(0.6, 0.4))

train <- z[ind == 1,]
test <- z[ind == 2,]

train$X <- as.factor(train$X)
test$X <- as.factor(test$X)

#RANDOM FOREST
#Accuracy : 0.6744 - All Classes
rfog <- randomForest(X~., data = train)
p <- predict(rfog, test)
confusionMatrix(p, test$X)

#Accuracy : 0.7674 - Confirmed + Tentative Classes
rf62 <- randomForest(X ~ X463 + X464 + X465 + X472 + X474 + X475 + X479 + X484 + X485 + 
                                                   X486 + X487 + X488 + X489 + X490 + X696 + X773 + X840 + X842 + 
                                                   X850 + X852 + X853 + X854 + X855 + X1000 + X1001 + X1002 + 
                                                   X1650 + X1652 + X1656 + X1657 + X1661 + X1669 + X1671 + X1678 + 
                                                   X1679 + X1680 + X1684 + X1687 + X1688 + X2057 + X2105 + X2106 + 
                                                   X2109 + X2113 + X2114 + X2124 + X2125 + X2127 + X2155 + X2184 + 
                                                   X2185 + X2186 + X2207 + X2211 + X2212 + X2214 + X2215 + X2220 + 
                                                   X2221 + X2222 + X2224 + X2225, data = train)

p62 <- predict(rf62, test)
confusionMatrix(p62, test$X)

#Accuracy : 0.7907 - Confirmed Classes only
rf51 <- randomForest(X ~ X463 + X464 + X465 + X475 + X479 + X485 + X486 + X487 + X488 + 
                       X489 + X490 + X696 + X773 + X840 + X850 + X852 + X853 + X854 + 
                       X855 + X1001 + X1002 + X1650 + X1652 + X1656 + X1657 + X1661 + 
                       X1669 + X1671 + X1678 + X1679 + X1680 + X1684 + X1687 + X1688 + 
                       X2057 + X2105 + X2113 + X2114 + X2125 + X2184 + X2185 + X2186 + 
                       X2207 + X2211 + X2212 + X2214 + X2215 + X2220 + X2221 + X2222 + 
                       X2224, data = train)
p51 <- predict(rf51, test)
confusionMatrix(p51, test$X)

#KNN
#set.seed(1234)
trControl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(222)

#Accuracy : 0.4884 - All Classes
fit <- train(X~., data = train, method = 'knn', trControl = trControl, preProc = c("center", "scale"))
pknn <- predict(fit, test)
confusionMatrix(pknn, test$X)

set.seed(222)
#Accuracy : 0.6744 - Confirmed + Tentative Classes
fit62 <- train(X ~ X463 + X464 + X465 + X472 + X474 + X475 + X479 + X484 + X485 + 
                 X486 + X487 + X488 + X489 + X490 + X696 + X773 + X840 + X842 + 
                 X850 + X852 + X853 + X854 + X855 + X1000 + X1001 + X1002 + 
                 X1650 + X1652 + X1656 + X1657 + X1661 + X1669 + X1671 + X1678 + 
                 X1679 + X1680 + X1684 + X1687 + X1688 + X2057 + X2105 + X2106 + 
                 X2109 + X2113 + X2114 + X2124 + X2125 + X2127 + X2155 + X2184 + 
                 X2185 + X2186 + X2207 + X2211 + X2212 + X2214 + X2215 + X2220 + 
                 X2221 + X2222 + X2224 + X2225, data = train, method = 'knn', trControl = trControl, preProc = c("center", "scale"))
pknn62 <- predict(fit62, test)
confusionMatrix(pknn62, test$X)

set.seed(222)
#Accuracy : 0.6744 - Confirmed Classes
fit51 <- train(X ~ X463 + X464 + X465 + X475 + X479 + X485 + X486 + X487 + X488 + 
                 X489 + X490 + X696 + X773 + X840 + X850 + X852 + X853 + X854 + 
                 X855 + X1001 + X1002 + X1650 + X1652 + X1656 + X1657 + X1661 + 
                 X1669 + X1671 + X1678 + X1679 + X1680 + X1684 + X1687 + X1688 + 
                 X2057 + X2105 + X2113 + X2114 + X2125 + X2184 + X2185 + X2186 + 
                 X2207 + X2211 + X2212 + X2214 + X2215 + X2220 + X2221 + X2222 + 
                 X2224, data = train, method = 'knn', trControl = trControl, preProc = c("center", "scale"))
pknn51 <- predict(fit51, test)
confusionMatrix(pknn51, test$X)

#Naive Bayes 
set.seed(100)
#Accuracy : 0.5349 - All Classes
model <- naive_bayes(X~. , data= train)
pnb <- predict(model, test)
confusionMatrix(pnb, test$X)

set.seed(100)
#Accuracy : 0.7442 - Confirmed + Tentative Classes
model62 <- naive_bayes(X ~ X463 + X464 + X465 + X472 + X474 + X475 + X479 + X484 + X485 + 
                       X486 + X487 + X488 + X489 + X490 + X696 + X773 + X840 + X842 + 
                       X850 + X852 + X853 + X854 + X855 + X1000 + X1001 + X1002 + 
                       X1650 + X1652 + X1656 + X1657 + X1661 + X1669 + X1671 + X1678 + 
                       X1679 + X1680 + X1684 + X1687 + X1688 + X2057 + X2105 + X2106 + 
                       X2109 + X2113 + X2114 + X2124 + X2125 + X2127 + X2155 + X2184 + 
                       X2185 + X2186 + X2207 + X2211 + X2212 + X2214 + X2215 + X2220 + 
                       X2221 + X2222 + X2224 + X2225, data= train)
pnb62 <- predict(model62, test)
confusionMatrix(pnb62, test$X)


set.seed(100)
#Accuracy : 0.6977 - Confirmed 
model51 <- naive_bayes(X ~ X463 + X464 + X465 + X475 + X479 + X485 + X486 + X487 + X488 + 
                         X489 + X490 + X696 + X773 + X840 + X850 + X852 + X853 + X854 + 
                         X855 + X1001 + X1002 + X1650 + X1652 + X1656 + X1657 + X1661 + 
                         X1669 + X1671 + X1678 + X1679 + X1680 + X1684 + X1687 + X1688 + 
                         X2057 + X2105 + X2113 + X2114 + X2125 + X2184 + X2185 + X2186 + 
                         X2207 + X2211 + X2212 + X2214 + X2215 + X2220 + X2221 + X2222 + 
                         X2224, data= train)
pnb51 <- predict(model51, test)
confusionMatrix(pnb51, test$X)

#Neural Network
set.seed(444)

#Accuracy : 0.6279 - All Classes
n <- neuralnet(X~. , data = train, hidden = c(200, 50, 9), err.fct = "ce", linear.output = FALSE)
a <- subset(test, select = -c(X))
output <- compute(n,a)
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
a <- subset(test, select = -c(X))
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
a <- subset(test, select = -c(X))
output51 <- compute(n62,a)
pnn51 <- output51$net.result
predpnn51 <- apply(pnn51,1,function(x) which(x==max(x)))
predpnn51 <- as.factor(predpnn51)
confusionMatrix(predpnn51, test$X)


#hyp tuning
runs <- tuning_run("exp.R", flags = list(dense_units = c(32,64)))


