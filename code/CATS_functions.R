########## GLOBAL PARAMETERS
NN_MAX_ITERATIONS <- 100
TRAINING_SIZE_OUTER <- 0.7

########## LIBRARIES
library(tidyr)
library(caret)
library(dplyr)
library(here)
library(ggplot2)

########## FUNCTIONS
predict.neural_network <- function(aCGH_training, aCGH_test, decay = '', size = ''){
  # Handle optional hyperparameters
  if (decay == '' | size == ''){
    nnetGrid <- expand.grid(decay = c(0.1, 0.01, 0.001),  size = c(5:7))
  } else {
    nnetGrid <- expand.grid(decay = decay, size = size)
  }
  
  # Select data from the training and test sets
  aCGH_training.data <- data.frame(subset(aCGH_training, select = -c(Subgroup)))
  aCGH_test.data <- data.frame(subset(aCGH_test, select = -c(Subgroup)))
  
  # Set other neural network parameters
  maxSize <- max(nnetGrid$size)
  numWts <-  10 * (maxSize * (length(aCGH_training.data) + 1) + maxSize + 1)
  
  # Train the neural network
  # print(paste('Training neural network, at:', Sys.time()))
  nn_model <- train(aCGH_training.data, aCGH_training$Subgroup,
                    method = 'nnet',                             # train neural network using 'nnet'
                    tuneGrid = nnetGrid,                         # grid search hyperparameters
                    trace = FALSE,                               # hide the training trace
                    MaxNWts = numWts*10,                         # maximum number of weights of the network
                    maxit = NN_MAX_ITERATIONS                    # maximum iterations
  )
  
  # Make predictions
  nn_model.predictions <- predict(nn_model, newdata = aCGH_test.data) 
  
  return(list('predictions' = nn_model.predictions, 'test' = aCGH_test.data, 'model' = nn_model))
}


predict.knn <- function(aCGH_training, aCGH_test, k = ''){
  # Handle optional hyperparameters
  if (k == ''){
    knnGrid <- expand.grid(k = c(1:10))
  } else {
    knnGrid <- expand.grid(k = k)
  }
  
  # Select data from the training and test sets
  aCGH_training.data <- data.frame(subset(aCGH_training, select = -c(Subgroup)))
  aCGH_test.data <- data.frame(subset(aCGH_test, select = -c(Subgroup)))
  
  # Train the K-nearest neighbors classifier
  # print(paste('Training K-nearest neighbors classifier, at:', Sys.time()))
  knn_model <- train(aCGH_training.data, aCGH_training$Subgroup,
                     method = 'knn',
                     tuneGrid = knnGrid
  )
  
  # Make predictions
  knn_model.predictions <- predict(knn_model, newdata = aCGH_test.data) 
  
  return(list('predictions' = knn_model.predictions, 'test' = aCGH_test.data, 'model' = knn_model))
}


predict.nsc <- function(aCGH_training, aCGH_test, threshold = ''){
  # Handle optional hyperparameters
  if (threshold == ''){
    nscGrid <- expand.grid(threshold = c(1:8))
  } else {
    nscGrid <- expand.grid(threshold = threshold)
  }
  
  # Select data from the training and test sets
  aCGH_training.data <- data.frame(subset(aCGH_training, select = -c(Subgroup)))
  aCGH_test.data <- data.frame(subset(aCGH_test, select = -c(Subgroup)))
  
  # Train the nearest shrunken centroids classifier
  # print(paste('Training nearest shrunken centroids classifier, at:', Sys.time()))
  nsc_model <- train(aCGH_training.data, aCGH_training$Subgroup,
                     method = 'pam',
                     tuneGrid = nscGrid
  )
  
  # Make predictions
  nsc_model.predictions <- predict(nsc_model, newdata = aCGH_test.data) 
  
  return(list('predictions' = nsc_model.predictions, 'test' = aCGH_test.data, 'model' = nsc_model))
} 


split_aCGH <- function(aCGH, training_size){
  aCGH_training.index <- createDataPartition(aCGH$Subgroup, p = training_size)
  
  # split out the training set
  resample.training <- rownames(aCGH)[aCGH_training.index$Resample1] 
  aCGH_training <- filter(aCGH, rownames(aCGH) %in% resample.training) 
  rownames(aCGH_training) <- resample.training
  
  # split out the test set
  resample.test <- rownames(aCGH)[-aCGH_training.index$Resample1] 
  aCGH_test <- filter(aCGH, rownames(aCGH) %in% resample.test) 
  rownames(aCGH_test) <- resample.test
  
  return(list(aCGH_training, aCGH_test))
}