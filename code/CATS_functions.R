########## LIBRARIES
# install.packages(c('readxl', 'tidyverse', 'ggplot2', 'kableExtra'), dependencies=TRUE)

#### Load packages
library('tidyr')
library('dplyr')
library('here')
library('class')
library('gmodels')
library('caret')
library('Boruta')
library('randomForest')
library('ggplot2')
library('tibble')
library('plotly')
library('ggdendro')
library('reshape2')
library('grid')
library('cluster')
library('fossil')
library('kohonen')
library('heatmaply')
library('ggcorrplot')
library('biomaRt')
library('ggpubr')
library('plyr')


########## GLOBAL PARAMETERS
NN_MAX_ITERATIONS <- 100 # maximum number of neural network iterations
TRAINING_SIZE_OUTER <- 0.7 # training set size of outer loop
CV_FOLDS_OUTER <- 10 # number of outer loop folds
TRAINING_SIZE_INNER <- 0.7 # training set size of inner loop
CV_FOLDS_INNER <- 10 # # number of innter loop folds
VARIMP_THRESHOLD <- 0.7 # Variable importance  threshold
DEFAULT_SIGNIFICANCE_THRESHOLD <- 0.05 # Statistical testing threshold
BORUTA_MAX_RUNS <- 1000 # Maximum number of Boruta runs
N_PERMUTATIONS <- 1000 # Number of permutations
models.SUBINDEX_MODEL <- 1 # Subindex of models in 'model' object
models.SUBINDEX_TYPE <- 2 # Subindex of model type in 'model' object


########## FUNCTIONS
attach_clinical_labels <- function(df_call, df_clinical){
  # Adding clinical labels to other dataframe
  tmp_df <- merge(df_call, df_clinical, by = 'row.names')
  row.names(tmp_df) <- tmp_df$Row.names
  tmp_df$Row.names <- NULL
  return(tmp_df)
}


detach_clinical_labels <- function(df){
  # Removing clinical labels from dataframe
  df$Subgroup <- NULL
  return(df)
}


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


feature_selection.VarImp <- function(aCGH, clinical_outcome){
  print(paste('Running filterVarImp feature selection, at:', Sys.time()))
  filter_variance_importances <- filterVarImp(aCGH, clinical_outcome, nonpara = TRUE)
  filter_variance_importances$mean <- apply(filter_variance_importances, 1, mean)
  filter_variance_importances <- filter_variance_importances[order(-filter_variance_importances$mean),]
  filter_variance_importances$features <- rownames(filter_variance_importances)
  features <- as.numeric((filter(filter_variance_importances, filter_variance_importances$mean > VARIMP_THRESHOLD))$features)
  return(features)
}


feature_selection.chi_squared <- function(aCGH, clinical_outcome){
  # Select features using Bonferroni-significant results of chi-squared statistical test
  print(paste('Running Pearson chi-squared test feature selection, at:', Sys.time()))
  chi_squared <- list()
  chi_squared_simulated <- list()
  
  for (feature in colnames(aCGH)){
    # Extract feature data
    feature_data <- aCGH[feature]
    colnames(feature_data) <- 'copy_number'
    
    # Attach clinical labels
    feature_data <- attach_clinical_labels(feature_data, clinical_outcome)

      # tally counts of all copy number-subgroup combinations
    feature_tally <- group_by(feature_data, copy_number, Subgroup) %>% tally(name = 'count')
    
    # Construct contigency table
    feature_contingency_table <- xtabs(count~Subgroup+copy_number, data = feature_tally)
  
    # Perform Pearson's chi-squared test
    chi_squared[[feature]] <- chisq.test(feature_contingency_table, simulate.p.value = TRUE)
  }
  # Set bonferroni significance threshold
  bonferroni_sign_threshold <- DEFAULT_SIGNIFICANCE_THRESHOLD / length(aCGH)
  
  # Perform chi-squared statistical test with Monte Carlo simulations
  chi_squared.results <- data.frame(P_value = sapply(chi_squared, function(x) x$p.value), P_value_bonferroni = sapply(chi_squared, function(x) x$p.value / length(aCGH)))
  chi_squared.significance <- data.frame(default = sapply(chi_squared.results$P_value, function(x) x < DEFAULT_SIGNIFICANCE_THRESHOLD), bonferroni = sapply(chi_squared.results$P_value_bonferroni, function(x) x < bonferroni_sign_threshold))
  
  # Selecting the features below the bonferroni significance threshold
  chi_squared.significance$features <- rownames(chi_squared.significance)
  features <- as.numeric((filter(chi_squared.significance, chi_squared.significance$bonferroni)$features))
  return(features)
}


feature_selection.boruta <- function(aCGH){
  # Select features using the Boruta algorithm
  print(paste('Running Borata feature selection, at:', Sys.time()))
  boruta <- Boruta(Subgroup~ .,
                   data = aCGH,
                   doTrace = 1,
                   maxRuns = BORUTA_MAX_RUNS)
  boruta.tentative_fixed <- TentativeRoughFix(boruta)
  return(boruta.tentative_fixed)
}


get_hyperparameters <- function(model){
  # Get the hyperparameters of a model
  hyperparameters <- models[[model]][[models.SUBINDEX_MODEL]]$bestTune
  return(hyperparameters)
}


get_n_features <- function(model, model_type){
  # Get the number of features for a specific model type
  if(model_type == 'neural network'){
    n_features <- length(models[[model]][[models.SUBINDEX_MODEL]]$finalModel$coefnames)
  }
  if(model_type == 'k-nearest neighbors'){
    n_features <- length(models[[model]][[models.SUBINDEX_MODEL]]$finalModel$xNames)
  }
  if(m.type == 'nearest shrunken centroid'){
    n_features <- length(models[[model]][[models.SUBINDEX_MODEL]]$finalModel$xNames)
  }
  return(n_features)
}