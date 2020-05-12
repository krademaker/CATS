---
  title: "CATS-Cross-Validation"
output: html_notebook
---
  
  #The dataset consist of an aCGH (array comparative genomic hybridization) with 100 samples of three different breast cancer subtypes. The chromosomal copy number abberations have been measured by comparison to regular DNA and from this the specific type of 

#```{r}
#### Install packages
library(tidyr)
library(dplyr)
root_dir = getwd()
#### Loading the aCGH data
#Load all the data (log2-scaled of copy number aberations between breast cancer subtypes)
train_call <- read.delim(file=file.path(root_dir, 'data', 'train_call.tsv'), header = TRUE, sep='\t')
train_call$ID <- seq.int(nrow(train_call))
genomic_regions <- train_call[c('ID','Chromosome','Start','End','Nclone')]
train_call <- subset(train_call, select = -c(Chromosome,Start,End,Nclone))
train_call <- as.data.frame(t(subset(train_call, select = -c(ID))))
#### Clinical data
train_clinical <- read.delim(file=file.path(root_dir, 'data', 'train_clinical.tsv'), header=TRUE, sep="\t")
row.names(train_clinical) <- train_clinical$Sample
train_clinical <- train_clinical[,-(1:1), drop=FALSE]
#### Merged clinical & array data
train_data <- merge(train_call, train_clinical, by='row.names')
rownames(train_data) <- train_data$Row.names
train_data$Row.names <- NULL
train_data <- train_data %>% select(Subgroup, everything())
#```

#For this part I will provide a cross-validation scheme. We will first split the data into a test and a training set. The training set will be used to do feature selection and create the model with various classification algorithms.

#For the scheme, we will split the data set into a test (for validation of final predictor) and training set (from which we select features and build the classification model). This outer loop will be 10-fold cross-validated


#```{r}
############## Cross-Validation & Model
### Load packages
library('class')
library('gmodels')
library('caret')
library(Boruta)
library(randomForest)
### Parameters
#Outer loop CV
cv_fold_outer <- 10
training_size_outer <- 0.7
validation_size_outer <- 1 - training_size_outer
#Inner loop CV
cv_fold_features <- 10
feature_size <- 0.5
#KNN
k_number = 10
#Variable importance
threshold.VarImp <- 0.7
#Statistical-testing
default_sign_threshold <- 0.05
### Functions
#Adding clinical labels to other dataframe
attach_clinical_labels <- function(df_call, df_clinical){
  tmp_df <- merge(df_call, df_clinical, by = 'row.names')
  row.names(tmp_df) <- tmp_df$Row.names
  tmp_df$Row.names <- NULL
  return(tmp_df)
}
#Removing clinical labels from dataframe
detach_clinical_labels <- function(df){
  df$Subgroup <- NULL
  return(df)
}
# Split the input data
split <- function(train, trainingsize){
  trainIndex <- createDataPartition(train$Subgroup, p = trainingsize)
  #Save the sampleIDs with index from partition & filter out the corresponding samples for training
  resample.train <- rownames(train)[trainIndex$Resample1] 
  train_set <- filter(train, rownames(train) %in% resample.train) 
  rownames(train_set) <- resample.train 
  #Save the sampleIDs with index from partition & filter out the corresponding samples for validation
  resample.validation <- rownames(train)[-trainIndex$Resample1] 
  validation_set <- filter(train, rownames(train) %in% resample.validation) 
  rownames(validation_set) <- resample.validation 
  return(list(train_set, validation_set))
}
accuracy <- function(predictions) {
  predictions.freq <- as.data.frame(predictions[1])
  #select correct predictions (diagonal)
  correct <- sum(predictions.freq[c(1,5,9),3])
  #percentage out of total observations
  accuracy <- (correct / sum(predictions.freq$t.Freq))*100
  return(accuracy)
}
features.VarImp <- function(data, clinical_outcome){
  filter_variance_importances <- filterVarImp(data, clinical_outcome, nonpara = TRUE)
  filter_variance_importances$mean <- apply(filter_variance_importances, 1, mean)
  filter_variance_importances <- filter_variance_importances[order(-filter_variance_importances$mean),]
  filter_variance_importances$features <- rownames(filter_variance_importances)
  features <- as.numeric((filter(filter_variance_importances, filter_variance_importances$mean > threshold.VarImp))$features)
  return(features)
}
features.chi_squared <- function(data, clinical_outcome){
  chi_squared <- list()
  chi_squared_simulated <- list()
  for (feature in colnames(data)){
    # Extract feature data
    feature_data <- data[feature]
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
  bonferroni_sign_threshold <- default_sign_threshold / length(data)
  # Perform chi-squared statistical test with Monte Carlo simulations
  chi_squared.results <- data.frame(P_value = sapply(chi_squared, function(x) x$p.value), P_value_bonferroni = sapply(chi_squared, function(x) x$p.value / length(data)))
  chi_squared.significance <- data.frame(default = sapply(chi_squared.results$P_value, function(x) x < default_sign_threshold), bonferroni = sapply(chi_squared.results$P_value_bonferroni, function(x) x < bonferroni_sign_threshold))
  # Selecting the features below the bonferroni significance threshold
  chi_squared.significance$features <- rownames(chi_squared.significance)
  features <- as.numeric((filter(chi_squared.significance, chi_squared.significance$bonferroni)$features))
  return (features)
}

features.boruta <- function(z){
  print('Running Boruta..')
  boruta <- Boruta(Subgroup~ ., data = z, doTrace = 1, maxRuns = 1000)
  Tboruta <- TentativeRoughFix(boruta)
  #k <- z[Tboruta$finalDecision == "Confirmed"]
  #k$Subgroup <- z$Subgroup
  #k$Tboruta <- Tboruta
  return(Tboruta)
}
####Example with KNN
prediction <- function(training, test, features) {
  #Select columns with selected features only
  
  selected.train <- select(training, features)
  selected.test <- select(test, features)
  cl <- training$Subgroup
  
  #Make predictions
  prediction <- knn(selected.train[,-c(1)], selected.test[,-c(1)], cl, k=k_number, prob=TRUE)
  return(list('predictions'=prediction, 'test'=selected.test))
}
pred_nn <- function(z, test){
  
  set.seed(444)
  ctrl <- trainControl(method = "cv",  # corss-validation
                       number = 3,  # 10 folds
                       classProbs = FALSE , # report class probability)
  )
  nnetGrid <- expand.grid(decay = c(0.1),  size = c(5))
  maxSize <- max(nnetGrid$size)
  xx <- data.frame(subset(z, select = -c(Subgroup)))
  xxt <- data.frame(subset(test, select = -c(Subgroup)))
  
  mml <- as.list("312", "328", "347", "2778", "2779", "2782", "2783", "2784", "2785", "2786")
  for (c in mml){
    xx[,c] <- NULL
    xxt[,c] <- NULL}
  
  numWts <- 1*(maxSize * (length(xx) + 1) + maxSize + 1)
  print("NN working")
  n <- train(xx, z$Subgroup,
             method = "nnet", # train neural network using `nnet` 
             tuneGrid = nnetGrid, # tuning grid
             #preProc = c("scale"), # standardize data
             trace = FALSE,  # hide the training trace
             #trControl = ctrl,
             MaxNWts = 2658,
             maxit = 100 # maximum iteration)
  )
  print(n)
  op <- predict(n, newdata = xxt) 
    
  return(list('predictions'=op, 'test'=test))
             
  
}
#####Start program

final_accuracy_nn <- c()
### Most outer loop - divide the training & validation set
#For a 10-fold CV
for (k in 1:cv_fold_outer) {
  set.seed(k)
  #Use partition function to split data set into a training & validation set
  training_test_set <- split(train_data, training_size_outer)
  
  training_set <- training_test_set[[1]]
  test_set <- training_test_set[[2]]

  
  average_list_nn = list()
  feature_list_nn = list(feature=list())
  
  ####Inner loop - training
  for (i in 1:10) {
    set.seed(10+i)
    
    #New training sets for the feature selection from the training set
    training_test_set.feature <- split(training_set, feature_size)
    training_set.feature <- training_test_set.feature[[1]]
    test_set.feature <- training_test_set.feature[[2]]
    
    ####Method of feature selection
    #Define control (random selection)
    #Insert feature selection method
    feature_selection.results <- features.chi_squared(subset(training_set.feature, select = -c(Subgroup)), subset(training_set.feature, select = c(Subgroup)))
    feature_selection.chi_train <- select(training_set.feature, feature_selection.results)
    feature_selection.chi_test <- select(test_set.feature, feature_selection.results)
    feature_selection.chi_train$Subgroup <- training_set.feature$Subgroup
    feature_selection.chi_test$Subgroup <- test_set.feature$Subgroup
    
    
    #Tboruta <- features.boruta(training_set.feature)
    #feature_selection.boruta_train <- training_set.feature[Tboruta$finalDecision == "Confirmed"]
    #feature_selection.boruta_train$Subgroup <- training_set.feature$Subgroup
    
    #feature_selection.boruta_test <- test_set.feature[Tboruta$finalDecision == "Confirmed"]
    #feature_selection.boruta_test$Subgroup <- test_set.feature$Subgroup

    #Accuracy values of fecature selection
    total_accuracy_nn = c()
    
    for (j in 1:5){
      print(k,i,j)
      nn_model = pred_nn(feature_selection.chi_train, feature_selection.chi_test)
      print("NN finished")
      #Evaluate performance Neural Network
      score_nn = CrossTable(test_set.feature$Subgroup, nn_model$predictions, prop.chisq = FALSE)
      acc_nn = accuracy(score_nn)
      

      #Append to the vector that represents the total accuracy Neural Network
      total_accuracy_nn <- c(total_accuracy_nn, acc_nn)
    }
    #Mean of the accuracy
    mean.acc_nn <- mean(total_accuracy_nn, na.rm = TRUE) 
    
    #Append the average of accuracy and corresponding features to lists
    average_list[i] <- mean.acc_nn
    #feature_list$feature[[i]] <- feature_selection
    
  }
  avgs_vector <- unlist(average_list)
  top_index <- which.max(avgs_vector)
  top_features <- feature_list$feature[[top_index]]
  

  nn_model = pred_nn(training_set, test_set)

  
  score_nn = CrossTable(test_set$Subgroup, nn_model$predictions, prop.chisq = FALSE)
  acc_nn = accuracy(score_nn)
  final_accuracy_nn <- c(final_accuracy_nn, acc_nn)
}
mean(final_accuracy_nn, na.rm = TRUE)
```

