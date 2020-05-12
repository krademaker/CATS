# Author: Chao (Cico) Zhang
# Date: 31 Mar 2017
# Usage: Rscript ~/CATS/data/predictions_deliverable/group03/model/run_model.R -i ~/CATS/data/Validation_call.txt -m ~/CATS/predictions_deliverable/group03/model/model.rds -o ~/CATS/predictions_deliverable/group03/results/predictions.txt
  # Set up R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

# Import required libraries
# You might need to load other packages here.
suppressPackageStartupMessages({
  library('getopt')
  library('caret')
  library('tibble')
  library('dplyr')
})

# Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Get options using the spec as defined by the enclosed list
# Read the options from the default: commandArgs(TRUE)
option_specification <- matrix(c(
  'input', 'i', 2, 'character',
  'model', 'm', 2, 'character',
  'output', 'o', 2, 'character'
), byrow=TRUE, ncol=4);

# Parse options
options <- getopt(option_specification);

# LOAD MODEL FROM MODEL .RDS FILE (options$model)
model <- readRDS(options$model)

# LOAD PREDICTION DATA FROM INPUT FILE (options$input)
validation_set <- read.delim(options$input, header = TRUE, sep = '\t')
validation_set <- subset(validation_set, select = -c(Chromosome,Start,End,Nclone))
validation_set.data <- as.data.frame(t(validation_set))
colnames(validation_set.data) <- paste('X', colnames(validation_set.data), sep = '')
validation_set.samples <- row.names(validation_set.data)
validation_set.sample_size <- nrow(validation_set.data)

# APPLY MODEL TO VALIDATION SET TO DO PREDICTIONS
validation_set.predictions <- predict(model, newdata = validation_set.data)

# WRITE PREDICTIONS TO OUTPUT FILE (options$output)
sink(options$output)
cat('"Sample"\t"Subgroup"', sep = '\n', file = options$output)
for (s in 1:validation_set.sample_size){
  cat(paste('"', validation_set.samples[s], '"\t"', validation_set.predictions[s], '"', sep = ''), sep = '\n', file = options$output, append=TRUE)
}
sink()

message ("Done!")