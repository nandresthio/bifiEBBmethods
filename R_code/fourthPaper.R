source("R_code/dataProcessor.R")
data <- combineArrayResults("experimentalRunSurrogateModelWithGivenBudgetKriging", 1, 4676)
features <- read.table("data/features/sampleAndRealFeaturesClean.txt", header = TRUE, sep = " ")




# Need a way to process performance for surrogate model with budget
augmentSurrogateModelWithBudgetData <- function(data, features, relativeBudgets, costRatios, models, acquisitionFunctions, initialSampleStrategies){
  # I think the first thing is to match dimension to the data
  # Then get the budget ratio, and use that
  split <- strsplit(features$instances, ",")
  instances <- gsub("[(]", "", sapply(split, "[[", 1))
  order <- match(data$functionName, instances)
  data$dimension <- features[order, "feature_dimension"]
  data$relativeBudget <- data$budget / data$dimension
  data$usedRelativeBudget <- data$usedBudget / data$dimension
  seeds <- unique(data$seed)
  i <- 0
  for(instance in unique(data$functionName)){
    temp <- data[data$functionName == instance, ]
    for(seed in seeds){
      temp1 <- temp[temp$seed == seed, ]
      for(relativeBudget in relativeBudgets){
        for(costRatio in costRatios){
          instanceName <- paste0(instance, "_B", relativeBudget, "_Cr", costRatio)
          performances <- c()
          performanceNames <- c()
          for(model in models){
            for(acquisitionFunction in acquisitionFunctions){
              for(initialSampleStrategy in initialSampleStrategies){
                if(initialSampleStrategy == "all" & acquisitionFunction != "globalVariance"){next}
                techniqueName <- paste0(model, "_", acquisitionFunction, "_", initialSampleStrategy)
                temp2 <- temp1[temp1$method == techniqueName, ]
                if(initialSampleStrategy == "all" & acquisitionFunction == "globalVariance"){
                  techniqueName <- paste0(model, "_", initialSampleStrategy)
                }
                if(nrow(temp2) == 0){
                  print(paste0("Something weird with instance ", instanceName, " and technique ", techniqueName, ", have no rows!"))
                }
                # Look at different sampling options to get the performance
                if(initialSampleStrategy == "all"){
                  temp3 <- temp2[temp2$relativeBudget == relativeBudget, ]
                  # Should only have one row
                  if(nrow(temp3) != 1){
                    print(paste0("Something weird with instance ", instanceName, " and technique ", techniqueName, ", should have 1 row but have ", nrow(temp2), "!"))
                  }
                }else if(initialSampleStrategy == "half"){
                  temp3 <- temp2[temp2$relativeBudget == relativeBudget, ]
                  temp3 <- temp3[nrow(temp3), ]
                }else if(initialSampleStrategy == "small"){
                  # Here need to extract it!
                  temp3 <- temp2[temp2$usedRelativeBudget <= relativeBudget, ]
                  temp3 <- temp3[nrow(temp3), ]
                }
                performances <- c(performances, temp3[c("modelError", "modelCorrelation", "time")])
                performanceNames <- c(performanceNames, paste0(techniqueName, "_", c("modelError", "modelCorrelation", "time")))
              }
            }
          }
          row <- c(instanceName, seed, performances)
          names(row) <- c("instance", "seed", performanceNames)
          i <- i + 1
          print(i)
          if(i == 1){
            augmentedData <- as.data.frame(row)
          }else{
            augmentedData <- rbind(augmentedData, as.data.frame(row))
          }
        }
      }
    }
  }
  
  return(augmentedData)
}

condenseSurrogateModelWithBudgetData <- function(augmentedData, methods){
  i <- 0
  for(instance in unique(augmentedData$instance)){
    print(instance)
    temp <- augmentedData[augmentedData$instance == instance, ]
    performances <- c()
    performanceNames <- c()
    for(method in methods){
      # Really want to do a statistical comparison, but for now just take the 
      # mean and median performance
      performances <- c(performances,
                        colMeans(temp[paste0(method, c("_modelError", "_modelCorrelation", "_time"))]), 
                        median(as.matrix(temp[paste0(method, "_modelError")])),
                        median(as.matrix(temp[paste0(method, "_modelCorrelation")])),
                        median(as.matrix(temp[paste0(method, "_time")])))
      performanceNames <- c(performanceNames,
                            paste0(method, "_mean", c("_modelError", "_modelCorrelation", "_time")),
                            paste0(method, "_median", c("_modelError", "_modelCorrelation", "_time")))
    }
    row <- c(instance, performances)
    row <- matrix(row, nrow = 1)
    i <- i + 1
    print(i)
    if(i == 1){
      condensedData <- as.data.frame(row)
      colnames(condensedData) <- c("instance", performanceNames)
    }else{
      newRow <- as.data.frame(row)
      colnames(newRow) <- c("instance", performanceNames)
      condensedData <- rbind(condensedData, newRow)
    }
  }
  return(condensedData)
}

test <- augmentSurrogateModelWithBudgetData(data,
                                            features,
                                            c(5, 10, 15, 20), c(0),
                                            c("kriging"),
                                            c("variance", "globalVariance"),
                                            c("small", "half", "all"))


test2 <- condenseSurrogateModelWithBudgetData(test, 
                                              c("kriging_variance_small",
                                                "kriging_variance_half",
                                                "kriging_globalVariance_small",
                                                "kriging_globalVariance_half",
                                                "kriging_all"))
# Finally plot this!
# Need to change it into something understandable
library(ggplot2)
library(reshape2)

test2$budget <- as.numeric(gsub("B", "", sapply(strsplit(test2$instance, "_"), "[[", 2)))
test2$name <- sapply(strsplit(test2$instance, "_"), "[[", 1)
split <- strsplit(features$instances, ",")
instances <- gsub("[(]", "", sapply(split, "[[", 1))
order <- match(test2$name, instances)
test2$dimension <- features[order, "feature_dimension"]
test2$relativeBudget <- test2$budget / test2$dimension

names <- colnames(test2[str_which(colnames(test2), "mean_modelError")])


for(name in names){
  aggData <- aggregate(as.numeric(unlist(test2[name])), list(test2$budget), FUN=mean)
  aggData$method <- name
  if(name == names[[1]]){
    finalAggData <- aggData
  }else{
    finalAggData <- rbind(finalAggData, aggData)
  }
}

colnames(finalAggData) <- c("RelativeBudget", "Error", "Method")

ggplot(finalAggData, aes(RelativeBudget, Error, col=Method)) + 
  geom_line()

# Save
