source("R_code/dataProcessor.R")
library(ggplot2)
library(stringr)

data <- combineArrayResults("experimentalRunSurrogateModelWithGivenBudgetSmallTest", 1, 1000)
for(i in 1:37){
  tempData <- combineArrayResults("experimentalRunSurrogateModelWithGivenBudgetSmallTest", 1 + i*1000, (1 + i)*1000)
  data <- rbind(data, tempData)
}
tempData <- combineArrayResults("experimentalRunSurrogateModelWithGivenBudgetSmallTest",38001, 38016)
data <- rbind(data, tempData)

features <- read.table("data/features/sampleAndRealFeaturesClean.txt", header = TRUE, sep = " ")


methods <- unique(data$method)

intervals <- seq(0,1,0.05)

generatePlottableData <- function(data, features, methods, intervals){
  # I think the first thing is to match dimension to the data
  # Then get the budget ratio, and use that
  split <- strsplit(features$instances, ",")
  instances <- gsub("[(]", "", sapply(split, "[[", 1))
  order <- match(data$functionName, instances)
  data$dimension <- features[order, "feature_dimension"]
  data$relativeBudget <- data$budget / data$dimension
  data$usedRelativeBudget <- data$usedBudget / data$dimension
  seeds <- unique(data$seed)
  costRatios <- unique(data$costRatio)
  relativeBudgets <- unique(data$relativeBudget)

  plottingData <- as.data.frame(matrix(nrow = 0, ncol = 7))
  colnames(plottingData) <- c("method",
                              "usedBudget",
                              "costRatio",
                              "relativeBudget",
                              "correlation",
                              "error",
                              "time")

  data$proportionUsedBudget <- data$usedBudget / data$budget
  # Ok let's just start with one with all instances
  for(method in methods){
    temp <- data[data$method == method, ]
    meanCorrs <- c()
    meanErrors <- c()
    meanTimes <- c()
    print(method)
    # Want to distinguish each run
    for(val in intervals){
      print(val)
      corrs <- c()
      errors <- c()
      times <- c()
      for(instance in unique(temp$functionName)){
        for(seed in seeds){
          for(costRatio in costRatios){
            for(relativeBudget in relativeBudgets){
              # This should give me a single run
              temp2 <- temp[temp$functionName == instance &
                              temp$seed == seed &
                              temp$costRatio == costRatio &
                              temp$relativeBudget == relativeBudget, ]

              if(nrow(temp2) == 0){
                print("Problem!!")
                print(method)
                print(seed)
                print(costRatio)
                print(relativeBudget)
              }
              # Go through each of the interval values and find the correct value.
              # If no value exists, take the smallest value
              temp3 <- temp2[temp2$proportionUsedBudget <= val, ]
              if(nrow(temp3) == 0){
                temp4 <- temp2[1, ]
              }else{
                temp4 <- temp3[nrow(temp3), ]
              }
              corrs <- c(corrs, temp4$modelCorrelation)
              errors <- c(errors, temp4$modelError)
              times <- c(times, temp4$time)
            }
          }
        }
      }
      meanCorrs <- c(meanCorrs, mean(corrs))
      meanErrors <- c(meanErrors, mean(errors))
      meanTimes <- c(meanTimes, mean(times))
    }
    tempPlottingData <- as.data.frame(matrix(nrow = length(intervals), ncol = 7))
    colnames(tempPlottingData) <- c("method",
                                    "usedBudget",
                                    "costRatio",
                                    "relativeBudget",
                                    "correlation",
                                    "error",
                                    "time")
    tempPlottingData$usedBudget <- intervals
    tempPlottingData$method <- method
    tempPlottingData$costRatio <- "all"
    tempPlottingData$relativeBudget <- "all"
    tempPlottingData$correlation <- meanCorrs
    tempPlottingData$error <- meanErrors
    tempPlottingData$time <- meanTimes
    plottingData <- rbind(plottingData, tempPlottingData)
  }
  
  
  # Going to repeat, first split by cost ratio, then by total budget, then both
  for(method in methods){
    for(costRatio in costRatios){
      print(paste0(method, " - ", costRatio))
      temp <- data[data$method == method &
                     data$costRatio == costRatio, ]
      meanCorrs <- c()
      meanErrors <- c()
      meanTimes <- c()
      # Want to distinguish each run
      for(val in intervals){
        print(val)
        corrs <- c()
        errors <- c()
        times <- c()
        for(instance in unique(temp$functionName)){
          for(seed in seeds){
            for(relativeBudget in relativeBudgets){
              # This should give me a single run
              temp2 <- temp[temp$functionName == instance &
                              temp$seed == seed &
                              temp$relativeBudget == relativeBudget, ]

              if(nrow(temp2) == 0){
                print("Problem!!")
                print(method)
                print(seed)
                print(relativeBudget)
              }
              # Go through each of the interval values and find the correct value.
              # If no value exists, take the smallest value
              temp3 <- temp2[temp2$proportionUsedBudget <= val, ]
              if(nrow(temp3) == 0){
                temp4 <- temp2[1, ]
              }else{
                temp4 <- temp3[nrow(temp3), ]
              }
              corrs <- c(corrs, temp4$modelCorrelation)
              errors <- c(errors, temp4$modelError)
              times <- c(times, temp4$time)
            }
          }
        }
        meanCorrs <- c(meanCorrs, mean(corrs))
        meanErrors <- c(meanErrors, mean(errors))
        meanTimes <- c(meanTimes, mean(times))
      }
      tempPlottingData <- as.data.frame(matrix(nrow = length(intervals), ncol = 7))
      colnames(tempPlottingData) <- c("method",
                                      "usedBudget",
                                      "costRatio",
                                      "relativeBudget",
                                      "correlation",
                                      "error",
                                      "time")
      tempPlottingData$usedBudget <- intervals
      tempPlottingData$method <- method
      tempPlottingData$costRatio <- costRatio
      tempPlottingData$relativeBudget <- "all"
      tempPlottingData$correlation <- meanCorrs
      tempPlottingData$error <- meanErrors
      tempPlottingData$time <- meanTimes
      plottingData <- rbind(plottingData, tempPlottingData)
    }
  }
  
  for(method in methods){
    for(relativeBudget in relativeBudgets){
      print(paste0(method, " - ", relativeBudget))
      temp <- data[data$method == method &
                     data$relativeBudget == relativeBudget, ]
      meanCorrs <- c()
      meanErrors <- c()
      meanTimes <- c()


      # Want to distinguish each run
      for(val in intervals){
        print(val)
        corrs <- c()
        errors <- c()
        times <- c()
        for(instance in unique(temp$functionName)){
          for(seed in seeds){
            for(costRatio in costRatios){
              # This should give me a single run
              temp2 <- temp[temp$functionName == instance &
                              temp$seed == seed &
                              temp$costRatio == costRatio, ]

              if(nrow(temp2) == 0){
                print("Problem!!")
                print(method)
                print(seed)
                print(costRatio)
              }
              # Go through each of the interval values and find the correct value.
              # If no value exists, take the smallest value
              temp3 <- temp2[temp2$proportionUsedBudget <= val, ]
              if(nrow(temp3) == 0){
                temp4 <- temp2[1, ]
              }else{
                temp4 <- temp3[nrow(temp3), ]
              }
              corrs <- c(corrs, temp4$modelCorrelation)
              errors <- c(errors, temp4$modelError)
              times <- c(times, temp4$time)
            }
          }
        }
        meanCorrs <- c(meanCorrs, mean(corrs))
        meanErrors <- c(meanErrors, mean(errors))
        meanTimes <- c(meanTimes, mean(times))
      }
      tempPlottingData <- as.data.frame(matrix(nrow = length(intervals), ncol = 7))
      colnames(tempPlottingData) <- c("method",
                                      "usedBudget",
                                      "costRatio",
                                      "relativeBudget",
                                      "correlation",
                                      "error",
                                      "time")
      tempPlottingData$usedBudget <- intervals
      tempPlottingData$method <- method
      tempPlottingData$costRatio <- "all"
      tempPlottingData$relativeBudget <- relativeBudget
      tempPlottingData$correlation <- meanCorrs
      tempPlottingData$error <- meanErrors
      tempPlottingData$time <- meanTimes
      plottingData <- rbind(plottingData, tempPlottingData)
    }
  }

  for(method in methods){
    for(relativeBudget in relativeBudgets){
      for(costRatio in costRatios){
        print(paste0(method, " - ", relativeBudget, " - ", costRatio))
        temp <- data[data$method == method &
                       data$relativeBudget == relativeBudget &
                       data$costRatio == costRatio, ]
        meanCorrs <- c()
        meanErrors <- c()
        meanTimes <- c()
        # Want to distinguish each run
        for(val in intervals){
          print(val)
          corrs <- c()
          errors <- c()
          times <- c()
          for(instance in unique(temp$functionName)){
            for(seed in seeds){
              # This should give me a single run
              temp2 <- temp[temp$functionName == instance &
                              temp$seed == seed, ]

              if(nrow(temp2) == 0){
                print("Problem!!")
                print(method)
                print(seed)
                print(costRatio)
              }
              # Go through each of the interval values and find the correct value.
              # If no value exists, take the smallest value
              temp3 <- temp2[temp2$proportionUsedBudget <= val, ]
              if(nrow(temp3) == 0){
                temp4 <- temp2[1, ]
              }else{
                temp4 <- temp3[nrow(temp3), ]
              }
              corrs <- c(corrs, temp4$modelCorrelation)
              errors <- c(errors, temp4$modelError)
              times <- c(times, temp4$time)
            }
          }
          meanCorrs <- c(meanCorrs, mean(corrs))
          meanErrors <- c(meanErrors, mean(errors))
          meanTimes <- c(meanTimes, mean(times))
        }
        tempPlottingData <- as.data.frame(matrix(nrow = length(intervals), ncol = 7))
        colnames(tempPlottingData) <- c("method",
                                        "usedBudget",
                                        "costRatio",
                                        "relativeBudget",
                                        "correlation",
                                        "error",
                                        "time")
        tempPlottingData$usedBudget <- intervals
        tempPlottingData$method <- method
        tempPlottingData$costRatio <- costRatio
        tempPlottingData$relativeBudget <- relativeBudget
        tempPlottingData$correlation <- meanCorrs
        tempPlottingData$error <- meanErrors
        tempPlottingData$time <- meanTimes
        plottingData <- rbind(plottingData, tempPlottingData)
      }
    }
  }
  
  
  
  return(plottingData)
}
  
test <- generatePlottableData(data, features, methods, intervals)

test$usedBudget <- as.numeric(test$usedBudget)
test$method <- as.factor(test$method)
test$costRatio <- as.factor(test$costRatio)
test$relativeBudget <- as.factor(test$relativeBudget)


temp <- test
temp <- temp[str_which(temp$method, "globalVariance_small"), ]
# temp <- temp[str_which(temp$method, "globalVariance"), ]


ggplot(temp, aes(x = usedBudget, y = correlation, color = method, pch = relativeBudget)) + 
  geom_line() + 
  geom_point()

ggplot(temp, aes(x = usedBudget, y = error, color = method, pch = relativeBudget)) + 
  geom_line() + 
  geom_point()

ggplot(test, aes(x = usedBudget, y = time, color = method, pch = method)) + 
  geom_line() + 
  geom_point()




augmentSurrogateModelWithBudgetData <- function(data, features, methods){
  split <- strsplit(features$instances, ",")
  instances <- gsub("[(]", "", sapply(split, "[[", 1))
  order <- match(data$functionName, instances)
  data$dimension <- features[order, "feature_dimension"]
  data$relativeBudget <- data$budget / data$dimension
  data$usedRelativeBudget <- data$usedBudget / data$dimension
  seeds <- unique(data$seed)
  costRatios <- unique(data$costRatio)
  relativeBudgets <- unique(data$relativeBudget)
  data$proportionUsedBudget <- data$usedBudget / data$budget
  index <- 1
  for(instance in unique(data$functionName)){
    for(relativeBudget in relativeBudgets){
      for(costRatio in costRatios){
        temp <- data[data$functionName == instance &
                       data$relativeBudget == relativeBudget &
                       data$costRatio == costRatio, ]
        
        instanceName <- paste0(instance, "_B", relativeBudget, "_Cr", costRatio)
        print(instanceName)
        
        if(nrow(temp) == 0){
          print("Didn't find data!!")
          next
        }
        # Now want to get the final corr, error and time for each seed and for
        # each method
        repData <- as.data.frame(matrix(nrow = length(seeds), ncol = 1+2*length(methods)))
        colnames(repData) <- c("seeds", paste0(methods, "_corr"), paste0(methods, "_err"))
        
        # Now need to populate
        repData$seeds <- seeds
        # Now populate every other entry
        for(method in methods){
          for(seed in seeds){
            temp2 <- temp[temp$method == method &
                            temp$seed == seed, ]
            
            if(nrow(temp2) == 0){
              print("Didn't find data!!")
              print(method)
              print(seed)
              next
            }
            
            repData[repData$seed == seed, c(paste0(method, "_corr"), paste0(method, "_err"))] <- temp2[nrow(temp2), c("modelCorrelation", "modelError")]
          }
        }
        tempData <- as.data.frame(matrix(nrow = 1, ncol = 3))
        colnames(tempData) <- c("instance", "relativeBudget", "costRatio")
        tempData$instance <- instanceName
        tempData$relativeBudget <- relativeBudget
        tempData$costRatio <- costRatio
        
        # Now that have all the information, should get the median, mean,
        # and wilcoxon tests
        for(method in methods){
          corrs <- unlist(repData[paste0(method, "_corr")])
          errs <- unlist(repData[paste0(method, "_err")])
          performances <- c(mean(corrs),
                            median(corrs),
                            mean(errs),
                            median(errs))
          minWilcoxonErr0 <- 1
          minWilcoxonErr0.001 <- 1
          minWilcoxonErr0.005 <- 1
          minWilcoxonCorr0 <- 1
          minWilcoxonCorr0.001 <- 1
          minWilcoxonCorr0.005 <- 1
          
          for(compMethod in methods){
            if(compMethod == method){next}
            compCorrs <- unlist(repData[paste0(compMethod, "_corr")])
            compErrs <- unlist(repData[paste0(compMethod, "_err")])
            # Now want to use the wilcoxon test
            hypSame <- wilcox.test(corrs, compCorrs, mu = 0, paired = FALSE, alternative = "l")
            minWilcoxonCorr0 <- min(minWilcoxonCorr0, hypSame$p.value)
            hypSame <- wilcox.test(corrs, compCorrs, mu = 0.001, paired = FALSE, alternative = "l")
            minWilcoxonCorr0.001 <- min(minWilcoxonCorr0.001, hypSame$p.value)
            hypSame <- wilcox.test(corrs, compCorrs, mu = 0.005, paired = FALSE, alternative = "l")
            minWilcoxonCorr0.005 <- min(minWilcoxonCorr0.005, hypSame$p.value)
            
            hypSame <- wilcox.test(errs, compErrs, mu = 0, paired = FALSE, alternative = "g")
            minWilcoxonErr0 <- min(minWilcoxonErr0, hypSame$p.value)
            hypSame <- wilcox.test(errs, compErrs, mu = 0.001, paired = FALSE, alternative = "g")
            minWilcoxonErr0.001 <- min(minWilcoxonErr0.001, hypSame$p.value)
            hypSame <- wilcox.test(errs, compErrs, mu = 0.005, paired = FALSE, alternative = "g")
            minWilcoxonErr0.005 <- min(minWilcoxonErr0.005, hypSame$p.value)
            
          }
          tempData[paste0(method, c("_corrMean", "_corrMedian", "_errMean", "_errMedian", 
                                    "_corrWilcoxon0", "_corrWilcoxon0.001", "_corrWilcoxon0.005",
                                    "_errWilcoxon0", "_errWilcoxon0.001", "_errWilcoxon0.005"))] <- c(performances,
                                                                                                      minWilcoxonCorr0, minWilcoxonCorr0.001, minWilcoxonCorr0.005,
                                                                                                      minWilcoxonErr0, minWilcoxonErr0.001, minWilcoxonErr0.005)
          
                            
        }
        if(index == 1){
          processedData <- tempData
        }else{
          processedData <- rbind(processedData, tempData)
        }
        index <- index + 1
      }
    }
  }
  return(processedData)
}


augmented <- augmentSurrogateModelWithBudgetData(data, features, methods)

testing <- augmented
testing <- augmented[augmented$relativeBudget == 5, ]

chosenMethods <- c("kriging_variance_half", "kriging_globalVariance_half",
                   "cokriging_variance_half", "cokriging_globalVariance_half", "cokriging_globalVarianceWithChoice_half",
                   "adaptiveCokriging_variance_half", "adaptiveCokriging_globalVariance_half", "adaptiveCokriging_globalVarianceWithChoice_half")
# Get some quick statistics
for(method in methods){
  error0.005 <- round(sum(testing[paste0(method, "_errMean")] - 0.005 <= apply(testing[str_which(colnames(testing), "errMean")], 1, FUN = min)) / nrow(testing), 3)
  error0.001 <- round(sum(testing[paste0(method, "_errMean")] - 0.001 <= apply(testing[str_which(colnames(testing), "errMean")], 1, FUN = min)) / nrow(testing), 3)
  error0.0005 <- round(sum(testing[paste0(method, "_errMean")] - 0.0005 <= apply(testing[str_which(colnames(testing), "errMean")], 1, FUN = min)) / nrow(testing), 3)
  corr0.001 <- round(sum(testing[paste0(method, "_corrMean")] + 0.001 >= apply(testing[str_which(colnames(testing), "corrMean")], 1, FUN = max)) / nrow(testing), 3)
  corr0.005 <- round(sum(testing[paste0(method, "_corrMean")] + 0.005 >= apply(testing[str_which(colnames(testing), "corrMean")], 1, FUN = max)) / nrow(testing), 3)
  corr0.01 <- round(sum(testing[paste0(method, "_corrMean")] + 0.01 >= apply(testing[str_which(colnames(testing), "corrMean")], 1, FUN = max)) / nrow(testing), 3)
  
  minTimesCorr <- round(sum(testing[paste0(method, "_corrMean")] == apply(testing[str_which(colnames(testing), "corrMean")], 1, FUN = min)) / nrow(testing), 3)
  minTimesErr <- round(sum(testing[paste0(method, "_errMean")] == apply(testing[str_which(colnames(testing), "errMean")], 1, FUN = min)) / nrow(testing), 3)
  
  print(paste0("Method ", method, " is best ", minTimesCorr, " of the time, and within 0.01, 0.005 and 0.001 of the best correlation ", corr0.01, ", ", corr0.005, " and ", corr0.001, " of the time."))
  print(paste0("Method ", method, " is best ", minTimesErr, " of the time, and within 0.001, 0.005 and 0.0001 of the best error ", error0.005, ", ", error0.001, " and ", error0.0005, " of the time."))
  
}

miniTest <- data[data$method == "adaptiveCokriging_globalVarianceWithChoice_all", ]
miniTest2 <- data[data$method == "adaptiveCokriging_globalVariance_all", ]

miniTest3 <- miniTest[miniTest$modelCorrelation != miniTest2$modelCorrelation, ]
miniTest4 <- miniTest2[miniTest$modelCorrelation != miniTest2$modelCorrelation, ]

  
  


# The next thing would be to correlate with features
instances <- sapply(strsplit(testing$instance, "_"), "[[", 1)
budgets <- as.numeric(gsub("B", "", sapply(strsplit(testing$instance, "_"), "[[", 2)))
costRatios <- as.numeric(gsub("Cr", "", sapply(strsplit(testing$instance, "_"), "[[", 3)))

# Will also want a feature subset
processedFeatures <- read.table("data/features/sampleAndRealFeaturesCleanStandarised.txt", header = TRUE, sep = " ")
processedFeatures$instanceName <- gsub("[(]", "", sapply(strsplit(processedFeatures$instances, ","), "[[", 1))
order <- match(instances, processedFeatures$instanceName)
featureSubset <- processedFeatures[order, str_which(colnames(processedFeatures), "_real_", negate = FALSE)]


tempPerfVals <- as.data.frame(matrix(nrow = nrow(testing, ncol = 0)))
tempPerfVals$instances <- testing$instance
for(method in methods){
  tempPerfVals[paste0(method, "_err")] <- testing[paste0(method, "_errMean")] - apply(testing[str_which(colnames(testing), "errMean")], 1, FUN = min)
  tempPerfVals[paste0(method, "_corr")] <- apply(testing[str_which(colnames(testing), "corrMean")], 1, FUN = max) - testing[paste0(method, "_corrMean")]
  
}

for(method in methods){
  print(method)
  if(abs(cor(tempPerfVals[paste0(method, "_err")], budgets)) > 0.3){
    print(cor(tempPerfVals[paste0(method, "_err")], budgets))
  }
  if(abs(cor(tempPerfVals[paste0(method, "_corr")], budgets)) > 0.3){
    print(cor(tempPerfVals[paste0(method, "_corr")], budgets))
  }
  if(abs(cor(tempPerfVals[paste0(method, "_err")], costRatios)) > 0.3){
    print(cor(tempPerfVals[paste0(method, "_err")], costRatios))
  }
  if(abs(cor(tempPerfVals[paste0(method, "_corr")], costRatios)) > 0.3){
    print(cor(tempPerfVals[paste0(method, "_corr")], costRatios))
  }
  for(name in colnames(featureSubset)){
    if(abs(cor(tempPerfVals[paste0(method, "_err")], featureSubset[name])) > 0.3){
      print(name)
      print(cor(tempPerfVals[paste0(method, "_err")], featureSubset[name]))
    }
  }
}



# Save
