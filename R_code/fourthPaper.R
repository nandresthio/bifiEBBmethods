source("R_code/dataProcessor.R")
library(ggplot2)
library(stringr)
library(ggpubr)
library(grid)

data <- combineArrayResults("experimentalRunSurrogateModelWithGivenBudgetHalfOnly", 1, 100)
for(i in 1:99){
  tempData <- combineArrayResults("experimentalRunSurrogateModelWithGivenBudgetHalfOnly", 1 + i*100, (1 + i)*100)
  data <- rbind(data, tempData)
}
data <- data[data$seed <= 5, ]

for(i in 100:246){
  tempData <- combineArrayResults("experimentalRunSurrogateModelWithGivenBudgetHalfOnlySmall", i*100, i*100 + 99)
  data <- rbind(data, tempData)
}
tempData <- combineArrayResults("experimentalRunSurrogateModelWithGivenBudgetHalfOnlySmall", 24700, 24752)
data <- rbind(data, tempData)



features <- read.table("data/features/sampleAndRealFeaturesClean.txt", header = TRUE, sep = " ")
processedFeatures <- read.table("data/features/sampleAndRealFeaturesCleanStandarised.txt", header = TRUE, sep = " ")




augmentSurrogateModelWithBudgetData <- function(givenData, givenFeatures, givenMethods){
  split <- strsplit(givenFeatures$instances, ",")
  instances <- gsub("[(]", "", sapply(split, "[[", 1))
  order <- match(givenData$functionName, instances)
  givenData$dimension <- givenFeatures[order, "feature_dimension"]
  givenData$relativeBudget <- givenData$budget / givenData$dimension
  givenData$usedRelativeBudget <- givenData$usedBudget / givenData$dimension
  seeds <- unique(givenData$seed)
  costRatios <- unique(givenData$costRatio)
  relativeBudgets <- unique(givenData$relativeBudget)
  givenData$proportionUsedBudget <- givenData$usedBudget / givenData$budget
  
  # Going to remove entries where went past available budget
  givenData <- givenData[givenData$proportionUsedBudget <= 1, ]
  index <- 1
  for(instance in unique(givenData$functionName)){
    print(instance)
    temp_1 <- givenData[givenData$functionName == instance, ]
    for(relativeBudget in relativeBudgets){
      temp_2 <- temp_1[temp_1$relativeBudget == relativeBudget, ]
      for(costRatio in costRatios){
        temp <- temp_2[temp_2$costRatio == costRatio, ]
        
        instanceName <- paste0(instance, "_B", relativeBudget, "_Cr", costRatio)
        # print(instanceName)
        
        if(nrow(temp) == 0){
          print("Didn't find data!!")
          print(instanceName)
          next
        }
        # Now want to get the final corr, error and time for each seed and for
        # each method
        repData <- as.data.frame(matrix(nrow = length(seeds), ncol = 1+2*length(givenMethods)))
        colnames(repData) <- c("seeds", paste0(givenMethods, "_corr"), paste0(givenMethods, "_err"))
        
        # Now need to populate
        repData$seeds <- seeds
        # Now populate every other entry
        for(method in givenMethods){
          for(seed in seeds){
            temp2 <- temp[temp$method == method &
                            temp$seed == seed, ]
            
            if(nrow(temp2) == 0){
              print("Didn't find data!!")
              print(method)
              print(seed)
              next
            }
            
            # Add a check here to see if all the budget was used up!
            if(temp2[nrow(temp2), "budget"] - temp2[nrow(temp2), "usedBudget"] > 1){
              print("Didn't use full budget!!")
              print(method)
              print(seed)
              next
            }
            repData[repData$seed == seed, c(paste0(method, "_corr"), paste0(method, "_err"))] <- temp2[nrow(temp2), c("modelCorrelation", "modelError")]
          }
        }
        for(seed in seeds){
          repData[repData$seed == seed, "bestCorr"] <- max(repData[repData$seed == seed, str_which(colnames(repData), "_corr")])
          repData[repData$seed == seed, "bestErr"] <- min(repData[repData$seed == seed, str_which(colnames(repData), "_err")])
          
        }
        # repData$bestErr <- max(repData[str_which(colnames(repData), "_corr")])
        # repData$bestErr <- min(repData[str_which(colnames(repData), "_err")])
        for(method in givenMethods){
          repData[paste0(method, "_corrDiff")] <- repData$bestCorr - repData[paste0(method, "_corr")]
          repData[paste0(method, "_errDiff")] <- repData$bestErr - repData[paste0(method, "_err")]
          
        }
        # print(repData)
        tempData <- as.data.frame(matrix(nrow = 1, ncol = 4))
        colnames(tempData) <- c("instance", "dimension", "relativeBudget", "costRatio")
        tempData$instance <- instanceName
        tempData$dimension <- temp[1, "dimension"]
        tempData$relativeBudget <- relativeBudget
        tempData$costRatio <- costRatio
        
        # Now that have all the information, should get the median, mean,
        # and wilcoxon tests
        for(method in givenMethods){
          corrs <- unlist(repData[paste0(method, "_corr")])
          errs <- unlist(repData[paste0(method, "_err")])
          corrsDiff <- unlist(repData[paste0(method, "_corrDiff")])
          errsDiff <- unlist(repData[paste0(method, "_errDiff")])
          performances <- c(mean(corrs),
                            median(corrs),
                            mean(errs),
                            median(errs),
                            mean(corrsDiff),
                            median(corrsDiff),
                            mean(errsDiff),
                            median(errsDiff))
          minWilcoxonErr0 <- 1
          minWilcoxonErr0.001 <- 1
          minWilcoxonErr0.005 <- 1
          minWilcoxonCorr0 <- 1
          minWilcoxonCorr0.005 <- 1
          minWilcoxonCorr0.01 <- 1
          
          for(compMethod in givenMethods){
            if(compMethod == method){next}
            compCorrs <- unlist(repData[paste0(compMethod, "_corr")])
            compErrs <- unlist(repData[paste0(compMethod, "_err")])
            # Now want to use the wilcoxon test
            hypSame <- wilcox.test(corrs, compCorrs, mu = 0, paired = FALSE, alternative = "l")
            minWilcoxonCorr0 <- min(minWilcoxonCorr0, hypSame$p.value)
            hypSame <- wilcox.test(corrs, compCorrs, mu = -0.005, paired = FALSE, alternative = "l")
            minWilcoxonCorr0.005 <- min(minWilcoxonCorr0.005, hypSame$p.value)
            hypSame <- wilcox.test(corrs, compCorrs, mu = -0.01, paired = FALSE, alternative = "l")
            minWilcoxonCorr0.01 <- min(minWilcoxonCorr0.01, hypSame$p.value)
            
            hypSame <- wilcox.test(errs, compErrs, mu = 0, paired = FALSE, alternative = "g")
            minWilcoxonErr0 <- min(minWilcoxonErr0, hypSame$p.value)
            hypSame <- wilcox.test(errs, compErrs, mu = 0.001, paired = FALSE, alternative = "g")
            minWilcoxonErr0.001 <- min(minWilcoxonErr0.001, hypSame$p.value)
            hypSame <- wilcox.test(errs, compErrs, mu = 0.005, paired = FALSE, alternative = "g")
            minWilcoxonErr0.005 <- min(minWilcoxonErr0.005, hypSame$p.value)
            
          }
          tempData[paste0(method, c("_corrMean", "_corrMedian", "_errMean", "_errMedian",
                                    "_corrDiffMean", "_corrDiffMedian", "_errDiffMean", "_errDiffMedian",
                                    "_corrWilcoxon0", "_corrWilcoxon0.005", "_corrWilcoxon0.01",
                                    "_errWilcoxon0", "_errWilcoxon0.001", "_errWilcoxon0.005"))] <- c(performances,
                                                                                                      minWilcoxonCorr0, minWilcoxonCorr0.005, minWilcoxonCorr0.01,
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


workingData <- data[str_which(data$functionName, "SOLAR", negate = TRUE), ]
workingData <- data[str_which(data$functionName, "SongToalForretal", negate = FALSE), ]

chosenMethods <- c("kriging_globalVariance_half"
                   ,"cokriging_globalVariance_half"
                   , "cokriging_globalVarianceWithChoice_half"
                   , "adaptiveCokriging_globalVariance_half"
                   , "adaptiveCokriging_globalVarianceWithChoice_half"
                   , "adaptiveCokrigingAdvanced_globalVariance_half"
                   , "adaptiveCokrigingAdvanced_globalVarianceWithChoice_half"
)

augmented <- augmentSurrogateModelWithBudgetData(workingData, features, chosenMethods)

chosenAugmented <- augmented

split <- strsplit(processedFeatures$instances, ",")
functionNames <- gsub("[(]", "", sapply(split, "[[", 1))
augFunctionNames <- sapply(strsplit(chosenAugmented$instance, "_"), "[[", 1)

order <- match(augFunctionNames, functionNames)
featNames <- colnames(processedFeatures[str_which(colnames(processedFeatures), "feature_real_")])
chosenAugmented[featNames] <- processedFeatures[order, featNames]

chosenAugmented["feature_real_costRatio"] <- log(chosenAugmented$costRatio)
chosenAugmented["feature_real_costRatio"] <- (chosenAugmented["feature_real_costRatio"] - min(chosenAugmented["feature_real_costRatio"])) / (max(chosenAugmented["feature_real_costRatio"]) - min(chosenAugmented["feature_real_costRatio"]))
chosenAugmented["feature_real_costRatio"] <- chosenAugmented["feature_real_costRatio"] * 4 - 2

chosenAugmented["feature_real_budget"] <- chosenAugmented$relativeBudget
chosenAugmented["feature_real_budget"] <- (chosenAugmented["feature_real_budget"] - min(chosenAugmented["feature_real_budget"])) / (max(chosenAugmented["feature_real_budget"]) - min(chosenAugmented["feature_real_budget"]))
chosenAugmented["feature_real_budget"] <- chosenAugmented["feature_real_budget"] * 4 - 2



featNames <- colnames(chosenAugmented[str_which(colnames(chosenAugmented), "feature_real_")])
# featNames <- c('feature_real_costRatio', "feature_real_budget")

for(feat in featNames){
  for(method in chosenMethods){
    # print(method)
    if(abs(cor(chosenAugmented[feat], chosenAugmented[paste0(method, "_corrDiffMean")])) > 0.3){
      print(abs(cor(chosenAugmented[feat], chosenAugmented[paste0(method, "_corrDiffMean")])))
    }
    if(abs(cor(chosenAugmented[feat], chosenAugmented[paste0(method, "_errDiffMean")])) > 0.3){
      print(abs(cor(chosenAugmented[feat], chosenAugmented[paste0(method, "_errMean")])))
    }
  }
}
    
    if(abs(cor(chosenAugmented[paste0(method, "_errMean")], budgets)) > 0.3){
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
}








testing <- augmentedNoChoice
testing <- testing[rowSums(is.na(testing)) == 0, ]
testing <- testing[c(1:3,
                     str_which(colnames(testing), "globalVariance_"))]
# testing <- testing[c(1:3, str_which(colnames(testing), "globalVariance_half"))]
testing <- testing[c(1:3, str_which(colnames(testing), "corrWilcoxon0.01"))]
# testing <- testing[c(1:3, str_which(colnames(testing), "half"))]
# testing$Best <- apply(testing[str_which(colnames(testing), "_corrMedian")], 1, FUN = max)
testing$Best <- apply(testing[str_which(colnames(testing), "_corrWilcoxon0.01")], 1, FUN = max)

testing <- testing[testing$adaptiveCokriging_globalVariance_half_corrMean < testing$cokriging_globalVariance_half_corrMean &
                     testing$adaptiveCokriging_globalVariance_half_corrMean < testing$kriging_globalVariance_half_corrMean, ]
colSums(testing[4:(ncol(testing)-1)] >= 0.5)

testing[4:(ncol(testing)-1)] <- testing$Best - testing[4:(ncol(testing)-1)]



secondTesting <- augmentedNoChoice[testing$adaptiveCokriging_globalVariance_half_corrWilcoxon0.01 < 0.5, ]
secondTesting <- secondTesting[c(1:3,
                                 str_which(colnames(secondTesting), "corrMean"))]

secondTesting <- augmentedNoChoice[testing$kriging_globalVariance_half_corrWilcoxon0.01 >= 0.5 & 
                                     testing$cokriging_globalVariance_half_corrWilcoxon0.01 < 0.5 &
                                     testing$adaptiveCokriging_globalVariance_half_corrWilcoxon0.01 >= 0.5, ]

secondTesting <- augmentedNoChoice[augmentedNoChoice$kriging_globalVariance_half_corrMean > augmentedNoChoice$adaptiveCokriging_globalVariance_half_corrMean & 
                                     augmentedNoChoice$cokriging_globalVariance_half_corrMean > augmentedNoChoice$adaptiveCokriging_globalVariance_half_corrMean, ]

secondTesting <- secondTesting[c(1:3,
                                 str_which(colnames(secondTesting), "corrMean"))]


testing[4:(ncol(testing)-1)] <- testing$Best - testing[4:(ncol(testing)-1)]
colSums(testing[4:(ncol(testing)-1)] < 0.01)


DisturbanceBasedFunction2-seed1-dists2-centres5-radius0.05-freq5-amp0.1_B20_Cr0.1

testingData <- data
testingData <- testingData[testingData$functionName == "DisturbanceBasedFunction2-seed1-dists2-centres5-radius0.05-freq5-amp0.1", ]
testingData <- testingData[str_which(testingData$method, "_globalVariance_"), ]
testingData <- testingData[testingData$costRatio == 0.1, ]
testingData <- testingData[testingData$budget == 20, ]



testingData <- testingData[testingData$usedBudget >= 19, ]

testingData <- testingData[testingData$seed == 9, ]
testingData <- testingData[testingData$usedBudget >= 29, ]


# Get some quick statistics
for(method in chosenMethods){
  error0.005 <- round(sum(testing[paste0(method, "_errMean")] - 0.005 <= apply(testing[str_which(colnames(testing), "errMean")], 1, FUN = min)) / nrow(testing), 3)
  error0.001 <- round(sum(testing[paste0(method, "_errMean")] - 0.001 <= apply(testing[str_which(colnames(testing), "errMean")], 1, FUN = min)) / nrow(testing), 3)
  error0.0005 <- round(sum(testing[paste0(method, "_errMean")] - 0.0005 <= apply(testing[str_which(colnames(testing), "errMean")], 1, FUN = min)) / nrow(testing), 3)
  corr0.001 <- round(sum(testing[paste0(method, "_corrMean")] + 0.001 >= apply(testing[str_which(colnames(testing), "corrMean")], 1, FUN = max)) / nrow(testing), 3)
  corr0.005 <- round(sum(testing[paste0(method, "_corrMean")] + 0.005 >= apply(testing[str_which(colnames(testing), "corrMean")], 1, FUN = max)) / nrow(testing), 3)
  corr0.01 <- round(sum(testing[paste0(method, "_corrMean")] + 0.01 >= apply(testing[str_which(colnames(testing), "corrMean")], 1, FUN = max)) / nrow(testing), 3)
  
  minTimesCorr <- round(sum(testing[paste0(method, "_corrMean")] == apply(testing[str_which(colnames(testing), "corrMean")], 1, FUN = max)) / nrow(testing), 3)
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
































methods <- unique(data$method)



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

  plottingData <- as.data.frame(matrix(nrow = 0, ncol = 9))
  colnames(plottingData) <- c("method",
                              "usedBudget",
                              "costRatio",
                              "relativeBudget",
                              "correlationMean",
                              "correlationMedian",
                              "errorMean",
                              "errorMedian",
                              "time")

  data$proportionUsedBudget <- data$usedBudget / data$budget
  # Ok let's just start with one with all instances
  for(method in methods){
    temp <- data[data$method == method, ]
    meanCorrs <- c()
    meanErrors <- c()
    medianCorrs <- c()
    medianErrors <- c()
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
      
      medianCorrs <- c(medianCorrs, median(corrs))
      medianErrors <- c(medianErrors, median(errors))
      
      meanTimes <- c(meanTimes, mean(times))
    }
    tempPlottingData <- as.data.frame(matrix(nrow = length(intervals), ncol = 9))
    colnames(tempPlottingData) <- c("method",
                                    "usedBudget",
                                    "costRatio",
                                    "relativeBudget",
                                    "correlationMean",
                                    "correlationMedian",
                                    "errorMean",
                                    "errorMedian",
                                    "time")
    tempPlottingData$usedBudget <- intervals
    tempPlottingData$method <- method
    tempPlottingData$costRatio <- "all"
    tempPlottingData$relativeBudget <- "all"
    tempPlottingData$correlationMean <- meanCorrs
    tempPlottingData$errorMean <- meanErrors
    tempPlottingData$correlationMedian <- medianCorrs
    tempPlottingData$errorMedian <- medianErrors
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
      medianCorrs <- c()
      medianErrors <- c()
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
        medianCorrs <- c(medianCorrs, median(corrs))
        medianErrors <- c(medianErrors, median(errors))
        meanTimes <- c(meanTimes, mean(times))
      }
      tempPlottingData <- as.data.frame(matrix(nrow = length(intervals), ncol = 9))
      colnames(tempPlottingData) <- c("method",
                                      "usedBudget",
                                      "costRatio",
                                      "relativeBudget",
                                      "correlationMean",
                                      "correlationMedian",
                                      "errorMean",
                                      "errorMedian",
                                      "time")
      tempPlottingData$usedBudget <- intervals
      tempPlottingData$method <- method
      tempPlottingData$costRatio <- costRatio
      tempPlottingData$relativeBudget <- "all"
      tempPlottingData$correlationMean <- meanCorrs
      tempPlottingData$errorMean <- meanErrors
      tempPlottingData$correlationMedian <- medianCorrs
      tempPlottingData$errorMedian <- medianErrors
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
      medianCorrs <- c()
      medianErrors <- c()
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
        medianCorrs <- c(medianCorrs, median(corrs))
        medianErrors <- c(medianErrors, median(errors))
        meanTimes <- c(meanTimes, mean(times))
      }
      tempPlottingData <- as.data.frame(matrix(nrow = length(intervals), ncol = 9))
      colnames(tempPlottingData) <- c("method",
                                      "usedBudget",
                                      "costRatio",
                                      "relativeBudget",
                                      "correlationMean",
                                      "correlationMedian",
                                      "errorMean",
                                      "errorMedian",
                                      "time")
      tempPlottingData$usedBudget <- intervals
      tempPlottingData$method <- method
      tempPlottingData$costRatio <- "all"
      tempPlottingData$relativeBudget <- relativeBudget
      tempPlottingData$correlationMean <- meanCorrs
      tempPlottingData$errorMean <- meanErrors
      tempPlottingData$correlationMedian <- medianCorrs
      tempPlottingData$errorMedian <- medianErrors
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
        medianCorrs <- c()
        medianErrors <- c()
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
          medianCorrs <- c(medianCorrs, median(corrs))
          medianErrors <- c(medianErrors, median(errors))
          meanTimes <- c(meanTimes, mean(times))
        }
        tempPlottingData <- as.data.frame(matrix(nrow = length(intervals), ncol = 9))
        colnames(tempPlottingData) <- c("method",
                                        "usedBudget",
                                        "costRatio",
                                        "relativeBudget",
                                        "correlationMean",
                                        "correlationMedian",
                                        "errorMean",
                                        "errorMedian",
                                        "time")
        tempPlottingData$usedBudget <- intervals
        tempPlottingData$method <- method
        tempPlottingData$costRatio <- costRatio
        tempPlottingData$relativeBudget <- relativeBudget
        tempPlottingData$correlationMean <- meanCorrs
        tempPlottingData$correlationMedian <- medianCorrs
        tempPlottingData$errorMean <- meanErrors
        tempPlottingData$errorMedian <- medianErrors
        tempPlottingData$time <- meanTimes
        plottingData <- rbind(plottingData, tempPlottingData)
      }
    }
  }
  
  
  
  return(plottingData)
}
 
intervals <- seq(0.5,1,0.05)
chosenMethodsNoChoice <- c("kriging_globalVariance_half",
                           "cokriging_globalVariance_half",
                           "cokriging_globalVarianceWithChoice_half",
                           "adaptiveCokriging_globalVariance_half",
                           "adaptiveCokriging_globalVarianceWithChoice_half",
                           "adaptiveCokrigingAdvanced_globalVariance_half",
                           "adaptiveCokrigingAdvanced_globalVarianceWithChoice_half")
 

smallPlottableData <- generatePlottableData(data[data$functionName == "RajnarayanWoods", ], features, chosenMethodsNoChoice, intervals)

plottableData <- generatePlottableData(data, features, chosenMethodsNoChoice, intervals)

plottableData$usedBudget <- as.numeric(plottableData$usedBudget)
plottableData$method <- as.factor(plottableData$method)
plottableData$costRatio <- as.factor(plottableData$costRatio)
plottableData$relativeBudget <- as.factor(plottableData$relativeBudget)

temp <- plottableData[str_which(plottableData$method, "adaptiveCokriging_", negate = TRUE), ]

test <- temp[!complete.cases(temp), ]


ggplot(temp[temp$costRatio == 0.01 & temp$relativeBudget == 5, ], aes(x = usedBudget, y = correlationMedian, color = method, pch = relativeBudget)) + 
  geom_line() + 
  geom_point()

ggplot(temp[temp$costRatio == 0.1 & temp$relativeBudget == 10, ], aes(x = usedBudget, y = errorMedian, color = method, pch = relativeBudget)) + 
  geom_line() + 
  geom_point()

ggplot(test, aes(x = usedBudget, y = time, color = method, pch = method)) + 
  geom_line() + 
  geom_point()

testingData <- data
# testingData <- testingData[testingData$functionName == "SongToalForretal0.90", ]
testingData <- testingData[str_which(testingData$method, "cokriging_globalVariance_"), ]
testingData <- testingData[testingData$costRatio == 0.01, ]
testingData <- testingData[testingData$budget == 20, ]
testingData <- testingData[testingData$modelCorrelation < 0.1, ]
testingData <- testingData[testingData$usedBudget - testingData$budget > -1, ]

testingData <- testingData[testingData$seed == 1, ]





augmentSurrogateModelWithBudgetDataPerSeed <- function(givenData, givenFeatures, givenMethods){
  split <- strsplit(givenFeatures$instances, ",")
  instances <- gsub("[(]", "", sapply(split, "[[", 1))
  order <- match(givenData$functionName, instances)
  givenData$dimension <- givenFeatures[order, "feature_dimension"]
  givenData$relativeBudget <- givenData$budget / givenData$dimension
  givenData$usedRelativeBudget <- givenData$usedBudget / givenData$dimension
  seeds <- unique(givenData$seed)
  costRatios <- unique(givenData$costRatio)
  relativeBudgets <- unique(givenData$relativeBudget)
  givenData$proportionUsedBudget <- givenData$usedBudget / givenData$budget
  
  # Going to remove entries where went past available budget
  givenData <- givenData[givenData$proportionUsedBudget <= 1, ]
  
  
  index <- 1
  for(instance in unique(givenData$functionName)){
    print(instance)
    temp_1 <- givenData[givenData$functionName == instance, ]
    for(relativeBudget in relativeBudgets){
      temp_2 <- temp_1[temp_1$relativeBudget == relativeBudget, ]
      for(costRatio in costRatios){
        temp_3 <- temp_2[temp_2$costRatio == costRatio, ]
        for(seed in seeds){
          temp <- temp_3[temp_3$seed == seed, ]
        
          instanceName <- paste0(instance, "_B", relativeBudget, "_Cr", costRatio, "_S", seed)
          # print(instanceName)
        
          if(nrow(temp) == 0){
            print("Didn't find data!!")
            print(instanceName)
            next
          }
          # Now want to get the final corr, error and time for each seed and for
          # each method
          vals <- c()
          names <- c()
          
          # Now populate every other entry
          for(method in givenMethods){
            temp2 <- temp[temp$method == method, ]
              
            if(nrow(temp2) == 0){
              print("Didn't find data!!")
              print(method)
              print(seed)
              next
            }
            
            # Add a check here to see if all the budget was used up!
            if(temp2[nrow(temp2), "budget"] - temp2[nrow(temp2), "usedBudget"] > 1){
              print("Didn't use full budget!!")
              print(method)
              print(seed)
              next
            }
            vals <- c(vals, temp2[nrow(temp2), c("modelCorrelation", "modelError")])
            names <- c(names, paste0(method, "_corr"), paste0(method, "_err"))
          
          }
          
          tempData <- as.data.frame(matrix(nrow = 1, ncol = 5))
          colnames(tempData) <- c("instance", "dimension", "relativeBudget", "costRatio", "seed")
          tempData$instance <- instanceName
          tempData$dimension <- temp[1, "dimension"]
          tempData$relativeBudget <- relativeBudget
          tempData$costRatio <- costRatio
          tempData$seed <- seed
          
          tempData[names] <- vals
          
          if(index == 1){
            processedData <- tempData
          }else{
            processedData <- rbind(processedData, tempData)
          }
          index <- index + 1
        }
      }
    }
  }
  return(processedData)
}

chosenMethods <- c("kriging_globalVariance_half",
                   "cokriging_globalVariance_half",
                   "cokriging_globalVarianceWithChoice_half",
                   "adaptiveCokriging_globalVariance_half",
                   "adaptiveCokriging_globalVarianceWithChoice_half",
                   "adaptiveCokrigingAdvanced_globalVariance_half",
                   "adaptiveCokrigingAdvanced_globalVarianceWithChoice_half")

workingData <- data[str_which(data$functionName, "SOLAR", negate = TRUE), ]

augmentedNoChoicePerSeed <- augmentSurrogateModelWithBudgetDataPerSeed(workingData, features, chosenMethods)

# Let's try a box plot
# Will want to do y is corr, x is cost ratio, group by technique
plotting <- augmentedNoChoicePerSeed
plotting <- plotting[plotting$relativeBudget == 5, ]
plotting <- plotting[plotting$costRatio == 5, ]

plotting$costRatio <- as.factor(plotting$costRatio)
ggplot(plotting, aes(x=costRatio, y=correlation, fill=method)) +
  geom_boxplot()





test <- augmentedNoChoicePerSeed
test <- test[c(1:5, str_which(colnames(test), "globalVariance_half_corr"))]

chosen <- test[test$kriging_globalVariance_half_corr > test$adaptiveCokriging_globalVariance_half_corr &
               test$cokriging_globalVariance_half_corr > test$adaptiveCokriging_globalVariance_half_corr, ]

nrow(chosen) / nrow(test)

chosen <- test[test$kriging_globalVariance_half_corr > test$adaptiveCokrigingAdvanced_globalVariance_half_corr &
               test$cokriging_globalVariance_half_corr > test$adaptiveCokrigingAdvanced_globalVariance_half_corr, ]


nrow(test[test$kriging_globalVariance_half_corr > test$adaptiveCokriging_globalVariance_half_corr &
                 test$cokriging_globalVariance_half_corr > test$adaptiveCokriging_globalVariance_half_corr, ]) / nrow(test)

nrow(test[test$kriging_globalVariance_half_corr > test$adaptiveCokrigingAdvanced_globalVariance_half_corr &
            test$cokriging_globalVariance_half_corr > test$adaptiveCokrigingAdvanced_globalVariance_half_corr, ]) / nrow(test)


nrow(test[test$kriging_globalVariance_half_corr > test$adaptiveCokriging_globalVarianceWithChoice_half_corr &
            test$cokriging_globalVarianceWithChoice_half_corr > test$adaptiveCokriging_globalVarianceWithChoice_half_corr, ]) / nrow(test)

nrow(test[test$kriging_globalVariance_half_corr > test$adaptiveCokrigingAdvanced_globalVarianceWithChoice_half_corr &
            test$cokriging_globalVarianceWithChoice_half_corr > test$adaptiveCokrigingAdvanced_globalVarianceWithChoice_half_corr, ]) / nrow(test)


# Let's see how often each is best, and within a certain distance of the best
test <- augmentedNoChoicePerSeed
test <- test[c(1:5, str_which(colnames(test), "half_corr"))]
test <- test[test$relativeBudget == 10, ]

test$Best <- apply(test[str_which(colnames(test), "_corr")], 1, FUN = max)

colSums(test$Best - test[6:(ncol(test)-1)] <= 0.05)



log(unique(test$costRatio))



# Need to start looking into ISA, seems to be 

chosen <- test[test$instance == 'ShiStyblinskiTang_B10_Cr0.025_S11', ]

testingData <- data
testingData <- testingData[testingData$functionName == "ShiStyblinskiTang", ]
testingData <- testingData[str_which(testingData$method, "_globalVariance_"), ]
testingData <- testingData[testingData$costRatio == 0.025, ]
testingData <- testingData[testingData$budget == 20, ]
testingData <- testingData[testingData$seed == 11, ]



testingData <- testingData[testingData$usedBudget >= 19, ]

testingData <- testingData[testingData$usedBudget >= 29, ]









createBoxPlotData <- function(givenData, givenFeatures, givenMethods){
  split <- strsplit(givenFeatures$instances, ",")
  instances <- gsub("[(]", "", sapply(split, "[[", 1))
  order <- match(givenData$functionName, instances)
  givenData$dimension <- givenFeatures[order, "feature_dimension"]
  givenData$relativeBudget <- givenData$budget / givenData$dimension
  givenData$usedRelativeBudget <- givenData$usedBudget / givenData$dimension
  seeds <- unique(givenData$seed)
  costRatios <- unique(givenData$costRatio)
  relativeBudgets <- unique(givenData$relativeBudget)
  givenData$proportionUsedBudget <- givenData$usedBudget / givenData$budget
  
  # Going to remove entries where went past available budget
  givenData <- givenData[givenData$proportionUsedBudget <= 1, ]
  
  index <- 1
  for(instance in unique(givenData$functionName)){
    print(instance)
    temp_1 <- givenData[givenData$functionName == instance, ]
    for(relativeBudget in relativeBudgets){
      temp_2 <- temp_1[temp_1$relativeBudget == relativeBudget, ]
      for(costRatio in costRatios){
        temp_3 <- temp_2[temp_2$costRatio == costRatio, ]
        for(seed in seeds){
          temp <- temp_3[temp_3$seed == seed, ]
          
          instanceName <- paste0(instance, "_B", relativeBudget, "_Cr", costRatio, "_S", seed)
          # print(instanceName)
          
          if(nrow(temp) == 0){
            print("Didn't find data!!")
            print(instanceName)
            next
          }
          # Now want to get the final corr, error and time for each seed and for
          # each method
          vals <- c()
          names <- c()
          
          tempData <- as.data.frame(matrix(nrow = length(givenMethods), ncol = 8))
          colnames(tempData) <- c("instance", "dimension", "relativeBudget", "costRatio", "seed", "method", "correlation", "error")
          tempData$instance <- instanceName
          tempData$dimension <- temp[1, "dimension"]
          tempData$relativeBudget <- relativeBudget
          tempData$costRatio <- costRatio
          tempData$seed <- seed
          
          localIndex <- 0
          
          # Now populate every other entry
          for(method in givenMethods){
            localIndex <- localIndex + 1
            temp2 <- temp[temp$method == method, ]
            
            if(nrow(temp2) == 0){
              print("Didn't find data!!")
              print(method)
              print(seed)
              next
            }
            
            # Add a check here to see if all the budget was used up!
            if(temp2[nrow(temp2), "budget"] - temp2[nrow(temp2), "usedBudget"] > 1){
              print("Didn't use full budget!!")
              print(method)
              print(seed)
              next
            }
            
            tempData[localIndex, c("method", "correlation", "error")] <- c(method, temp2[nrow(temp2), c("modelCorrelation", "modelError")])
          }
          tempData$bestCorrelation <- max(tempData$correlation)
          tempData$bestError <- min(tempData$error)
          tempData$correlationDiff <- tempData$bestCorrelation - tempData$correlation
          tempData$errorDiff <-  tempData$error - tempData$bestError
          
          
          
          if(index == 1){
            processedData <- tempData
          }else{
            processedData <- rbind(processedData, tempData)
          }
          index <- index + 1
        }
      }
    }
  }
  return(processedData)
}


# To do

# Look into grants
# Look into study away
# Discuss seminar
# Maybe set next run of experiments
# Look at crashing of SOLAR and see what happened
# Finish appendix
# Apply for tutoring
# Renew working with children check



chosenMethods <- c("kriging_globalVariance_half",
                   "cokriging_globalVariance_half",
                   "cokriging_globalVarianceWithChoice_half",
                   "adaptiveCokriging_globalVariance_half",
                   "adaptiveCokriging_globalVarianceWithChoice_half",
                   "adaptiveCokrigingAdvanced_globalVariance_half",
                   "adaptiveCokrigingAdvanced_globalVarianceWithChoice_half")

workingData <- data[str_which(data$functionName, "SOLAR", negate = TRUE), ]

boxPlotData <- createBoxPlotData(workingData, features, chosenMethods)

# Let's try a box plot
# Will want to do y is corr, x is cost ratio, group by technique
plotting <- boxPlotData
useMethods <- c("cokriging_globalVariance_half",
                # "adaptiveCokriging_globalVariance_half",
                "adaptiveCokrigingAdvanced_globalVariance_half",
                "cokriging_globalVarianceWithChoice_half",
                # "adaptiveCokriging_globalVarianceWithChoice_half",
                "adaptiveCokrigingAdvanced_globalVarianceWithChoice_half",
                "kriging_globalVariance_half")
useMethodsNames <- c("Co-Kriging",
                     "Adaptive Co-Kriging",
                     # "Adaptive Co-Kriging (2.0)",
                     "Co-Kriging with choice",
                     "Adaptive Co-Kriging with choice",
                     # "Adaptive Co-Kriging with choice (2.0)",
                     "Kriging")

plotting <- plotting[plotting$method %in% useMethods, ]
for(i in 1:length(useMethods)){
  plotting[plotting$method == useMethods[[i]], "method"] <- useMethodsNames[[i]]
}

plotting$costRatio <- as.factor(plotting$costRatio)
plotting$relativeBudget <- as.factor(plotting$relativeBudget)
plotting$method <- factor(plotting$method, levels = useMethodsNames)


range <- c(0, 0.25)
plot5 <- ggplot(plotting[plotting$relativeBudget == 5, ], aes(x=costRatio, y=errorDiff, fill=method)) +
            geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
            coord_cartesian(ylim = range) + 
            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
            ggtitle("B = 5d") + 
            theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))
plot10 <- ggplot(plotting[plotting$relativeBudget == 10, ], aes(x=costRatio, y=errorDiff, fill=method)) +
            geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
            coord_cartesian(ylim = range) + 
            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
            ggtitle("B = 10d") + 
            theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))
plot15 <- ggplot(plotting[plotting$relativeBudget == 15, ], aes(x=costRatio, y=errorDiff, fill=method)) +
            geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
            coord_cartesian(ylim = range) + 
            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
            ggtitle("B = 15d") + 
            theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))
plot20 <- ggplot(plotting[plotting$relativeBudget == 20, ], aes(x=costRatio, y=errorDiff, fill=method)) +
            geom_boxplot(outlier.shape = NA)+
            coord_cartesian(ylim = range) + 
            theme(legend.title.align=0.5,
                  legend.position = c(.8, .7),
                  legend.margin = margin(6, 6, 6, 6),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank()) +
            scale_fill_discrete(name = "Method")  + 
            ggtitle("B = 20d") + 
            theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))

plot <- ggarrange(plot5, plot10, plot15, plot20, ncol = 1, nrow = 4, common.legend = TRUE, legend="right")
plot <- annotate_figure(plot, left = textGrob("Difference in RRMSE from best", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                bottom = textGrob("Cost ratio", gp = gpar(cex = 1.3), hjust = 1.5))

factor <- 1.5
ggsave("errorByBudget.pdf", plot, device = "pdf", width = factor*210, height = factor*297, units = "mm", dpi  = 1200)

plot5 <- ggplot(plotting[plotting$costRatio == 0.01, ], aes(x=relativeBudget, y=errorDiff, fill=method)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
  coord_cartesian(ylim = range) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ggtitle("Cr = 0.01") + 
  theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))
plot10 <- ggplot(plotting[plotting$costRatio == 0.025, ], aes(x=relativeBudget, y=errorDiff, fill=method)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
  coord_cartesian(ylim = range) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ggtitle("Cr = 0.025") + 
  theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))
plot15 <- ggplot(plotting[plotting$costRatio == 0.1, ], aes(x=relativeBudget, y=errorDiff, fill=method)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
  coord_cartesian(ylim = range) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ggtitle("Cr = 0.1") + 
  theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))
plot20 <- ggplot(plotting[plotting$costRatio == 0.5, ], aes(x=relativeBudget, y=errorDiff, fill=method)) +
  geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim = range) + 
  theme(legend.title.align=0.5,
        legend.position = c(.8, .7),
        legend.margin = margin(6, 6, 6, 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_discrete(name = "Method")  + 
  ggtitle("Cr = 0.5") + 
  theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))

plot <- ggarrange(plot5, plot10, plot15, plot20, ncol = 1, nrow = 4, common.legend = TRUE, legend="right")
plot <- annotate_figure(plot, left = textGrob("Difference in RRMSE from best", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                        bottom = textGrob("Relative budget", gp = gpar(cex = 1.3), hjust = 1.5))

factor <- 1.5
ggsave("errorByCostRatio.pdf", plot, device = "pdf", width = factor*210, height = factor*297, units = "mm", dpi  = 1200)




range <- c(0, 0.9)
plot5 <- ggplot(plotting[plotting$relativeBudget == 5, ], aes(x=costRatio, y=correlationDiff, fill=method)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
  coord_cartesian(ylim = range) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ggtitle("B = 5d") + 
  theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))
plot10 <- ggplot(plotting[plotting$relativeBudget == 10, ], aes(x=costRatio, y=correlationDiff, fill=method)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
  coord_cartesian(ylim = range) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ggtitle("B = 10d") + 
  theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))
plot15 <- ggplot(plotting[plotting$relativeBudget == 15, ], aes(x=costRatio, y=correlationDiff, fill=method)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
  coord_cartesian(ylim = range) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ggtitle("B = 15d") + 
  theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))
plot20 <- ggplot(plotting[plotting$relativeBudget == 20, ], aes(x=costRatio, y=correlationDiff, fill=method)) +
  geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim = range) + 
  theme(legend.title.align=0.5,
        legend.position = c(.8, .7),
        legend.margin = margin(6, 6, 6, 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_discrete(name = "Method")  + 
  ggtitle("B = 20d") + 
  theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))

plot <- ggarrange(plot5, plot10, plot15, plot20, ncol = 1, nrow = 4, common.legend = TRUE, legend="right")
plot <- annotate_figure(plot, left = textGrob("Difference in CC from best", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                        bottom = textGrob("Cost ratio", gp = gpar(cex = 1.3), hjust = 1.5))

factor <- 1.5
ggsave("correlationByBudget.pdf", plot, device = "pdf", width = factor*210, height = factor*297, units = "mm", dpi  = 1200)

plot5 <- ggplot(plotting[plotting$costRatio == 0.01, ], aes(x=relativeBudget, y=correlationDiff, fill=method)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
  coord_cartesian(ylim = range) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ggtitle("Cr = 0.01") + 
  theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))
plot10 <- ggplot(plotting[plotting$costRatio == 0.025, ], aes(x=relativeBudget, y=correlationDiff, fill=method)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
  coord_cartesian(ylim = range) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ggtitle("Cr = 0.025") + 
  theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))
plot15 <- ggplot(plotting[plotting$costRatio == 0.1, ], aes(x=relativeBudget, y=correlationDiff, fill=method)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE)+
  coord_cartesian(ylim = range) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ggtitle("Cr = 0.1") + 
  theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))
plot20 <- ggplot(plotting[plotting$costRatio == 0.5, ], aes(x=relativeBudget, y=correlationDiff, fill=method)) +
  geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim = range) + 
  theme(legend.title.align=0.5,
        legend.position = c(.8, .7),
        legend.margin = margin(6, 6, 6, 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_discrete(name = "Method")  + 
  ggtitle("Cr = 0.5") + 
  theme(plot.title=element_text(hjust=0.5, vjust=1, margin=margin(t=40,b=-30)))

plot <- ggarrange(plot5, plot10, plot15, plot20, ncol = 1, nrow = 4, common.legend = TRUE, legend="right")
plot <- annotate_figure(plot, left = textGrob("Difference in CC from best", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                        bottom = textGrob("Relative budget", gp = gpar(cex = 1.3), hjust = 1.5))

factor <- 1.5
ggsave("correlationByCostRatio.pdf", plot, device = "pdf", width = factor*210, height = factor*297, units = "mm", dpi  = 1200)






verifier <- function(givenData, givenFeatures, givenMethods){
  split <- strsplit(givenFeatures$instances, ",")
  instances <- gsub("[(]", "", sapply(split, "[[", 1))
  order <- match(givenData$functionName, instances)
  givenData$dimension <- givenFeatures[order, "feature_dimension"]
  givenData$relativeBudget <- givenData$budget / givenData$dimension
  givenData$usedRelativeBudget <- givenData$usedBudget / givenData$dimension
  seeds <- unique(givenData$seed)
  costRatios <- unique(givenData$costRatio)
  relativeBudgets <- unique(givenData$relativeBudget)
  givenData$proportionUsedBudget <- givenData$usedBudget / givenData$budget
  
  # Going to remove entries where went past available budget
  givenData <- givenData[givenData$proportionUsedBudget <= 1, ]
  
  reruns <- as.data.frame(matrix(nrow = 0, ncol = 3))
  colnames(reruns) <- c('problemType', 'problem', 'method')
  
  index <- 1
  for(instance in unique(givenData$functionName)){
    print(instance)
    temp <- givenData[givenData$functionName == instance, ]
    for(relativeBudget in relativeBudgets){
      temp2 <- temp[temp$relativeBudget == relativeBudget, ]
      for(costRatio in costRatios){
        temp3 <- temp2[temp2$costRatio == costRatio, ]
        for(method in givenMethods){
          temp4 <- temp3[temp3$method == method, ]
          for(seed in seeds){
            temp5 <- temp4[temp4$seed == seed, ]
            if(nrow(temp5) == 0){
              print("Didn't find data!!")
              print(instance)
              print(relativeBudget)
              print(costRatio)
              print(seed)
              print(method)
              reruns[index, ] <- c("surrogateModelWithBudget", 
                                   paste0("(", instance, ",", relativeBudget, ",", costRatio, ",", seed, "-", seed, ")"),
                                   method)
              index <- index + 1
              next
            }
            if(temp5[nrow(temp5), "budget"] - temp5[nrow(temp5), "usedBudget"] > 1){
              print("Didn't use full budget!!")
              print(instance)
              print(relativeBudget)
              print(costRatio)
              print(seed)
              print(method)
              reruns[index, ] <- c("surrogateModelWithBudget", 
                                   paste0("(", instance, ",", relativeBudget, ",", costRatio, ",", seed, "-", seed, ")"),
                                   method)
              index <- index + 1
            }
          }
        }
      }
    }
  }
  return(reruns)
}

reruns <- verifier(data, features, chosenMethods)







