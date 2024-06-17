######################################################################################

# Copyright 2023, Nicolau Andres-Thio

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

######################################################################################

library(plyr)
library(stringr)
# Functions which take the experimental results and process them so that they can be plotted
# and analysed

# Basic functions which combines all of the experimental results given into a single
# dataframe. No processing is performed, only a storing of the results.
combineArrayResults <- function(runName, arrayStart, arrayEnd, jobsPerArray = 1, printInfo = TRUE){
  options(scipen=999)
  if(jobsPerArray == 1){combinedData <- read.table(paste0("data/clusterResults/", runName, "_arrayJob", arrayStart, ".txt"), header = TRUE, sep = " ", fill = TRUE)}
  else{combinedData <- read.table(paste0("data/clusterResults/", runName, "_arrayJob", (arrayStart - 1) * jobsPerArray + 1, "-", arrayStart * jobsPerArray, ".txt"), header = TRUE, sep = " ", fill = TRUE)}
  if(arrayEnd == arrayStart){return(combinedData)}
  for(i in (arrayStart + 1):arrayEnd){
    if(printInfo){cat(paste0("\rCombining files, done ", i, "/", arrayEnd))}
    if(jobsPerArray == 1){filename <- paste0("data/clusterResults/", runName, "_arrayJob", i, ".txt")}
    else{filename <- paste0("data/clusterResults/", runName, "_arrayJob", (i - 1) * jobsPerArray + 1, "-", i * jobsPerArray, ".txt")}
    
    if(!file.exists(filename)){
      print(paste0("Skipping file ", filename, " as it does not exist!"))
      next
    }
    newData <- read.table(filename, header = TRUE, sep = " ", fill = TRUE)
    combinedData <- rbind.fill(combinedData, newData)
  }
  if(printInfo){cat(" - done.\n")}
  options(scipen=0)
  return(combinedData)
}

# Function which combines algorithm performance and features into a single data frame
augmentData <- function(experimentalData, featuresData, algorithms){
  # Each instance is defined by the function pair, the high and low budgets, and the random seed
  instances <- unique(experimentalData$instances)
  numInstances <- length(instances)
  finalData <- data.frame(instances = instances)
  # Genius move, simply sort featuresData based on instances to match experimental data!!
  for(algorithm in algorithms){
    temp <- experimentalData[experimentalData$method == algorithm, ]
    order <- match(instances, temp$instances)
    temp <- temp[order, ]
    if(nrow(temp) != numInstances){
      print(paste0("Could not find matching algorithm performance for every instance for algorithm ", algorithm, "!"))
      return(NA)
    }
    finalData[, paste0(algorithm, c("_modelError", "_modelCorrelation", "_time"))] <- temp[, c("modelError", "modelCorrelation", "time")]
  }
  order <- match(instances, featuresData$instances)
  relevantSampleFeatureData <- featuresData[order, -which(names(temp) %in% c("instances"))]
  if(sum(is.na(order)) > 0){
    print(instances[!(instances %in% featuresData$instances)])
    print("Could not find matching feature data for every instance!")
    print(instances[is.na(order)])
    return(NA)
  }
  finalData[, colnames(relevantSampleFeatureData)] <- relevantSampleFeatureData
  # Finally add source
  finalData$Source <- "Literature"
  finalData[c(str_which(finalData$instances, "COCO", negate = FALSE),
              str_which(finalData$instances, "Disturbance", negate = FALSE)), "Source"] = "Disturbance-based"
  finalData[str_which(finalData$instances, "SOLAR", negate = FALSE), "Source"] = "SOLAR"
  return(finalData)
}

# Function which "condenses" information from the same instances with different
# seeds into a single row in the final data frame. Algorithm performance
# is aggregated using the mean and the median, and compared to other algorithms
# using the wilcoxon test.
condenseData <- function(data, algorithms){
  split <- strsplit(data$instance, ",")
  instances <- unique(paste0(gsub("[(]", "", sapply(split, "[[", 1)), ",", sapply(split, "[[", 2), ",", sapply(split, "[[", 3)))
  
  featNames <- colnames(data[, str_which(colnames(data), "feature_")])
  colnames <- c("instances", "Source", 
                "cohenCorr",
                paste0(algorithms, "_corrWilcoxon0"), paste0(algorithms, "_corrWilcoxon0.001"), paste0(algorithms, "_corrWilcoxon0.0025"), paste0(algorithms, "_corrWilcoxon0.005"), paste0(algorithms, "_corrWilcoxon0.01"),
                paste0(algorithms, "_corrWilcoxon0Bad"), paste0(algorithms, "_corrWilcoxon0.001Bad"), paste0(algorithms, "_corrWilcoxon0.0025Bad"), paste0(algorithms, "_corrWilcoxon0.005Bad"), paste0(algorithms, "_corrWilcoxon0.01Bad"),
                paste0(algorithms, "_corrMean"), paste0(algorithms, "_corrMedian"),
                paste0(algorithms, "_errorWilcoxon0"), paste0(algorithms, "_errorWilcoxon0.001"), paste0(algorithms, "_errorWilcoxon0.0025"), paste0(algorithms, "_errorWilcoxon0.005"), paste0(algorithms, "_errorWilcoxon0.01"),
                paste0(algorithms, "_errorWilcoxon0Bad"), paste0(algorithms, "_errorWilcoxon0.001Bad"), paste0(algorithms, "_errorWilcoxon0.0025Bad"), paste0(algorithms, "_errorWilcoxon0.005Bad"), paste0(algorithms, "_errorWilcoxon0.01Bad"),
                paste0(algorithms, "_errorMean"), paste0(algorithms, "_errorMedian"), featNames)
  output <- setNames(data.frame(matrix(ncol = length(colnames), nrow = 0)), colnames)
  row <- 0
  for(instance in unique(instances)){
    row <- row + 1
    cat(paste0("\rWorking on row ", row))
    tempData <- data[str_which(data$instance, instance), ]
    if(nrow(tempData) != 40){
      print(paste0("Something weird with instance ", instance, ", have ", nrow(tempData), " rows instead of 40!"))
    }
    output[row, "instances"] <- instance
    output[row, "Source"] <- tempData[1, "Source"]
    output[row, featNames] <- colMeans(tempData[, featNames])
    for(algorithm in algorithms){
      output[row, c(paste0(algorithm, "_corrMean"), paste0(algorithm, "_corrMedian"), paste0(algorithm, "_errorMean"), paste0(algorithm, "_errorMedian"))] <- 
        c(mean(tempData[, paste0(algorithm, "_modelCorrelation")]), median(tempData[, paste0(algorithm, "_modelCorrelation")]), mean(tempData[, paste0(algorithm, "_modelError")]), median(tempData[, paste0(algorithm, "_modelError")]))
    }
    # Formula for Cohen's
    meanKrig <- mean(tempData$kriging_modelCorrelation)
    meanCoKrig <- mean(tempData$cokriging_modelCorrelation)
    varKrig <- var(tempData$kriging_modelCorrelation)
    varCoKrig <- var(tempData$cokriging_modelCorrelation)
    if(varKrig == 0 & varCoKrig == 0){
      output[row, "cohenCorr"] <- 0
    }else{
      output[row, "cohenCorr"] <- (meanKrig - meanCoKrig) / sqrt((varKrig^2 + varCoKrig^2) / 2)
    }
    
    hypSame <- wilcox.test(tempData$kriging_modelCorrelation, tempData$cokriging_modelCorrelation, paired = TRUE, alternative = "less")
    output[row, "kriging_corrWilcoxon0"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelCorrelation, tempData$kriging_modelCorrelation, paired = TRUE, alternative = "less")
    output[row, "cokriging_corrWilcoxon0"] <- hypSame$p.value
    
    hypSame <- wilcox.test(tempData$kriging_modelCorrelation, tempData$cokriging_modelCorrelation, mu = -0.001, paired = TRUE, alternative = "less")
    output[row, "kriging_corrWilcoxon0.001"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelCorrelation, tempData$kriging_modelCorrelation, mu = -0.001, paired = TRUE, alternative = "less")
    output[row, "cokriging_corrWilcoxon0.001"] <- hypSame$p.value
    
    hypSame <- wilcox.test(tempData$kriging_modelCorrelation, tempData$cokriging_modelCorrelation, mu = -0.0025, paired = TRUE, alternative = "less")
    output[row, "kriging_corrWilcoxon0.0025"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelCorrelation, tempData$kriging_modelCorrelation, mu = -0.0025, paired = TRUE, alternative = "less")
    output[row, "cokriging_corrWilcoxon0.0025"] <- hypSame$p.value
    
    hypSame <- wilcox.test(tempData$kriging_modelCorrelation, tempData$cokriging_modelCorrelation, mu = -0.005, paired = TRUE, alternative = "less")
    output[row, "kriging_corrWilcoxon0.005"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelCorrelation, tempData$kriging_modelCorrelation, mu = -0.005, paired = TRUE, alternative = "less")
    output[row, "cokriging_corrWilcoxon0.005"] <- hypSame$p.value
    
    hypSame <- wilcox.test(tempData$kriging_modelCorrelation, tempData$cokriging_modelCorrelation, mu = -0.01, paired = TRUE, alternative = "less")
    output[row, "kriging_corrWilcoxon0.01"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelCorrelation, tempData$kriging_modelCorrelation, mu = -0.01, paired = TRUE, alternative = "less")
    output[row, "cokriging_corrWilcoxon0.01"] <- hypSame$p.value
    
    
    # Here the null hypothesis is that the method is bad
    hypSame <- wilcox.test(tempData$kriging_modelCorrelation, tempData$cokriging_modelCorrelation, paired = TRUE, alternative = "g")
    output[row, "kriging_corrWilcoxon0Bad"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelCorrelation, tempData$kriging_modelCorrelation, paired = TRUE, alternative = "g")
    output[row, "cokriging_corrWilcoxon0Bad"] <- hypSame$p.value
    
    hypSame <- wilcox.test(tempData$kriging_modelCorrelation, tempData$cokriging_modelCorrelation, mu = -0.001, paired = TRUE, alternative = "g")
    output[row, "kriging_corrWilcoxon0.001Bad"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelCorrelation, tempData$kriging_modelCorrelation, mu = -0.001, paired = TRUE, alternative = "g")
    output[row, "cokriging_corrWilcoxon0.001Bad"] <- hypSame$p.value
    
    hypSame <- wilcox.test(tempData$kriging_modelCorrelation, tempData$cokriging_modelCorrelation, mu = -0.0025, paired = TRUE, alternative = "g")
    output[row, "kriging_corrWilcoxon0.0025Bad"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelCorrelation, tempData$kriging_modelCorrelation, mu = -0.0025, paired = TRUE, alternative = "g")
    output[row, "cokriging_corrWilcoxon0.0025Bad"] <- hypSame$p.value
    
    hypSame <- wilcox.test(tempData$kriging_modelCorrelation, tempData$cokriging_modelCorrelation, mu = -0.005, paired = TRUE, alternative = "g")
    output[row, "kriging_corrWilcoxon0.005Bad"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelCorrelation, tempData$kriging_modelCorrelation, mu = -0.005, paired = TRUE, alternative = "g")
    output[row, "cokriging_corrWilcoxon0.005Bad"] <- hypSame$p.value
    
    hypSame <- wilcox.test(tempData$kriging_modelCorrelation, tempData$cokriging_modelCorrelation, mu = -0.01, paired = TRUE, alternative = "g")
    output[row, "kriging_corrWilcoxon0.01Bad"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelCorrelation, tempData$kriging_modelCorrelation, mu = -0.01, paired = TRUE, alternative = "g")
    output[row, "cokriging_corrWilcoxon0.01Bad"] <- hypSame$p.value
    
    
    hypSame <- wilcox.test(tempData$kriging_modelError, tempData$cokriging_modelError, paired = TRUE, alternative = "g")
    output[row, "kriging_errorWilcoxon0"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelError, tempData$kriging_modelError, paired = TRUE, alternative = "g")
    output[row, "cokriging_errorWilcoxon0"] <- hypSame$p.value

    hypSame <- wilcox.test(tempData$kriging_modelError, tempData$cokriging_modelError, mu = 0.001, paired = TRUE, alternative = "g")
    output[row, "kriging_errorWilcoxon0.001"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelError, tempData$kriging_modelError, mu = 0.001, paired = TRUE, alternative = "g")
    output[row, "cokriging_errorWilcoxon0.001"] <- hypSame$p.value

    hypSame <- wilcox.test(tempData$kriging_modelError, tempData$cokriging_modelError, mu = 0.0025, paired = TRUE, alternative = "g")
    output[row, "kriging_errorWilcoxon0.0025"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelError, tempData$kriging_modelError, mu = 0.0025, paired = TRUE, alternative = "g")
    output[row, "cokriging_errorWilcoxon0.0025"] <- hypSame$p.value

    hypSame <- wilcox.test(tempData$kriging_modelError, tempData$cokriging_modelError, mu = 0.005, paired = TRUE, alternative = "g")
    output[row, "kriging_errorWilcoxon0.005"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelError, tempData$kriging_modelError, mu = 0.005, paired = TRUE, alternative = "g")
    output[row, "cokriging_errorWilcoxon0.005"] <- hypSame$p.value

    hypSame <- wilcox.test(tempData$kriging_modelError, tempData$cokriging_modelError, mu = 0.01, paired = TRUE, alternative = "g")
    output[row, "kriging_errorWilcoxon0.01"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelError, tempData$kriging_modelError, mu = 0.01, paired = TRUE, alternative = "g")
    output[row, "cokriging_errorWilcoxon0.01"] <- hypSame$p.value
    
    
    hypSame <- wilcox.test(tempData$kriging_modelError, tempData$cokriging_modelError, paired = TRUE, alternative = "l")
    output[row, "kriging_errorWilcoxon0Bad"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelError, tempData$kriging_modelError, paired = TRUE, alternative = "l")
    output[row, "cokriging_errorWilcoxon0Bad"] <- hypSame$p.value
    
    hypSame <- wilcox.test(tempData$kriging_modelError, tempData$cokriging_modelError, mu = 0.001, paired = TRUE, alternative = "l")
    output[row, "kriging_errorWilcoxon0.001Bad"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelError, tempData$kriging_modelError, mu = 0.001, paired = TRUE, alternative = "l")
    output[row, "cokriging_errorWilcoxon0.001Bad"] <- hypSame$p.value
    
    hypSame <- wilcox.test(tempData$kriging_modelError, tempData$cokriging_modelError, mu = 0.0025, paired = TRUE, alternative = "l")
    output[row, "kriging_errorWilcoxon0.0025Bad"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelError, tempData$kriging_modelError, mu = 0.0025, paired = TRUE, alternative = "l")
    output[row, "cokriging_errorWilcoxon0.0025Bad"] <- hypSame$p.value
    
    hypSame <- wilcox.test(tempData$kriging_modelError, tempData$cokriging_modelError, mu = 0.005, paired = TRUE, alternative = "l")
    output[row, "kriging_errorWilcoxon0.005Bad"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelError, tempData$kriging_modelError, mu = 0.005, paired = TRUE, alternative = "l")
    output[row, "cokriging_errorWilcoxon0.005Bad"] <- hypSame$p.value
    
    hypSame <- wilcox.test(tempData$kriging_modelError, tempData$cokriging_modelError, mu = 0.01, paired = TRUE, alternative = "l")
    output[row, "kriging_errorWilcoxon0.01Bad"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelError, tempData$kriging_modelError, mu = 0.01, paired = TRUE, alternative = "l")
    output[row, "cokriging_errorWilcoxon0.01Bad"] <- hypSame$p.value
  }
  return(output)
}






