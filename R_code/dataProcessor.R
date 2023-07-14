source("R_code/libraries.R")
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
    # if(nrow(newData) != 200 & nrow(newData) != 10){
    #   print(paste0("Weird, for this file have ", nrow(newData), " rows instead of 10 or 200!"))
    #   print(filename)
    # }
    combinedData <- rbind.fill(combinedData, newData)
  }
  if(printInfo){cat(" - done.\n")}
  options(scipen=0)
  return(combinedData)
}

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
  finalData$Source <- "Fixed"
  finalData[str_which(finalData$instances, "Wang", negate = FALSE), "Source"] = "Error-based"
  finalData[str_which(finalData$instances, "COCO", negate = FALSE), "Source"] = "COCO-Disturbance-based"
  finalData[str_which(finalData$instances, "Disturbance", negate = FALSE), "Source"] = "Lit-Disturbance-based"
  finalData[str_which(finalData$instances, "Toal", negate = FALSE), "Source"] = "Parameter-based"
  finalData[str_which(finalData$instances, "SOLAR", negate = FALSE), "Source"] = "SOLAR"
  return(finalData)
}

condenseData <- function(data, algorithms){
  
  split <- strsplit(data$instance, ",")
  instances <- unique(paste0(gsub("[(]", "", sapply(split, "[[", 1)), ",", sapply(split, "[[", 2), ",", sapply(split, "[[", 3)))
  
  featNames <- colnames(data[, str_which(colnames(data), "feature_")])
  colnames <- c("instances", "Source", paste0(algorithms, "_corrWilcoxon"), paste0(algorithms, "_corrMean"), paste0(algorithms, "_corrMedian"),
                            paste0(algorithms, "_errorWilcoxon"), paste0(algorithms, "_errorMean"), paste0(algorithms, "_errorMedian"), featNames)
  output <- setNames(data.frame(matrix(ncol = length(colnames), nrow = 0)), colnames)
  row <- 0
  for(instance in unique(instances)){
    row <- row + 1
    cat(paste0("\rWorking on row ", row))
    tempData <- data[str_which(data$instance, instance), ]
    if(nrow(tempData) != 40){
      print(paste0("Something weird with instance ", instance, ", have ", nrow(tempData), " rows instead of 80!"))
    }
    output[row, "instances"] <- instance
    output[row, "Source"] <- tempData[1, "Source"]
    output[row, featNames] <- colMeans(tempData[, featNames])
    for(algorithm in algorithms){
      output[row, c(paste0(algorithm, "_corrMean"), paste0(algorithm, "_corrMedian"), paste0(algorithm, "_errorMean"), paste0(algorithm, "_errorMedian"))] <- 
        c(mean(tempData[, paste0(algorithm, "_modelCorrelation")]), median(tempData[, paste0(algorithm, "_modelCorrelation")]), mean(tempData[, paste0(algorithm, "_modelError")]), median(tempData[, paste0(algorithm, "_modelError")]))
    }
    
    hypSame <- wilcox.test(tempData$kriging_modelCorrelation, tempData$cokriging_modelCorrelation, paired = TRUE, alternative = "less")
    output[row, "kriging_corrWilcoxon"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelCorrelation, tempData$kriging_modelCorrelation, paired = TRUE, alternative = "less")
    output[row, "cokriging_corrWilcoxon"] <- hypSame$p.value
    
    hypSame <- wilcox.test(tempData$kriging_modelError, tempData$cokriging_modelError, paired = TRUE, alternative = "g")
    output[row, "kriging_errorWilcoxon"] <- hypSame$p.value
    hypSame <- wilcox.test(tempData$cokriging_modelError, tempData$kriging_modelError, paired = TRUE, alternative = "g")
    output[row, "cokriging_errorWilcoxon"] <- hypSame$p.value
  }
  return(output)
}




