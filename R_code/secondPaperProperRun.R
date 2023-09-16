source("R_code/dataProcessor.R")

# Get the information from the runs, and the sample and "perfect" features, then combine into a single data set
rawResultsSynth <- combineArrayResults("experimentalRunSurrogateModelWithFixedSample", 1, 5780, 10)
rawResultsSOLAR <- combineArrayResults("experimentalRunSurrogateModelWithFixedSampleSOLAR", 1, 2160)
rawResults <- rbind(rawResultsSOLAR, rawResultsSynth)
features <- read.table("data/features/sampleAndRealFeaturesCleanStandarised.txt", header = TRUE, sep = " ")

rawResults$instances <- paste0("(", rawResults$functionName, 
                               ",", rawResults$nh,
                               ",", rawResults$nl,
                               ",", rawResults$seed, ")")
# features$feature_sample_seed <- as.numeric(gsub('[)]', '', sapply(strsplit(features$instances, ","), "[[", 4)))

augmentedResults <- augmentData(rawResults, features, c("kriging", "cokriging"))
# Remove feature seed which should not be there
# augmentedResults <- augmentedResults[-str_which(colnames(augmentedResults), "seed")]


write.table(augmentedResults, "data/isaMetadata/surrogateModelWithFixedSampleAggregatedData.txt", quote = FALSE, row.names = FALSE)

# Round performance to 3 decimal places so that Wilcoxon test does not see a difference
# in model correlation of 0.9951 and 0.9952
# augmentedResults$kriging_modelCorrelation <- round(augmentedResults$kriging_modelCorrelation / 0.001) * 0.001
# augmentedResults$cokriging_modelCorrelation <- round(augmentedResults$cokriging_modelCorrelation / 0.001) * 0.001
# augmentedResults$kriging_modelError <- round(augmentedResults$kriging_modelError / 0.001)*0.001
# augmentedResults$cokriging_modelError <- round(augmentedResults$cokriging_modelError / 0.001)*0.001
condensedData <- condenseData(augmentedResults, c("kriging", "cokriging"))
# Sadly cannot use , for instance names, need to change it to underscore
condensedData$instances <- gsub(',', '_', condensedData$instances)
write.table(condensedData, "data/isaMetadata/surrogateModelWithFixedSampleMetadata.txt", quote = FALSE, row.names = FALSE)


# Instance filtering algorithm
purifyInstances <- function(distMatrix, algorithmLabels, eps){
  # First create a dataframe in which to store all the information
  labels <- data.frame(matrix(ncol = 0, nrow = nrow(distMatrix)))
  labels$preclude <- FALSE
  labels$visa <- FALSE
  labels$dissimlar <- TRUE
  if(eps == 0){return(labels)}
  for(i in 1:(nrow(labels)-1)){
    cat(paste0("\rWorking on row ", i, "/", nrow(labels)))
    if(labels[i, "preclude"]){next}
    for(j in (i+1):nrow(labels)){
      if(labels[j, "preclude"]){next}
      # dist <- sqrt(sum((featureData[i,] - featureData[j,])^2))
      dist <- distMatrix[i,j]
      if(dist > eps){next}
      labels[j, "dissimlar"] <- FALSE
      if(sum(algorithmLabels[i, ] == algorithmLabels[j, ]) == ncol(algorithmLabels)){
        labels[j, "preclude"] <- TRUE
        labels[j, "visa"] <- FALSE
      }else{
        labels[j, "visa"] <- TRUE
      }
    }
  }
  cat(paste0("\rWorking on row ", i, "/", nrow(labels), " - done\n"))
  return(labels)
}

calculateCVNND <- function(subsetDistMatrix){
  labels <- data.frame(matrix(ncol = 0, nrow = nrow(subsetDistMatrix)))
  labels$minDist <- 0
  if(nrow(labels) == 1){return(c(0,1))}
  for(i in 1:nrow(labels)){
    cat(paste0("\rWorking on row ", i, "/", nrow(labels)))
    temp <- subsetDistMatrix[-i,i]
    # temp <- featureData[-i, ]
    # temp <- sweep(x = temp, MARGIN = 2, STATS = featureData[i, ], FUN = "-")
    # temp <- temp^2
    # temp <- rowSums(temp)
    # temp <- sqrt(temp)
    labels[i, "minDist"] <- min(temp)
  }
  CV = sd(labels$minDist)/mean(labels$minDist)
  Uniformity = 1 - CV
  cat(paste0("\rWorking on row ", i, "/", nrow(labels), " - done\n"))
  return(c(CV, Uniformity))
}

instancePurificationProcedure <- function(epsilons, instancesData, name){
  results <- data.frame(matrix(nrow = 0, ncol = 0))
  numResults <- 0
  
  # Really should do the instance distances once, save me a lot of time
  Xfeat <- instancesData[, str_which(colnames(instancesData), "feature_")]
  Xfeat <- data.matrix(Xfeat)
  distMatrix <- matrix(nrow = nrow(instancesData), ncol = nrow(instancesData))
  for(i in 1:(nrow(instancesData)-1)){
    cat(paste0("\rWorking on distance matrix row ", i, "/", nrow(instancesData)))
    for(j in i:nrow(instancesData)){
      distMatrix[i,j] <- sqrt(sum((Xfeat[i,] - Xfeat[j,])^2))
      distMatrix[j,i] <- distMatrix[i,j]
    }
  }
  
  for(eps in epsilons){
    print(eps)
    numResults <- numResults + 1
    results[numResults, "epsilon"] <- eps
    # Xfeat <- instancesData[, str_which(colnames(instancesData), "feature_")]
    # Xfeat <- data.matrix(Xfeat)
    Ybin <- instancesData[, c('algo_kriging', 'algo_cokriging')]
    Ybin$algo_kriging <- as.numeric(instancesData$algo_kriging >= 0.5)
    Ybin$algo_cokriging <- as.numeric(instancesData$algo_cokriging >= 0.5)
    Ybin <- data.matrix(Ybin)
    
    # Have the data, first purify
    labels <- purifyInstances(distMatrix, Ybin, eps)
    results[numResults, "VisaRatio"] <- sum(labels$visa) / sum(!labels$preclude)
    results[numResults, "InstancesRatio"] <- sum(!labels$preclude) / nrow(labels)
    # Xfeat <- Xfeat[labels$dissimlar, ]
    # Ybin <- Ybin[labels$dissimlar, ]
    # Now calculate uniformity
    subsetDistMatrix <- distMatrix[labels$dissimlar,labels$dissimlar]
    vals <- calculateCVNND(subsetDistMatrix)
    results[numResults, c("CV", "UniformVal")] <- vals
    print(results)
    if(0 %in% results$epsilon){
      maxCV <- results[results$epsilon == 0, "CV"]
      results$UniformValStandarised <- (results$UniformVal - (1 - maxCV)) / (1 - (1 - maxCV))
    }
    # Save results
    write.table(labels, paste0("data/isaMetadata/instanceFiltering/", name, "eps", eps, ".txt"), quote = FALSE, row.names = FALSE)
    # Also save the data
    write.table(results, paste0("data/isaMetadata/instancePurification", name, ".txt"), quote = FALSE, row.names = FALSE)
    # Plot spread of data to know which epsilon value to use
    if(nrow(results) > 1){
      plot(c(), c(), ylim = c(0, 1), xlim = c(min(epsilons), max(epsilons)))
      lines(results$epsilon, results$UniformVal, col="red", lty = 1)
      lines(results$epsilon, results$VisaRatio, col="green", lty = 1)
      lines(results$epsilon, results$InstancesRatio, col="blue", lty = 1)
    }
  }
  
  # If have done it with eps = 0, calculate the relative uniform value
  # if(0 %in% epsilons){
  #   maxCV <- results[results$epsilon == 0, "CV"]
  #   results$UniformValStandarised <- (results$UniformVal - (1 - maxCV)) / (1 - (1 - maxCV))
  # }
  # write.table(results, paste0("data/isaMetadata/instancePurification", name, ".txt"), quote = FALSE, row.names = FALSE)
  
}


findCorrs <- function(data){
  corrs <- as.data.frame(matrix(ncol = 3, nrow = 0))
  colnames(corrs) <- c("feature", "kriging", "cokriging")
  for(i in 1:length(str_which(colnames(data), "feature_"))){
    feat <- colnames(data[str_which(colnames(data), "feature_")])[[i]]
    corrs[i, ] <- c(feat, abs(cor(data$algo_kriging, data[feat])), abs(cor(data$algo_cokriging, data[feat])))
  }
  corrs$maxCorr <- pmax(corrs$kriging, corrs$cokriging)
  return(corrs)
}



# Get data subset for feature selection using instance filtering
condensedData <- read.table("data/isaMetadata/surrogateModelWithFixedSampleMetadata.txt", header = TRUE, sep = " ")

# Will want filtering for different features, performance matrix and tolerance!
# First order instances
instancesSolar <- condensedData[str_which(condensedData$instances, "SOLAR"), ]
instancesDist <- condensedData[c(str_which(condensedData$instances, "Disturbance"),
                                 str_which(condensedData$instances, "COCO")), ]

instancesLit <- condensedData[-c(str_which(condensedData$instances, "SOLAR"),
                                 str_which(condensedData$instances, "Disturbance"),
                                 str_which(condensedData$instances, "COCO")), ]
set.seed(1)
instancesSolar <- instancesSolar[sample(1:nrow(instancesSolar)), ]
instancesDist <- instancesDist[sample(1:nrow(instancesDist)), ]
instancesLit <- instancesLit[sample(1:nrow(instancesLit)), ]
# Put them together again and save them
instancesOrdered <- rbind(instancesSolar, instancesLit, instancesDist)
write.table(instancesOrdered, "data/isaMetadata/surrogateModelWithFixedSampleMetadataRandomOrder.txt", quote = FALSE, row.names = FALSE)

# Now choose what the performance is, and do the filtering!
tempInstancesAllFeatures <- instancesOrdered[c("instances", "Source", colnames(condensedData[str_which(colnames(condensedData), "feature_")]))]
tempInstancesSampleFeatures <- instancesOrdered[c("instances", "Source", colnames(condensedData[str_which(colnames(condensedData), "feature_sample_")]))]
tempInstancesSampleFeatures$feature_sample_dimension <- tempInstancesAllFeatures$feature_real_dimension

# Now assign performance and do the filtering!
tempInstances <- tempInstancesAllFeatures

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0
instancePurificationProcedure(5:0, tempInstances, "allFeaturesCorrTol0")

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.001
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.001
instancePurificationProcedure(5:0, tempInstances, "allFeaturesCorrTol0.001")

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.0025
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.0025
instancePurificationProcedure(5:0, tempInstances, "allFeaturesCorrTol0.0025")

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.005
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.005
instancePurificationProcedure(5:0, tempInstances, "allFeaturesCorrTol0.005")
                                                                               
tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0
instancePurificationProcedure(5:0, tempInstances, "allFeaturesErrorTol0")

tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0.001
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0.001
instancePurificationProcedure(5:0, tempInstances, "allFeaturesErrorTol0.001")

tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0.0025
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0.0025
instancePurificationProcedure(5:0, tempInstances, "allFeaturesErrorTol0.0025")

tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0.005
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0.005
instancePurificationProcedure(5:0, tempInstances, "allFeaturesErrorTol0.005")

# Repeat with sample features!
tempInstances <- tempInstancesSampleFeatures

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0
instancePurificationProcedure(5:0, tempInstances, "sampleFeaturesCorrTol0")

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.001
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.001
instancePurificationProcedure(5:0, tempInstances, "sampleFeaturesCorrTol0.001")

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.0025
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.0025
instancePurificationProcedure(5:0, tempInstances, "sampleFeaturesCorrTol0.0025")

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.005
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.005
instancePurificationProcedure(5:0, tempInstances, "sampleFeaturesCorrTol0.005")

tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0
instancePurificationProcedure(5:0, tempInstances, "sampleFeaturesErrorTol0")

tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0.001
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0.001
instancePurificationProcedure(5:0, tempInstances, "sampleFeaturesErrorTol0.001")

tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0.0025
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0.0025
instancePurificationProcedure(5:0, tempInstances, "sampleFeaturesErrorTol0.0025")

tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0.005
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0.005
instancePurificationProcedure(5:0, tempInstances, "sampleFeaturesErrorTol0.005")


# Ok, have filtered everything, now is the time to choose features
# Use hash to keep track of all chosen features
# library(hash)
# chosenFeatures <- hash()

name <- "sampleFeaturesCorrTol0.001"
data <- tempInstancesSampleFeatures
# Need to add algorithm performance
data$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.001
data$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.001

results <- read.table(paste0("data/isaMetadata/instancePurification", name, ".txt"), header = TRUE, sep = " ")
# Don't really need to plot, just take the smallest epsilon for which the standarised uniform value is above 0.5
results <- results[results$UniformValStandarised >= 0.5, ]
eps <- min(results$epsilon)
eps
labels <- read.table(paste0("data/isaMetadata/instanceFiltering/", name, "eps", eps, ".txt"), header = TRUE, sep = " ")
filteredData <- data[!labels$preclude, ]
for(feat in colnames(filteredData[str_which(colnames(filteredData), "feature")])){
  filteredData[feat] <- (filteredData[feat] - mean(unlist(filteredData[feat])))/ sd(unlist(filteredData[feat]))
}
corrs <- findCorrs(filteredData)
write.csv(filteredData, "matlab_code/featureSelection/metadata.csv", quote = FALSE, row.names = FALSE)
chosenFeatures[[name]] <- paste0("feature_", read.table("matlab_code/featureSelection/features.txt", header = FALSE, sep = ",")[1,])
# chosenFeatures[[name]] <- c(chosenFeatures[[name]], "feature_sample_highFiBudgetRatio")

test <- corrs
test <- test[str_which(test$feature, "Budget"), ]

# NEED TO SAVE THIS!!
# Want to pad with nas
chosenFeaturesDataframe <- as.data.frame(matrix(nrow = length(keys(chosenFeatures)), ncol = 20))
rownames(chosenFeaturesDataframe) <- keys(chosenFeatures)
for(key in keys(chosenFeatures)){
  line <- rep(NA, 20)
  line[1:length(chosenFeatures[[key]])] <- chosenFeatures[[key]]
  chosenFeaturesDataframe[key, ] <- line
  
}
write.table(chosenFeaturesDataframe, "data/isaMetadata/chosenFeatures.txt", quote = FALSE, row.names = TRUE, col.names = TRUE)



# Ok so now have all of the chosen features, need to re run the filtering with the selected features
library(hash)
chosenFeatures <- hash()
chosenFeaturesDataframe <- read.table("data/isaMetadata/chosenFeatures.txt", header = TRUE, sep = " ")
for(rowName in rownames(chosenFeaturesDataframe)){
  row <- chosenFeaturesDataframe[rowName, ]
  row <- row[!is.na(row)]
  chosenFeatures[[rowName]] <- row
}


tempInstances <- tempInstancesAllFeatures

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["allFeaturesCorrTol0"]])
instancePurificationProcedure(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0), tempInstances[colnames], "instanceFilteringAllFeaturesCorrTol0")

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.001
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.001
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["allFeaturesCorrTol0.001"]])
instancePurificationProcedure(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0), tempInstances[colnames], "instanceFilteringAllFeaturesCorrTol0.001")

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.0025
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.0025
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["allFeaturesCorrTol0.0025"]])
instancePurificationProcedure(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0), tempInstances[colnames], "instanceFilteringAllFeaturesCorrTol0.0025")

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.005
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.005
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["allFeaturesCorrTol0.005"]])
instancePurificationProcedure(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0), tempInstances[colnames], "instanceFilteringAllFeaturesCorrTol0.005")

tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["allFeaturesErrorTol0"]])
instancePurificationProcedure(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0), tempInstances[colnames], "instanceFilteringAllFeaturesErrorTol0")

tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0.001
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0.001
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["allFeaturesErrorTol0.001"]])
instancePurificationProcedure(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0), tempInstances[colnames], "instanceFilteringAllFeaturesErrorTol0.001")

tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0.0025
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0.0025
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["allFeaturesErrorTol0.0025"]])
instancePurificationProcedure(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0), tempInstances[colnames], "instanceFilteringAllFeaturesErrorTol0.0025")

tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0.005
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0.005
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["allFeaturesErrorTol0.005"]])
instancePurificationProcedure(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0), tempInstances[colnames], "instanceFilteringAllFeaturesErrorTol0.005")

# Repeat with sample features!
tempInstances <- tempInstancesSampleFeatures

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["sampleFeaturesCorrTol0"]])
instancePurificationProcedure(c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), tempInstances[colnames], "instanceFilteringSampleFeaturesCorrTol0")

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.001
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.001
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["sampleFeaturesCorrTol0.001"]])
instancePurificationProcedure(c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), tempInstances[colnames], "instanceFilteringSampleFeaturesCorrTol0.001")

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.0025
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.0025
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["sampleFeaturesCorrTol0.0025"]])
instancePurificationProcedure(c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), tempInstances[colnames], "instanceFilteringSampleFeaturesCorrTol0.0025")

tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.005
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.005
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["sampleFeaturesCorrTol0.005"]])
instancePurificationProcedure(c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), tempInstances[colnames], "instanceFilteringSampleFeaturesCorrTol0.005")

tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["sampleFeaturesErrorTol0"]])
instancePurificationProcedure(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0), tempInstances[colnames], "instanceFilteringSampleFeaturesErrorTol0")

tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0.001
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0.001
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["sampleFeaturesErrorTol0.001"]])
instancePurificationProcedure(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0), tempInstances[colnames], "instanceFilteringSampleFeaturesErrorTol0.001")

tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0.0025
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0.0025
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["sampleFeaturesErrorTol0.0025"]])
instancePurificationProcedure(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0), tempInstances[colnames], "instanceFilteringSampleFeaturesErrorTol0.0025")

tempInstances$algo_kriging <- instancesOrdered$kriging_errorWilcoxon0.005
tempInstances$algo_cokriging <- instancesOrdered$cokriging_errorWilcoxon0.005
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["sampleFeaturesErrorTol0.005"]])
instancePurificationProcedure(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0), tempInstances[colnames], "instanceFilteringSampleFeaturesErrorTol0.005")





# This deals with all the features
isaName <- "sampleFeaturesCorrTol0.001"
name <- "instanceFilteringSampleFeaturesCorrTol0.001"
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[[isaName]])
data <- tempInstancesAllFeatures
# Need to add algorithm performance
data$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.001
data$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.001
results <- read.table(paste0("data/isaMetadata/instancePurification", name, ".txt"), header = TRUE, sep = " ")
# Don't really need to plot, just take the smallest epsilon for which the standarised uniform value is above 0.5
results <- results[results$UniformValStandarised >= 0.5, ]
eps <- min(results$epsilon)
eps
labels <- read.table(paste0("data/isaMetadata/instanceFiltering/", name, "eps", eps, ".txt"), header = TRUE, sep = " ")
filteredData <- data[!labels$preclude, colnames]
filteredData <- data[!labels$preclude, ]

corrs <- findCorrs(filteredData)
dir.create(paste0("matlab_code/ISA/", isaName, "/"))
write.csv(filteredData, paste0("matlab_code/ISA/", isaName, "/metadata.csv"), quote = FALSE, row.names = FALSE)







# Now plot some information on feature compared to performance
# Ah, here we have the problem where I actually want real feature values
# Get the "real" feature value
preProcessedFeatures <- read.table("data/features/sampleAndRealFeaturesClean.txt", header = TRUE, sep = " ")
preProcessedFeatures$instanceName <- paste0(gsub('[(]', '', sapply(strsplit(preProcessedFeatures$instances, ","), "[[", 1)), "_",
                                            sapply(strsplit(preProcessedFeatures$instances, ","), "[[", 2), "_",
                                            sapply(strsplit(preProcessedFeatures$instances, ","), "[[", 3))
                                            
projection <- read.table(paste0("matlab_code/ISA/sampleFeaturesCorrTol0.001/projection_matrix.csv"), header = TRUE, sep = ",")

featuresOfInterest <- paste0("feature_", colnames(projection[-1]))
# Add highFiBudgetRatio
featuresOfInterest <- c(featuresOfInterest, "feature_sample_highFiBudgetRatio", "feature_sample_CC")

dataOfInterest <- filteredData[c("instances", "Source", "algo_kriging", "algo_cokriging", featuresOfInterest)]

for(feat in featuresOfInterest){
  vals <- c()
  i <- 0
  print(feat)
  for(instance in dataOfInterest$instances){
    i <- i + 1
    temp <- preProcessedFeatures[preProcessedFeatures$instanceName == instance, feat]
    vals <- c(vals, mean(temp))
  }
  dataOfInterest[paste0("original_", feat)] <- vals
}

dataOfInterest[dataOfInterest$original_feature_sample_mid_ela_meta_lin_simple_adj_r2 < 0, 'original_feature_sample_mid_ela_meta_lin_simple_adj_r2'] <- 0
dataOfInterest[dataOfInterest$original_feature_sample_mid_ela_meta_lin_w_interact_adj_r2 < 0, 'original_feature_sample_mid_ela_meta_lin_w_interact_adj_r2'] <- 0


feats <- c("original_feature_sample_highFiBudgetRatio", 
           "original_feature_sample_budgetRatio",
           "original_feature_sample_LCCrel_0_4",
           "original_feature_sample_LCCrel_0_95",
           'original_feature_sample_mid_ela_meta_lin_w_interact_adj_r2',
           'original_feature_sample_mid_ela_meta_lin_simple_adj_r2',
           "original_feature_sample_RRMSE",
           "original_feature_sample_CC")

plottingData <- as.data.frame(matrix(nrow = 0, ncol = 5))           
colnames(plottingData) <- c("proportion", "number", "xVals", "Method", "feature")
for(featName in paste0("original_", featuresOfInterest)){
  print(featName)
  intervals <- seq(min(unlist(dataOfInterest[featName])),max(unlist(dataOfInterest[featName])), 
                   (max(unlist(dataOfInterest[featName])) - min(unlist(dataOfInterest[featName]))) / 10)
  
  proportionKrigGood <- c()
  proportionCoKrigGood <- c()
  proportionKrigBad <- c()
  proportionCoKrigBad <- c()
  numsKrigGood <- c()
  numsCoKrigGood <- c()
  xVals <- c()
  for(i in 1:(length(intervals) - 1)){
    xStart <- intervals[[i]]
    xEnd <- intervals[[i+1]]
    # xVals <- c(xVals, xStart + (xEnd - xStart)/2)
    xVals <- c(xVals, xStart, xEnd)
    
    
    numKrigGood <- nrow(dataOfInterest[dataOfInterest$algo_kriging >= 0.5 &
                                         dataOfInterest[featName] >= (xStart - 0.001) &
                                         dataOfInterest[featName] <= (xEnd + 0.001), ])
    numCoKrigGood <- nrow(filteredData[dataOfInterest$algo_cokriging >= 0.5 &
                                         dataOfInterest[featName] >= (xStart - 0.001) &
                                         dataOfInterest[featName] <= (xEnd + 0.001), ])
    
    numKrigBad <- nrow(dataOfInterest[dataOfInterest$algo_kriging < 0.5 &
                                        dataOfInterest[featName] >= (xStart - 0.001) &
                                        dataOfInterest[featName] <= (xEnd + 0.001), ])
    numCoKrigBad <- nrow(filteredData[dataOfInterest$algo_cokriging < 0.5 &
                                        dataOfInterest[featName] >= (xStart - 0.001) &
                                        dataOfInterest[featName] <= (xEnd + 0.001), ])
    
    # numsKrigGood <- c(numsKrigGood, numKrigGood)
    # numsCoKrigGood <- c(numsCoKrigGood, numCoKrigGood)
    # 
    # proportionKrigGood <- c(proportionKrigGood, numKrigGood / (numKrigGood + numKrigBad))
    # proportionCoKrigGood <- c(proportionCoKrigGood, numCoKrigGood / (numCoKrigGood + numCoKrigBad))
    # proportionKrigBad <- c(proportionKrigBad, numKrigBad / (numKrigGood + numKrigBad))
    # proportionCoKrigBad <- c(proportionCoKrigBad, numCoKrigBad / (numCoKrigGood + numCoKrigBad))
    
    numsKrigGood <- c(numsKrigGood, numKrigGood, numKrigGood)
    numsCoKrigGood <- c(numsCoKrigGood, numCoKrigGood, numCoKrigGood)
    
    proportionKrigGood <- c(proportionKrigGood, numKrigGood / (numKrigGood + numKrigBad), numKrigGood / (numKrigGood + numKrigBad))
    proportionCoKrigGood <- c(proportionCoKrigGood, numCoKrigGood / (numCoKrigGood + numCoKrigBad), numCoKrigGood / (numCoKrigGood + numCoKrigBad))
    proportionKrigBad <- c(proportionKrigBad, numKrigBad / (numKrigGood + numKrigBad), numKrigBad / (numKrigGood + numKrigBad))
    proportionCoKrigBad <- c(proportionCoKrigBad, numCoKrigBad / (numCoKrigGood + numCoKrigBad), numCoKrigBad / (numCoKrigGood + numCoKrigBad))
    
  }
  currIndex <- nrow(plottingData)
  range <- (currIndex + 1) : (currIndex + 2*length(proportionKrigGood))
  plottingData[range, "proportion"] <- c(proportionKrigGood, proportionCoKrigGood)
  plottingData[range, "number"] <- c(numsKrigGood, numsCoKrigGood)
  plottingData[range, "xVals"] <- c(xVals, xVals)
  plottingData[range, "Method"] <- c(rep("Kriging", length(proportionKrigGood)), rep("Co-Kriging", length(proportionCoKrigGood)))
  plottingData[range, "feature"] <- featName
}
library(ggplot2)

ggplot(plottingData[plottingData$feature == "original_feature_sample_budgetRatio", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(B^{r}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_proportionGoodBudgetRatio.png", width = 3.5, height = 3.5)

ggplot(plottingData[plottingData$feature == "original_feature_sample_highFiBudgetRatio", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(B[h]^{r}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_proportionGoodHighFiBudgetRatio.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_LCCrel_0_4", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(LCC[0.4]^{0.2^{1/d}}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_proportionGoodLCC4.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_LCCrel_0_95", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(LCC[0.95]^{0.2^{1/d}}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_proportionGoodLCC95.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_LCC_sd", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(LCC[sd]^{0.2}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_proportionGoodLCCsd.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_mid_ela_meta_lin_w_interact_adj_r2", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(R[LI]))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_proportionGoodR2LI.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_mid_ela_meta_lin_simple_adj_r2", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(R[L]))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_proportionGoodR2L.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_CC", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(CC))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_proportionGoodCC.png", width = 4, height = 4)


ggplot(plottingData[plottingData$feature == "original_feature_sample_high_ela_level_mmce_lda_50", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(MMCE[lda]^{0.5}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_proportionGoodMMCElda.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_RRMSE", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(RRMSE))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_proportionGoodRRMSE.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_high_ic_m0", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(M[0]))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_proportionGoodM0.png", width = 4, height = 4)



# ACTUAL NUMBER
ggplot(plottingData[plottingData$feature == "original_feature_sample_budgetRatio", ], aes(x=xVals, y=number, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  xlab(expression(B^{r}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_numberGoodBudgetRatio.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_highFiBudgetRatio", ], aes(x=xVals, y=number, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  xlab(expression(B[h]^{r}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_numberGoodHighFiBudgetRatio.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_LCCrel_0_4", ], aes(x=xVals, y=number, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  xlab(expression(LCC[0.4]^{0.2^{1/d}}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_numberGoodLCC4.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_LCCrel_0_95", ], aes(x=xVals, y=number, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  xlab(expression(LCC[0.95]^{0.2^{1/d}}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_numberGoodLCC95.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_LCC_sd", ], aes(x=xVals, y=number, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  xlab(expression(LCC[sd]^{0.2}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_numberGoodLCCsd.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_mid_ela_meta_lin_w_interact_adj_r2", ], aes(x=xVals, y=number, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  xlab(expression(R[LI]))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_numberGoodR2LI.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_mid_ela_meta_lin_simple_adj_r2", ], aes(x=xVals, y=number, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  xlab(expression(R[L]))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_numberGoodR2L.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_CC", ], aes(x=xVals, y=number, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  xlab(expression(CC))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_numberGoodCC.png", width = 4, height = 4)


ggplot(plottingData[plottingData$feature == "original_feature_sample_high_ela_level_mmce_lda_50", ], aes(x=xVals, y=number, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  xlab(expression(MMCE[lda]^{0.5}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_numberGoodMMCElda.png", width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_RRMSE", ], aes(x=xVals, y=number, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  xlab(expression(RRMSE))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_numberGoodRRMSE.png", width = 4, height = 4)


ggplot(plottingData[plottingData$feature == "original_feature_sample_high_ic_m0", ], aes(x=xVals, y=number, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  xlab(expression(M[0]))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave("matlab_code/ISA/sampleFeaturesCorrTol0.001/custom_numberGoodM0.png", width = 4, height = 4)








# Get some statistics on how accurate the guidelines are
test <- dataOfInterest




vals <- seq(0,1,0.05)
accuracies <- c()
for(val in vals){
  chooseKrig <- test[test$original_feature_sample_CC < val, ]
  chooseCoKrig <- test[!(test$original_feature_sample_CC < val), ]
  accuracies <- c(accuracies, (sum(chooseKrig$algo_kriging >= 0.5) + sum(chooseCoKrig$algo_cokriging >= 0.5)) / nrow(test))
}

plot(vals, accuracies)


test <- dataOfInterest

# Rule based on high fi budget ratio
chooseKrig <- test[test$original_feature_sample_highFiBudgetRatio >= 18, ]
undecided <- test[test$original_feature_sample_highFiBudgetRatio < 18, ]

# Rule based on budget ratio
chooseKrig <- rbind(chooseKrig, undecided[undecided$original_feature_sample_budgetRatio == 1, ])
undecided <- undecided[undecided$original_feature_sample_budgetRatio < 1, ]

# Rule based on LCC 0.4
chooseKrig <- rbind(chooseKrig, undecided[undecided$original_feature_sample_LCCrel_0_4 <= 0.7, ])
undecided <- undecided[undecided$original_feature_sample_LCCrel_0_4 > 0.7, ]

# Rule based on LCC 0.95
chooseCoKrig <- undecided[undecided$original_feature_sample_LCCrel_0_95 >= 0.5 | 
                            undecided$original_feature_sample_mid_ela_meta_lin_simple_adj_r2 >= 0.4,]

undecided <- undecided[undecided$original_feature_sample_LCCrel_0_95 < 0.5 & 
                                          undecided$original_feature_sample_mid_ela_meta_lin_simple_adj_r2 < 0.4,]


# Sanity check that these are all distinct
intersect(chooseKrig$instances, undecided$instances)
intersect(chooseKrig$instances, chooseCoKrig$instances)
intersect(chooseCoKrig$instances, undecided$instances)

nrow(chooseKrig)
nrow(chooseCoKrig)
nrow(undecided)


sum(chooseKrig$algo_kriging >= 0.5) / nrow(chooseKrig)
sum(chooseCoKrig$algo_cokriging >= 0.5) / nrow(chooseCoKrig)


# Now check what the remainder of the set looks like
sum(undecided$algo_kriging >= 0.5)
sum(undecided$algo_cokriging >= 0.5)

# Might be nice to get a plot of these instances to know where they lie
write.csv(undecided, paste0("matlab_code/ISA/sampleFeaturesCorrTol0.001/undecidedInstances.csv"), quote = FALSE, row.names = FALSE)

# See about separating what is left using budget ratio or high fi budget ratio
vals <- seq(0,1,0.05)
accuracies <- c()
for(val in vals){
  subsetChooseKrig <- undecided[undecided$original_feature_sample_budgetRatio >= val, ]
  subsetChooseCoKrig <- undecided[!(undecided$original_feature_sample_budgetRatio >= val), ]
  accuracies <- c(accuracies, (sum(subsetChooseKrig$algo_kriging >= 0.5) + sum(subsetChooseCoKrig$algo_cokriging >= 0.5)) / nrow(undecided))
}
plot(vals, accuracies)

vals <- seq(0,20,1)
accuracies <- c()
for(val in vals){
  subsetChooseKrig <- undecided[undecided$original_feature_sample_highFiBudgetRatio > val, ]
  subsetChooseCoKrig <- undecided[!(undecided$original_feature_sample_highFiBudgetRatio > val), ]
  accuracies <- c(accuracies, (sum(subsetChooseKrig$algo_kriging >= 0.5) + sum(subsetChooseCoKrig$algo_cokriging >= 0.5)) / nrow(undecided))
}
plot(vals, accuracies)


# Seems better to split based on relative number of high fidelity samples, with a split of 5d
chooseKrig <- rbind(chooseKrig, undecided[undecided$original_feature_sample_highFiBudgetRatio >= 5, ])
chooseCoKrig <- rbind(chooseCoKrig, undecided[undecided$original_feature_sample_highFiBudgetRatio < 5, ])

nrow(chooseKrig)
nrow(chooseCoKrig)


(sum(chooseKrig$algo_kriging >= chooseKrig$algo_cokriging) + sum(chooseCoKrig$algo_cokriging >= chooseCoKrig$algo_kriging)) / nrow(test)


# Something weird is going on with the predictor accuracy, I want to see
# how often the chosen algorithm is labelled good






# Get metadata with all features in order to do some extra plotting
isaName <- "sampleFeaturesCorrTol0.001"
name <- "instanceFilteringSampleFeaturesCorrTol0.001"
data <- tempInstancesAllFeatures
# Need to add algorithm performance
data$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.001
data$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.001
results <- read.table(paste0("data/isaMetadata/instancePurification", name, ".txt"), header = TRUE, sep = " ")
# Don't really need to plot, just take the smallest epsilon for which the standarised uniform value is above 0.5
results <- results[results$UniformValStandarised >= 0.5, ]
eps <- min(results$epsilon)
eps
labels <- read.table(paste0("data/isaMetadata/instanceFiltering/", name, "eps", eps, ".txt"), header = TRUE, sep = " ")
filteredData <- data[!labels$preclude, ]
for(feat in colnames(filteredData[str_which(colnames(filteredData), "feature")])){
  filteredData[feat] <- (filteredData[feat] - mean(unlist(filteredData[feat])))/ sd(unlist(filteredData[feat]))
}


write.csv(filteredData, paste0("matlab_code/ISA/", isaName, "/metadataAllFeatures.csv"), quote = FALSE, row.names = FALSE)

processedFeatures <- read.table(paste0("matlab_code/ISA/", isaName, "/feature_process.csv"), header = TRUE, sep = ",")






# Ok, time to plot all of these things!
# First focus on the sources
folderName <- "sampleFeaturesCorrTol0.001"
projection <- read.table(paste0("matlab_code/ISA/", folderName, "/projection_matrix.csv"), header = TRUE, sep = ",")
z1Vals <- projection[projection$Row == "Z_{1}", ]
z1Vals <- z1Vals[, -1]
z2Vals <- projection[projection$Row == "Z_{2}", ]
z2Vals <- z2Vals[, -1]

eps <- 0.75
labels <- read.table(paste0("data/isaMetadata/instanceFiltering/instanceFilteringSampleFeaturesCorrTol0.001eps0.3.txt"), header = TRUE, sep = " ")
vals <- instancesOrdered[!labels$preclude, c("instances", "Source", "algo_kriging", "algo_cokriging", featuresKeptAll)]

vals$z1 <- 0
vals$z2 <- 0
for(name in colnames(projection[-1])){
  realName <- paste0("feature_", name)
  vals$z1 <- vals$z1 + z1Vals[1, name] * vals[, realName]
  vals$z2 <- vals$z2 + z2Vals[1, name] * vals[, realName]
}


# vals$Source <- factor(vals$Source, levels = c('Disturbance based instances', "Parameter-based", "Error-based", 'Fixed'))
vals$Source <- factor(vals$Source, levels = c("Disturbance-based", "Literature", "SOLAR"))
vals <- vals[order(vals$Source), ]
mult <- 1

p <- ggplot(data = vals, aes(x = z1, y = z2, col = Source)) +
  geom_point(size = 1) +
  scale_color_manual(values=c("deepskyblue1", "yellow3", "brown3"), drop = FALSE) +
  xlim(c(-4,4))+
  ylim(c(-6, 8)) +
  theme_bw() +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Sources")

ggsave(paste0("matlab_code/ISA/", folderName, "/plots/sourcesAll.png"), plot = egg::set_panel_size(p=p, width=unit(mult*9, "cm"), height=unit(mult*9, "cm")),
       width = mult*17, height = mult*12, units = "cm")

p <- ggplot(data = vals[vals$Source != "Disturbance-based", ], aes(x = z1, y = z2, col = Source)) +
  geom_point(size = 1) +
  scale_color_manual(values=c("deepskyblue1", "yellow3", "brown3"), drop = FALSE) +
  xlim(c(-4,4))+
  ylim(c(-6,8)) +
  theme_bw() +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Sources")

ggsave(paste0("matlab_code/ISA/", folderName, "/plots/sourcesSOLARandLit.png"), plot = egg::set_panel_size(p=p, width=unit(mult*9, "cm"), height=unit(mult*9, "cm")),
       width = mult*17, height = mult*12, units = "cm")

p <- ggplot(data = vals[vals$Source == "SOLAR", ], aes(x = z1, y = z2, col = Source)) +
  geom_point(size = 1) +
  scale_color_manual(values=c("deepskyblue1", "yellow3", "brown3"), drop = FALSE) +
  xlim(c(-4,4))+
  ylim(c(-6,8)) +
  theme_bw() +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Sources")

ggsave(paste0("matlab_code/ISA/", folderName, "/plots/sourcesSOLAR.png"), plot = egg::set_panel_size(p=p, width=unit(mult*9, "cm"), height=unit(mult*9, "cm")),
       width = mult*17, height = mult*12, units = "cm")


# Will probably want the original... not sure how to go about this.
# Maybe start with processing 

p <- ggplot(data = vals, aes(x = z1, y = z2, col = feature_sample_LCC_0_5)) +
  geom_point(size = 0.5) +
  # scale_color_gradient2(low = "blue1", mid = "orange", high = "red2") + 
  theme_classic() + 
  xlim(c(-3.5,3.5))+
  ylim(c(-6,8)) +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  ggtitle(expression(paste("Sample ", LCC[0.5]^{0.2}))) + 
  labs(color = NULL)
ggsave(paste0("matlab_code/ISA/", folderName, "/plots/featuresSampleLCC05.png"), plot = egg::set_panel_size(p=p, width=unit(mult*9, "cm"), height=unit(mult*9, "cm")),
       width = mult*17, height = mult*12, units = "cm")


p <- ggplot(data = vals, aes(x = z1, y = z2, col = feature_sample_LCCrel_0_9)) +
  geom_point(size = 0.5) +
  theme_classic() + 
  xlim(c(-3.5,3.5))+
  ylim(c(-6,8)) +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  ggtitle(expression(paste("Sample ", LCC[0.5]^{0.2^{1/d}}))) + 
  labs(color = NULL)
ggsave(paste0("matlab_code/ISA/", folderName, "/plots/featuresSampleLCCrel09.png"), plot = egg::set_panel_size(p=p, width=unit(mult*9, "cm"), height=unit(mult*9, "cm")),
       width = mult*17, height = mult*12, units = "cm")

p <- ggplot(data = vals, aes(x = z1, y = z2, col = feature_sample_mid_ela_meta_lin_simple_adj_r2)) +
  geom_point(size = 0.5) +
  theme_classic() + 
  xlim(c(-3.5,3.5))+
  ylim(c(-6,8)) +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  ggtitle(expression(paste("Sample ", f[h] - f[l], " ", bar(R)[I]^2))) + 
  labs(color = NULL)
ggsave(paste0("matlab_code/ISA/", folderName, "/plots/featuresSamplemid_ela_meta_lin_simple_adj_r2.png"), plot = egg::set_panel_size(p=p, width=unit(mult*9, "cm"), height=unit(mult*9, "cm")),
       width = mult*17, height = mult*12, units = "cm")















# END OF USED CODE









# This deals with sample features
name <- "instanceFilteringSampleFeaturesCorrTol0"
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures[["sampleFeaturesCorrTol0"]])
data <- tempInstancesSampleFeatures
# Need to add algorithm performance
data$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0
data$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0
results <- read.table(paste0("data/isaMetadata/instancePurification", name, ".txt"), header = TRUE, sep = " ")
# Don't really need to plot, just take the smallest epsilon for which the standarised uniform value is above 0.5
results <- results[results$UniformValStandarised >= 0.5, ]
eps <- min(results$epsilon)
labels <- read.table(paste0("data/isaMetadata/instanceFiltering/", name, "eps", eps, ".txt"), header = TRUE, sep = " ")
filteredData <- data[!labels$preclude, colnames]
write.csv(filteredData, "matlab_code/ISA/isaSampleFeatures/metadata.csv", quote = FALSE, row.names = FALSE)





for(feat in colnames(filteredData[str_which(colnames(filteredData), "feature")])){
  filteredData[feat] <- (filteredData[feat] - mean(unlist(filteredData[feat])))/ sd(unlist(filteredData[feat]))
}
corrs <- findCorrs(filteredData)
write.csv(filteredData, "matlab_code/featureSelection/metadata.csv", quote = FALSE, row.names = FALSE)
chosenFeatures[[name]] <- paste0("feature_", read.table("matlab_code/featureSelection/features.txt", header = FALSE, sep = ",")[1,])
# chosenFeatures[[name]] <- c(chosenFeatures[[name]], "feature_sample_highFiBudgetRatio")











# Look at actual correlation difference to see what it looks like
# condensedData$algo_kriging <- condensedData$kriging_corrMean - condensedData$cokriging_corrMean
# condensedData$algo_cokriging <- condensedData$cokriging_corrMean - condensedData$kriging_corrMean
# condensedData <- condensedData[c("instances", "Source", "algo_kriging", "algo_cokriging",
#                                  colnames(condensedData[str_which(colnames(condensedData), "feature_")]))]




# # # Will only want instances, sources, single algorithm measure and features
# condensedData <- condensedData[c("instances", "Source", "kriging_corrWilcoxon", "cokriging_corrWilcoxon",
#                                  colnames(condensedData[str_which(colnames(condensedData), "feature_")]))]
# colnames(condensedData)[c(3,4)] <- c("algo_kriging", "algo_cokriging")

# # Will only want instances, sources, single algorithm measure and features
condensedData <- condensedData[c("instances", "Source", "kriging_corrWilcoxon0.005", "cokriging_corrWilcoxon0.005",
                                 colnames(condensedData[str_which(colnames(condensedData), "feature_")]))]
colnames(condensedData)[c(3,4)] <- c("algo_kriging", "algo_cokriging")

# Want to "scramble" instances so to even out what gets filtered out
instancesSolar <- condensedData[str_which(condensedData$instances, "SOLAR"), ]
instancesDist <- condensedData[c(str_which(condensedData$instances, "Disturbance"),
                                 str_which(condensedData$instances, "COCO")), ]

instancesLit <- condensedData[-c(str_which(condensedData$instances, "SOLAR"),
                                 str_which(condensedData$instances, "Disturbance"),
                                 str_which(condensedData$instances, "COCO")), ]
set.seed(1)
instancesSolar <- instancesSolar[sample(1:nrow(instancesSolar)), ]
instancesDist <- instancesDist[sample(1:nrow(instancesDist)), ]
instancesLit <- instancesLit[sample(1:nrow(instancesLit)), ]
# Put them together again and save them
instancesOrdered <- rbind(instancesSolar, instancesLit, instancesDist)

write.table(instancesOrdered, "data/isaMetadata/surrogateModelWithFixedSampleMetadataRandomOrder.txt", quote = FALSE, row.names = FALSE)


# Do a filtering with all of the features
instancePurificationProcedure(10:0, instancesOrdered, "allFeatures")

instancesOrderedSampleFeatures <- instancesOrdered[c(1:4, str_which(colnames(instancesOrdered), "feature_sample"))]
instancePurificationProcedure(10:0, instancesOrderedSampleFeatures, "sampleFeatures")


# Next, for each of the two plot the filtering results, choose the right epsilon,
# and get the features using ISA feature selection
results <- read.table("data/isaMetadata/instancePurificationallFeatures.txt", header = TRUE, sep = " ")
epsRange <- results$epsilon
plot(c(), c(), ylim = c(0, 1), xlim = c(min(epsRange), max(epsRange)))
lines(results$epsilon, results$UniformValStandarised, col="red", lty = 1)
lines(results$epsilon, results$VisaRatio, col="green", lty = 1)
lines(results$epsilon, results$InstancesRatio, col="blue", lty = 1)
lines(c(min(epsRange), max(epsRange)), c(0.5, 0.5))
lines(c(3, 3), c(0, 1))
# Choose epsilon and filter data
eps <- 3
labels <- read.table(paste0("data/isaMetadata/instanceFiltering/allFeatureseps", eps, ".txt"), header = TRUE, sep = " ")
filteredData <- instancesOrdered
filteredData <- instancesOrdered[!labels$preclude, ]
# Scale data to have mean 0 and sd 1 as this helps ISA in general
for(feat in colnames(filteredData[str_which(colnames(filteredData), "feature")])){
  filteredData[feat] <- (filteredData[feat] - mean(unlist(filteredData[feat])))/ sd(unlist(filteredData[feat]))
}
write.csv(filteredData, "matlab_code/featureSelection/metadata.csv", quote = FALSE, row.names = FALSE)


# Currently running k = 10, k = 6 does not give nice results, going with k = 15 and adding highFiBudget feature
featuresKeptAll <- paste0("feature_", read.table("matlab_code/featureSelection/features.txt", header = FALSE, sep = ",")[1,])
featuresKeptAll <- c(featuresKeptAll, "feature_sample_highFiBudgetRatio")



corrs <- findCorrs(filteredData)
# # Remove LCC coeff as it is correlated with RRMSE feature, and LCC mean as it is
# # correlated with CC feature (want these as they are "traditional" features)
# filteredData <- filteredData[str_which(colnames(filteredData), "LCC_coeff", negate = TRUE)]
# filteredData <- filteredData[str_which(colnames(filteredData), "LCCrel_coeff", negate = TRUE)]
# filteredData <- filteredData[str_which(colnames(filteredData), "LCC_mean", negate = TRUE)]
# filteredData <- filteredData[str_which(colnames(filteredData), "LCCrel_mean", negate = TRUE)]
# write.csv(filteredData, "matlab_code/featureSelection/metadata.csv", quote = FALSE, row.names = FALSE)
# # Run the feature choice (note that this should probably be run from MATLAB in order
# # to see how many features matlab recommends)
# system("matlab -nodisplay -r \"run('matlab_code/featureSelection/featureSelection.m'); exit\"")
# featuresKeptAll <- paste0("feature_", read.table("matlab_code/featureSelection/features.txt", header = FALSE, sep = ",")[1,])


# Very nice! Repeat process with sample features
results <- read.table("data/isaMetadata/instancePurificationsampleFeatures.txt", header = TRUE, sep = " ")
epsRange <- results$epsilon
plot(c(), c(), ylim = c(0, 1), xlim = c(min(epsRange), max(epsRange)))
lines(results$epsilon, results$UniformValStandarised, col="red", lty = 1)
lines(results$epsilon, results$VisaRatio, col="green", lty = 1)
lines(results$epsilon, results$InstancesRatio, col="blue", lty = 1)
lines(c(min(epsRange), max(epsRange)), c(0.5, 0.5))
lines(c(2, 2), c(0, 1))
# Choose epsilon and filter data
eps <- 2
labels <- read.table(paste0("data/isaMetadata/instanceFiltering/sampleFeatureseps", eps, ".txt"), header = TRUE, sep = " ")
filteredData <- instancesOrdered[, c(1:4, str_which(colnames(instancesOrdered), "feature_sample"))]
filteredData <- instancesOrdered[!labels$preclude, c(1:4, str_which(colnames(instancesOrdered), "feature_sample"))]
# Scale data to have mean 0 and sd 1 as this helps ISA in general
for(feat in colnames(filteredData[str_which(colnames(filteredData), "feature")])){
  filteredData[feat] <- (filteredData[feat] - mean(unlist(filteredData[feat])))/ sd(unlist(filteredData[feat]))
}
corrs <- findCorrs(filteredData)
# Currently running asking for 15 features
write.csv(filteredData, "matlab_code/featureSelection/metadata.csv", quote = FALSE, row.names = FALSE)
featuresKeptSample <- paste0("feature_", read.table("matlab_code/featureSelection/features.txt", header = FALSE, sep = ",")[1,])


# filteredData <- filteredData[str_which(colnames(filteredData), "LCC_coeff", negate = TRUE)]
# filteredData <- filteredData[str_which(colnames(filteredData), "LCCrel_coeff", negate = TRUE)]
# filteredData <- filteredData[str_which(colnames(filteredData), "LCC_mean", negate = TRUE)]
# filteredData <- filteredData[str_which(colnames(filteredData), "LCCrel_mean", negate = TRUE)]
# write.csv(filteredData, "matlab_code/featureSelection/metadata.csv", quote = FALSE, row.names = FALSE)
# system("matlab -nodisplay -r \"run('matlab_code/featureSelection/featureSelection.m'); exit\"")
# featuresKeptSample <- paste0("feature_", read.table("matlab_code/featureSelection/features.txt", header = FALSE, sep = ",")[1,])





# Rerunning the process leads to different features being chosen, here is a set that produced nice results
# featuresKeptAll <- c("feature_sample_LCC_0_5",
#                      "feature_sample_LCCrel_0_9",
#                      "feature_sample_mid_ela_meta_lin_simple_adj_r2",
#                      "feature_sample_mid_ela_meta_lin_w_interact_adj_r2",
#                      "feature_sample_mid_ela_meta_quad_w_interact_adj_r2",
#                      "feature_sample_high_nbc_nb_fitness_cor",
#                      "feature_sample_RRMSE",
#                      "feature_sample_mid_nbc_dist_ratio_coeff_var",
#                      "feature_sample_highFiBudgetRatio",
#                      "feature_sample_budgetRatio",
#                      "feature_real_CC",
#                      "feature_diff_high_ela_meta_lin_w_interact_adj_r2",
#                      "feature_diff_mid_ela_meta_lin_simple_adj_r2",
#                      "feature_diff_high_ela_distr_kurtosis",
#                      "feature_diff_mid_ela_distr_kurtosis")



# One more round of instance filtering with the now chosen features
instancePurificationProcedure(c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0), instancesOrdered[c("instances", "Source", "algo_kriging", "algo_cokriging", featuresKeptAll)], "chosenFeatures")
instancePurificationProcedure(c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0), instancesOrdered[c("instances", "Source", "algo_kriging", "algo_cokriging", featuresKeptSample)], "chosenSampleFeatures")

# Get instances for ISA
results <- read.table("data/isaMetadata/instancePurificationchosenFeatures.txt", header = TRUE, sep = " ")
epsRange <- results$epsilon
plot(c(), c(), ylim = c(0, 1), xlim = c(min(epsRange), max(epsRange)))
lines(results$epsilon, results$UniformValStandarised, col="red", lty = 1)
lines(results$epsilon, results$VisaRatio, col="green", lty = 1)
lines(results$epsilon, results$InstancesRatio, col="blue", lty = 1)
lines(c(min(epsRange), max(epsRange)), c(0.5, 0.5))
lines(c(0.5, 0.5), c(0, 1))
# Choose epsilon and filter data
eps <- 0.5
labels <- read.table(paste0("data/isaMetadata/instanceFiltering/chosenFeatureseps", eps, ".txt"), header = TRUE, sep = " ")
filteredData <- instancesOrdered[!labels$preclude, c("instances", "Source", "algo_kriging", "algo_cokriging", featuresKeptAll)]

filteredData <- instancesOrdered[, c("instances", "Source", "algo_kriging", "algo_cokriging", featuresKeptAll)]

# One last round of standarisation
for(feat in colnames(filteredData[str_which(colnames(filteredData), "feature")])){
  filteredData[feat] <- (filteredData[feat] - mean(unlist(filteredData[feat])))/ sd(unlist(filteredData[feat]))
}
write.csv(filteredData, "matlab_code/ISA/isaAllFeatures/metadata.csv", quote = FALSE, row.names = FALSE)



# Repeat for sample only features
results <- read.table("data/isaMetadata/instancePurificationchosenSampleFeatures.txt", header = TRUE, sep = " ")
epsRange <- results$epsilon
plot(c(), c(), ylim = c(0, 1), xlim = c(min(epsRange), max(epsRange)))
lines(results$epsilon, results$UniformValStandarised, col="red", lty = 1)
lines(results$epsilon, results$VisaRatio, col="green", lty = 1)
lines(results$epsilon, results$InstancesRatio, col="blue", lty = 1)
lines(c(min(epsRange), max(epsRange)), c(0.5, 0.5))
lines(c(0.3, 0.3), c(0, 1))
# Choose epsilon and filter data
eps <- 0.3
labels <- read.table(paste0("data/isaMetadata/instanceFiltering/chosenSampleFeatureseps", eps, ".txt"), header = TRUE, sep = " ")
filteredData <- instancesOrdered[!labels$preclude, c("instances", "Source", "algo_kriging", "algo_cokriging", featuresKeptSample)]
for(feat in colnames(filteredData[str_which(colnames(filteredData), "feature")])){
  filteredData[feat] <- (filteredData[feat] - mean(unlist(filteredData[feat])))/ sd(unlist(filteredData[feat]))
}
write.csv(filteredData, "matlab_code/ISA/isaSampleFeatures/metadata.csv", quote = FALSE, row.names = FALSE)







# Ignore the plotting, let's look at playing with p-value stuff
dataPreprocessed <- read.table(paste0("data/isaMetadata/surrogateModelWithFixedSampleAggregatedData.txt"), header = TRUE, sep = " ")
featPreprocessed <- 
  
  
  # Round performance to 3 decimal places so that Wilcoxon test does not see a difference
  # in model correlation of 0.9951 and 0.9952
  augmentedResults$kriging_modelCorrelation <- round(augmentedResults$kriging_modelCorrelation / 0.001) * 0.001
augmentedResults$cokriging_modelCorrelation <- round(augmentedResults$cokriging_modelCorrelation / 0.001) * 0.001
augmentedResults$kriging_modelError <- round(augmentedResults$kriging_modelError / 0.001)*0.001
augmentedResults$cokriging_modelError <- round(augmentedResults$cokriging_modelError / 0.001)*0.001
condensedData <- condenseData(augmentedResults, c("kriging", "cokriging"))
# Sadly cannot use , for instance names, need to change it to underscore
condensedData$instances <- gsub(',', '_', condensedData$instances)
write.table(condensedData, "data/isaMetadata/surrogateModelWithFixedSampleMetadata.txt", quote = FALSE, row.names = FALSE)





eps <- 0.75
labels <- read.table(paste0("data/isaMetadata/instanceFiltering/chosenFeatureseps", eps, ".txt"), header = TRUE, sep = " ")
filteredData <- instancesOrdered[!labels$preclude, c("instances", "Source", "algo_kriging", "algo_cokriging", featuresKeptAll)]
criticalData <- instancesOrdered[labels$visa, c("instances", "Source", "algo_kriging", "algo_cokriging", featuresKeptAll)]

predictions <- read.csv("matlab_code/ISA/isaAllFeatures/algorithm_bin.csv")
krigGood <- filteredData[predictions$kriging == 1, ]
krigBad <- filteredData[predictions$kriging == 0, ]



# Here 


featName <- "feature_sample_highFiBudgetRatio"
featName <- "feature_sample_budgetRatio"


intervals <- seq(min(unlist(filteredData[featName])),max(unlist(filteredData[featName])), 
                 (max(unlist(filteredData[featName])) - min(unlist(filteredData[featName]))) / 10)

proportionKrigGood <- c()
proportionCoKrigGood <- c()
proportionKrigBad <- c()
proportionCoKrigBad <- c()
xVals <- c()
for(i in 1:(length(intervals) - 1)){
  xStart <- intervals[[i]]
  xEnd <- intervals[[i+1]]
  xVals <- c(xVals, xStart + (xEnd - xStart)/2)
  # numKrig <- nrow(filteredData[predictions$kriging == 1 & 
  #                           filteredData[featName] >= xStart &
  #                           filteredData[featName] <= xEnd, ])
  # numCoKrig <- nrow(filteredData[predictions$cokriging == 1 & 
  #                           filteredData[featName] >= xStart &
  #                           filteredData[featName] <= xEnd, ])
  # numKrig <- nrow(filteredData[predictions$kriging == 1 & 
  #                                filteredData[featName] >= xStart, ])
  # numCoKrig <- nrow(filteredData[predictions$cokriging == 1 & 
  #                                  filteredData[featName] >= xStart, ])
  
  numKrigGood <- nrow(filteredData[predictions$kriging == 1 &
                                     filteredData[featName] >= xStart &
                                     filteredData[featName] <= xEnd, ])
  numCoKrigGood <- nrow(filteredData[predictions$cokriging == 1 &
                                       filteredData[featName] >= xStart &
                                       filteredData[featName] <= xEnd, ])
  
  numKrigBad <- nrow(filteredData[predictions$kriging == 0 &
                                    filteredData[featName] >= xStart &
                                    filteredData[featName] <= xEnd, ])
  numCoKrigBad <- nrow(filteredData[predictions$cokriging == 0 &
                                      filteredData[featName] >= xStart &
                                      filteredData[featName] <= xEnd, ])
  
  
  proportionKrigGood <- c(proportionKrigGood, numKrigGood / (numKrigGood + numKrigBad))
  proportionCoKrigGood <- c(proportionCoKrigGood, numCoKrigGood / (numCoKrigGood + numCoKrigBad))
  proportionKrigBad <- c(proportionKrigBad, numKrigBad / (numKrigGood + numKrigBad))
  proportionCoKrigBad <- c(proportionCoKrigBad, numCoKrigBad / (numCoKrigGood + numCoKrigBad))
  # proportionKrig <- c(proportionKrig, numKrig)
  # proportionCoKrig <- c(proportionCoKrig, numCoKrig)
  
  
}

# Here should be able to plot
plot(xVals, proportionKrigGood, type = 'l', col = 'blue', ylim = c(0,1))
# lines(xVals, proportionKrigBad, type = 'l', col = 'red')
lines(xVals, proportionCoKrigGood, lty = 2, col = 'blue', ylim = c(0,1))
# lines(xVals, proportionCoKrigBad, lty = 2, col = 'red')



test <- read.table("data/isaMetadata/surrogateModelWithFixedSampleMetadata.txt", header = TRUE, sep = " ")
test$algo_kriging <- test$kriging_errorMedian - test$cokriging_errorMedian
test$algo_cokriging <- test$cokriging_errorMedian - test$kriging_errorMedian
corrs <- findCorrs(test)













wilcox.test(unlist(krigGood[featName]), unlist(krigBad[featName]))

mean <- mean(unlist(krigGood[featName]))
sd <- sd(unlist(krigGood[featName]))
n <- nrow(krigGood)
error <- qnorm(0.975)*sd/sqrt(n)
mean + error
mean - error

quantile(unlist(krigGood[featName]), probs = c(0.05, 0.95))

boxplot(krigGood[featName])

boxplot(krigBad$feature_sample_highFiBudgetRatio)


dev.print(png, file = "myplot.png", width = 1024, height = 768)
jpeg("good.jpg")


krigGood$val <- 1
krigBad$val <- 2

krig <- rbind(krigGood, krigBad)

ggplot(krig, aes(x = val, y = feature_sample_highFiBudgetRatio)) +
  boxplot(x = val, y = feature_sample_highFiBudgetRatio)

boxplot(krig$feature_sample_highFiBudgetRatio)
boxplot(krigBad$feature_sample_highFiBudgetRatio)






# Ok, time to plot all of these things!
# First focus on the sources
folderName <- "isaAllFeatures"
projection <- read.table(paste0("matlab_code/ISA/", folderName, "/projection_matrix.csv"), header = TRUE, sep = ",")
z1Vals <- projection[projection$Row == "Z_{1}", ]
z1Vals <- z1Vals[, -1]
z2Vals <- projection[projection$Row == "Z_{2}", ]
z2Vals <- z2Vals[, -1]

eps <- 0.75
labels <- read.table(paste0("data/isaMetadata/instanceFiltering/chosenFeatureseps", eps, ".txt"), header = TRUE, sep = " ")
vals <- instancesOrdered[!labels$preclude, c("instances", "Source", "algo_kriging", "algo_cokriging", featuresKeptAll)]

vals$z1 <- 0
vals$z2 <- 0
for(name in colnames(projection[-1])){
  realName <- paste0("feature_", name)
  vals$z1 <- vals$z1 + z1Vals[1, name] * vals[, realName]
  vals$z2 <- vals$z2 + z2Vals[1, name] * vals[, realName]
}


# vals$Source <- factor(vals$Source, levels = c('Disturbance based instances', "Parameter-based", "Error-based", 'Fixed'))
vals$Source <- factor(vals$Source, levels = c("Disturbance-based", "Literature", "SOLAR"))
vals <- vals[order(vals$Source), ]
mult <- 1

p <- ggplot(data = vals, aes(x = z1, y = z2, col = Source)) +
  geom_point(size = 1) +
  scale_color_manual(values=c("deepskyblue1", "yellow3", "brown3"), drop = FALSE) +
  xlim(c(-4,4))+
  ylim(c(-6, 8)) +
  theme_bw() +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Sources")

ggsave(paste0("matlab_code/ISA/", folderName, "/plots/sourcesAll.png"), plot = egg::set_panel_size(p=p, width=unit(mult*9, "cm"), height=unit(mult*9, "cm")),
       width = mult*17, height = mult*12, units = "cm")

p <- ggplot(data = vals[vals$Source != "Disturbance-based", ], aes(x = z1, y = z2, col = Source)) +
  geom_point(size = 1) +
  scale_color_manual(values=c("deepskyblue1", "yellow3", "brown3"), drop = FALSE) +
  xlim(c(-4,4))+
  ylim(c(-6,8)) +
  theme_bw() +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Sources")

ggsave(paste0("matlab_code/ISA/", folderName, "/plots/sourcesSOLARandLit.png"), plot = egg::set_panel_size(p=p, width=unit(mult*9, "cm"), height=unit(mult*9, "cm")),
       width = mult*17, height = mult*12, units = "cm")

p <- ggplot(data = vals[vals$Source == "SOLAR", ], aes(x = z1, y = z2, col = Source)) +
  geom_point(size = 1) +
  scale_color_manual(values=c("deepskyblue1", "yellow3", "brown3"), drop = FALSE) +
  xlim(c(-4,4))+
  ylim(c(-6,8)) +
  theme_bw() +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Sources")

ggsave(paste0("matlab_code/ISA/", folderName, "/plots/sourcesSOLAR.png"), plot = egg::set_panel_size(p=p, width=unit(mult*9, "cm"), height=unit(mult*9, "cm")),
       width = mult*17, height = mult*12, units = "cm")


# Will probably want the original... not sure how to go about this.
# Maybe start with processing 

p <- ggplot(data = vals, aes(x = z1, y = z2, col = feature_sample_LCC_0_5)) +
  geom_point(size = 0.5) +
  # scale_color_gradient2(low = "blue1", mid = "orange", high = "red2") + 
  theme_classic() + 
  xlim(c(-3.5,3.5))+
  ylim(c(-6,8)) +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  ggtitle(expression(paste("Sample ", LCC[0.5]^{0.2}))) + 
  labs(color = NULL)
ggsave(paste0("matlab_code/ISA/", folderName, "/plots/featuresSampleLCC05.png"), plot = egg::set_panel_size(p=p, width=unit(mult*9, "cm"), height=unit(mult*9, "cm")),
       width = mult*17, height = mult*12, units = "cm")


p <- ggplot(data = vals, aes(x = z1, y = z2, col = feature_sample_LCCrel_0_9)) +
  geom_point(size = 0.5) +
  theme_classic() + 
  xlim(c(-3.5,3.5))+
  ylim(c(-6,8)) +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  ggtitle(expression(paste("Sample ", LCC[0.5]^{0.2^{1/d}}))) + 
  labs(color = NULL)
ggsave(paste0("matlab_code/ISA/", folderName, "/plots/featuresSampleLCCrel09.png"), plot = egg::set_panel_size(p=p, width=unit(mult*9, "cm"), height=unit(mult*9, "cm")),
       width = mult*17, height = mult*12, units = "cm")

p <- ggplot(data = vals, aes(x = z1, y = z2, col = feature_sample_mid_ela_meta_lin_simple_adj_r2)) +
  geom_point(size = 0.5) +
  theme_classic() + 
  xlim(c(-3.5,3.5))+
  ylim(c(-6,8)) +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  ggtitle(expression(paste("Sample ", f[h] - f[l], " ", bar(R)[I]^2))) + 
  labs(color = NULL)
ggsave(paste0("matlab_code/ISA/", folderName, "/plots/featuresSamplemid_ela_meta_lin_simple_adj_r2.png"), plot = egg::set_panel_size(p=p, width=unit(mult*9, "cm"), height=unit(mult*9, "cm")),
       width = mult*17, height = mult*12, units = "cm")




























































# Ok at this point I am all done with ISA, the interesting thing would be to look
# and the prediction model


library(R.matlab)
Matlab$startServer()
matlab <- Matlab()
isOpen <- open(matlab)

name <- 'SpecifiedFeatures'

chosenData <- read.table(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrCondensedFinalWilcoxon", name, "/metadata.csv"), header = TRUE, sep = ",")
featProcessed <- featProcessed <- read.table(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrCondensedFinalWilcoxon", name, "/feature_process.csv"), header = TRUE, sep = ",")
chosenData[paste0("feature_", colnames(featProcessed[-1]))] <- featProcessed[-1]
coords <- read.table(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrCondensedFinalWilcoxon", name, "/coordinates.csv"), header = TRUE, sep = ",")

# projection <- read.table(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrCondensedFinalWilcoxon", name, "/projection_matrix.csv"), header = TRUE, sep = ",")
# features <- paste0("feature_", colnames(projection[-1]))
# chosenData <- dataNormalised[colnames(dataNormalised) %in% features]

# Load up the model
evaluate(matlab, paste0('model = load("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrCondensedFinalWilcoxon', name, '/model.mat");'))
evaluate(matlab, 'disp(model.pythia.selection0);')

algos <- data.matrix(chosenData[str_which(colnames(chosenData), "algo")])
feats <- data.matrix(chosenData[str_which(colnames(chosenData), "feature")])
setVariable(matlab, X = feats)
# setVariable(matlab, X = data.matrix(chosenData))
evaluate(matlab, "Z = X*model.pilot.A';")
evaluate(matlab, 'disp(model.pilot.A);')
evaluate(matlab, 'disp(model.pythia.mu);')

evaluate(matlab, "Z = (Z-model.pythia.mu)./model.pythia.sigma;")
evaluate(matlab, "[Yhat,aux] = model.pythia.svm{1}.predict(Z);")
krigZvals <- getVariable(matlab, "Z")[[1]]
krigPredictions <- getVariable(matlab, "Yhat")[[1]]
test <- getVariable(matlab, "Pr0hat")[[1]]
evaluate(matlab, "[Yhat,aux] = model.pythia.svm{2}.predict(Z);")
cokrigZvals <- getVariable(matlab, "Z")[[1]]
cokrigPredictions <- getVariable(matlab, "Yhat")[[1]]
# close(matlab)



#
# p <- ggplot(data = vals, aes(x = z1, y = z2, col = Source)) +
#   # geom_point(shape = ".") +
#   geom_point(size = 0.5) +
#   # geom_point() +
#   # scale_x_continuous(limits = xVals) +
#   # scale_y_continuous(limits = yVals) +
#   scale_x_continuous(limits = c(-10, 10)) +
#   scale_y_continuous(limits = c(-6, 6)) +
#   # scale_colour_gradientn(colours = colourGradient) +
#   scale_color_manual(values=c("chartreuse4", "darkgoldenrod1", "darkorchid", "brown3"), drop = FALSE) +
#   theme_bw() +
#   ggtitle("Sources")
#
# mult <- 1
#
# ggsave("3-sourcesAll.png", plot = egg::set_panel_size(p=p, width=unit(mult*10, "cm"), height=unit(mult*10, "cm")),
#        width = mult*17, height = mult*12, units = "cm")
#








# # p <- ggplot(data = vals[vals$Source %in% c('Disturbance based instances', "Parameter-based", "Error-based", 'Fixed')], aes(x = z1, y = z2, col = Source)) +
# p <- ggplot(data = vals[vals$Source %in% c('Fixed'),], aes(x = z1, y = z2, col = Source)) +
#   # geom_point(shape = ".") +
#   geom_point(size = 0.5) +
#   # geom_point() +
#   # scale_x_continuous(limits = xVals) +
#   # scale_y_continuous(limits = yVals) +
#   scale_x_continuous(limits = c(-10, 10)) +
#   scale_y_continuous(limits = c(-6, 6)) +
#   # scale_colour_gradientn(colours = colourGradient) +
#   scale_color_manual(values=c("chartreuse4", "darkgoldenrod1", "darkorchid", "brown3"), drop = FALSE) +
#   theme_bw() +
#   ggtitle("Sources")
# 
# mult <- 1
# 
# ggsave("1-sourcesFixed.png", plot = egg::set_panel_size(p=p, width=unit(mult*10, "cm"), height=unit(mult*10, "cm")),
#        width = mult*17, height = mult*12, units = "cm")

# p <- ggplot(data = vals[vals$Source %in% c('Fixed', "Parameter-based", "Error-based"),], aes(x = z1, y = z2, col = Source)) +
#   # geom_point(shape = ".") +
#   geom_point(size = 0.5) +
#   # geom_point() +
#   # scale_x_continuous(limits = xVals) +
#   # scale_y_continuous(limits = yVals) +
#   scale_x_continuous(limits = c(-10, 10)) +
#   scale_y_continuous(limits = c(-6, 6)) +
#   # scale_colour_gradientn(colours = colourGradient) +
#   scale_color_manual(values=c("chartreuse4", "darkgoldenrod1", "darkorchid", "brown3"), drop = FALSE) +
#   theme_bw() +
#   ggtitle("Sources")
# 
# mult <- 1
# 
# ggsave("2-sourcesLit.png", plot = egg::set_panel_size(p=p, width=unit(mult*10, "cm"), height=unit(mult*10, "cm")),
#        width = mult*17, height = mult*12, units = "cm")
# 
# p <- ggplot(data = vals, aes(x = z1, y = z2, col = Source)) +
#   # geom_point(shape = ".") +
#   geom_point(size = 0.5) +
#   # geom_point() +
#   # scale_x_continuous(limits = xVals) +
#   # scale_y_continuous(limits = yVals) +
#   scale_x_continuous(limits = c(-10, 10)) +
#   scale_y_continuous(limits = c(-6, 6)) +
#   # scale_colour_gradientn(colours = colourGradient) +
#   scale_color_manual(values=c("chartreuse4", "darkgoldenrod1", "darkorchid", "brown3"), drop = FALSE) +
#   theme_bw() +
#   ggtitle("Sources")
# 
# mult <- 1
# 
# ggsave("3-sourcesAll.png", plot = egg::set_panel_size(p=p, width=unit(mult*10, "cm"), height=unit(mult*10, "cm")),
#        width = mult*17, height = mult*12, units = "cm")
# 





sum(condensedData$kriging_corrWilcoxon > 0.5) / nrow(condensedData)
sum(condensedData$cokriging_corrWilcoxon > 0.5) / nrow(condensedData)





# 
# 
# # At this point, look at data from filtering and choose appropriate epsilon for each of the 4 cases
# labels <- read.table("data/combinedData/instanceFiltering/chosenFeatures_condensedErrWilcoxoneps0.3.txt", header = TRUE, sep = " ")
# chosenData <- condensedData[, 1:14]
# projection <- read.table("C:/Users/nandr/Documents/InstanceSpace-master/trial/errorCondensedFeatureSelectionWilcoxon/projection_matrix.csv", header = TRUE, sep = ",")
# featNames <- paste0("feature_", colnames(projection[, -1]))
# chosenData$algo_kriging <- chosenData$kriging_errorWilcoxon
# chosenData$algo_cokriging <- chosenData$cokriging_errorWilcoxon
# chosenData <- chosenData[, str_which(names(chosenData), "kriging_", negate = TRUE)]
# chosenData[, featNames] <- condensedData[, featNames]
# chosenData <- chosenData[!labels$preclude, ]
# # Ready to run feature selection!
# dir.create(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/errorCondensedSpecifiedFeaturesWilcoxon/"))
# write.csv(chosenData, paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/errorCondensedSpecifiedFeaturesWilcoxon/metadata.csv"), quote = FALSE, row.names = FALSE)
# 
# 
# # At this point, look at data from filtering and choose appropriate epsilon for each of the 4 cases
# labels <- read.table("data/combinedData/instanceFiltering/chosenFeatures_condensedCorrMedianeps0.4.txt", header = TRUE, sep = " ")
# chosenData <- condensedData[, 1:14]
# projection <- read.table("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrCondensedFeatureSelectionMedian/projection_matrix.csv", header = TRUE, sep = ",")
# featNames <- paste0("feature_", colnames(projection[, -1]))
# chosenData$algo_kriging <- chosenData$kriging_corrMedian + 1
# chosenData$algo_cokriging <- chosenData$cokriging_corrMedian + 1
# chosenData[, featNames] <- condensedData[, featNames]
# chosenData <- chosenData[!labels$preclude, ]
# # Ready to run feature selection!
# dir.create(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrCondensedSpecifiedFeaturesMedian/"))
# write.csv(chosenData, paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrCondensedSpecifiedFeaturesMedian/metadata.csv"), quote = FALSE, row.names = FALSE)
# 
# 
# # At this point, look at data from filtering and choose appropriate epsilon for each of the 4 cases
# labels <- read.table("data/combinedData/instanceFiltering/chosenFeatures_condensedErrMedianeps0.3.txt", header = TRUE, sep = " ")
# chosenData <- condensedData[, 1:14]
# projection <- read.table("C:/Users/nandr/Documents/InstanceSpace-master/trial/errorCondensedFeatureSelectionMedian/projection_matrix.csv", header = TRUE, sep = ",")
# featNames <- paste0("feature_", colnames(projection[, -1]))
# chosenData$algo_kriging <- chosenData$kriging_errorMedian + 1
# chosenData$algo_cokriging <- chosenData$cokriging_errorMedian + 1
# chosenData <- chosenData[, str_which(names(chosenData), "kriging_", negate = TRUE)]
# chosenData[, featNames] <- condensedData[, featNames]
# chosenData <- chosenData[!labels$preclude, ]
# # Ready to run feature selection!
# dir.create(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/errorCondensedSpecifiedFeaturesMedian/"))
# write.csv(chosenData, paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/errorCondensedSpecifiedFeaturesMedian/metadata.csv"), quote = FALSE, row.names = FALSE)
# 













# Time for some plotting!
library(ggpubr)

plottingData <- condensedData
plottingData$instance <- sapply(strsplit(plottingData$instances, "_"), "[[", 1)
plottingData$feature_sample_highFiBudget <- as.numeric(gsub('H', '', sapply(strsplit(plottingData$instances, "_"), "[[", 2)))
plottingData$feature_sample_lowFiBudget <- as.numeric(gsub('L', '', sapply(strsplit(plottingData$instances, "_"), "[[", 3)))

for(i in 1:length(unique(plottingData$instance))){
  name <- unique(plottingData$instance)[[i]]
  plotName <- paste0(i, "-", name)
  plotWilcoxon(name, plotName)
}
name <- "LiuPedagogical"

plotWilcoxon <- function(name, plotName){
  test <- plottingData
  test <- test[test$instance == name, ]
  # Add in missing points
  for(highFi in unique(test$feature_sample_highFiBudget)){
    for(lowFi in unique(test$feature_sample_lowFiBudget)){
      if(nrow(test[test$feature_sample_highFiBudget == highFi & 
                   test$feature_sample_lowFiBudget == lowFi, ]) == 0){
        test[nrow(test) + 1, c("feature_sample_highFiBudget", 'feature_sample_lowFiBudget',
                               'kriging_corrWilcoxon', 'cokriging_corrWilcoxon',
                               'kriging_corrMedian', 'cokriging_corrMedian',
                               'k 
                               riging_corrMean', 'cokriging_corrMean')] <- c(highFi, lowFi, 2, 2, 0, 0, 0, 0)
      }
    }
  }
  
  test$diffCorrMedian <- test$cokriging_corrMedian - test$kriging_corrMedian
  test$diffCorrMean <- test$cokriging_corrMean - test$kriging_corrMean
  test$krigingMedianPercent <- 100*(1 - test$kriging_corrMedian / pmax(test$kriging_corrMedian, test$cokriging_corrMedian))
  test$cokrigingMedianPercent <- 100*(1 - test$cokriging_corrMedian / pmax(test$kriging_corrMedian, test$cokriging_corrMedian))
  test[is.na(test$instances), c('krigingMedianPercent', 'cokrigingMedianPercent')] <- 0
  
  
  colors<-colorRampPalette(c("red", "darkblue"))(8)
  # mybreaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.00001, 2.00001)
  mybreaks <- c(0, 0.25, 0.5, 0.75, 1.00001, 1.25, 1.5, 1.75, 2.00001)
  
  krigWilcoxon <- ggplot(test, aes(feature_sample_lowFiBudget, feature_sample_highFiBudget, z = kriging_corrWilcoxon)) +
    geom_contour_filled(breaks = mybreaks) +
    scale_colour_manual(aesthetics = 'fill', drop = FALSE, values = colors) +
    labs(x = "", y = "High fidelity samples") +
    geom_point(data = test[!is.na(test$instances), ])
  
  cokrigWilcoxon <- ggplot(test, aes(feature_sample_lowFiBudget, feature_sample_highFiBudget, z = cokriging_corrWilcoxon)) +
    geom_contour_filled(breaks = mybreaks) +
    scale_colour_manual(aesthetics = 'fill', drop = FALSE, values = colors) +
    labs(x = "", y = "") +
    geom_point(data = test[!is.na(test$instances), ])
  
  
  colors<-colorRampPalette(c("red", "darkblue"))(11)
  mybreaks <- c(-2, -0.05, -0.025, -0.01, -0.0025, -0.001, 0.001, 0.0025, 0.01, 0.025, 0.05, 2)
  # mybreaks <- c(-2, -0.3, -0.25, -0.2, -0.15, -0.1, -0.075, -0.05, -0.025, 0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 2)
  
  corrDiffMedian <- ggplot(test, aes(feature_sample_lowFiBudget, feature_sample_highFiBudget, z = diffCorrMedian)) + 
    geom_contour_filled(breaks = mybreaks) +
    scale_colour_manual(aesthetics = 'fill', drop = FALSE, values = colors) +
    labs(x = "", y = "High fidelity samples") + 
    geom_point(data = test[!is.na(test$instances), ])
  
  corrDiffMean <- ggplot(test, aes(feature_sample_lowFiBudget, feature_sample_highFiBudget, z = diffCorrMean)) + 
    geom_contour_filled(breaks = mybreaks) +
    scale_colour_manual(aesthetics = 'fill', drop = FALSE, values = colors) +
    labs(x = "", y = "") + 
    geom_point(data = test[!is.na(test$instances), ])
  
  
  colors<-colorRampPalette(c("red", "darkblue"))(10)
  mybreaks <- c(-1, -0.5, -0.25, 0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1.001)
  
  corrPerfCoKrig <- ggplot(test, aes(feature_sample_lowFiBudget, feature_sample_highFiBudget, z = cokriging_corrMedian)) + 
    geom_contour_filled(breaks = mybreaks) +
    scale_colour_manual(aesthetics = 'fill', drop = FALSE, values = colors) +
    labs(x = "Low fidelity samples", y = "High fidelity samples") + 
    geom_point(data = test[!is.na(test$instances), ])
  
  corrPerfKrig <- ggplot(test, aes(feature_sample_lowFiBudget, feature_sample_highFiBudget, z = kriging_corrMedian)) + 
    geom_contour_filled(breaks = mybreaks) +
    scale_colour_manual(aesthetics = 'fill', drop = FALSE, values = colors) +
    labs(x = "Low fidelity samples", y = "") + 
    geom_point(data = test[!is.na(test$instances), ])
  
  
  
  colors<-colorRampPalette(c("red", "darkblue"))(10)
  mybreaks <- c(50, 40, 30, 20, 10, 5, 1, 0.5, 0.25, 0.1, -0.0001)
  
  percentKrig <- ggplot(test, aes(feature_sample_lowFiBudget, feature_sample_highFiBudget, z = krigingMedianPercent)) + 
    geom_contour_filled(breaks = mybreaks) +
    scale_colour_manual(aesthetics = 'fill', drop = FALSE, values = colors) +
    labs(x = "Low fidelity samples", y = "High fidelity samples") + 
    geom_point(data = test[!is.na(test$instances), ])
  
  percentCoKrig <- ggplot(test, aes(feature_sample_lowFiBudget, feature_sample_highFiBudget, z = cokrigingMedianPercent)) + 
    geom_contour_filled(breaks = mybreaks) +
    scale_colour_manual(aesthetics = 'fill', drop = FALSE, values = colors) +
    labs(x = "Low fidelity samples", y = "") + 
    geom_point(data = test[!is.na(test$instances), ])
  
  
  
  # plot <- ggarrange(krigWilcoxon, cokrigWilcoxon, 
  #           labels = c("Krig", "CoKrig"),
  #           ncol = 2, nrow = 1)
  # ggsave(paste0("data/plots/wilcoxonContourPlots/", plotName, ".png"), plot = plot, width = 30, height = 15, units = "cm")
  
  # plot <- ggarrange(krigWilcoxon, cokrigWilcoxon, 
  #                   corrDiffMean, corrDiffMedian,
  #                   corrPerfKrig, corrPerfCoKrig, 
  #                   labels = c("Krig", "CoKrig", "DiffMean", "DiffMedian", "KrigMedian", "CoKrigMedian"),
  #                   ncol = 2, nrow = 3)
  # ggsave(paste0("data/plots/wilcoxonContourPlots/", plotName, "All.png"), plot = plot, width = 40, height = 45, units = "cm")
  
  plot <- ggarrange(krigWilcoxon, cokrigWilcoxon,
                    corrDiffMean, corrDiffMedian,
                    percentKrig, percentCoKrig,
                    labels = c("Krig", "CoKrig", "DiffMean", "DiffMedian", "KrigPerc", "CoKrigPerc"),
                    ncol = 2, nrow = 3)
  ggsave(paste0("data/plots/wilcoxonContourPlots/", plotName, "New.png"), plot = plot, width = 40, height = 45, units = "cm")
  
}



# Going to do a bit of a play around
dataTieBreaker <- combineArrayResults("smallerExpandedTesting", 1, 890, 10)
dataNoTieBreaker <- combineArrayResults("smallerExpandedTestingRepeat", 1, 890, 10)

chosenData <- condensedData[1:14]
chosenData$algo_kriging <- chosenData$kriging_corrWilcoxon
chosenData$algo_cokriging <- chosenData$cokriging_corrWilcoxon
chosenData <- chosenData[, str_which(names(chosenData), "kriging_", negate = TRUE)]
# Only give LCC and budget features
chosenFeat <- paste0("feature_", read.csv("C:/Users/nandr/Documents/InstanceSpace-master/test.csv")[1,])


# chosenFeat <- c(colnames(condensedData)[str_which(colnames(condensedData), "LCC")])#,
#colnames(condensedData)[str_which(colnames(condensedData), "budget")],
#colnames(condensedData)[str_which(colnames(condensedData), "Budget")])

pcaData <- prcomp(condensedDataNormalised[chosenFeat], center = TRUE,scale. = TRUE)
summary(pcaData)
plot(cumsum(pcaData$sdev^2 / sum(pcaData$sdev^2)), type="b")
sum(pcaData$sdev^2 > 1)
relevantPcaData <- pcaData$x[, 1:10]

chosenData[paste0("feature_", colnames(relevantPcaData))] <- relevantPcaData

dir.create(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/nicoPlayAroundAgain/"))
write.csv(chosenData, paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/nicoPlayAroundAgain/metadata.csv"), quote = FALSE, row.names = FALSE)


results <- data.frame(matrix(nrow = 0, ncol = 0))
numResults <- 0
epsRange <- c(5, 4, 3, 2, 1.5, 1)

for(eps in epsRange){
  print(eps)
  numResults <- numResults + 1
  results[numResults, "epsilon"] <- eps
  featNames <- colnames(chosenData)[str_which(colnames(chosenData), "feature_")]
  Xfeat <- chosenData[, featNames]
  Ybin <- chosenData[, c('algo_kriging', 'algo_cokriging')]
  Ybin$algo_kriging <- as.numeric(chosenData$algo_kriging >= 0.5)
  Ybin$algo_cokriging <- as.numeric(chosenData$algo_cokriging >= 0.5)
  name <- 'condensedCorr'
  Xfeat <- data.matrix(Xfeat)
  Ybin <- data.matrix(Ybin)
  # Have the data, first purify
  labels <- purifyInstances(Xfeat, Ybin, eps)
  results[numResults, paste0(name, "VisaRatio")] <- sum(labels$visa) / sum(!labels$preclude)
  results[numResults, paste0(name, "InstancesRatio")] <- sum(!labels$preclude) / nrow(labels)
  Xfeat <- Xfeat[!labels$preclude, ]
  Ybin <- Ybin[!labels$preclude, ]
  # Now calculate uniformity
  vals <- calculateCVNND(Xfeat)
  results[numResults, c(paste0(name, "CV"), paste0(name, "UniformVal"))] <- vals
  print(results)
  # Here should be able to plot what is going on, maybe do it if have at least 2 rows
  if(nrow(results) > 1){
    plot(c(), c(), ylim = c(0, 1), xlim = c(min(epsRange), max(epsRange)))
    lines(results$epsilon, results$condensedCorrUniformVal, col="red", lty = 1)
    lines(results$epsilon, results$condensedCorrVisaRatio, col="green", lty = 1)
    lines(results$epsilon, results$condensedCorrInstancesRatio, col="blue", lty = 1)
    
    legend("topleft", legend=c("condensedCorrUniformVal", "condensedCorrVisaRatio", "condensedCorrInstancesRatio"),
           col = c("red", "green", "blue"),
           lty = c(1, 1, 1), cex = 0.5)
  }
}
secondChosenData <- chosenData
secondChosenData <- secondChosenData[!labels$preclude, ]
dir.create(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/nicoPlayAroundTwo/"))
write.csv(secondChosenData, paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/nicoPlayAroundTwo/metadata.csv"), quote = FALSE, row.names = FALSE)




# Next, have another go at using decision trees for prediction
dataset <- condensedData
dataset$algoGood <- dataset$kriging_corrWilcoxon >= 0.5
prefix <- "correlation"
# Get all features
allFeatures <- colnames(dataset)[str_which(colnames(dataset), "feature", negate = FALSE)]
# Get features with limited sample
limitedFeatures <- colnames(dataset)[str_which(colnames(dataset), "feature_sample", negate = FALSE)]
# Get perfect features
perfectFeatures <- allFeatures
for(feature in limitedFeatures){perfectFeatures <- perfectFeatures[perfectFeatures != feature]}

performDecisionTreePrediction(dataset, 'cokriging', 'majorityRule', 5, paste0(prefix, 'majority'), 0.05, 0)
performDecisionTreePrediction(dataset, 'cokriging', limitedFeatures, 5, paste0(prefix, "SampleFeatures"), 0.05, 0)
performDecisionTreePrediction(dataset, 'cokriging', perfectFeatures, 5, paste0(prefix, "PerfectFeatures"), 0.05, 0)
performDecisionTreePrediction(dataset, 'cokriging', allFeatures, 5, paste0(prefix, 'AllFeatures'), 0.05, 0)
performDecisionTreePrediction(dataset, 'cokriging', customFeatures, 5, paste0(prefix, 'CustomFeatures'), 0.05, 0)







# Selected instances

# Next, save them and run ISA. Would be interesting to see if supplying performance measure 
# changes anything (just to make sure all is running well). Will check it out.















# # Next is to get rid of some features based on correlation, to speed up instance selection
# keptFeat <- 0
# featCorrData <- data.frame(matrix(nrow = 0, ncol = 9))
# colnames(featCorrData) <- c("feature", "krigCorr", "krigAggCorr", "cokrigCorr", "cokrigAggCorr", "krigErr", "krigAggErr", "cokrigErr", "cokrigAggErr")
# for(feat in colnames(dataNormalised[, str_which(colnames(dataNormalised), "feature_")])){
#   keptFeat <- keptFeat + 1
#   featCorrData[keptFeat, ] <- c(feat, as.numeric(abs(cor(dataNormalised$kriging_modelCorrelation, dataNormalised[, feat]))),
#                                 as.numeric(abs(cor(condensedDataNormalised$kriging_corrMedian, condensedDataNormalised[, feat]))),
#                                 as.numeric(abs(cor(dataNormalised$cokriging_modelCorrelation, dataNormalised[, feat]))),
#                                 as.numeric(abs(cor(condensedDataNormalised$cokriging_corrMedian, condensedDataNormalised[, feat]))),
#                                 as.numeric(abs(cor(dataNormalised$kriging_modelError, dataNormalised[, feat]))),
#                                 as.numeric(abs(cor(condensedDataNormalised$kriging_errorMedian, condensedDataNormalised[, feat]))),
#                                 as.numeric(abs(cor(dataNormalised$cokriging_modelError, dataNormalised[, feat]))),
#                                 as.numeric(abs(cor(condensedDataNormalised$cokriging_errorMedian, condensedDataNormalised[, feat]))))
# 
# }
# 
# featCorrData$min <- apply(featCorrData[, -1], 1, FUN = min)
# featCorrData$max <- apply(featCorrData[, -1], 1, FUN = max)













# Excellent! Can now run instance filtering





















featCorrData <- featCorrData[rowSums(is.na(featCorrData)) == 0, ]

featCorrData$diff <- as.numeric(featCorrData$cokrigCorr) - as.numeric(featCorrData$cokrigAggCorr)

keptFeatures <- c('feature_sample_LCC_0_975',
                  'feature_error_diff_LCC_0_5',
                  'feature_error_diff_high_ela_distr_kurtosis',
                  'feature_error_perc_high_ela_meta_quad_w_interact_adj_r2',
                  'feature_sample_lowFiBudgetRatio',
                  'feature_sample_highFiBudgetRatio')
















# # HERE CODE CHOOSES A SUBSET OF THE DATA BY CHOOSING ABOUT 252 instances 
# # AND RUNS WITH SEED < 20. THIS REDUCES DATASET FROM 785160 TO 22680 (ROUGHLY 30 TIMES SMALLER)
# cleanedData <- read.table("data/combinedData/cleanedDataSet.txt", header = TRUE, sep = " ")
# augmentedResultsNormalised <- read.table("data/combinedData/normalisedFeatures.txt", header = TRUE, sep = ",")
# data <- augmentedResultsNormalised
# 
# # Going to remove instances by doing instance filtering only on the features calculated with prefect features
# data$instanceName <- sapply(strsplit(data$instance, "_"), "[[", 1)
# instanceNames <- unique(data$instanceName)
# singleInstanceData <- data[match(instanceNames, data$instanceName), ]
# singleInstanceData <- singleInstanceData[, str_which(colnames(singleInstanceData), "feature_sample", negate = TRUE)]
# singleInstanceData <- singleInstanceData[, c(1, str_which(colnames(singleInstanceData), "feature_"))]
# 
# # Here should have instances with only the large sample features, time to filter them
# name <- "initialFiltering"
# write.csv(singleInstanceData, paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/data/", name, ".csv"), quote = FALSE, row.names = FALSE)
# 
# epsilonValues <- seq(10.5, 4, -0.25)
# originalUniformValues <- c()
# originalVisaRatio <- c()
# originalRatioInstances <- c()
# 
# for(eps in epsilonValues){
#   print(eps)
#   originalUniformValues <- c(originalUniformValues, read.table(paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/cvnnd/CVNND_Pur_Ftr_", name, "_Dist_", sprintf(eps, fmt = '%#.3f'), ".csv"), header = TRUE, sep = ",")$Uniformity_All)
#   originalRatioInstances <- c(originalRatioInstances, nrow(read.table(paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/purified/Pur_Ftr_", name, "_Dist_", sprintf(eps, fmt = '%#.3f'), ".csv"), header = TRUE, sep = ",")) / nrow(singleInstanceData))
# }
# plot(c(), c(), ylim = c(0, 1), xlim = c(min(epsilonValues), max(epsilonValues)))
# lines(epsilonValues, originalUniformValues, col="red")
# lines(epsilonValues, originalVisaRatio, col="green")
# lines(epsilonValues, originalRatioInstances, col="blue")
# lines(c(2, 2), c(0,3))
# 
# eps <- 10.5
# chosenInstances <- read.table(paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/purified/Pur_Ftr_", name, "_Dist_", sprintf(eps, fmt = '%#.3f'), ".csv"), header = TRUE, sep = ",")
# nrow(chosenInstances)
# chosenInstances <- sapply(strsplit(chosenInstances$instance, "_"), "[[", 1)
# write(chosenInstances, "data/availableFunctions/chosenSubsubsetAllFunctions.txt")
# 
# data <- cleanedData
# data$instanceName <- sapply(strsplit(data$instance, "_"), "[[", 1)
# dataSubset <- data[data$instanceName == chosenInstances[[1]], ]
# for(instance in chosenInstances[2:length(chosenInstances)]){
#   print(instance)
#   dataSubset <- rbind(dataSubset, data[data$instanceName == instance, ])
# }
# dataSubset <- dataSubset[as.numeric(gsub('S', '', sapply(strsplit(dataSubset$instances, "_"), "[[", 4))) <= 5, ]
# dataSubset <- dataSubset[, -428]
# 
# # Can now save this!!
# write.table(dataSubset, "data/combinedData/cleanedDataSubset.txt", quote = FALSE, row.names = FALSE)




# 
# 
# # Let's start with literature everything
# cleanedData <- read.table("data/combinedData/cleanedDataSetReduced.txt", header = TRUE, sep = " ")
# augmentedResultsNormalised <- read.table("data/combinedData/normalisedFeaturesSubset.txt", header = TRUE, sep = ",")
# currData <- cleanedData
# 
# currData$cokriging_modelError <- round(currData$cokriging_modelError, 4)
# currData$kriging_modelError <- round(currData$kriging_modelError, 4)
# 
# 
# 
# # Here will do instance filtering with sample features, to hopefully 
# # reduce bias! Have removed data normalisation and algorithm labelling, so do that here first
# # Start with correlation
# for(name in c("correlationPerformance", "errorPerformance", "relErrorPerformance")){
#   print(name)
#   dataPerf <- data.frame(instances = cleanedData$instances)
#   if(name == "correlationPerformance"){
#     dataPerf$algo_cokriging <- as.numeric(cleanedData$cokriging_modelCorrelation >= 0.99*cleanedData$kriging_modelCorrelation)
#     dataPerf$algo_kriging <- as.numeric(cleanedData$kriging_modelCorrelation >= 0.99*cleanedData$cokriging_modelCorrelation)
#   }else if(name == "errorPerformance"){
#     dataPerf$algo_cokriging <- as.numeric(cleanedData$cokriging_modelError <= 1.01*cleanedData$kriging_modelError)
#     dataPerf$algo_kriging <- as.numeric(cleanedData$kriging_modelError <= 1.01*cleanedData$cokriging_modelError)
#   }else if(name == "relErrorPerformance"){
#     dataPerf$algo_cokriging <- (cleanedData$cokriging_modelError - cleanedData$kriging_modelError) / pmax(cleanedData$cokriging_modelError, cleanedData$kriging_modelError)
#     dataPerf$algo_kriging <- (cleanedData$kriging_modelError - cleanedData$cokriging_modelError) / pmax(cleanedData$cokriging_modelError, cleanedData$kriging_modelError)
#     bothZero <- (cleanedData$cokriging_modelError == 0) & (cleanedData$kriging_modelError == 0)
#     dataPerf[bothZero, 'algo_cokriging'] <- 0
#     dataPerf[bothZero, 'algo_kriging'] <- 0
#     dataPerf$algo_cokriging <- as.numeric(dataPerf$algo_cokriging <= 0.01)
#     dataPerf$algo_kriging <- as.numeric(dataPerf$algo_kriging <= 0.01)
#   }
#   # sampleFeat <- colnames(augmentedResultsNormalised[, str_which(colnames(augmentedResultsNormalised), "feature_sample")])
#   # dataPerf[, sampleFeat] <- augmentedResultsNormalised[match(dataPerf$instances, augmentedResultsNormalised$instances), sampleFeat]
#   feat <- colnames(augmentedResultsNormalised[, str_which(colnames(augmentedResultsNormalised), "feature_")])
#   dataPerf[, feat] <- augmentedResultsNormalised[match(dataPerf$instances, augmentedResultsNormalised$instances), feat]
#   write.csv(dataPerf, paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/data/", name, ".csv"), quote = FALSE, row.names = FALSE)
# }
# 
# name <- c("correlationPerformance", "errorPerformance", "relErrorPerformance")[[3]]
# 
# epsilonValues <- seq(15, 5.75, -0.25)
# originalUniformValues <- c()
# originalVisaRatio <- c()
# originalRatioInstances <- c()
# 
# for(eps in epsilonValues){
#   print(eps)
#   originalUniformValues <- c(originalUniformValues, read.table(paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/cvnnd/CVNND_Pur_Ftr&Good_", name, "_Dist_", sprintf(eps, fmt = '%#.3f'), ".csv"), header = TRUE, sep = ",")$Uniformity_All)
#   originalVisaRatio <- c(originalVisaRatio, nrow(read.table(paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/visa/ViSA_Ftr&Good_", name, "_Dist_", sprintf(eps, fmt = '%#.3f'), ".csv"), header = TRUE, sep = ",")) / nrow(read.table(paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/purified/Pur_Ftr&Good_", name, "_Dist_", sprintf(eps, fmt = '%#.3f'), ".csv"), header = TRUE, sep = ",")))
#   originalRatioInstances <- c(originalRatioInstances, nrow(read.table(paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/purified/Pur_Ftr&Good_", name, "_Dist_", sprintf(eps, fmt = '%#.3f'), ".csv"), header = TRUE, sep = ",")) / nrow(dataPerf))
# }
# plot(c(), c(), ylim = c(0, 1), xlim = c(min(epsilonValues), max(epsilonValues)))
# lines(epsilonValues, originalUniformValues, col="red")
# lines(epsilonValues, originalVisaRatio, col="green")
# lines(epsilonValues, originalRatioInstances, col="blue")
# lines(c(7, 7), c(0,1))
# 
# 
# 
# epsCorr <- 7
# epsErr <- 7
# epsErrRel <- 6.75
# 
# chosenInstancesCorr <- read.table(paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/purified/Pur_Ftr&Good_correlationPerformance_Dist_", sprintf(epsCorr, fmt = '%#.3f'), ".csv"), header = TRUE, sep = ",")$instances
# chosenInstancesErr <- read.table(paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/purified/Pur_Ftr&Good_errorPerformance_Dist_", sprintf(epsErr, fmt = '%#.3f'), ".csv"), header = TRUE, sep = ",")$instances
# chosenInstancesErrRel <- read.table(paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/purified/Pur_Ftr&Good_relErrorPerformance_Dist_", sprintf(epsErrRel, fmt = '%#.3f'), ".csv"), header = TRUE, sep = ",")$instances
# 
# 
# dataCorr <- cleanedData[match(chosenInstancesCorr, cleanedData$instances), ]
# dataError <- cleanedData[match(chosenInstancesErr, cleanedData$instances), ]
# dataErrorRel <- cleanedData[match(chosenInstancesErrRel, cleanedData$instances), ]
# 
# 
# dataError$cokriging_modelError <- round(dataError$cokriging_modelError, 4)
# dataError$kriging_modelError <- round(dataError$kriging_modelError, 4)
# dataErrorRel$cokriging_modelError <- round(dataErrorRel$cokriging_modelError, 4)
# dataErrorRel$kriging_modelError <- round(dataErrorRel$kriging_modelError, 4)
# 
# # 
# # # Now have data to work with! First take into separate datasets
# dataCorr <- cleanedData
# dataError <- cleanedData
# dataErrorRel <- cleanedData
# 
# # Get performance
# dataCorr$algo_kriging <- dataCorr$kriging_modelCorrelation
# dataCorr$algo_cokriging <- dataCorr$cokriging_modelCorrelation
# 
# dataError$algo_kriging <- dataError$kriging_modelError
# dataError$algo_cokriging <- dataError$cokriging_modelError
# 
# dataErrorRel$algo_cokriging <- (dataErrorRel$cokriging_modelError - dataErrorRel$kriging_modelError) / pmax(dataErrorRel$cokriging_modelError, dataErrorRel$kriging_modelError)
# dataErrorRel$algo_kriging <- (dataErrorRel$kriging_modelError - dataErrorRel$cokriging_modelError) / pmax(dataErrorRel$cokriging_modelError, dataErrorRel$kriging_modelError)
# bothZero <- (dataErrorRel$cokriging_modelError == 0) & (dataErrorRel$kriging_modelError == 0)
# dataErrorRel[bothZero, 'algo_cokriging'] <- 0
# dataErrorRel[bothZero, 'algo_kriging'] <- 0
# 
# dataCorr <- dataCorr[, -which(names(dataCorr) %in% c('kriging_modelError', 'kriging_modelCorrelation', 'kriging_time', 'cokriging_modelError', 'cokriging_modelCorrelation', 'cokriging_time'))]
# dataError <- dataError[, -which(names(dataError) %in% c('kriging_modelError', 'kriging_modelCorrelation', 'kriging_time', 'cokriging_modelError', 'cokriging_modelCorrelation', 'cokriging_time'))]
# dataErrorRel <- dataErrorRel[, -which(names(dataErrorRel) %in% c('kriging_modelError', 'kriging_modelCorrelation', 'kriging_time', 'cokriging_modelError', 'cokriging_modelCorrelation', 'cokriging_time'))]
# 
# corrPerc <- 0.05
# dataCorr <- dataCorr[dataCorr$feature_sample_dimension == 1, ]
# # dataCorr <- dataCorr[dataCorr$algo_cokriging < 0.95 | dataCorr$algo_kriging < 0.95, ]
# dataCorr <- dataCorr[dataCorr$algo_cokriging < (1-corrPerc)* dataCorr$algo_kriging | dataCorr$algo_kriging < (1-corrPerc)* dataCorr$algo_cokriging, ]
# dataCorr <- dataCorr[dataCorr$algo_cokriging < 0.95 | dataCorr$algo_kriging < 0.95, ]
# 
# # Removing basic_lower and basic_upper as these should not have an impact 
# # dataCorr <- dataCorr[, -c(str_which(colnames(dataCorr), "basic_lower"), str_which(colnames(dataCorr), "basic_upper"))]
# dataCorr <- dataCorr[, -str_which(colnames(dataCorr), "basic")]
# 
# 
# customFeatures <- c("feature_sample_highFiBudgetRatio",
#                     "feature_sample_CC",
#                     "feature_sample_LCC_0_5",
#                     "feature_sample_LCC_sd",
#                     "feature_sample_high_disp_diff_mean_10")
# keptRows <- c("instances", "Source", "algo_kriging", "algo_cokriging", customFeatures)
# 
# 
# 
# 
# for(set in c("corr", "err", "relErr")){
#   if(set == "corr"){
#     print("Correlation based")
#     dataset <- dataCorr
#     dataset$algoGood <- dataset$algo_cokriging >= (1-corrPerc)*dataset$algo_kriging
#     prefix <- "correlation"
#    
#   }
#   if(set == "err"){
#     print("Error based")
#     dataset <- dataError
#     dataset$algoGood <- dataset$algo_cokriging <= 1.01*dataset$algo_kriging
#     prefix <- "error"
#   }
#   if(set == "relErr"){
#     print("Relative error based")
#     dataset <- dataError
#     dataset$algoGood <- dataset$algo_cokriging <= 0.01
#     prefix <- "errorRel"
#   }
#   
#   # Get all features
#   allFeatures <- colnames(dataset)[str_which(colnames(dataset), "feature", negate = FALSE)]
#   # Get features with limited sample
#   limitedFeatures <- colnames(dataset)[str_which(colnames(dataset), "feature_sample", negate = FALSE)]
#   # Get perfect features
#   perfectFeatures <- allFeatures
#   for(feature in limitedFeatures){perfectFeatures <- perfectFeatures[perfectFeatures != feature]}
#   
#   performDecisionTreePrediction(dataset, 'cokriging', 'majorityRule', 5, paste0(prefix, 'majority'), 0.05, 0)
#   performDecisionTreePrediction(dataset, 'cokriging', limitedFeatures, 5, paste0(prefix, "SampleFeatures"), 0.05, 0)
#   performDecisionTreePrediction(dataset, 'cokriging', perfectFeatures, 5, paste0(prefix, "PerfectFeatures"), 0.05, 0)
#   performDecisionTreePrediction(dataset, 'cokriging', allFeatures, 5, paste0(prefix, 'AllFeatures'), 0.05, 0)
#   performDecisionTreePrediction(dataset, 'cokriging', customFeatures, 5, paste0(prefix, 'CustomFeatures'), 0.05, 0)
# }
# 
# 
# krigTest <- combineArrayResults("smallerExpandedTesting", 1, 756, 10)
# # write.table(krigTest, "data/combinedData/smallerExpandedTestingKrigWithoutTieBreak.txt", quote = FALSE, row.names = FALSE)
# test <- krigTest
# test <- test[test$det > 0.999, ]
# test <- test[test$modelCorrelation < 0.2, ]
# test$seeds <- paste0(test$seed, "-", test$seed)
# test <- test[c("instance", "technique", "highFiBudget", "lowFiBudget", "seeds")]
# test <- test[str_which(test$instance, "DisturbanceBasedFunction36", negate = TRUE), ]
# test <- test[str_which(test$instance, "DisturbanceBasedFunction37", negate = TRUE), ]
# test <- test[str_which(test$instance, "DisturbanceBasedFunction38", negate = TRUE), ]
# 
# write.table(test, "data/runScripts/miniKrigTestNoTieBreaker.txt", quote = FALSE, row.names = FALSE)
# write.table(test, "data/runScripts/miniKrigTestWithTieBreaker.txt", quote = FALSE, row.names = FALSE)
# 
# testNoTieBreaker <- combineArrayResults("miniKrigTestNoTieBreaker", 1, 1, 954)
# testWithTieBreaker <- combineArrayResults("miniKrigTestWithTieBreaker", 1, 1, 954)
# 
# processed <- min(nrow(testNoTieBreaker), nrow(testWithTieBreaker))
# testNoTieBreaker <- testNoTieBreaker[1:processed, ]
# testWithTieBreaker <- testWithTieBreaker[1:processed, ]
# 
# combined <- testNoTieBreaker[c("instance", "highFiBudget", "lowFiBudget", "seed")]
# combined$likelihoodBefore <- testNoTieBreaker$likelihood
# combined$likelihoodAfter <- testWithTieBreaker$likelihood
# # combined$detBefore <- testNoTieBreaker$det
# # combined$detAfter <- testWithTieBreaker$det
# combined$corrBefore <- testNoTieBreaker$modelCorrelation
# combined$corrAfter <- testWithTieBreaker$modelCorrelation
# 
# combined$likelihoodChange <- combined$likelihoodAfter - combined$likelihoodBefore
# # combined$detChange <- combined$detAfter - combined$detBefore
# combined$corrChange <- combined$corrAfter - combined$corrBefore
# 
# # Only care about negative loglikelihood change, otherwise the tie breaker just helped with improvement of the log likelihood
# improved <- combined[combined$likelihoodChange > 0, ]
# combined <- combined[combined$likelihoodChange < 0, ]
# 
# noDimension1 <- combined
# noDimension1 <- noDimension1[str_which(noDimension1$instance, "dim1-", negate = TRUE), ]
# noDimension1 <- noDimension1[str_which(noDimension1$instance, "DisturbanceBasedFunction1-", negate = TRUE), ]
# noDimension1 <- noDimension1[str_which(noDimension1$instance, "DisturbanceBasedFunction2-", negate = TRUE), ]
# noDimension1 <- noDimension1[str_which(noDimension1$instance, "DisturbanceBasedFunction3-", negate = TRUE), ]
# noDimension1 <- noDimension1[str_which(noDimension1$instance, "DisturbanceBasedFunction4-", negate = TRUE), ]
# noDimension1 <- noDimension1[str_which(noDimension1$instance, "DisturbanceBasedFunction5-", negate = TRUE), ]
# noDimension1 <- noDimension1[noDimension1$instance != "LiuPedagogical", ]
# noDimension1 <- noDimension1[noDimension1$instance != "ShiGramacyLee", ]
# noDimension1 <- noDimension1[noDimension1$instance != "ShiCurrinSin", ]
# noDimension1 <- noDimension1[noDimension1$instance != "ShiHolsclaw", ]
# noDimension1 <- noDimension1[noDimension1$instance != "ShiSantner", ]
# noDimension1 <- noDimension1[str_which(noDimension1$instance, "SongToalForretal", negate = TRUE), ]
# noDimension1 <- noDimension1[rownames(unique( noDimension1[ , -1 ] )), ]
# 
# 
# yesDimension1 <- combined
# yesDimension1 <- yesDimension1[!(yesDimension1$instance %in% noDimension1$instance), ]
# yesDimension1 <- yesDimension1[rownames(unique( yesDimension1[ , -1 ] )), ]
# 
# 
# 
# # Repeat work with cokrig
# cokrigTest <- combineArrayResults("smallerExpandedTestingCoKrig", 1, 7560)
# 
# test <- read.table("data/combinedData/smallerExpandedTestingCoKrig.txt", sep = " ")
# colnames(test) <- test[1, ]
# test <- test[-1, ]
# test <- test[1:(992*5), ]
# 
# # write.table(cokrigTest, "data/combinedData/smallerExpandedTestingCoKrig.txt", quote = FALSE, row.names = FALSE)
# test <- cokrigTest
# test <- test[!is.na(test$likelihoodLow), ]
# noDimension1 <- test
# noDimension1 <- noDimension1[str_which(noDimension1$instance, "dim1-", negate = TRUE), ]
# noDimension1 <- noDimension1[str_which(noDimension1$instance, "DisturbanceBasedFunction1-", negate = TRUE), ]
# noDimension1 <- noDimension1[str_which(noDimension1$instance, "DisturbanceBasedFunction2-", negate = TRUE), ]
# noDimension1 <- noDimension1[str_which(noDimension1$instance, "DisturbanceBasedFunction3-", negate = TRUE), ]
# noDimension1 <- noDimension1[str_which(noDimension1$instance, "DisturbanceBasedFunction4-", negate = TRUE), ]
# noDimension1 <- noDimension1[str_which(noDimension1$instance, "DisturbanceBasedFunction5-", negate = TRUE), ]
# noDimension1 <- noDimension1[noDimension1$instance != "LiuPedagogical", ]
# noDimension1 <- noDimension1[noDimension1$instance != "ShiGramacyLee", ]
# noDimension1 <- noDimension1[noDimension1$instance != "ShiCurrinSin", ]
# noDimension1 <- noDimension1[noDimension1$instance != "ShiHolsclaw", ]
# noDimension1 <- noDimension1[noDimension1$instance != "ShiSantner", ]
# noDimension1 <- noDimension1[str_which(noDimension1$instance, "SongToalForretal", negate = TRUE), ]
# yesDimension1 <- test
# yesDimension1 <- yesDimension1[!(yesDimension1$instance %in% noDimension1$instance), ]
# noDimension1 <- noDimension1[rownames(unique( noDimension1[ , -1 ] )), ]
# yesDimension1 <- yesDimension1[rownames(unique( yesDimension1[ , -1 ] )), ]
# 
# 
# # As a continuation of this, currently looking at models which give all function change in only one direction
# Kriging2d <- combineArrayResults("specialTesting", 1, 1, 1494)
# Kriging2dBad <- combineArrayResults("specialTestingBad", 1, 1, 1494)
# 
# processed <- min(nrow(Kriging2d), nrow(Kriging2dBad))
# Kriging2d <- Kriging2d[1:processed, ]
# Kriging2dBad <- Kriging2dBad[1:processed, ]
# 
# combined <- Kriging2d[c("instance", "highFiBudget", "lowFiBudget", "seed")]
# combined$likelihoodBefore <- Kriging2dBad$likelihood
# combined$likelihoodAfter <- Kriging2d$likelihood
# # combined$detBefore <- testNoTieBreaker$det
# # combined$detAfter <- testWithTieBreaker$det
# combined$corrBefore <- Kriging2dBad$modelCorrelation
# combined$corrAfter <- Kriging2d$modelCorrelation
# combined$errorBefore <- Kriging2dBad$modelError
# combined$errorAfter <- Kriging2d$modelError
# combined$distBefore <- Kriging2dBad$minDist
# combined$distAfter <- Kriging2d$minDist
# 
# 
# combined$likelihoodChange <- combined$likelihoodAfter - combined$likelihoodBefore
# # combined$detChange <- combined$detAfter - combined$detBefore
# combined$corrChange <- combined$corrAfter - combined$corrBefore
# combined$errorChange <- combined$errorAfter - combined$errorBefore
# 
# combined$distChange <- combined$distAfter - combined$distBefore
# 
# 
# 
# 
# test <- combineArrayResults("specialTesting", 1, 82, 5)
# test <- test[test$seed == 2 & test$instance == 'ShiGramacyLee', ]
# 








# I believe the code has now been improved; moving back to data while new data runs,
# trying features which are ratios and differences (see if they work better than analysis on the intermediate function)
cleanedData <- read.table("data/combinedData/cleanedDataSetReduced.txt", header = TRUE, sep = " ")
cleanedDataCondensed <- read.table("data/combinedData/cleanedDataSetReducedCondensed.txt", header = TRUE, sep = " ")
dataCorr <- cleanedData
dataCorrCondensed <- cleanedDataCondensed
colnames(dataCorr) <- c("instances", colnames(dataCorr)[2:length(colnames(dataCorr))])
colnames(dataCorrCondensed) <- c("instances", colnames(dataCorrCondensed)[2:length(colnames(dataCorrCondensed))])

# Get performance
dataCorr$algo_kriging <- dataCorr$kriging_modelCorrelation + 1
dataCorr$algo_cokriging <- dataCorr$cokriging_modelCorrelation + 1
dataCorr$labelled_algo_kriging <- as.numeric(dataCorr$algo_kriging >= 0.99*dataCorr$algo_cokriging)
dataCorr$labelled_algo_cokriging <- as.numeric(dataCorr$algo_cokriging >= 0.99*dataCorr$algo_kriging)

dataCorrCondensed$algo_kriging <- dataCorrCondensed$kriging_corrMedian + 1
dataCorrCondensed$algo_cokriging <- dataCorrCondensed$cokriging_corrMedian + 1
dataCorrCondensed$labelled_algo_kriging <- as.numeric(dataCorrCondensed$kriging_corrWilcoxon >= 0.5)
dataCorrCondensed$labelled_algo_cokriging <- as.numeric(dataCorrCondensed$cokriging_corrWilcoxon >= 0.5)


dataCorr <- dataCorr[, -which(names(dataCorr) %in% c('kriging_modelError', 'kriging_modelCorrelation', 'kriging_time', 'cokriging_modelError', 'cokriging_modelCorrelation', 'cokriging_time'))]
dataCorrCondensed <- dataCorrCondensed[, -which(names(dataCorrCondensed) %in% c('kriging_mean', 'kriging_median', 'kriging_wilcoxon', 'cokriging_mean', 'cokriging_median', 'cokriging_wilcoxon'))]

# Add feature which compares sample feature with "real" feature value
allColNames <- colnames(dataCorr)
featNames <- allColNames[str_which(allColNames, "feature_")]
featNames <- featNames[str_which(featNames, "mid", negate = TRUE)]
featNames <- featNames[c(str_which(featNames, "_low_"), str_which(featNames, "_high_"))]
featNames <- gsub('feature_','', featNames)
sampleFeats <- featNames[str_which(featNames, "sample_")]
featNames <- featNames[str_which(featNames, "sample_", negate = TRUE)]
sampleFeatsHigh <- sampleFeats[str_which(sampleFeats, "high_")]
sampleFeatsLow <- sampleFeats[str_which(sampleFeats, "low_")]
sampleFeatsHigh <- gsub('sample_high_','', sampleFeatsHigh)
sampleFeatsLow <- gsub('sample_low_','', sampleFeatsLow)

featNamesSampleAll <- allColNames[str_which(allColNames, "feature_sample")]
featNamesSampleAll <- gsub('feature_sample_','', featNamesSampleAll)


for(i in c(1, 2)){
  if(i == 1){
    data <- dataCorr
  }else{
    data <- dataCorrCondensed
  }
  print(i)
  # for(feat in sampleFeatsHigh){
  #   print(paste0("Working on ", feat))
  #   if(!(feat %in% sampleFeatsLow)){print(feat)}
  #   data[, paste0("feature_sample_diff_", feat)] <- data[, paste0("feature_sample_low_", feat)] - data[, paste0("feature_sample_high_", feat)]
  #   data[, paste0("feature_sample_ratio_", feat)] <- data[, paste0("feature_sample_low_", feat)] / data[, paste0("feature_sample_high_", feat)]
  #   data[, paste0("feature_sample_perc_", feat)] <- 100*(data[, paste0("feature_sample_low_", feat)] - data[, paste0("feature_sample_high_", feat)]) / data[, paste0("feature_sample_high_", feat)]
  #   
  #   indices <- (data[, paste0("feature_sample_low_", feat)] == 0) & (data[, paste0("feature_sample_high_", feat)] == 0)
  #   if(sum(indices) > 0){
  #     print(feat)
  #   }
  #   data[indices, paste0("feature_sample_ratio_", feat)] <- 1
  #   data[indices, paste0("feature_sample_perc_", feat)] <- 0
  # }
  
  for(feat in featNamesSampleAll){
    print(paste0("Working on ", feat))
    if(!(paste0("feature_", feat) %in% colnames(data))){next}
    data[, paste0("feature_error_diff_", feat)] <- data[, paste0("feature_sample_", feat)] - data[, paste0("feature_", feat)]
    data[, paste0("feature_error_ratio_", feat)] <- data[, paste0("feature_sample_", feat)] / data[, paste0("feature_", feat)]
    data[, paste0("feature_error_perc_", feat)] <- 100*(data[, paste0("feature_sample_", feat)] - data[, paste0("feature_", feat)]) / data[, paste0("feature_", feat)]
    
    indices <- (data[, paste0("feature_", feat)] == 0) & (data[, paste0("feature_sample_", feat)] == 0)
    if(sum(indices) > 0){
      print(feat)
    }
    data[indices, paste0("feature_error_ratio_", feat)] <- 1
    data[indices, paste0("feature_error_perc_", feat)] <- 0
  }
  
  
  
  for(name in colnames(data[, str_which(colnames(data), "feature_", negate = FALSE)])){
    if(sum(is.na(data[, name])) > 0 | !is.finite(sum(data[, name]))){
      print(paste0("Removing ", name))
      data <- data[colnames(data) != name]    
    }
  }
  data <- data[str_which(colnames(data), "_basic_", negate = TRUE)]
  if(i == 1){
    dataCorr <- data
  }else{
    dataCorrCondensed <- data
  }
}

dataCorrSubset <- dataCorr
# for(i in 20:20){
#   dataCorrSubset <- dataCorrSubset[str_which(dataCorrSubset$instances, paste0("S", i), negate = TRUE), ]
# }

dir.create(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrTest/"))
dir.create(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrTestCondensed/"))
dir.create(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrTestCondensedSuppliedPerformanceLabelNotPoly/"))
dir.create(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrTestCondensedSuppliedPerformanceLabelAndFeatures/"))
dir.create(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrTestCondensedSuppliedFeatures/"))


write.csv(dataCorrCondensed, paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrTestCondensed/metadata.csv"), quote = FALSE, row.names = FALSE)
write.csv(dataCorrCondensed, paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrTestCondensedSuppliedPerformanceLabelNotPoly/metadata.csv"), quote = FALSE, row.names = FALSE)
write.csv(dataCorrCondensed[, c("instances", "algo_kriging", "algo_cokriging", "labelled_algo_kriging", "labelled_algo_cokriging", keptFeatures)], paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrTestCondensedSuppliedPerformanceLabelAndFeatures/metadata.csv"), quote = FALSE, row.names = FALSE)
write.csv(dataCorrCondensed[, c("instances", "algo_kriging", "algo_cokriging", "labelled_algo_kriging", "labelled_algo_cokriging", keptFeatures)], paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrTestCondensedSuppliedFeatures/metadata.csv"), quote = FALSE, row.names = FALSE)


write.csv(dataCorr, paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrTest/metadata.csv"), quote = FALSE, row.names = FALSE)





write.csv(dataCorr, paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrTest/metadata.csv"), quote = FALSE, row.names = FALSE)

dataCorrCondensedFeatChosen <- dataCorrCondensed
dataCorrCondensedFeatChosen <- dataCorrCondensedFeatChosen[, c("instances", "algo_kriging", "algo_cokriging", "labelled_algo_kriging", 'labelled_algo_cokriging',
                                                               "feature_error_diff_LCC_0_7", "feature_LCC_0_7",
                                                               "feature_sample_highFiBudget", "feature_sample_highFiBudgetRatio",
                                                               "feature_sample_lowFiBudget", "feature_sample_lowFiBudgetRatio")]

write.csv(dataCorrCondensedFeatChosen, paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrTestCondensed/metadata.csv"), quote = FALSE, row.names = FALSE)


test <- dataCorrCondensed[1:100, str_which(colnames(dataCorrCondensed), "labelled", negate = TRUE)]







# Find correlation myself!
keptFeat <- 0
featCorrData <- data.frame(matrix(nrow = 0, ncol = 5))
colnames(featCorrData) <- c("feature", "krigCorr", "krigAggCorr", "cokrigCorr", "cokrigAggCorr")
for(feat in colnames(dataCorr[, str_which(colnames(dataCorr), "feature_")])){
  keptFeat <- keptFeat + 1
  featCorrData[keptFeat, ] <- c(feat, as.numeric(abs(cor(dataCorr$algo_kriging, dataCorr[, feat]))), 
                                as.numeric(abs(cor(dataCorrCondensed$algo_kriging, dataCorrCondensed[, feat]))),
                                as.numeric(abs(cor(dataCorr$algo_cokriging, dataCorr[, feat]))),
                                as.numeric(abs(cor(dataCorrCondensed$algo_cokriging, dataCorrCondensed[, feat]))))
  print(paste0(feat, " corr with krig ", cor(dataCorr$algo_kriging, dataCorr[, feat]), "(", 
               cor(dataCorrCondensed$algo_kriging, dataCorrCondensed[, feat]), ") cokrig ", 
               cor(dataCorr$algo_cokriging, dataCorr[, feat]), "(", 
               cor(dataCorrCondensed$algo_cokriging, dataCorrCondensed[, feat]), ")"))
}

featCorrData <- featCorrData[rowSums(is.na(featCorrData)) == 0, ]

featCorrData$diff <- as.numeric(featCorrData$cokrigCorr) - as.numeric(featCorrData$cokrigAggCorr)

keptFeatures <- c('feature_sample_LCC_0_975',
                  'feature_error_diff_LCC_0_5',
                  'feature_error_diff_high_ela_distr_kurtosis',
                  'feature_error_perc_high_ela_meta_quad_w_interact_adj_r2',
                  'feature_sample_lowFiBudgetRatio',
                  'feature_sample_highFiBudgetRatio')


test <- featCorrData[str_which(featCorrData$feature, "LC"), ]




# dataCorr <- dataCorr[dataCorr$feature_sample_dimension == 1, ]
# dataCorr <- dataCorr[abs(dataCorr$algo_cokriging - dataCorr$algo_kriging) > 0.1, ]

print("Correlation based")
dataset <- dataCorr
dataset$algoGood <- dataset$algo_cokriging  - dataset$algo_kriging >= 0.1
prefix <- "correlation"

# Get all features
allFeatures <- colnames(dataset)[str_which(colnames(dataset), "feature", negate = FALSE)]
# Get features with limited sample
limitedFeatures <- colnames(dataset)[str_which(colnames(dataset), "feature_sample", negate = FALSE)]
# Get perfect features
perfectFeatures <- allFeatures
for(feature in limitedFeatures){perfectFeatures <- perfectFeatures[perfectFeatures != feature]}

performDecisionTreePrediction(dataset, 'cokriging', 'majorityRule', 5, paste0(prefix, 'majority'), 0.05, 0)
performDecisionTreePrediction(dataset, 'cokriging', limitedFeatures, 5, paste0(prefix, "SampleFeatures"), 0.025, 0)
performDecisionTreePrediction(dataset, 'cokriging', perfectFeatures, 5, paste0(prefix, "PerfectFeatures"), 0.05, 0)
performDecisionTreePrediction(dataset, 'cokriging', allFeatures, 5, paste0(prefix, 'AllFeatures'), 0.05, 0)



dataCorr[, "labelled_algo_kriging"] <- as.numeric(dataCorr$algo_kriging - dataCorr$algo_cokriging >= -0.05)
dataCorr[, "labelled_algo_cokriging"] <- as.numeric(dataCorr$algo_cokriging - dataCorr$algo_kriging >= -0.05)


subsetData <- dataCorr
subsetData <- subsetData[subsetData$feature_sample_dimension == 2, ]
write(unique(sapply(strsplit(subsetData$instance, "_"), "[[", 1)), "data/availableFunctions/chosenSubsetAllFunctions2d.txt") 

subsetData <- subsetData[subsetData$feature_sample_dimension == 2 & 
                           subsetData$feature_sample_highFiBudget == 10 & 
                           subsetData$feature_sample_lowFiBudget == 25, ]

dir.create(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrTest/"))
write.csv(subsetData, paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrTest/metadata.csv"), quote = FALSE, row.names = FALSE)









projection <- read.table("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrCondensedOriginalWilcoxon/projection_matrix.csv", header = TRUE, sep = ",")
z1Vals <- projection[projection$Row == "Z_{1}", ]
z1Vals <- z1Vals[, -1]
z2Vals <- projection[projection$Row == "Z_{2}", ]
z2Vals <- z2Vals[, -1]

paste0("feature_", colnames(projection[-1]))

vals <- condensedDataNormalised[c(paste0("feature_", colnames(projection[-1])))]
# processedFeatures <- read.table("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrTest/feature_process.csv", header = TRUE, sep = ",")
# order <- match(subsetData$instances, processedFeatures$Row)
# vals <- subsetData
# vals <- vals[str_which(colnames(vals), "feature_", negate = TRUE)]
# vals[colnames(processedFeatures[-1])] <- processedFeatures[order,-1]
vals <- condensedDataNormalised
vals$z1 <- 0
vals$z2 <- 0
for(name in colnames(projection[-1])){
  realName <- paste0("feature_", name)
  vals$z1 <- vals$z1 + z1Vals[1, name] * vals[, realName]
  vals$z2 <- vals$z2 + z2Vals[1, name] * vals[, realName]
}

# vals$z1 <- vals$z1 + 2 * vals$feature_sample_mid_ela_meta_quad_w_interact_adj_r2


xVals <- c(min(vals$z1), max(vals$z1))
yVals <- c(min(vals$z2), max(vals$z2))

vals <- vals[order(-vals$feature_sample_lowFiBudgetRatio), ]

ggplot(data = vals, aes(x = z1, y = z2, col = feature_sample_lowFiBudgetRatio)) +
  # geom_point(shape = ".") +
  geom_point(size = 0.5) +
  # geom_point() +
  scale_x_continuous(limits = xVals) +
  scale_y_continuous(limits = yVals) +
  # scale_colour_gradientn(colours = colourGradient) +
  theme_bw() +
  ggtitle("Y vs. X with Z-Color")



vals <- vals[vals$z1 > 1.5 & 
               vals$z1 < 1.75 &
               vals$z2 > -1.25 &
               vals$z2 < -0.75, ]



test <- dataCorr[!dataCorr$feature_sample_high_disp_diff_median_25 >= -1.08, ]
test <- test[, c("instances", "algo_kriging", "algo_cokriging")]
test$diff <- test$algo_cokriging - test$algo_kriging





test <- subsetData[subsetData$instances == 'DisturbanceBasedFunction9-seed1-disth2-height0-radius0.025-freq2-amp1.5_H10_L25_S1' | 
                     subsetData$instances == 'DisturbanceBasedFunction9-seed1-disth2-height0-radius0.025-freq2-amp1.5_H10_L25_S3' |
                     subsetData$instances == 'DisturbanceBasedFunction9-seed1-disth2-height0-radius0.025-freq2-amp1.5_H10_L25_S4', ]























# Clean up to remove all rows which really are doubles
for(i in 1:nrow(noDimension1)){
  if(i > nrow(noDimension1)){break}
  noDimension1 <- noDimension1[noDimension1[-1] != noDimension1[i, -1]]
}



testNoTieBreaker$likelihoodImproved <- testWithTieBreaker$likelihood
testNoTieBreaker$likelihoodDiff <- testNoTieBreaker$likelihoodImproved - testNoTieBreaker$likelihood





test1 <- test[test$lowFiBudget == 25, ]
test2 <- test[test$lowFiBudget == 50, ]


test <- condensedDataNormalised
test <- test[, c("instances", "feature_sample_dimension", "feature_sample_highFiBudgetRatio", "algo_kriging", "algo_cokriging")]

test$diff <- test$algo_cokriging - test$algo_kriging
test$abs_krig <- abs(test$algo_kriging)
test$abs_cokrig <- abs(test$algo_cokriging)

test <- test[str_which(test$instances, "LiuAckley10"), ]

test <- test[test$feature_sample_dimension == 1, ]






# See if there is any correlation between features and this difference
# NICOCORRS
labels <- read.table("data/combinedData/instanceFiltering/condensedCorrWilcoxoneps8.txt", header = TRUE, sep = " ")
test <- condensedDataNormalised
test$kriging_corrWilcoxon <- condensedData$kriging_corrWilcoxon
test$cokriging_corrWilcoxon <- condensedData$cokriging_corrWilcoxon

test <- test[!labels$preclude, ]
correlations <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(correlations) <- c("feature", "krigCorr", "cokrigCorr")
# max <- 0.3
for(feat in colnames(test[, str_which(colnames(test), "feature_")])){
  correlations[nrow(correlations) + 1, ] <- c(feat, 
                                              abs(cor(test[, feat], test$cokriging_corrWilcoxon)), 
                                              abs(cor(test[, feat], test$kriging_corrWilcoxon)))
  # if(is.na(cor(test[, feat], test$cokriging_corrWilcoxon))){next}
  # if(cor(test[, feat], test$cokriging_corrWilcoxon) > max){
  #   # max <- cor(test[, feat], test$kriging_corrWilcoxon)
  #   print(paste0(feat, " - ", cor(test[, feat], test$cokriging_corrWilcoxon)))
  # }
}

nicoFeats <- c()
# High corr features to use
nicoFeats <- c('feature_error_diff_high_ela_meta_lin_w_interact_adj_r2',
               'feature_error_diff_high_pca_expl_var_PC1_cor_init',
               'feature_high_ela_meta_quad_simple_adj_r2',
               'feature_sample_high_nbc_nb_fitness_cor')
# Low corr features to use
# None - correlation below 0.2 for both Krig and CoKrig
nicoFeats <- c(nicoFeats,
               'feature_sample_mid_pca_expl_var_cor_init',
               'feature_sample_mid_nbc_dist_ratio_coeff_var',
               'feature_error_diff_mid_pca_expl_var_cor_init')
# Now sample features
nicoFeats <- c(nicoFeats, 'feature_sample_budgetRatio', 'feature_sample_highFiBudgetRatio')

# Finally, look at relationship features (i.e. LCC, CC, RRMSE, etc...)
nicoFeats <- c(nicoFeats,
               'feature_LCC_0_5',
               'feature_sample_LCC_coeff',
               'feature_sample_LCC_0_5',
               'feature_error_diff_LCC_0_5')

secondTest <- correlations
secondTest <- secondTest[c(7:20, 146:160, 328:341), ]

secondTest <- secondTest[c(str_which(secondTest$feature, "budget"),
                           str_which(secondTest$feature, "Budget")), ]



test <- test[test$feature_sample_high_ela_meta_quad_w_interact_adj_r2 >= 0.985, ]


test <- test[test$feature_sample_high_ela_meta_lin_simple_adj_r2 >= -0.102, ]
test <- test[test$feature_sample_high_ela_meta_lin_simple_coef_min < 2.86, ]
test <- test[test$feature_sample_low_ela_distr_kurtosis >= -0.263, ]

nrow(test) / nrow(dataCorr) * 100
test <- test[, c("instances", "algo_kriging", "algo_cokriging", "diff")]
test$cokrigingGood <- test$algo_cokriging >= 0.99* test$algo_kriging


# 
# dataCorr <- dataCorr[dataCorr$feature_sample_dimension == 1, ]
# dataCorr <- dataCorr[dataCorr$algo_cokriging < 0.95 | dataCorr$algo_kriging < 0.95, ]
# 


print("Percentage cokriging is good")
sum(dataCorr$algo_cokriging >= 0.99*dataCorr$algo_kriging) / nrow(dataCorr) * 100
print("Percentage kriging is good")
sum(dataCorr$algo_kriging >= 0.99*dataCorr$algo_cokriging) / nrow(dataCorr) * 100
print("Percentage only cokriging is good")
sum((dataCorr$algo_cokriging >= 0.99*dataCorr$algo_kriging) & !(dataCorr$algo_kriging >= 0.99*dataCorr$algo_cokriging)) / nrow(dataCorr) * 100
print("Percentage only kriging is good")
sum(!(dataCorr$algo_cokriging >= 0.99*dataCorr$algo_kriging) & (dataCorr$algo_kriging >= 0.99*dataCorr$algo_cokriging)) / nrow(dataCorr) * 100


dataCorr <- dataCorr[, str_which(colnames(dataCorr), "feature_", negate = TRUE)]
dataCorr[, sampleFeatureNames] <- augmentedResultsNormalised[match(dataCorr$instances, augmentedResultsNormalised$instances), sampleFeatureNames]
# randomSubset <- dataCorr[sample(1:nrow(dataCorr), nrow(dataCorr)/100, replace = FALSE), ]
dir.create(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrDataFeatureSelection/"))
write.csv(dataCorr, paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrDataFeatureSelection/metadata.csv"), quote = FALSE, row.names = FALSE)


# Should now be able to extract data and plot it; look for instances that are similar,
# plot them, and see why performance is different.
projection <- read.table("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrDataFeatureSelection/projection_matrix.csv", header = TRUE, sep = ",")
z1Vals <- projection[projection$Row == "Z_{1}", ]
z1Vals <- z1Vals[, -1]
z2Vals <- projection[projection$Row == "Z_{2}", ]
z2Vals <- z2Vals[, -1]

vals <- dataCorr
vals$z1 <- 0
vals$z2 <- 0
for(name in colnames(z1Vals)){
  realName <- paste0("feature_", name)
  vals$z1 <- vals$z1 + z1Vals[1, name] * vals[, realName]
  vals$z2 <- vals$z2 + z2Vals[1, name] * vals[, realName]
}

xVals <- c(min(vals$z1), max(vals$z1))
yVals <- c(min(vals$z2), max(vals$z2))

vals <- vals[order(vals$algo_kriging), ]

ggplot(data = vals, aes(x = z1, y = z2, col = algo_kriging)) +
  # geom_point(shape = ".") +
  geom_point(size = 2) +
  # geom_point() +
  scale_x_continuous(limits = xVals) +
  scale_y_continuous(limits = yVals) +
  # scale_colour_gradientn(colours = colourGradient) +
  theme_bw() +
  ggtitle("Y vs. X with Z-Color")



vals <- vals[vals$z1 > -0.5 & 
               vals$z1 < 0 &
               vals$z2 > -1.25 &
               vals$z2 < -0.75, ]










# Here want to save it and run instance space analysis

# To do so, use the already normalised data
sampleFeatureNames <- colnames(augmentedResultsNormalised[, str_which(colnames(augmentedResultsNormalised), "feature_sample", negate = FALSE)])
# Drop all feature information, then extract sample feature information
dataError <- dataError[, str_which(colnames(dataError), "feature_", negate = TRUE)]
dataError[, sampleFeatureNames] <- augmentedResultsNormalised[match(dataError$instances, augmentedResultsNormalised$instances), sampleFeatureNames]

dataErrorRel <- dataErrorRel[, str_which(colnames(dataErrorRel), "feature_", negate = TRUE)]
dataErrorRel[, sampleFeatureNames] <- augmentedResultsNormalised[match(dataErrorRel$instances, augmentedResultsNormalised$instances), sampleFeatureNames]

dataCorr <- dataCorr[, str_which(colnames(dataCorr), "feature_", negate = TRUE)]
dataCorr[, sampleFeatureNames] <- augmentedResultsNormalised[match(dataCorr$instances, augmentedResultsNormalised$instances), sampleFeatureNames]


# randomSubset <- dataError[sample(1:nrow(dataError), nrow(dataError)/100, replace = FALSE), ]
dir.create(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/errorDataFeatureSelection/"))
write.csv(dataError, paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/errorDataFeatureSelection/metadata.csv"), quote = FALSE, row.names = FALSE)

dir.create(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/errorCompDataFeatureSelection/"))
write.csv(dataErrorRel, paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/errorCompDataFeatureSelection/metadata.csv"), quote = FALSE, row.names = FALSE)

dataCorr <- dataCorr[dataCorr$feature_sample_dimension == 1, ]
dataCorr <- dataCorr[, str_which(colnames(dataCorr), "feature_", negate = TRUE)]
dataCorr[, sampleFeatureNames] <- augmentedResultsNormalised[match(dataCorr$instances, augmentedResultsNormalised$instances), sampleFeatureNames]
# randomSubset <- dataCorr[sample(1:nrow(dataCorr), nrow(dataCorr)/100, replace = FALSE), ]
dir.create(paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrDataFeatureSelection/"))
write.csv(dataCorr, paste0("C:/Users/nandr/Documents/InstanceSpace-master/trial/corrDataFeatureSelection/metadata.csv"), quote = FALSE, row.names = FALSE)



for(instance in unique(dataError$instances)){
  if(nrow(dataError[dataError$instances == instance, ]) > 1){
    print(instance)
  }
}

test <- dataError[dataError$instances == "ShiGramacyLee_H25_L25_S1", ]

test <- rawResults[rawResults$instance == "ShiGramacyLee" & rawResults$highFiBudget == 25 & rawResults$lowFiBudget == 25 & rawResults$seed == 1, ]







# Here have a previous step where we are going to try to get rid of some instances using instance filtering
# Combine the features with 0,1 performance
dataPerf <- data.frame(instances = augmentedResults$instances)
# Start with correlation
dataPerf$algo_cokriging <- as.numeric(augmentedResults$cokriging_modelCorrelation >= 0.99*augmentedResults$kriging_modelCorrelation)
dataPerf$algo_kriging <- as.numeric(augmentedResults$kriging_modelCorrelation >= 0.99*augmentedResults$cokriging_modelCorrelation)

sampleFeat <- colnames(augmentedResultsNormalised[, str_which(colnames(augmentedResultsNormalised), "feature_sample")])
name <- "correlationPerformance"
dataPerf[, sampleFeat] <- augmentedResultsNormalised[match(dataPerf$instances, augmentedResultsNormalised$instances), sampleFeat]
write.csv(dataPerf, paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/data/", name, ".csv"), quote = FALSE, row.names = FALSE)


epsilonValues <- seq(5, 0.25, -0.25)
originalUniformValues <- c()
originalVisaRatio <- c()
originalRatioInstances <- c()

for(eps in epsilonValues){
  print(eps)
  originalUniformValues <- c(originalUniformValues, read.table(paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/cvnnd/CVNND_Pur_Ftr&Good_", name, "_Dist_", sprintf(eps, fmt = '%#.3f'), ".csv"), header = TRUE, sep = ",")$Uniformity_All)
  originalVisaRatio <- c(originalVisaRatio, nrow(read.table(paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/visa/ViSA_Ftr&Good_", name, "_Dist_", sprintf(eps, fmt = '%#.3f'), ".csv"), header = TRUE, sep = ",")) / nrow(read.table(paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/purified/Pur_Ftr&Good_", name, "_Dist_", sprintf(eps, fmt = '%#.3f'), ".csv"), header = TRUE, sep = ",")))
  originalRatioInstances <- c(originalRatioInstances, nrow(read.table(paste0("C:/Users/nandr/Documents/repos/InstanceFilteringNico/purified/Pur_Ftr&Good_", name, "_Dist_", sprintf(eps, fmt = '%#.3f'), ".csv"), header = TRUE, sep = ",")) / nrow(dataPerf))
}
plot(c(), c(), ylim = c(0, 1), xlim = c(0, 3))
lines(epsilonValues, originalUniformValues, col="red")
lines(epsilonValues, originalVisaRatio, col="green")
lines(epsilonValues, originalRatioInstances, col="blue")
lines(c(2, 2), c(0,3))




# Here start work on the pupeline suggested by andres
# Should have normalised the feature data, read it in and play with it
augmentedResultsNormalised <- read.table("data/combinedData/normalisedFeatures.txt", header = TRUE, sep = ",")

# First will try to use pca on all features, but work only on the ones that are sample based
featureSampleData <- augmentedResultsNormalised[, str_which(colnames(augmentedResultsNormalised), "feature_sample", negate = FALSE)]

pcaData <- prcomp(featureSampleData, center = TRUE,scale. = TRUE)
summary(pcaData)
plot(cumsum(pcaData$sdev^2 / sum(pcaData$sdev^2)), type="b")
sum(pcaData$sdev^2 > 1)
relevantPcaData <- pcaData$x[, 1:26]

tsne(relevantPcaData[])


colnames(relevantPcaDataHigh) <- paste0("feature_high", colnames(relevantPcaDataHigh))








pcaDataLow <- prcomp(featureLowData, center = TRUE,scale. = TRUE)
summary(pcaDataLow)
plot(cumsum(pcaDataLow$sdev^2 / sum(pcaDataLow$sdev^2)), type="b")
sum(pcaDataLow$sdev^2 > 1)
relevantPcaDataLow <- pcaDataLow$x[, 1:6]
colnames(relevantPcaDataLow) <- paste0("feature_low", colnames(relevantPcaDataLow))

pcaDataMid <- prcomp(featureMidData, center = TRUE,scale. = TRUE)
summary(pcaDataMid)
plot(cumsum(pcaDataMid$sdev^2 / sum(pcaDataMid$sdev^2)), type="b")
sum(pcaDataMid$sdev^2 > 1)
relevantPcaDataMid <- pcaDataMid$x[, 1:5]
colnames(relevantPcaDataMid) <- paste0("feature_mid", colnames(relevantPcaDataMid))

ispSubset <- ispSubset[, str_which(colnames(ispSubset), "feature_", negate = TRUE)]
ispSubset[, colnames(relevantPcaDataHigh)] <- relevantPcaDataHigh
ispSubset[, colnames(relevantPcaDataLow)] <- relevantPcaDataLow
ispSubset[, colnames(relevantPcaDataMid)] <- relevantPcaDataMid







chosenInstances <- unique(sapply(strsplit(cleanedData$instances, "_"), "[[", 1))







testSaved <- combineArrayResults("smallerExpandedTesting", 1, 890, 10)
testOriginal <- read.table("data/combinedData/smallerExpandedTestingCoKrig.txt", header = TRUE, sep = " ")


test <- testSaved

testOriginal <- testOriginal[str_which(testOriginal$instance, "COCOfunction16", negate = TRUE), ]
test <- test[test$seed <= 5, ]


test$combined <- paste0(test$instance, "_",
                        test$highFiBudget, "_",
                        test$lowFiBudget, "_",
                        test$seed, "_",
                        test$technique)

testOriginal$combined <- paste0(testOriginal$instance, "_",
                                testOriginal$highFiBudget, "_",
                                testOriginal$lowFiBudget, "_",
                                testOriginal$seed, "_",
                                testOriginal$technique)


order <- match(testOriginal$combined, test$combined)
test <- test[order, ]
order <- match(testOriginal$combined, test$combined)
test$combined2 <- testOriginal[order, "combined"]
test$oldVal <- testOriginal[order, "modelCorrelation"]
test$diff <- test$modelCorrelation - test$oldVal

test <- test[, c("combined", "combined2", "highFiBudget", "technique", 'modelCorrelation', 'oldVal', 'diff')]

test <- test[!is.na(test$combined), ]

sum(test$diff > 0.5)
sum(test$diff < -0.5)


for(i in 1:nrow(test)){
  if(test[i, "instance"] != testOriginal[i, "instance"] | 
     test[i, "seed"] != testOriginal[i, "seed"] | 
     test[i, "highFiBudget"] != testOriginal[i, "highFiBudget"] | 
     test[i, "lowFiBudget"] != testOriginal[i, "lowFiBudget"]){
    print(i)
    break
  }
}



test <- test[1:nrow(testOriginal), ]

test$diff <- test$modelCorrelation - testOriginal$modelCorrelation



test1 <- read.table("data/combinedData/smallerExpandedTestingCoKrig.txt", header = TRUE, sep = " ")
test1 <- test1[test1$lowFiBudget <= 100, c("instance", "technique", "highFiBudget", "lowFiBudget", "seed", "modelError", "modelCorrelation", "time")]
test1 <- test1[str_which(test1$instance, "COCOfunction16", negate = TRUE), ]
test2 <- read.table("data/combinedData/smallerExpandedTestingKrigWithoutTieBreak.txt", header = TRUE, sep = " ")
test2 <- test2[test2$lowFiBudget <= 100, c("instance", "technique", "highFiBudget", "lowFiBudget", "seed", "modelError", "modelCorrelation", "time")]
test2 <- test2[str_which(test2$instance, "COCOfunction16", negate = TRUE), ]

test <- rbind(test1, test2)

test$instances <- paste0(test$instance, "_H", test$highFiBudget,
                         "_L", test$lowFiBudget,
                         "_S", test$seed)

augmentedResultsReducedOld <- augmentData(test, limitedInformationFeatures, perfectInformationFeatures, c("kriging", "cokriging"))
augmentedResultsReducedOld <- augmentedResultsReducedOld[!is.na(augmentedResultsReducedOld$cokriging_modelCorrelation), ]


augmentedResultsReducedOld$labelled_algo_kriging <- as.numeric(augmentedResultsReducedOld$kriging_modelCorrelation >= 0.99*augmentedResultsReducedOld$cokriging_modelCorrelation)
augmentedResultsReducedOld$labelled_algo_cokriging <- as.numeric(augmentedResultsReducedOld$cokriging_modelCorrelation >= 0.99*augmentedResultsReducedOld$kriging_modelCorrelation)





temp <- condensedData
temp$diff <- temp$feature_real_CC - temp$feature_real_LCCrel_0_5




test <- condensedData[condensedData$feature_sample_highFiBudgetRatio == 2, c(1:30, 135)]
colnames(test)

sum(condensedData$kriging_corrWilcoxon0.001 > 0.5)
sum(condensedData$cokriging_corrWilcoxon0.001 > 0.5)
sum(condensedData$kriging_errorWilcoxon0.001 > 0.5)
sum(condensedData$cokriging_errorWilcoxon0.001 > 0.5)

