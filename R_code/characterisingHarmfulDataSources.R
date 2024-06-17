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

# The code that follows analyses the data presented in the paper
# "Characterising harmful data sources when constructing multi-fidelity surrogate models",
# It is divided into two sections, each of which can be rerun
# for reproducibility purposes:

# - The first combines the files from the cluster into a single file. It combines 
#   the algorithm performance and feature values into a single data frame,
#   as well as condensing multiple seeds for the same instance by using the wilcoxon
#   test. The resulting files can be found in data/isaMetadata, so new users 
#   should not require to run this.

# - The second filters the instances using all the features,
#   uses the automatic feature selection of Instance Space Analysis to
#   choose a set of features, uses the chosen features to filter the original
#   instances, and finally constructs the instance space.

source("R_code/dataProcessor.R")


# # COMBINE FEATURE VALUES AND MODEL PERFORMANCE INTO A SINGLE FILE
# # Get the information from the runs, and the sample and "perfect" features, then combine into a single data set
# rawResultsSynth <- combineArrayResults("experimentalRunSurrogateModelWithFixedSample", 1, 5780, 10)
# rawResultsSOLAR <- combineArrayResults("experimentalRunSurrogateModelWithFixedSampleSOLAR", 1, 2160)
# rawResults <- rbind(rawResultsSOLAR, rawResultsSynth)
# features <- read.table("data/features/sampleAndRealFeaturesCleanStandarised.txt", header = TRUE, sep = " ")
# 
# rawResults$instances <- paste0("(", rawResults$functionName, 
#                                ",", rawResults$nh,
#                                ",", rawResults$nl,
#                                ",", rawResults$seed, ")")
# 
# augmentedResults <- augmentData(rawResults, features, c("kriging", "cokriging"))
# write.table(augmentedResults, "data/isaMetadata/surrogateModelWithFixedSampleAggregatedData.txt", quote = FALSE, row.names = FALSE)
# 
# augmentedResults <- read.table("data/isaMetadata/surrogateModelWithFixedSampleAggregatedData.txt", header = TRUE, sep = " ")
# condensedData <- condenseData(augmentedResults, c("kriging", "cokriging"))
# # Sadly cannot use "," for instance names, need to change it to underscore
# condensedData$instances <- gsub(',', '_', condensedData$instances)
# write.table(condensedData, "data/isaMetadata/surrogateModelWithFixedSampleMetadata.txt", quote = FALSE, row.names = FALSE)
# # COMPLETED COMBINING FEATURE FILES


# INSTANCE FILTERING, FEATURE SELECTION AND SPACE CONSTRUCTION
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
    labels[i, "minDist"] <- min(temp)
  }
  CV = sd(labels$minDist)/mean(labels$minDist)
  Uniformity = 1 - CV
  cat(paste0("\rWorking on row ", i, "/", nrow(labels), " - done\n"))
  return(c(CV, Uniformity))
}

instancePurificationProcedure <- function(epsilons, instancesData, name, Ybin, goodThreshold = 0.05){
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
    # Have the data, first purify
    labels <- purifyInstances(distMatrix, Ybin, eps)
    results[numResults, "VisaRatio"] <- sum(labels$visa) / sum(!labels$preclude)
    results[numResults, "InstancesRatio"] <- sum(!labels$preclude) / nrow(labels)
    # Now calculate uniformity
    subsetDistMatrix <- distMatrix[labels$dissimlar,labels$dissimlar]
    vals <- calculateCVNND(subsetDistMatrix)
    results[numResults, c("CV", "UniformVal")] <- vals
    if(0 %in% results$epsilon){
      maxCV <- results[results$epsilon == 0, "CV"]
      results$UniformValStandarised <- (results$UniformVal - (1 - maxCV)) / (1 - (1 - maxCV))
    }
    print(results)
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

# First order instances. Randomise order but make sure first we have SOLAR instances,
# then literature instances, then disturbance based instances. This ensures the
# order is the priority given to instances when filtering.
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

# Do an initial filtering using sample features i.e. features calculated using the
# same sample used to train the models.
tempInstances <- instancesOrdered[c("instances", "Source", colnames(condensedData[str_which(colnames(condensedData), "feature_sample_")]))]
# As specified in the paper, only use MMCE for lda (and not qda)
tempInstances <- tempInstances[str_which(colnames(tempInstances), "mmce_qda", negate = TRUE)]
tempInstances$feature_sample_dimension <- instancesOrdered$feature_real_dimension
tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.001Bad
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.001Bad
# kriging_corrWilcoxon0.001Bad gives the p-value of a one-tailed wilcoxon test
# of the probability of having observed the model accuracies assuming
# Kriging is worse (by a margin of 0.001) than Co-Kriging.
name <- "CorrTol0.001Bad"

# Need to provide the binary performance for the filtering
Ybin <- tempInstances[, c('algo_kriging', 'algo_cokriging')]
Ybin$algo_kriging <- as.numeric(tempInstances$algo_kriging <= 0.05)
Ybin$algo_cokriging <- as.numeric(tempInstances$algo_cokriging <= 0.05)
# For instances where no algorithm is good, need to decide based on mean
Ybin[tempInstances$algo_cokriging > 0.05 &
       tempInstances$algo_kriging > 0.05 &
       tempInstances$algo_kriging < tempInstances$algo_cokriging, "algo_kriging"] <- TRUE
Ybin[tempInstances$algo_cokriging > 0.05 &
       tempInstances$algo_kriging > 0.05 &
       tempInstances$algo_kriging >= tempInstances$algo_cokriging, "algo_cokriging"] <- TRUE

Ybin <- data.matrix(Ybin)

instancePurificationProcedure(0:5, tempInstances, paste0("features", name), Ybin)
results <- read.table(paste0("data/isaMetadata/instancePurificationfeatures", name, ".txt"), header = TRUE, sep = " ")
results <- results[results$UniformValStandarised >= 0.5, ]
# Don't really need to plot, just take the smallest epsilon for which the standarised uniform value is above 0.5
eps <- min(results$epsilon)
labels <- read.table(paste0("data/isaMetadata/instanceFiltering/features", name, "eps", eps, ".txt"), header = TRUE, sep = " ")
filteredData <- tempInstances[!labels$preclude, ]
# Standarise features to have mean 0 and sd 1 (helps with the feature selection)
for(feat in colnames(filteredData[str_which(colnames(filteredData), "feature")])){
  filteredData[feat] <- (filteredData[feat] - mean(unlist(filteredData[feat])))/ sd(unlist(filteredData[feat]))
}
corrs <- findCorrs(filteredData)

write.table(filteredData, "matlab_code/featureSelection/metadata.csv", sep=",", quote = FALSE, row.names = FALSE)
# system("matlab -nodisplay -r \"run('matlab_code/featureSelection/featureSelection.m'); exit\"")
# Option 1: Take the features with the lowest error, having looked using a search (heuristic)
# chosenFeatures <- paste0("feature_", read.table("matlab_code/featureSelection/features.txt", header = FALSE, sep = ",")[1,])

# Option 2: This is the best known solution with a value of 0.2715742
chosenFeatures <- c("feature_sample_LCCrel_0_5",
                    "feature_sample_high_ela_level_mmce_lda_50",
                    "feature_sample_RRMSE",
                    "feature_sample_high_ic_m0",
                    "feature_sample_budgetRatio",
                    "feature_sample_LCCrel_0_9",
                    "feature_sample_mid_nbc_nb_fitness_cor",
                    "feature_sample_mid_ela_meta_lin_w_interact_adj_r2",
                    "feature_sample_LCC_sd")


# Ok so now have all of the chosen features, need to re run the filtering with the selected features
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures)
instancePurificationProcedure(c(0, 0.1, 0.2, 0.3, 0.4, 0.5), tempInstances[colnames], paste0("instanceFiltering", name), Ybin)
results <- read.table(paste0("data/isaMetadata/instancePurificationinstanceFiltering", name, ".txt"), header = TRUE, sep = " ")
# Don't really need to plot, just take the smallest epsilon for which the standarised uniform value is above 0.5
results <- results[results$UniformValStandarised >= 0.5, ]
eps <- min(results$epsilon)
labels <- read.table(paste0("data/isaMetadata/instanceFiltering/instanceFiltering", name, "eps", eps, ".txt"), header = TRUE, sep = " ")
filteredData <- tempInstances[!labels$preclude, colnames]
dir.create(paste0("matlab_code/ISA/", name, "/"))
write.csv(filteredData, paste0("matlab_code/ISA/", name, "/metadata.csv"), quote = FALSE, row.names = FALSE)
write.csv(Ybin[!labels$preclude, ], paste0("matlab_code/ISA/", name, "/algorithm_binary_user.csv"), quote = FALSE, row.names = FALSE)

temp <- Ybin[!labels$preclude, ]
sum(temp[, 1] == 1)/nrow(temp)
sum(temp[, 2] == 1)/nrow(temp)

# Now run the ISA code. Note that if you changed the folder name
# you will need to edit the first line of the file runISA.m
# system("matlab -nodisplay -r \"run('matlab_code/ISA/runISA.m'); exit\"")
# The file extraPlots.m can also be run to recreate the plots used in the 
# paper. Some editing will be required for new instance spaces.
# The final instance space with the metadata, model, plots, etc... 
# has been saved to data/characterisingHarmfulDataSourcesISA/ to make sure 
# it is not overriden.







# Testing happens now, calculate precision and whatnot
name <- "CorrTol0.001Bad"
algorithmOriginalBin <- read.table(paste0("matlab_code/ISA/", name, "/algorithm_binary_user.csv"), header = TRUE, sep = ",")
algorithmRaw <- read.table(paste0("matlab_code/ISA/", name, "/algorithm_raw.csv"), header = TRUE, sep = ",")
algorithmBin <- read.table(paste0("matlab_code/ISA/", name, "/algorithm_bin.csv"), header = TRUE, sep = ",")
predictions <- read.table(paste0("matlab_code/ISA/", name, "/algorithm_svm.csv"), header = TRUE, sep = ",")
portfolio <- read.table(paste0("matlab_code/ISA/", name, "/portfolio.csv"), header = TRUE, sep = ",")
portfolioPrediction <- read.table(paste0("matlab_code/ISA/", name, "/portfolio_svm.csv"), header = TRUE, sep = ",")
# Kriging
print(paste0("Kriging, probability of good ", 100*sum(algorithmBin$kriging) / nrow(algorithmBin), "%"))
print(paste0("Kriging SVM, probability of good ", 100*sum(algorithmBin$kriging == predictions$kriging) / nrow(algorithmBin), "%"))
truePositives <- sum(algorithmBin$kriging == predictions$kriging & algorithmBin$kriging == 1)
trueNegatives <- sum(algorithmBin$kriging == predictions$kriging & algorithmBin$kriging == 0)
falsePositives <- sum(algorithmBin$kriging != predictions$kriging & algorithmBin$kriging == 0)
falseNegatives <- sum(algorithmBin$kriging != predictions$kriging & algorithmBin$kriging == 1)
print(paste0("Kriging SVM precision ",  truePositives / (truePositives + falsePositives)))
print(paste0("Kriging SVM 'inverse' precision ",  trueNegatives / (trueNegatives + falseNegatives)))
print(paste0("Kriging SVM recall ",  truePositives / (truePositives + falseNegatives)))
print(paste0("Kriging SVM 'inverse' recall ",  trueNegatives / (trueNegatives + falsePositives)))

# CoKriging
print(paste0("CoKriging, probability of good ", 100*sum(algorithmBin$cokriging) / nrow(algorithmBin), "%"))
print(paste0("CoKriging SVM, probability of good ", 100*sum(algorithmBin$cokriging == predictions$cokriging) / nrow(algorithmBin), "%"))
truePositives <- sum(algorithmBin$cokriging == predictions$cokriging & algorithmBin$cokriging == 1)
trueNegatives <- sum(algorithmBin$cokriging == predictions$cokriging & algorithmBin$cokriging == 0)
falsePositives <- sum(algorithmBin$cokriging != predictions$cokriging & algorithmBin$cokriging == 0)
falseNegatives <- sum(algorithmBin$cokriging != predictions$cokriging & algorithmBin$cokriging == 1)
print(paste0("CoKriging SVM precision ",  truePositives / (truePositives + falsePositives)))
print(paste0("CoKriging SVM 'inverse' precision ",  trueNegatives / (trueNegatives + falseNegatives)))
print(paste0("CoKriging SVM recall ",  truePositives / (truePositives + falseNegatives)))
print(paste0("CoKriging SVM 'inverse' recall ",  trueNegatives / (trueNegatives + falsePositives)))

# Accuracy of selector, where when the portfolio chooses none, pick Kriging
(sum(algorithmBin[portfolioPrediction$Best_Algorithm == 1, "kriging"] == 1) +
    sum(algorithmBin[portfolioPrediction$Best_Algorithm == 2, "cokriging"] == 1) +
    sum(algorithmBin[portfolioPrediction$Best_Algorithm == 0, "kriging"] == 1))/nrow(portfolio)


# Finally, in order to derive simple rules, create some plots which analyse when
# each of the techniques should be used.
# Get the "real" feature value
name <- "CorrTol0.001Bad"
preProcessedFeatures <- read.table("data/features/sampleAndRealFeaturesClean.txt", header = TRUE, sep = " ")
preProcessedFeatures$instanceName <- paste0(gsub('[(]', '', sapply(strsplit(preProcessedFeatures$instances, ","), "[[", 1)), "_",
                                            sapply(strsplit(preProcessedFeatures$instances, ","), "[[", 2), "_",
                                            sapply(strsplit(preProcessedFeatures$instances, ","), "[[", 3))

projection <- read.table(paste0("matlab_code/ISA/", name, "/projection_matrix.csv"), header = TRUE, sep = ",")

featuresOfInterest <- paste0("feature_", colnames(projection[-1]))
# Add highFiBudgetRatio and CC 
featuresOfInterest <- c(featuresOfInterest, 
                        "feature_sample_highFiBudgetRatio", 
                        "feature_sample_CC")

dataOfInterest <- tempInstances[!labels$preclude, c("instances", "Source", "algo_kriging", "algo_cokriging", featuresOfInterest)]
dataOfInterest$algo_kriging_binary <- Ybin[!labels$preclude,1]
dataOfInterest$algo_cokriging_binary <- Ybin[!labels$preclude,2]
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
# For adjusted R^2, whenever the feature is below 0 simply set it to 0
dataOfInterest[dataOfInterest$original_feature_sample_mid_ela_meta_lin_w_interact_adj_r2 < 0, 'original_feature_sample_mid_ela_meta_lin_w_interact_adj_r2'] <- 0
write.table(dataOfInterest, paste0("matlab_code/ISA/", name, "/metadataAllFeatures.csv"), sep=",", quote = FALSE, row.names = FALSE)
feats <- paste0("original_", featuresOfInterest)
# Now get data into a format that ggplot will use
# For each of the intervals, find the proportion of instances for which
# Kriging / Co-Kriging are good
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
    
    
    numKrigGood <- nrow(dataOfInterest[dataOfInterest$algo_kriging_binary == 1 &
                                         dataOfInterest[featName] >= (xStart - 0.001) &
                                         dataOfInterest[featName] <= (xEnd + 0.001), ])
    numCoKrigGood <- nrow(filteredData[dataOfInterest$algo_cokriging_binary == 1 &
                                         dataOfInterest[featName] >= (xStart - 0.001) &
                                         dataOfInterest[featName] <= (xEnd + 0.001), ])
    
    numKrigBad <- nrow(dataOfInterest[dataOfInterest$algo_kriging_binary == 0 &
                                        dataOfInterest[featName] >= (xStart - 0.001) &
                                        dataOfInterest[featName] <= (xEnd + 0.001), ])
    numCoKrigBad <- nrow(filteredData[dataOfInterest$algo_cokriging_binary == 0 &
                                        dataOfInterest[featName] >= (xStart - 0.001) &
                                        dataOfInterest[featName] <= (xEnd + 0.001), ])
    
    proportionKrigGood <- c(proportionKrigGood, numKrigGood / (numKrigGood + numKrigBad), numKrigGood / (numKrigGood + numKrigBad))
    proportionCoKrigGood <- c(proportionCoKrigGood, numCoKrigGood / (numCoKrigGood + numCoKrigBad), numCoKrigGood / (numCoKrigGood + numCoKrigBad))
    proportionKrigBad <- c(proportionKrigBad, numKrigBad / (numKrigGood + numKrigBad), numKrigBad / (numKrigGood + numKrigBad))
    proportionCoKrigBad <- c(proportionCoKrigBad, numCoKrigBad / (numCoKrigGood + numCoKrigBad), numCoKrigBad / (numCoKrigGood + numCoKrigBad))
    
  }
  currIndex <- nrow(plottingData)
  range <- (currIndex + 1) : (currIndex + 2*length(proportionKrigGood))
  plottingData[range, "proportion"] <- c(proportionKrigGood, proportionCoKrigGood)
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
ggsave(paste0("matlab_code/ISA/", name ,"/custom_proportionGoodBudgetRatio.png"), width = 3.5, height = 3.5)

ggplot(plottingData[plottingData$feature == "original_feature_sample_highFiBudgetRatio", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(B[h]^{r}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave(paste0("matlab_code/ISA/", name ,"/custom_proportionGoodHighFiBudgetRatio.png"), width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_LCCrel_0_5", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(LCC[0.5]^{0.2^{1/d}}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave(paste0("matlab_code/ISA/", name ,"/custom_proportionGoodLCC5.png"), width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_LCCrel_0_9", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(LCC[0.9]^{0.2^{1/d}}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.2))
ggsave(paste0("matlab_code/ISA/", name ,"/custom_proportionGoodLCC9.png"), width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_LCC_sd", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(LCC[coeff]^{0.2}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave(paste0("matlab_code/ISA/", name ,"/custom_proportionGoodLCCsd.png"), width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_mid_ela_meta_lin_w_interact_adj_r2", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(R[LI]))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave(paste0("matlab_code/ISA/", name ,"/custom_proportionGoodR2LI.png"), width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_mid_nbc_nb_fitness_cor", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(R[L]))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave(paste0("matlab_code/ISA/", name ,"/custom_proportionGoodNBC.png"), width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_CC", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(CC))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave(paste0("matlab_code/ISA/", name ,"/custom_proportionGoodCC.png"), width = 4, height = 4)


ggplot(plottingData[plottingData$feature == "original_feature_sample_high_ela_level_mmce_lda_50", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(MMCE[lda]^{0.5}))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave(paste0("matlab_code/ISA/", name ,"/custom_proportionGoodMMCElda.png"), width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_RRMSE", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(RRMSE))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave(paste0("matlab_code/ISA/", name ,"/custom_proportionGoodRRMSE.png"), width = 4, height = 4)

ggplot(plottingData[plottingData$feature == "original_feature_sample_high_ic_m0", ], aes(x=xVals, y=proportion, group=Method)) +
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method))+
  ylim(0,1)+
  xlab(expression(M[0]))+
  ylab("Proportion of instances with good performance") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
ggsave(paste0("matlab_code/ISA/", name ,"/custom_proportionGoodM0.png"), width = 4, height = 4)




# Get some statistics on how accurate the guidelines are
test <- dataOfInterest
# First check a rule based only on CC
vals <- seq(0,1,0.05)
accuracies <- c()
for(val in vals){
  chooseKrig <- test[test$original_feature_sample_CC < val, ]
  chooseCoKrig <- test[!(test$original_feature_sample_CC < val), ]
  accuracies <- c(accuracies, (sum(chooseKrig$algo_kriging_binary == 1) + sum(chooseCoKrig$algo_cokriging_binary == 1)) / nrow(test))
}
plot(vals, accuracies)
max(accuracies)

test <- dataOfInterest
# Rule based on high fi budget ratio
chooseKrig <- test[test$original_feature_sample_highFiBudgetRatio >= 20, ]
undecided <- test[test$original_feature_sample_highFiBudgetRatio < 20, ]

# Rule based on budget ratio
chooseKrig <- rbind(chooseKrig, undecided[undecided$original_feature_sample_budgetRatio >= 1, ])
undecided <- undecided[undecided$original_feature_sample_budgetRatio < 1, ]

# Rule based on LCC mean
chooseKrig <- rbind(chooseKrig, undecided[undecided$original_feature_sample_LCCrel_0_5 <= 0.7, ])
undecided <- undecided[undecided$original_feature_sample_LCCrel_0_5 > 0.7, ]

# Rule based on LCC 0.9
chooseCoKrig <- undecided[undecided$original_feature_sample_LCCrel_0_9 >= 0.8 | 
                            (undecided$original_feature_sample_mid_ela_meta_lin_w_interact_adj_r2 >= 0.6),]

undecided <- undecided[undecided$original_feature_sample_LCCrel_0_9 < 0.8 & 
                         (undecided$original_feature_sample_mid_ela_meta_lin_w_interact_adj_r2 < 0.6),]


# Sanity check that these are all distinct
intersect(chooseKrig$instances, undecided$instances)
intersect(chooseKrig$instances, chooseCoKrig$instances)
intersect(chooseCoKrig$instances, undecided$instances)
# Check how many instances went into each 
nrow(chooseKrig)
nrow(chooseCoKrig)
nrow(undecided)

# Accuracies of the first two classes
sum(chooseKrig$algo_kriging_binary == 1) / nrow(chooseKrig)
sum(chooseCoKrig$algo_cokriging_binary == 1) / nrow(chooseCoKrig)


# Now check what the remainder of the set looks like
sum(undecided$algo_kriging_binary == 1)
sum(undecided$algo_cokriging_binary == 1)


# Pretty evenly divided, separate based on budget
vals <- seq(0,1,0.05)
accuracies <- c()
for(val in vals){
  subsetChooseKrig <- undecided[undecided$original_feature_sample_budgetRatio >= val, ]
  subsetChooseCoKrig <- undecided[!(undecided$original_feature_sample_budgetRatio >= val), ]
  accuracies <- c(accuracies, (sum(subsetChooseKrig$algo_kriging_binary == 1) + sum(subsetChooseCoKrig$algo_cokriging_binary == 1)) / nrow(undecided))
}
plot(vals, accuracies)
max(accuracies)

vals <- seq(0,20,1)
accuracies <- c()
for(val in vals){
  subsetChooseKrig <- undecided[undecided$original_feature_sample_highFiBudgetRatio > val, ]
  subsetChooseCoKrig <- undecided[!(undecided$original_feature_sample_highFiBudgetRatio > val), ]
  accuracies <- c(accuracies, (sum(subsetChooseKrig$algo_kriging_binary == 1) + sum(subsetChooseCoKrig$algo_cokriging_binary == 1)) / nrow(undecided))
}
plot(vals, accuracies)
max(accuracies)


# Seems better to split based on relative number of high fidelity samples, with a split of 6d
chooseKrig <- rbind(chooseKrig, undecided[undecided$original_feature_sample_highFiBudgetRatio > 6, ])
chooseCoKrig <- rbind(chooseCoKrig, undecided[undecided$original_feature_sample_highFiBudgetRatio <= 6, ])

nrow(chooseKrig)
nrow(chooseCoKrig)

(sum(chooseKrig$algo_kriging_binary == 1) + sum(chooseCoKrig$algo_cokriging_binary == 1)) / nrow(dataOfInterest)

# Add them to matlab to plot them
portfolioPrediction <- read.table(paste0("matlab_code/ISA/", name, "/portfolio_svm.csv"), header = TRUE, sep = ",")
# indexes <- !(algorithmBin$kriging == 0 & algorithmBin$cokriging == 0)
rulesBasedPredictor <- portfolioPrediction
rulesBasedPredictor[rulesBasedPredictor$Row %in% chooseKrig$instances, "Best_Algorithm"] <- 1
rulesBasedPredictor[rulesBasedPredictor$Row %in% chooseCoKrig$instances, "Best_Algorithm"] <- 2
rulesBasedPredictor[!(rulesBasedPredictor$Row %in% chooseCoKrig$instances) & 
                      !(rulesBasedPredictor$Row %in% chooseKrig$instances), "Best_Algorithm"] <- 0

nrow(chooseKrig) + nrow(chooseCoKrig)
write.table(rulesBasedPredictor[2], paste0("matlab_code/ISA/", name, "/rulesSelector.csv"), sep=",", quote = FALSE, row.names = FALSE)


