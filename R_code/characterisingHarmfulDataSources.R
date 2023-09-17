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
# condensedData <- condenseData(augmentedResults, c("kriging", "cokriging"))
# # Sadly cannot use , for instance names, need to change it to underscore
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

# Now choose what the performance is, and do the filtering!
tempInstancesSampleFeatures <- instancesOrdered[c("instances", "Source", colnames(condensedData[str_which(colnames(condensedData), "feature_sample_")]))]
tempInstancesSampleFeatures$feature_sample_dimension <- instancesOrdered$feature_real_dimension
tempInstances <- tempInstancesSampleFeatures
# Here the performance is a wilcoxon comparison that a model has a correlation
# with f_h to less than epsilon = 0.001 than other models. Other available
# epsilons include 0, 0.0025 and 0.005. It is also possible to use
# median performance, mean performance, and wilcoxon comparison of the error.
tempInstances$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.001
tempInstances$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.001
instancePurificationProcedure(5:0, tempInstances, "sampleFeaturesCorrTol0.001")

# Choose a filtered set
name <- "sampleFeaturesCorrTol0.001"
data <- tempInstancesSampleFeatures
# Need to add algorithm performance
data$algo_kriging <- instancesOrdered$kriging_corrWilcoxon0.001
data$algo_cokriging <- instancesOrdered$cokriging_corrWilcoxon0.001

results <- read.table(paste0("data/isaMetadata/instancePurification", name, ".txt"), header = TRUE, sep = " ")
# Don't really need to plot, just take the smallest epsilon for which the standarised uniform value is above 0.5
results <- results[results$UniformValStandarised >= 0.5, ]
eps <- min(results$epsilon)
labels <- read.table(paste0("data/isaMetadata/instanceFiltering/", name, "eps", eps, ".txt"), header = TRUE, sep = " ")
filteredData <- data[!labels$preclude, ]
# Standarise features to have mean 0 and sd 1 (helps with the feature selection)
for(feat in colnames(filteredData[str_which(colnames(filteredData), "feature")])){
  filteredData[feat] <- (filteredData[feat] - mean(unlist(filteredData[feat])))/ sd(unlist(filteredData[feat]))
}
write.csv(filteredData, "matlab_code/featureSelection/metadata.csv", quote = FALSE, row.names = FALSE)
system("matlab -nodisplay -r \"run('matlab_code/featureSelection/featureSelection.m'); exit\"")
chosenFeatures <- paste0("feature_", read.table("matlab_code/featureSelection/features.txt", header = FALSE, sep = ",")[1,])

# Ok so now have all of the chosen features, need to re run the filtering with the selected features
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures)
instancePurificationProcedure(c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), tempInstances[colnames], "instanceFilteringSampleFeaturesCorrTol0.001")


# This deals with all the features
isaName <- "sampleFeaturesCorrTol0.001"
name <- "instanceFilteringSampleFeaturesCorrTol0.001"
colnames <- c("instances", "Source", "algo_kriging", "algo_cokriging", chosenFeatures)
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
filteredData <- data[!labels$preclude, colnames]
dir.create(paste0("matlab_code/ISA/", isaName, "/"))
write.csv(filteredData, paste0("matlab_code/ISA/", isaName, "/metadata.csv"), quote = FALSE, row.names = FALSE)
# Now run the ISA code. Note that if you changed the folder name
# you will need to edit the first line of the file runISA.m
system("matlab -nodisplay -r \"run('matlab_code/ISA/runISA.m'); exit\"")
# The file extraPlots.m can also be run to recreate the plots used in the 
# paper. Some editing will be required for new instance spaces.
# The final instance space with the metadata, model, plots, etc... 
# has been saved to data/characterisingHarmfulDataSourcesISA/ to make sure 
# it is not overriden.


# Finally, in order to derive simple rules, create some plots which analyse when
# each of the techniques should be used.
# Get the "real" feature value
preProcessedFeatures <- read.table("data/features/sampleAndRealFeaturesClean.txt", header = TRUE, sep = " ")
preProcessedFeatures$instanceName <- paste0(gsub('[(]', '', sapply(strsplit(preProcessedFeatures$instances, ","), "[[", 1)), "_",
                                            sapply(strsplit(preProcessedFeatures$instances, ","), "[[", 2), "_",
                                            sapply(strsplit(preProcessedFeatures$instances, ","), "[[", 3))
                                            
projection <- read.table(paste0("matlab_code/ISA/sampleFeaturesCorrTol0.001/projection_matrix.csv"), header = TRUE, sep = ",")

featuresOfInterest <- paste0("feature_", colnames(projection[-1]))
# Add highFiBudgetRatio and CC 
featuresOfInterest <- c(featuresOfInterest, "feature_sample_highFiBudgetRatio", "feature_sample_CC")

dataOfInterest <- data[!labels$preclude, c("instances", "Source", "algo_kriging", "algo_cokriging", featuresOfInterest)]

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
dataOfInterest[dataOfInterest$original_feature_sample_mid_ela_meta_lin_simple_adj_r2 < 0, 'original_feature_sample_mid_ela_meta_lin_simple_adj_r2'] <- 0
dataOfInterest[dataOfInterest$original_feature_sample_mid_ela_meta_lin_w_interact_adj_r2 < 0, 'original_feature_sample_mid_ela_meta_lin_w_interact_adj_r2'] <- 0

# Choose the feats to plot
feats <- c("original_feature_sample_highFiBudgetRatio", 
           "original_feature_sample_budgetRatio",
           "original_feature_sample_LCCrel_0_4",
           "original_feature_sample_LCCrel_0_95",
           'original_feature_sample_mid_ela_meta_lin_w_interact_adj_r2',
           'original_feature_sample_mid_ela_meta_lin_simple_adj_r2',
           "original_feature_sample_RRMSE",
           "original_feature_sample_CC")
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




# Get some statistics on how accurate the guidelines are
test <- dataOfInterest
# First check a rule based only on CC
vals <- seq(0,1,0.05)
accuracies <- c()
for(val in vals){
  chooseKrig <- test[test$original_feature_sample_CC < val, ]
  chooseCoKrig <- test[!(test$original_feature_sample_CC < val), ]
  accuracies <- c(accuracies, (sum(chooseKrig$algo_kriging >= 0.5) + sum(chooseCoKrig$algo_cokriging >= 0.5)) / nrow(test))
}
plot(vals, accuracies)
max(accuracies)

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
# Check how many instances went into each 
nrow(chooseKrig)
nrow(chooseCoKrig)
nrow(undecided)

# Accuracies of the first two classes
sum(chooseKrig$algo_kriging >= 0.5) / nrow(chooseKrig)
sum(chooseCoKrig$algo_cokriging >= 0.5) / nrow(chooseCoKrig)


# Now check what the remainder of the set looks like
sum(undecided$algo_kriging >= 0.5)
sum(undecided$algo_cokriging >= 0.5)
# Pretty evenly divided, separate based on budget
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

(sum(chooseKrig$algo_kriging >= 0.5) + sum(chooseCoKrig$algo_cokriging >= 0.5)) / nrow(dataOfInterest)


