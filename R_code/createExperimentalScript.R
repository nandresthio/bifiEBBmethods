source("R_code/libraries.R")

# This function eases the creation of the experimental specifications for 
# experiments with the same single source, high fi, and low fi budgets, and
# the same seeds per run.
createScriptStructure <- function(problemTypes, methods, functions, highFiBudget, lowFiBudget, seedsStart = 1, seedsEnd = 20 , seedsPerRun = 20, printInfo = FALSE){
  jobs <- 0
  runData <- data.frame(matrix(ncol = 3, nrow = 0))
  for(problemType in problemTypes){
    for(method in methods){
      for(func in functions){
        for(seedStart in seq(seedsStart, seedsEnd, seedsPerRun)){
          jobs <- jobs + 1
          if(printInfo){print(paste0("Working on job ", jobs, " for krig"))}
          data <- c(problemType, paste0("(", func, "," , highFiBudget, ",", lowFiBudget, ",", paste0(seedStart, "-", min(seedStart + seedsPerRun - 1, seedsEnd)), ")"), method)
          runData[jobs, ] <- data
        }
      }
    }
  }
  colnames(runData) <- c("problemType", "problem", "method")
  if(printInfo){print(paste0("Created script with ", jobs, " jobs"))}
  return(runData)
}



# Create a script for the creation of all the sampling
print("CREATING SAMPLING SCRIPT.")
dim <- 0
for(func in c('LiuPedagogical', 'LiuBranin', 'RajnarayanHartmannH3', 'RajnarayanWoods', 
              'LiuStyblinskiTang', 'RajnarayanHartmannH6',
              'COCOfunction1-dim7-seed1-disth1-height0-radius0.025-freq2-amp0', 'ShiDettePepelyshev',
              'COCOfunction1-dim9-seed1-disth1-height0-radius0.025-freq2-amp0', 'LiuAckley10', 'LiuEllipsoid')){
  dim <- dim + 1
  if(dim == 11){dim = 20}
  print(func)
  print(dim)
  for(size in c(5, 10, 25, 50, 75, 100, 250, 500)){
    seedsPerRun <- 40
    if(dim >= 100){seedsPerRun <- 1}
    runData <- createScriptStructure(c("sampleCreation"), c("morrisMitchellLHS"), c(func), size, 0, seedsStart = 1, seedsEnd = 40, seedsPerRun = seedsPerRun)
    if(size != 5 || func != 'LiuPedagogical'){
      existingRunData <- rbind(existingRunData, runData)
    }else{
      existingRunData <- runData
    }
  }
  for(size in c(20, 30, 40, 60, 70, 80, 90)){
    seedsPerRun <- 20
    runData <- createScriptStructure(c("sampleCreation"), c("morrisMitchellLHS"), c(func), size, 0, seedsStart = 1, seedsEnd = 20, seedsPerRun = seedsPerRun)
    existingRunData <- rbind(existingRunData, runData)
  }
  for(size in seq(2, 20, by = 2)){
    if(size*dim %in% c(5, 10, 25, 50, 75, 100, 250, 500)){next}
    seedsPerRun <- 40
    if(size*dim >= 100){seedsPerRun <- 1}
    runData <- createScriptStructure(c("sampleCreation"), c("morrisMitchellLHS"), c(func), size * dim, 0, seedsStart = 1, seedsEnd = 40, seedsPerRun = seedsPerRun)
    # if(size != 5 || dim != 1){
    existingRunData <- rbind(existingRunData, runData)
  }  
}
write.table(existingRunData, "data/runScripts/sampleCreationScript.txt", quote = FALSE, row.names = FALSE)
print("COMPLETED SAMPLING SCRIPT.")

print("CREATING SUBSET SAMPLING SCRIPT.")
dim <- 0
for(func in c('LiuPedagogical', 'LiuBranin', 'RajnarayanHartmannH3', 'RajnarayanWoods', 
              'LiuStyblinskiTang', 'RajnarayanHartmannH6',
              'COCOfunction1-dim7-seed1-disth1-height0-radius0.025-freq2-amp0', 'ShiDettePepelyshev',
              'COCOfunction1-dim9-seed1-disth1-height0-radius0.025-freq2-amp0', 'LiuAckley10', 'LiuEllipsoid')){
  dim <- dim + 1
  if(dim == 11){dim = 20}
  print(func)
  print(dim)
  for(lowFi in c(5, 10, 25, 50, 75, 100, 250, 500)){
    for(highFi in c(5, 10, 25, 50, 75, 100, 250, 500)){
      if(highFi > lowFi){next}
      runData <- createScriptStructure(c("sampleCreation"), c("morrisMitchellLHS"), c(func), highFi, lowFi, seedsStart = 1, seedsEnd = 40, seedsPerRun = 40)
      if(lowFi != 5 || func != 'LiuPedagogical' || highFi != 5){
        existingRunData <- rbind(existingRunData, runData)
      }else{
        existingRunData <- runData
      }
    }
  }
  for(lowFi in c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)){
    for(highFi in c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)){
      if(highFi > lowFi){next}
      runData <- createScriptStructure(c("sampleCreation"), c("morrisMitchellLHS"), c(func), highFi, lowFi, seedsStart = 1, seedsEnd = 20, seedsPerRun = 20)
      existingRunData <- rbind(existingRunData, runData)
      
    }
  }
  for(sizeLow in seq(2, 20, by = 2)){
    for(sizeHigh in seq(2, 20, by = 2)){
      if(sizeHigh > sizeLow){next}
      if((sizeLow * dim) %in% c(5, 10, 25, 50, 75, 100, 250, 500) &
         (sizeHigh * dim) %in% c(5, 10, 25, 50, 75, 100, 250, 500)){next}
      
      runData <- createScriptStructure(c("sampleCreation"), c("morrisMitchellLHS"), c(func), sizeHigh * dim, sizeLow * dim, seedsStart = 1, seedsEnd = 40, seedsPerRun = 40)
      existingRunData <- rbind(existingRunData, runData)
    }
  }
}
write.table(existingRunData, "data/runScripts/subsetSampleCreationScript.txt", quote = FALSE, row.names = FALSE)
print("COMPLETED SUBSET SAMPLING SCRIPT.")








# Now want a script for the actual experimental setup
functions <- read.table("cpp_code/bifiEBBbenchmarks/data/availableFunctions/chosenTestSuiteN207.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
features <- read.table("data/features/featuresClean.txt", header = TRUE, sep = " ")
order <- match(functions, features$instances)
dims <- features[order, 'feature_dimension']
for(i in 1:length(functions)){
  func <- functions[[i]]
  dim <- dims[[i]]
  print(func)
  print(dim)
  for(lowFi in c(4, 8, 12, 16, 20)){
    for(highFi in c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)){
      if(highFi > lowFi){next}
      #print(paste0("Low ", lowFi, " high ", highFi))
      seedsStart <- 1
      seedsEnd <- 40
      seedsPerRun <- 40
      if(highFi*dim >= 100){seedsPerRun <- 5}
      if(highFi*dim >= 250){seedsPerRun <- 1}
      runDataKrig <- createScriptStructure(c("surrogateModelWithFixedSample"), c("kriging"), c(func), highFi*dim, lowFi*dim, seedsStart = 1, seedsEnd = 40, seedsPerRun = seedsPerRun)
      
      seedsPerRun <- 40
      if(lowFi*dim >= 100){seedsPerRun <- 5}
      if(lowFi*dim >= 250){seedsPerRun <- 1}
      runDataCoKrig <- createScriptStructure(c("surrogateModelWithFixedSample"), c("cokriging"), c(func), highFi*dim, lowFi*dim, seedsStart = 1, seedsEnd = 40, seedsPerRun = seedsPerRun)
      if(lowFi != 4 || highFi != 2 || func != functions[[1]]){
        existingRunData <- rbind(existingRunData, runDataKrig, runDataCoKrig)
      }else{
        existingRunData <- rbind(runDataKrig, runDataCoKrig)
      }
    }
  }
}
write.table(existingRunData, "data/runScripts/experimentalRunSurrogateModelWithFixedSample.txt", quote = FALSE, row.names = FALSE)

# Do the same with the SOLAR features
functions <- paste0("SOLAR", c("0.10", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70", "0.80", "0.90"))
features <- read.table("data/features/featuresClean.txt", header = TRUE, sep = " ")
order <- match(functions, features$instances)
dims <- features[order, 'feature_dimension']
for(i in 1:length(functions)){
  func <- functions[[i]]
  dim <- dims[[i]]
  print(func)
  print(dim)
  for(lowFi in c(4, 8, 12, 16, 20)){
    for(highFi in c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)){
      if(highFi > lowFi){next}
      #print(paste0("Low ", lowFi, " high ", highFi))
      seedsStart <- 1
      seedsEnd <- 40
      seedsPerRun <- 10
      runDataKrig <- createScriptStructure(c("surrogateModelWithFixedSample"), c("kriging"), c(func), highFi*dim, lowFi*dim, seedsStart = 1, seedsEnd = 40, seedsPerRun = seedsPerRun)
      runDataCoKrig <- createScriptStructure(c("surrogateModelWithFixedSample"), c("cokriging"), c(func), highFi*dim, lowFi*dim, seedsStart = 1, seedsEnd = 40, seedsPerRun = seedsPerRun)
      if(lowFi != 4 || highFi != 2 || func != functions[[1]]){
        existingRunData <- rbind(existingRunData, runDataKrig, runDataCoKrig)
      }else{
        existingRunData <- rbind(runDataKrig, runDataCoKrig)
      }
    }
  }
}
write.table(existingRunData, "data/runScripts/experimentalRunSurrogateModelWithFixedSampleSOLAR.txt", quote = FALSE, row.names = FALSE)

































source("R_code/experimentalRunScriptCreator.R")
options(warn = -1)







# Current run, using a "diverse" set of instances
print("CREATING SMALLER EXPERIMENTAL SCRIPT.")
functions <- read.table("data/availableFunctions/chosenSubsubsetAllFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
# Take out function 16
functions <- functions[str_which(functions, "COCOfunction16", negate = TRUE)]
for(algo in c("krigingMaxVar", "krigingMaxImpr", "cokrigingMaxVar", "cokrigingMaxImpr")){
  for(budget in c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)){
    for(costRatio in c(0.5, 0.2, 0.1, 0.05, 0.02, 0.01)){
      seedsStart <- 1
      seedsEnd <- 20
      seedsPerRun <- 20
      runData <- createScriptStructure(algo, functions, budget, costRatio, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
      if(budget != 10 || costRatio != 0.5 || algo != "krigingMaxVar"){
        existingRunData <- read.table("data/runScripts/runWithBudget.txt", header = TRUE, sep = " ", fill = TRUE)
        runData <- rbind(existingRunData, runData)
      }
      write.table(runData, "data/runScripts/runWithBudget.txt", quote = FALSE, row.names = FALSE)
    }
  }
}
print("COMPLETED SMALLER EXPERIMENTAL SCRIPT.")













# print("CREATING LARGE EXPERIMENTAL SCRIPT.")
# litFunctions <- read.table("data/availableFunctions/biSourceAll.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
# # newFunctions <- read.table("data/availableFunctions/newFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
# # functions <- c(litFunctions, newFunctions)
# functions <- litFunctions
# for(lowFi in c(25, 50)){
#   for(highFi in c(10, 25, 50)){
#     if(highFi > lowFi){next}
#     print(paste0("Low ", lowFi, " high ", highFi))
#     seedsStart <- 1
#     seedsEnd <- 40
#     seedsPerRun <- 40
#     runData <- createScriptStructure(c('kriging', 'cokriging'), functions, highFi, lowFi, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
#     if(lowFi != 25 || highFi != 10){
#       existingRunData <- read.table("data/runScripts/wwsmallTesting.txt", header = TRUE, sep = " ", fill = TRUE)
#       runData <- rbind(existingRunData, runData)
#     }
#     write.table(runData, "data/runScripts/wwsmallTesting.txt", quote = FALSE, row.names = FALSE)
#   }
# }
# 
# print("COMPLETED LARGE EXPERIMENTAL SCRIPT.")






print("CREATING SUBSET SAMPLING SCRIPT.")
# for(lowFi in c(25, 50, 75, 100, 250, 500)){
#   for(dim in c(1:10, 20)){
dim <- 0
for(func in c('LiuPedagogical', 'LiuBranin', 'RajnarayanHartmannH3', 'RajnarayanWoods', 
              'LiuStyblinskiTang', 'RajnarayanHartmannH6',
              'COCOfunction1-dim7-seed1-disth1-height0-radius0.025-freq2-amp0', 'ShiDettePepelyshev',
              'COCOfunction1-dim9-seed1-disth1-height0-radius0.025-freq2-amp0', 'LiuAckley10', 'LiuEllipsoid')){
  dim <- dim + 1
  if(dim == 11){dim = 20}
  print(func)
  print(dim)
  for(lowFi in c(5, 10, 25, 50, 75, 100, 250, 500)){
    for(highFi in c(5, 10, 25, 50, 75, 100, 250, 500)){
      if(highFi > lowFi){next}
    
      runData <- createScriptStructure(c('generateSample'), c(func), highFi, lowFi, seedsStart = 1, seedsEnd = 40, seedsPerRun = 40)
      if(lowFi != 5 || func != 'LiuPedagogical' || highFi != 5){
        existingRunData <- read.table("data/runScripts/subsetSampleCreationScript.txt", header = TRUE, sep = " ", fill = TRUE)
        runData <- rbind(existingRunData, runData)
      }
      write.table(runData, "data/runScripts/subsetSampleCreationScript.txt", quote = FALSE, row.names = FALSE)
    }
  }
  for(sizeLow in seq(2, 20, by = 2)){
    for(sizeHigh in seq(2, 20, by = 2)){
      if(sizeHigh > sizeLow){next}
      if(sizeLow * dim %in% c(5, 10, 25, 50, 75, 100, 250, 500) &
         sizeHigh * dim %in% c(5, 10, 25, 50, 75, 100, 250, 500)){next}
      
      runData <- createScriptStructure(c('generateSample'), c(func), sizeHigh * dim, sizeLow * dim, seedsStart = 1, seedsEnd = 40, seedsPerRun = 40)
      existingRunData <- read.table("data/runScripts/subsetSampleCreationScript.txt", header = TRUE, sep = " ", fill = TRUE)
      runData <- rbind(existingRunData, runData)
      write.table(runData, "data/runScripts/subsetSampleCreationScript.txt", quote = FALSE, row.names = FALSE)
    }
  }
  
}
print("COMPLETED SUBSET SAMPLING SCRIPT.")



# print("CREATING LARGE EXPERIMENTAL SCRIPT.")
# 
# for(lowFi in c(25, 50, 75, 100)){
#   for(highFi in c(5, 10, 25, 50, 75, 100)){
#     if(highFi > lowFi){next}
#     print(paste0("Low ", lowFi, " high ", highFi))
#     seedsStart <- 1
#     seedsEnd <- 20
#     seedsPerRun <- 20
#     runData <- createScriptStructure(c('kriging', 'cokriging'), functions, highFi, lowFi, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
#     if(lowFi != 25 || highFi != 5){
#       existingRunData <- read.table("data/runScripts/expandedTesting.txt", header = TRUE, sep = " ", fill = TRUE)
#       runData <- rbind(existingRunData, runData)
#     }
#     write.table(runData, "data/runScripts/expandedTesting.txt", quote = FALSE, row.names = FALSE)
#   }
# }
# 
# print("COMPLETED LARGE EXPERIMENTAL SCRIPT.")




# Script which creates the experimental run for the literature test suite.
# Seeds per run decreases for larger dimensions as it takes longer to run,
# this is done so that each experimental specification takes a similar amount of
# time to run.
print("CREATING SMALLER EXPERIMENTAL SCRIPT.")
# newFunctions <- read.table("data/availableFunctions/newFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
functions <- read.table("data/availableFunctions/chosenSubsubsetAllFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
# Take out function 16
functions <- functions[str_which(functions, "COCOfunction16", negate = TRUE)]
for(lowFi in c(25, 50, 75, 100)){
  for(highFi in c(5, 10, 25, 50, 75, 100)){
    if(highFi > lowFi){next}
    print(paste0("Low ", lowFi, " high ", highFi))
    seedsStart <- 1
    seedsEnd <- 20
    seedsPerRun <- 20
    runData <- createScriptStructure(c('kriging', 'cokriging'), functions, highFi, lowFi, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
    if(lowFi != 25 || highFi != 5){
      existingRunData <- read.table("data/runScripts/smallerExpandedTesting.txt", header = TRUE, sep = " ", fill = TRUE)
      runData <- rbind(existingRunData, runData)
    }
    write.table(runData, "data/runScripts/smallerExpandedTesting.txt", quote = FALSE, row.names = FALSE)
  }
}
print("COMPLETED SMALLER EXPERIMENTAL SCRIPT.")





print("CREATING SMALLER EXPERIMENTAL SCRIPT.")
# newFunctions <- read.table("data/availableFunctions/newFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
functions <- read.table("data/availableFunctions/chosenSubsubsetAllFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
# Take out function 16
functions <- functions[str_which(functions, "COCOfunction16", negate = TRUE)]
for(lowFi in c(20, 40, 60, 80, 100)){
  for(highFi in c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)){
    if(highFi > lowFi){next}
    print(paste0("Low ", lowFi, " high ", highFi))
    seedsStart <- 1
    seedsEnd <- 20
    seedsPerRun <- 20
    runData <- createScriptStructure(c('kriging', 'cokriging'), functions, highFi, lowFi, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
    if(lowFi != 20 || highFi != 5){
      existingRunData <- read.table("data/runScripts/smallerExpandedTestingLargerSampleRange.txt", header = TRUE, sep = " ", fill = TRUE)
      runData <- rbind(existingRunData, runData)
    }
    write.table(runData, "data/runScripts/smallerExpandedTestingLargerSampleRange.txt", quote = FALSE, row.names = FALSE)
  }
}
print("COMPLETED SMALLER EXPERIMENTAL SCRIPT.")














functions <- read.table("data/availableFunctions/chosenSubsubsetAllFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
functions <- functions[str_which(functions, "COCOfunction16", negate = TRUE)]
limitedInformationFeatures <- read.table("data/combinedData/limitedInformationFeatures.txt", header = TRUE, sep = " ")
limitedInformationFeatures$instanceName <- sapply(strsplit(limitedInformationFeatures$instance, "_"), "[[", 1)
order <- match(functions, limitedInformationFeatures$instanceName)
dims <- limitedInformationFeatures[order, 'feature_sample_dimension']
for(i in 1:length(functions)){
  func <- functions[[i]]
  dim <- dims[[i]]
  print(func)
  print(dim)
  for(lowFi in c(4, 8, 12, 16, 20)){
    for(highFi in c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)){
      if(highFi > lowFi){next}
      #print(paste0("Low ", lowFi, " high ", highFi))
      seedsStart <- 1
      seedsEnd <- 20
      if(lowFi*dim > 100){seedsPerRun <- 1}
      else{seedsPerRun <- 20}
      runData <- createScriptStructure(c('kriging', 'cokriging'), c(func), highFi*dim, lowFi*dim, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
      if(lowFi != 4 || highFi != 2 || func != functions[[1]]){
        existingRunData <- read.table("data/runScripts/smallerExpandedTestingRelativeSampleSize.txt", header = TRUE, sep = " ", fill = TRUE)
        runData <- rbind(existingRunData, runData)
      }
      write.table(runData, "data/runScripts/smallerExpandedTestingRelativeSampleSize.txt", quote = FALSE, row.names = FALSE)
    }
  }
}
















functions <- read.table("data/availableFunctions/biSourceAll.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
for(i in 1:length(functions)){
  func <- functions[[i]]
  if(i <= 16){dim <- 1}
  else if(i <= 16){dim <- 1}
  else if(i <= 69){dim <- 2}
  else if(i <= 83){dim <- 3}
  else if(i <= 99){dim <- 4}
  else if(i <= 153){dim <- 5}
  else if(i <= 157){dim <- 6}
  else if(i <= 158){dim <- 8}
  else if(i <= 214){dim <- 10}
  else if(i <= 217){dim <- 20}
  print(func)
  print(dim)
  for(lowFi in c(4, 8, 12, 16, 20)){
    for(highFi in c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)){
      if(highFi > lowFi){next}
      seedsStart <- 1
      seedsEnd <- 20
      if(lowFi*dim > 100){seedsPerRun <- 1}
      else if(substr(func, 1, 5) == "SOLAR"){seedsPerRun <- 1}
      else{seedsPerRun <- 20}
      runData <- createScriptStructure(c('kriging', 'cokriging'), c(func), highFi*dim, lowFi*dim, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
      if(lowFi != 4 || highFi != 2 || func != functions[[1]]){
        existingRunData <- read.table("data/runScripts/LitRelativeSampleSize.txt", header = TRUE, sep = " ", fill = TRUE)
        runData <- rbind(existingRunData, runData)
      }
      write.table(runData, "data/runScripts/LitRelativeSampleSize.txt", quote = FALSE, row.names = FALSE)
    }
  }
}
# Need to add something to the solar name to make sure there is no conflict when running, this is a horrible fix
test <- runData
runData[str_which(runData$instance, "SOLAR"), "instance"] <- paste0(runData[str_which(runData$instance, "SOLAR"), "instance"], str_which(runData$instance, "SOLAR"))
write.table(runData, "data/runScripts/LitRelativeSampleSize.txt", quote = FALSE, row.names = FALSE)
















functions <- read.table("data/availableFunctions/chosenSubsubsetAllFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
functions <- functions[str_which(functions, "COCOfunction16", negate = TRUE)]
limitedInformationFeatures <- read.table("data/combinedData/limitedInformationFeatures.txt", header = TRUE, sep = " ")
limitedInformationFeatures$instanceName <- sapply(strsplit(limitedInformationFeatures$instance, "_"), "[[", 1)
order <- match(functions, limitedInformationFeatures$instanceName)
dims <- limitedInformationFeatures[order, 'feature_sample_dimension']
for(i in 1:length(functions)){
  func <- functions[[i]]
  dim <- dims[[i]]
  print(func)
  print(dim)
  for(lowFi in c(4, 8, 12, 16, 20)){
    for(highFi in c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)){
      if(highFi > lowFi){next}
      #print(paste0("Low ", lowFi, " high ", highFi))
      seedsStart <- 1
      seedsEnd <- 20
      if(lowFi*dim > 100){seedsPerRun <- 1}
      else{seedsPerRun <- 20}
      runData <- createScriptStructure(c('kriging', 'cokriging'), c(func), highFi*dim, lowFi*dim, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
      if(lowFi != 4 || highFi != 2 || func != functions[[1]]){
        existingRunData <- read.table("data/runScripts/smallerExpandedTestingRelativeSampleSize.txt", header = TRUE, sep = " ", fill = TRUE)
        runData <- rbind(existingRunData, runData)
      }
      write.table(runData, "data/runScripts/smallerExpandedTestingRelativeSampleSize.txt", quote = FALSE, row.names = FALSE)
    }
  }
}























# This is a script to get the features using the given sample for each of the instances
functions <- read.table("data/availableFunctions/chosenSubsubsetAllFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
# Take out function 16
functions <- functions[str_which(functions, "COCOfunction16", negate = TRUE)]
limitedInformationFeatures <- read.table("data/combinedData/limitedInformationFeatures.txt", header = TRUE, sep = " ")
limitedInformationFeatures$instanceName <- sapply(strsplit(limitedInformationFeatures$instance, "_"), "[[", 1)
order <- match(functions, limitedInformationFeatures$instanceName)
dims <- limitedInformationFeatures[order, 'feature_sample_dimension']
seedsStart <- 1
seedsEnd <- 20
seedsPerRun <- 20
for(i in 1:length(functions)){
  func <- functions[[i]]
  dim <- dims[[i]]
  print(paste0(i, " - ", func, " - ", dim))
  for(lowFi in c(20, 40, 60, 80, 100)){
    for(highFi in c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)){
      if(highFi > lowFi){next}
      newRunData <- createScriptStructure(c('sampleFeature'), func, highFi, lowFi, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
      if(lowFi != 20 || highFi != 5){
        funcRunData <- rbind(funcRunData, newRunData)
      }else{
        funcRunData <- newRunData
      }
    }
  }
  for(lowFi in c(4, 8, 12, 16, 20)){
    for(highFi in c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)){
      if(highFi > lowFi){next}
      if((lowFi*dim) %in% c(20, 40, 60, 80, 100) && (highFi*dim) %in% c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)){
        print(paste0("Skipped ", highFi, "(", highFi*dim, ") - ", lowFi, "(", lowFi*dim, ")"))
        next
      }
      newRunData <- createScriptStructure(c('sampleFeature'), func, highFi*dim, lowFi*dim, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
      funcRunData <- rbind(funcRunData, newRunData)
    }
  }
  if(func != functions[[1]]){
    existingRunData <- read.table("data/runScripts/sampleBasedFeatureCalculation.txt", header = TRUE, sep = " ", fill = TRUE)
    runData <- rbind(existingRunData, funcRunData)
  }else{
    runData <- funcRunData
  }
  write.table(runData, "data/runScripts/sampleBasedFeatureCalculation.txt", quote = FALSE, row.names = FALSE)
}
















print("CREATING SPECIAL EXPERIMENTAL SCRIPT.")
# newFunctions <- read.table("data/availableFunctions/newFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
functions <- paste0("NicoSongToalForretal0.00d", c(2, 5, 10))
functions <- c(functions, read.table("data/availableFunctions/singleSourceLit.txt", header = FALSE, sep = " ", fill = TRUE)[[1]])
# Really only care about high fi (as only working with kriging)
for(i in 1:24){
  for(j in c(1, 2, 5, 10)){
    if(j == 1 & i %in% c(8, 9, 17, 18, 19)){next}
    functions <- c(functions, paste0("COCOfunction", i, "-dim", j, "-seed1-disth2-height0-radius0.025-freq2-amp0.1"))
  }
}
for(basemethod in c('krigingARS', 'krigingHT', 'krigingRHT', 'krigingCMAES')){
  for(addedmethod in c('basic', 'local', 'localmaxmin')){
    method <- paste0(basemethod, addedmethod)
    print(method)
    for(lowFi in c(0)){
      for(highFi in seq(5, 100, by = 5)){
        # if(highFi > lowFi){next}
        print(paste0("Low ", lowFi, " high ", highFi))
        seedsStart <- 1
        seedsEnd <- 5
        seedsPerRun <- 5
        runData <- createScriptStructure(c(method), functions, highFi, lowFi, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
        if(lowFi != 0 || highFi != 5 || method != 'krigingARSbasic'){
          existingRunData <- read.table("data/runScripts/specialTesting.txt", header = TRUE, sep = " ", fill = TRUE)
          runData <- rbind(existingRunData, runData)
        }
        write.table(runData, "data/runScripts/specialTesting.txt", quote = FALSE, row.names = FALSE)
      }
    }
  }
}
print("COMPLETED SPECIAL EXPERIMENTAL SCRIPT.")






# Script which creates the experimental run for the literature test suite.
# Seeds per run decreases for larger dimensions as it takes longer to run,
# this is done so that each experimental specification takes a similar amount of
# time to run.
print("CREATING LARGE EXPERIMENTAL SCRIPT.")
litFunctions <- read.table("data/availableFunctions/biSourceAll.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
# newFunctions <- read.table("data/availableFunctions/newFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
newFunctions <- read.table("data/availableFunctions/chosenSubsetNewFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]



functions <- c(litFunctions, newFunctions)
for(lowFi in c(25, 50, 75, 100)){
  for(highFi in c(5, 10, 25, 50, 75, 100)){
    if(highFi > lowFi){next}
    print(paste0("Low ", lowFi, " high ", highFi))
    seedsStart <- 1
    seedsEnd <- 20
    seedsPerRun <- 20
    runData <- createScriptStructure(c('kriging', 'cokriging'), functions, highFi, lowFi, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
    if(lowFi != 25 || highFi != 5){
      existingRunData <- read.table("data/runScripts/expandedTesting.txt", header = TRUE, sep = " ", fill = TRUE)
      runData <- rbind(existingRunData, runData)
    }
    write.table(runData, "data/runScripts/expandedTesting.txt", quote = FALSE, row.names = FALSE)
  }
}

print("COMPLETED LARGE EXPERIMENTAL SCRIPT.")





# Script which creates the experimental run for the literature test suite.
# Seeds per run decreases for larger dimensions as it takes longer to run,
# this is done so that each experimental specification takes a similar amount of
# time to run.
print("CREATING SAMPLE FEATURE EXPERIMENTAL SCRIPT.")
litFunctions <- read.table("data/availableFunctions/biSourceAll.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
# newFunctions <- read.table("data/availableFunctions/newFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
newFunctions <- read.table("data/availableFunctions/chosenSubsetNewFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
functions <- c(litFunctions, newFunctions)
for(lowFi in c(25, 50, 75, 100, 250, 500)){
  for(highFi in c(5, 10, 25, 50, 75, 100)){
    if(highFi > lowFi){next}
    print(paste0("Low ", lowFi, " high ", highFi))
    seedsStart <- 1
    seedsEnd <- 40
    seedsPerRun <- 40
    runData <- createScriptStructure(c('sampleFeature'), functions, highFi, lowFi, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
    if(lowFi != 25 || highFi != 5){
      existingRunData <- read.table("data/runScripts/allExperimentalRunsNoMethod.txt", header = TRUE, sep = " ", fill = TRUE)
      runData <- rbind(existingRunData, runData)
    }
    write.table(runData, "data/runScripts/allExperimentalRunsNoMethod.txt", quote = FALSE, row.names = FALSE)
  }
}

print("COMPLETED FULL LARGE EXPERIMENTAL SCRIPT.")





print("CREATING LARGE EXPERIMENTAL SCRIPT 20 - 40 seeds.")
litFunctions <- read.table("data/availableFunctions/biSourceAll.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
# newFunctions <- read.table("data/availableFunctions/newFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
newFunctions <- read.table("data/availableFunctions/chosenSubsetNewFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
functions <- c(litFunctions, newFunctions)
for(lowFi in c(25, 50, 75, 100)){
  for(highFi in c(5, 10, 25, 50, 75, 100)){
    if(highFi > lowFi){next}
    print(paste0("Low ", lowFi, " high ", highFi))
    seedsStart <- 21
    seedsEnd <- 40
    seedsPerRun <- 40
    runData <- createScriptStructure(c('kriging', 'cokriging'), functions, highFi, lowFi, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
    if(lowFi != 25 || highFi != 5){
      existingRunData <- read.table("data/runScripts/expandedTesting20-40.txt", header = TRUE, sep = " ", fill = TRUE)
      runData <- rbind(existingRunData, runData)
    }
    write.table(runData, "data/runScripts/expandedTesting20-40.txt", quote = FALSE, row.names = FALSE)
  }
}

print("COMPLETED LARGE EXPERIMENTAL SCRIPT 20 - 40 seeds.")



















print("CREATING LARGE 250 EXPERIMENTAL SCRIPT.")
litFunctions <- read.table("data/availableFunctions/biSourceAll.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
newFunctions <- read.table("data/availableFunctions/newFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
functions <- c(litFunctions, newFunctions)
for(lowFi in c(250)){
  for(highFi in c(5, 10, 25, 50, 75, 100)){
    seedsStart <- 21
    seedsEnd <- 40
    seedsPerRun <- 40
    runData <- createScriptStructure(c('kriging', 'cokriging'), functions, highFi, lowFi, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
    if(highFi != 5){
      existingRunData <- read.table("data/runScripts/expandedTesting250-2.txt", header = TRUE, sep = " ", fill = TRUE)
      runData <- rbind(existingRunData, runData)
    }
    write.table(runData, "data/runScripts/expandedTesting250-2.txt", quote = FALSE, row.names = FALSE)
  }
}

print("COMPLETED LARGE 250 EXPERIMENTAL SCRIPT.")


print("CREATING LARGE 500 EXPERIMENTAL SCRIPT.")
litFunctions <- read.table("data/availableFunctions/biSourceAll.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
newFunctions <- read.table("data/availableFunctions/newFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
functions <- c(litFunctions, newFunctions)
for(lowFi in c(500)){
  for(highFi in c(5, 10, 25, 50, 75, 100)){
    seedsStart <- 21
    seedsEnd <- 40
    seedsPerRun <- 10
    runData <- createScriptStructure(c('kriging', 'cokriging'), functions, highFi, lowFi, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
    if(highFi != 5){
      existingRunData <- read.table("data/runScripts/expandedTesting500-2.txt", header = TRUE, sep = " ", fill = TRUE)
      runData <- rbind(existingRunData, runData)
    }
    write.table(runData, "data/runScripts/expandedTesting500-2.txt", quote = FALSE, row.names = FALSE)
  }
}

print("COMPLETED LARGE 500 EXPERIMENTAL SCRIPT.")




print("CREATING CHEAP COKRIG TESTING SCRIPT.")
litFunctions <- read.table("data/availableFunctions/biSourceAll.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
newFunctions <- read.table("data/availableFunctions/newFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
functions <- c(litFunctions, newFunctions)
for(lowFi in c(100)){
  for(highFi in c(5, 10, 25, 50, 75, 100)){
    seedsStart <- 1
    seedsEnd <- 40
    seedsPerRun <- 40
    runData <- createScriptStructure(c('cokrigingcheap-0', 'cokrigingcheap-s'), functions, highFi, lowFi, seedsStart = seedsStart, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
    if(highFi != 5){
      existingRunData <- read.table("data/runScripts/testingCheapLowFi-100.txt", header = TRUE, sep = " ", fill = TRUE)
      runData <- rbind(existingRunData, runData)
    }
    write.table(runData, "data/runScripts/testingCheapLowFi-100.txt", quote = FALSE, row.names = FALSE)
  }
}

print("COMPLETED CHEAP LOW FI EXPERIMENTAL SCRIPT.")


print("CREATING SAMPLE GATHERING SCRIPT.")
litFunctions <- read.table("data/availableFunctions/biSourceAll.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
newFunctions <- read.table("data/availableFunctions/newFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
functions <- c(litFunctions, newFunctions)
runData <- createScriptStructure(c("sample"), functions, 10000, 10000, 1, 1, TRUE, TRUE)
write.table(runData, "data/runScripts/sampleGathering.txt", quote = FALSE, row.names = FALSE)
print("COMPLETED SAMPLE GATHERING SCRIPT.")

















print("CREATING NEW INSTANCES FEATURE ANALYSIS SCRIPT.")
functions <- 1:38
centres <- c(1, 3, 5)
count <- 0
radii <- c(0.025, 0.05, 0.1, 0.2seq(0.05, 0.25, 0.05))
heights <- c(0, 0.25, 0.5, 0.75, 1)
amps <- c(0.1, 0.5, 1.0, 1.5, 2.0)
freqs <- c(2, 5, 10)
savedInstances <- c()
instances <- c()
# First create height based disturbances
for(height in heights){
  for(radius in radii){
    for(amp in amps){
      for(freq in freqs){
        for(func in functions){
          instances <- c(instances, paste0("DisturbanceBasedFunction", func, "-seed1-disth1-height", height, "-radius", radius, "-freq", freq, "-amp", amp))
          instances <- c(instances, paste0("DisturbanceBasedFunction", func, "-seed1-disth2-height", height, "-radius", radius, "-freq", freq, "-amp", amp))
          cat(paste0("\rCreated ", length(instances), " instances "))
        }
      }
    }
  }
}
savedInstances <- c(savedInstances, instances)
# Now source based disturbances
instances <- c()
for(centre in centres){
  for(radius in radii){
    for(amp in amps){
      for(freq in freqs){
        for(func in functions){
          instances <- c(instances, paste0("DisturbanceBasedFunction", func, "-seed1-dists1-centres", centre, "-radius", radius, "-freq", freq, "-amp", amp))
          instances <- c(instances, paste0("DisturbanceBasedFunction", func, "-seed1-dists2-centres", centre, "-radius", radius, "-freq", freq, "-amp", amp))
          cat(paste0("\rCreated ", length(instances), " instances "))
        }
      }
    }
  }
}
# Add on the coco instances
functions <- 1:24
dims <- c(1, 2, 5, 10)
# First create height based disturbances
savedInstances <- c(savedInstances, instances)
instances <- c()
for(height in heights){
  for(radius in radii){
    for(amp in amps){
      for(freq in freqs){
        for(func in functions){
          for(dim in dims){
            if(dim == 1 & func%in%c(8, 9, 17, 18, 19)){next}
            instances <- c(instances, paste0("COCOfunction", func, "-dim", dim, "-seed1-disth1-height", height, "-radius", radius, "-freq", freq, "-amp", amp))
            instances <- c(instances, paste0("COCOfunction", func, "-dim", dim, "-seed1-disth2-height", height, "-radius", radius, "-freq", freq, "-amp", amp))
            cat(paste0("\rCreated ", length(instances), " instances "))
          }
        }
      }
    }
  }
}
# Now source based disturbances
savedInstances <- c(savedInstances, instances)
instances <- c()
for(centre in centres){
  for(radius in radii){
    for(amp in amps){
      for(freq in freqs){
        for(func in functions){
          for(dim in dims){
            if(dim == 1 & func%in%c(8, 9, 17, 18, 19)){next}
            instances <- c(instances, paste0("COCOfunction", func, "-dim", dim, "-seed1-dists1-centres", centre, "-radius", radius, "-freq", freq, "-amp", amp))
            instances <- c(instances, paste0("COCOfunction", func, "-dim", dim, "-seed1-dists2-centres", centre, "-radius", radius, "-freq", freq, "-amp", amp))
            cat(paste0("\rCreated ", length(instances), " instances "))
          }
        }
      }
    }
  }
}
savedInstances <- c(savedInstances, instances)
cat(" - done.\n")

write(savedInstances, "data/availableFunctions/allNewFunctions.txt")


# # Create run, just give sample budget of 0 so that only feature is calculated
# runData <- data.frame("instance" = instances)
# runData$technique <- "featuresOnly"
# runData$highFiBudget <- 0
# runData$lowFiBudget <- 0
# runData$seeds <- "1-1"
# write.table(runData, paste0("data/runScripts/newInstancesFeatures.txt"), quote = FALSE, row.names = FALSE)

print("CREATED NEW INSTANCES FEATURE ANALYSIS SCRIPT.")










