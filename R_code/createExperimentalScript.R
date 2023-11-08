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

# This R code is used to generate the scripts which the c++ code will then read
# from. For details on what the scripts need to look like, consult
# the file cpp_code/run_experiments.hpp

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
          if(printInfo){print(paste0("Working on job ", jobs))}
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



# Create a script for the creation of all the sampling plans which are then stored
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

# Create a script for all of the subset sampling plans which are then stored
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
functions <- read.table("cpp_code/bifiEBBbenchmarks/data/availableFunctions/chosenTestSuiteN312.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
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




functions <- read.table("cpp_code/bifiEBBbenchmarks/data/availableFunctions/chosenTestSuiteN312.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
# Order them by dimension
features <- read.table("data/features/featuresClean.txt", header = TRUE, sep = " ")
order <- match(functions, features$instances)
dims <- features[order, 'feature_dimension']
functions <- functions[order(match(dims, 1:20))]
dims <- dims[order(match(dims, 1:20))]
index <- 0
for(i in 1:length(functions)){
  func <- functions[[i]]
  dim <- dims[[i]]
  if(dim != 2){next}
  print(func)
  print(dim)
  seedsEnd <- 5
  seedsPerRun <- 40
  if(dim > 5){seedsPerRun <- 5}
  if(dim > 10){seedsPerRun <- 1}
  for(model in c("kriging", "cokriging", "adaptiveCokriging")){
    for(acquisition in c("variance", "globalVariance", "globalVarianceWithChoice")){
      for(doe in c("small", "half", "all")){
        if(model == "kriging" & acquisition == "globalVarianceWithChoice"){next}
        method <- paste0(model, "_", acquisition, "_", doe)
        for(budget in c(5, 10, 15, 20)){
          # if(doe == "all" & acquisition != "globalVariance"){next}
          for(costRatio in c(0.5, 0.1, 0.025, 0.01)){
            if(model == "kriging" & costRatio != 0.5){next}
            index <- index + 1
            runData <- createScriptStructure(c("surrogateModelWithBudget"), c(method), c(func), budget, costRatio, seedsStart = 1, seedsEnd = seedsEnd, seedsPerRun = seedsPerRun)
            if(index != 1){
              existingRunData <- rbind(existingRunData, runData)
            }else{
              existingRunData <- runData
            }
          }
        }
      }
    }
  }
}
write.table(existingRunData, "data/runScripts/experimentalRunSurrogateModelWithGivenBudgetSmallTest.txt", quote = FALSE, row.names = FALSE)

# write.table(existingRunData, "data/runScripts/experimentalRunSurrogateModelWithGivenBudgetKrigingTest.txt", quote = FALSE, row.names = FALSE)




