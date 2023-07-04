#ifndef RUN_EXPERIMENTS_HPP
#define RUN_EXPERIMENTS_HPP

#include "libraries.hpp"
#include "aux_solvers.hpp"
#include "functions.hpp"
#include "input_output.hpp"
#include "surrogate_models.hpp"


void processExperiment(string outputFilename, string instructionLine, bool printInfo = true);

void generateAndSaveSample(string outputFilename, string problemType, string functionName, int highFiBudget, int lowFiBudget, int seed, string method, bool printInfo = true);

// Either reads in or generates initial sample. If the sample cannot be read in from the file, a random LHS sample is generated instead.
pair<vector<VectorXd>, vector<VectorXd> > readInOrGenerateInitialSample(BiFidelityFunction* function, int highFiBudget, int lowFiBudget, int seed, bool printInfo = true);


SurrogateModel* processModelName(string name, BiFidelityFunction* function, int seed);

void assessSurrogateModelWithFixedSample(string outputFilename, string problemType, string functionName, int highFiBudget, int lowFiBudget, int seed, string method, bool printInfo = true);



// // Executes the experiment as speficied by "processExperiment". 
// // double executeExperiment(string filename, string functionName, string technique, int highFiBudget, int lowFiBudget, int seed, double r, vector<double> pVals,
// // 							bool printSolverInfo = true, bool printAllInfo = true, int testSampleSize = 1000, int auxSolverIterations = 1000);
// vector<double> executeExperiment(string filename, string functionName, string technique, int problemType, int budget, double lowFiBudgetORcostRatio, int seed,
// 									bool printInfo = true, bool printSolverInfo = false, int testSampleSize = 5000, int auxSolverIterations = 1000);


// // Generates a high and low fidelity sample. A sample of size "lowFiBudget" (or "highFiBudget" if 'lowFiBudget' = 0) is chosen by finding a random LHS sample plan
// // and doing local swaps to find a morris-mitchell locally optimal sample.
// // If 'lowFiBudget' > 0, a subset of size 'highFiBudget' is chosen to also be morris-mitchell locally optimal.
// // Note if the sample plan already exists (i.e. a file with the right name exists) this work is skipped.
// void generateAndSaveSample(string functionName, int highFiBudget, int lowFiBudget, int seed, bool printInfo);




// void plotAndAnalysieFunction(string functionName, int highFiBudget, int lowFiBudget, int seed, bool printInfo, bool scaling = true);




#endif