// ######################################################################################

// # Copyright 2023, Nicolau Andrés-Thió

// # Permission is hereby granted, free of charge, to any person obtaining a copy
// # of this software and associated documentation files (the "Software"), to deal
// # in the Software without restriction, including without limitation the rights
// # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// # copies of the Software, and to permit persons to whom the Software is
// # furnished to do so, subject to the following conditions:

// # The above copyright notice and this permission notice shall be included in all
// # copies or substantial portions of the Software.

// # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// # SOFTWARE.

// ######################################################################################

#ifndef RUN_EXPERIMENTS_HPP
#define RUN_EXPERIMENTS_HPP

#include "libraries.hpp"
#include "aux_solvers.hpp"
#include "functions.hpp"
#include "input_output.hpp"
#include "surrogate_models.hpp"
#include "hyperparameter_solvers.hpp"

// Function which deciphers the experiment instructions given in instructionLine, and writes the output to outputFilename. 
// The line should contain the problemType, the problem specifications, and the method used to solve them.
// The current implemented instructions are the following:
// - sampleCreation: Creates a two samples for a specified input function, and saves the scaled (i.e. to lie in [0,1]^d) sample
//						in a file. The problem specification should follow the format (function,nH,nL,sFirst-sLast), where 
//						function should be a string understood by processFunctionName, nH is the size of the high-fidelity sample, 
//						nL is the size of the low-fidelity sample, and samples are gathered for seeds starting from sFirst and ending 
//						with sLast, inclusive. A sample of size nL is first gathered, and a subset of size nH is then chosen. If 
//						nL = 0, then only a sample of size nH is gathered. The two available methods are "randomLHS" which generate
//						both sets randomly based on Latin Hypercube Sampling (LHS), and "morrisMitchellLHS" which starts
//						with random samples and locally optimises them to increase the minimum distance between pairs of points.
//						
// - surrogateModelWithFixedSample: Assesses the accuracy of a surrogate model trained on a low- and high-fidelity sample.
//									The problem specification should follow the same format as for sampleCreation. The available 
//									techniques currently are kriging and cokriging. The performance is measured both in terms of
//									model error and model correlation with f_h, both of which are stored in data/clusterResults
//
// - surrogateModelWithBudget: 	Similar to the problem above, assesses the accuracy of a surrogate model approach given 
//								an overall budget and cost ratio. As detailed in the thesis, the decisions to be made include
//								how to split the overall budget between an initial sample and further sampling, the aquisition
//								function and the choice of model.
//
// - optimisationWithBudget: Natural follow up to the above but not implemented here.
void processExperiment(string outputFilename, string instructionLine, bool printInfo = true);

// Function which implementes the sampleCreation detailed above. Generates a sample of size lowFiBudget and a subset of size highFiBudget
// and saves the generated samples to a file in data\samplePlans
void generateAndSaveSample(string outputFilename, string problemType, string functionName, int highFiBudget, int lowFiBudget, int seed, string method, bool printInfo = true);

// Either reads in or generates initial sample. If the sample cannot be read in from the file, a random LHS sample is generated instead.
pair<vector<VectorXd>, vector<VectorXd> > readInOrGenerateInitialSample(BiFidelityFunction* function, int highFiBudget, int lowFiBudget, int seed, bool printInfo = true);

// Given a name, this function initialises a child of the surrogateModel class. Currently accepts kriging, cokriging and Rules-based Co-Kriging as the
// options for a surrogate model. 
SurrogateModel* processModelName(string name, BiFidelityFunction* function, int seed, bool printInfo = true, bool printSolverInfo = false);

// Function which implements the surrogateModelWithFixedSample problem type. Initialises a bifidelity function and a surrogate model,
// and reads in (or generates if the file is non-existent) the sampling plan to be used for this particular sample size and seed. 
void assessSurrogateModelWithFixedSample(string outputFilename, string problemType, string functionName, int highFiBudget, int lowFiBudget, int seed, string method, bool printInfo = true);

// Function which implementes the surrogateModelWithBudget problem type. Initialises a bifidelity function and a surrogate model
// with a specified acquisition function, reads in (or generates if the file is non-existence) a sampling plan, and iteratively
// samples the bifidelity function until the whole budget has been used up.\
// The options for the initial use of the budget include "small" (n_h = d + 1, n_l = 2(d+1)), "medium" (half the budget is used in
// the initial sample with n_l = 2 * n_h), and "all" (the whole budget is used in the initial spread out samping plan, with n_l = 2 * n_h).
// The options for the surrogate model include Kriging, Co-Kriging and Rules-based Co-Kriging (here specified as adaptiveCoKriging).
// The options for the acquisition function include max uncertainty ("variance"), max global uncertainty reduction ("globalVariance"),
// and a max global uncertainty reduction which compares potentially sampling the high- and low-fidelity sources and chooses only one
// as detailed in the thesis ("globalVarianceWithChoice").
// An example call for method is "adaptiveCokriging_globalVarianceWithChoice_half"; for further examples consult the folder data/runScripts
void assessSurrogateModelWithBudget(string outputFilename, string problemType, string functionName, int budget, double costRatio, int seed, string method, bool printInfo = true);


#endif