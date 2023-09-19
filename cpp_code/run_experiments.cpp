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

#ifndef RUN_EXPERIMENTS_CPP
#define RUN_EXPERIMENTS_CPP

#include "run_experiments.hpp"


void processExperiment(string outputFilename, string instructionLine, bool printInfo){
	// Extract information, it should be made up of problemType, problem, and method
	stringstream ss(instructionLine);
	string token;
	getline(ss, token, ' ');
	string problemType = token;
	getline(ss, token, ' ');
	string problem = token;
	getline(ss, token, ' ');
	string method = token;
	printf("Got specified problem type %s, with problem %s and method %s\n", problemType.c_str(), problem.c_str(), method.c_str());
	// Remove brackets from problem
	problem.erase(remove(problem.begin(), problem.end(), '('), problem.end());
	problem.erase(remove(problem.begin(), problem.end(), ')'), problem.end());
	// Define behaviour dependent on the problem type
	if(problemType.compare("sampleCreation") == 0 || problemType.compare("surrogateModelWithFixedSample") == 0){
		// Here need to extract function name, high and low budgets, and seeds
		stringstream ss2(problem);
		string token2;
		getline(ss2, token2, ',');
		string functionName = token2;
		getline(ss2, token2, ',');
		int nh = stoi(token2);
		getline(ss2, token2, ',');
		int nl = stoi(token2);
		getline(ss2, token2, ',');
		stringstream ss3(token2);
		string token3;
		// Get the seeds
		getline(ss3, token3, '-');
		int seedStart = stoi(token3);
		getline(ss3, token, '-');
		int seedEnd = stoi(token);
		for(int seed = seedStart; seed <= seedEnd; seed++){
			if(problemType.compare("sampleCreation") == 0){
				if(printInfo){printf("Running seed %d\n", seed);}
				generateAndSaveSample(outputFilename, problemType, functionName, nh, nl, seed, method);
				if(printInfo){printf("\n");}
			}else{
				// So here is the whole read in sample, initialise model and assess accuracy
				if(printInfo){printf("Running seed %d\n", seed);}
				assessSurrogateModelWithFixedSample(outputFilename, problemType, functionName, nh, nl, seed, method);
				if(printInfo){printf("\n");}
			}
		}
	}else if(problemType.compare("surrogateModelWithBudget") == 0 || problemType.compare("optimisationWithBudget") == 0){
		// Need to do something here, not yet implemented
		printf("Not yet implemented!\n");
	}else{
		printf("Problem, do not recognise problem type in run file! Must be one of sampleCreation, surrogateModelWithFixedSample, surrogateModelWithBudget or optimisationWithBudget. Stopping now...\n");
		exit(0);
	}
}


void generateAndSaveSample(string outputFilename, string problemType, string functionName, int highFiBudget, int lowFiBudget, int seed, string method, bool printInfo){
	clock_t tstart;
	tstart = clock();
	// Check that the specified method makes sense
	if(method.compare("randomLHS") != 0 && method.compare("morrisMitchellLHS") != 0){
		printf("Do not understand the method for creating the sample, must be one of randomLHS or morrisMitchellLHS. Stopping now...\n");
		exit(0);
	}
	// First initialise function, sample generator, solver and model
	BiFidelityFunction* function = processFunctionName(functionName);
	SampleGenerator* sampleGenerator = new SampleGenerator(function, seed, printInfo);
	// Note will not really use this, just needed when initialising a model
	AuxSolver* auxSolver = new ARSsolver(function);
	Kriging* noModel = new Kriging(function, auxSolver, 0, seed, printInfo);
	vector<VectorXd> sampledPointsLow;
	vector<VectorXd> sampledPoints;
	string samplePlanFilename;
	if(lowFiBudget > 0){
		if(printInfo){
			printf("Generate low fi data sample: ");
			sampleGenerator->prePrint_ = "Generate low fi data sample: ";
		}
		// If here, specifically want to spend the time to find "good" sampling plans, saving them, and then stopping.
		// Saved the points scaled from 0 to 1!
		// First check whether this exists or not.
		samplePlanFilename = "data/samplePlans/morrisMitchellLHS-dim" + to_string(function->d_) + "-n" + to_string(lowFiBudget) + "-s" + to_string(seed) + ".txt";
		sampledPointsLow = readPointsFromFile(samplePlanFilename, lowFiBudget, function->d_);
		// If this is empty, the points did not exist and need to create them
		if(sampledPointsLow.empty()){
			if(printInfo){printf("Did not read sample file, generating sample.\n");}
			sampledPointsLow = sampleGenerator->randomLHS(lowFiBudget);
			noModel->scalePoints(sampledPointsLow);
			// Check if we want locally optimal points or are happy with a random LHS plan
			if(method.compare("morrisMitchellLHS") == 0){
				sampleGenerator->morrisMitchellLocalOptimal(sampledPointsLow);
			}
			writePointsToFile(samplePlanFilename, sampledPointsLow, lowFiBudget, function->d_);
		}
		if(printInfo){
			printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/ 1000);
			printf("Generate high fi data sample: ");
			sampleGenerator->prePrint_ = "Generate high fi data sample: ";
		}
		// Choose a subset
		samplePlanFilename = "data/samplePlans/morrisMitchellSubset-dim" + to_string(function->d_) + "-nL" + to_string(lowFiBudget) + "-nH" + to_string(highFiBudget) + "-s" + to_string(seed) + ".txt";
		sampledPoints = readPointsFromFile(samplePlanFilename, highFiBudget, function->d_);
		// Again, if empty need to find them and store them
		// Always take a MorrisMitchell subset as this is relatively quick eben for large datasets
		if(sampledPoints.empty()){
			if(printInfo){printf("Did not read subset file, generating subset.\n");}
			sampledPoints = sampleGenerator->morrisMitchellSubset(sampledPointsLow, highFiBudget);
			writePointsToFile(samplePlanFilename, sampledPoints, highFiBudget, function->d_);
		}
	}else{
		if(printInfo){
			printf("Generate high fi data sample: ");
			sampleGenerator->prePrint_ = "Generate high fi data sample: ";
		}
		samplePlanFilename = "data/samplePlans/morrisMitchellLHS-dim" + to_string(function->d_) + "-n" + to_string(highFiBudget) + "-s" + to_string(seed) + ".txt";
		sampledPoints = readPointsFromFile(samplePlanFilename, highFiBudget, function->d_);
		if(sampledPoints.empty()){
			if(printInfo){printf("Did not read sample file, generating sample.\n");}
			sampledPoints = sampleGenerator->randomLHS(highFiBudget);
			noModel->scalePoints(sampledPoints);
		}
		// Gonna add something here to run a portion of the optimisation, then save it, then repeat until no improvement is found
		// First need to know the current "score" i.e. min distance between a pair of points.
		tuple<double, int, int> info = sampleGenerator->morrisMitchellCriterion(sampledPoints);
		double distance = get<0>(info);
		double nextDistance = distance;
		// Now a do while loop
		// Don't do any of this if happy with a random LHS though
		if(method.compare("morrisMitchellLHS") == 0){
			while(true){
				sampleGenerator->morrisMitchellLocalOptimal(sampledPoints, 100);
				info = sampleGenerator->morrisMitchellCriterion(sampledPoints);
				nextDistance = get<0>(info);
				writePointsToFile(samplePlanFilename, sampledPoints, highFiBudget, function->d_);
				if(abs(distance - nextDistance) > TOL){
					if(printInfo){printf("\nPausing sample optimisation, improved from %.4f to %.4f. Elapsed time since start: %.2f seconds. Saving to file and continuing...\n", distance, nextDistance, (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);}
					// Improved, want to rewrite file to save the current points, then continue
					writePointsToFile(samplePlanFilename, sampledPoints, highFiBudget, function->d_);
					distance = nextDistance;
				}else{
					// Otherwise optimisation is complete, get out of here!
					if(printInfo){printf("\nSample plan is optimised, final min point distance is %.4f.", nextDistance);}
					break;
				}
			}
		}else{
			writePointsToFile(samplePlanFilename, sampledPoints, highFiBudget, function->d_);
		}
			
		
	}
	if(printInfo){printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);}

	// Output to file
	ofstream outputFile;
	outputFile.open(outputFilename, ios_base::app);
	outputFile << problemType <<
					" " << functionName <<
					" " << to_string(highFiBudget) << 
					" " << to_string(lowFiBudget) << 
					" " << to_string(seed) << 
					" " << method <<
					" " << (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000 << 
					"\n";

	outputFile.close();

	delete sampleGenerator;
	delete noModel;
	delete auxSolver;
	delete function;

}



pair<vector<VectorXd>, vector<VectorXd> > readInOrGenerateInitialSample(BiFidelityFunction* function, int highFiBudget, int lowFiBudget, int seed, bool printInfo){
	clock_t tstart;
	tstart = clock();
	string samplePlanFilename;
	SampleGenerator* sampleGenerator = new SampleGenerator(function, seed, printInfo);
	// Note will not really use this, just needed when initialising a model
	AuxSolver* auxSolver = new ARSsolver(function);
	SurrogateModel* noModel = new SurrogateModel(function, auxSolver, 0, seed, printInfo);

	vector<VectorXd> sampledPoints = {};
	vector<VectorXd> sampledPointsLow = {};
	pair<vector<VectorXd>, vector<VectorXd> > points;
	
	// First generate the sample, read it if the file exists
	if(lowFiBudget > 0){
		// First try to read both sets
		samplePlanFilename = "data/samplePlans/morrisMitchellLHS-dim" + to_string(function->d_) + "-n" + to_string(lowFiBudget) + "-s" + to_string(seed) + ".txt";
		sampledPointsLow = readPointsFromFile(samplePlanFilename, lowFiBudget, function->d_);
		samplePlanFilename = "data/samplePlans/morrisMitchellSubset-dim" + to_string(function->d_) + "-nL" + to_string(lowFiBudget) + "-nH" + to_string(highFiBudget) + "-s" + to_string(seed) + ".txt";
		sampledPoints = readPointsFromFile(samplePlanFilename, highFiBudget, function->d_);
		
		// If there is no low fidelity sample, need to choose both samples
		if(sampledPointsLow.empty()){
			if(printInfo){
				printf("Generating low fi data sample, choosing random LHS sample, NOT Morris-Mitchell locally optimal.\n");
				printf("Generate low fi data sample: ");
				sampleGenerator->prePrint_ = "Generate low fi data sample: ";
			}
			sampledPointsLow = sampleGenerator->randomLHS(lowFiBudget);
			if(printInfo){printf("\n");}
			noModel->scalePoints(sampledPointsLow);
			if(printInfo){
				printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);
				printf("Generate high fi data sample: ");
				sampleGenerator->prePrint_ = "Generate high fi data sample: ";
			}
			sampledPoints = sampleGenerator->morrisMitchellSubset(sampledPointsLow, highFiBudget);
			if(printInfo){printf("\n");}
		}else if(sampledPoints.empty()){
			if(printInfo){
				printf("Read in sample locations for low fidelity sample only.");
				printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);
				printf("Generate high fi data sample: ");
				sampleGenerator->prePrint_ = "Generate high fi data sample: ";
			}
			sampledPoints = sampleGenerator->morrisMitchellSubset(sampledPointsLow, highFiBudget);
			if(printInfo){printf("\n");}

		}else{
			if(printInfo){
				printf("Read in sample locations for low and high fidelity samples.");
				printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);
			}
		}

		// Unscale both samples and return!
		noModel->unscalePoints(sampledPointsLow);
		noModel->unscalePoints(sampledPoints);
		
	
	}else{
		samplePlanFilename = "data/samplePlans/morrisMitchellLHS-dim" + to_string(function->d_) + "-n" + to_string(highFiBudget) + "-s" + to_string(seed) + ".txt";
		sampledPoints = readPointsFromFile(samplePlanFilename, highFiBudget, function->d_);
		if(sampledPoints.empty()){
			if(printInfo){
				printf("Generating high fi data sample, choosing random LHS sample, NOT Morris-Mitchell locally optimal.\n");
				printf("Generate high fi data sample: ");
				sampleGenerator->prePrint_ = "Generate high fi data sample: ";
			}
			sampledPoints = sampleGenerator->randomLHS(highFiBudget);
			if(printInfo){
				printf("Read in sample locations for low and high fidelity samples.");
				printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);
			}
		}else{
			noModel->unscalePoints(sampledPoints);
			if(printInfo){
				printf("Read in sample locations for high fidelity samples.");
				printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);
			}
		}
	}

	delete sampleGenerator;
	delete noModel;
	delete auxSolver;

	points = make_pair(sampledPoints, sampledPointsLow);

	return points;
}


SurrogateModel* processModelName(string name, BiFidelityFunction* function, int seed, bool printInfo, bool printSolverInfo){
	SurrogateModel* model;
	if(name.compare("kriging") == 0){
		AuxSolver* auxSolver = new RandomHeavyTuneSolver(function, false, seed, printSolverInfo);
		model = new Kriging(function, auxSolver, seed, printInfo, true);
		
	}else if(name.compare("cokriging") == 0){
		AuxSolver* auxSolver = new RandomHeavyTuneSolver(function, false, seed, printSolverInfo);
		model = new CoKriging(function, auxSolver, seed, printInfo, true);
	
	}else{
		printf("Problem: Could not match method name!! Stopping here...\n");							
		exit(0);
	}

	return model;
}




void assessSurrogateModelWithFixedSample(string outputFilename, string problemType, string functionName, int highFiBudget, int lowFiBudget, int seed, string method, bool printInfo){
	clock_t tstart;
	tstart = clock();
	// First initialise function generator and solver
	BiFidelityFunction* function;
	// Have an extra line here for SOLAR instances, to add a unique file name to the end of the instance name
	if(functionName.compare(0, 5, "SOLAR") == 0){
		string appendedName = functionName + to_string(highFiBudget) + to_string(lowFiBudget) + to_string(seed) + method;
		function = processFunctionName(appendedName);	
	}else{
		function = processFunctionName(functionName);
	}
	pair<vector<VectorXd>, vector<VectorXd> > points = readInOrGenerateInitialSample(function, highFiBudget, lowFiBudget, seed, printInfo);
	vector<VectorXd> sampledPoints = points.first;
	vector<VectorXd> sampledPointsLow = points.second;
	vector<double> sampledPointsValues = function->evaluateMany(sampledPoints);
	vector<double> sampledPointsValuesLow = function->evaluateManyLow(sampledPointsLow);
	SurrogateModel* model = processModelName(method, function, seed);

	model->saveSample(sampledPoints, sampledPointsLow, sampledPointsValues, sampledPointsValuesLow);
	model->trainModel();

	int testSampleSize = 1000 * function->d_;
	vector<VectorXd> samples;
	vector<double> trueVals;
	if(functionName.compare(0, 5, "SOLAR") == 0){
		// If we are here, read in the known values; taking a large sample will take forever
		string fileLocation = "cpp_code/bifiEBBbenchmarks/data/misc/solarPointsLocations.txt";
		samples = readPointsFromFile(fileLocation, 5000, 5);
		fileLocation = "cpp_code/bifiEBBbenchmarks/data/misc/solarPointsValues.txt";
		vector<VectorXd> trueVals2 = readPointsFromFile(fileLocation, 5000, 1);
		for(int i = 0; i < (int)trueVals2.size(); i++){trueVals.push_back(trueVals2[i](0));}
	}else{
		// Read in large sample plan for reproducibility
		string fileLocation = "data/samplePlans/LHS-dim" + to_string(function->d_) + "-n" + to_string(testSampleSize) + ".txt";
		samples = readPointsFromFile(fileLocation, testSampleSize, function->d_);
		unscalePoints(samples, function);
		// If the points were not read in, generate them
		if((int)samples.size() == 0){
			if(printInfo){printf("Could not read in pre-generated large sample, generating it instead!\n");}
			SampleGenerator* sampleGenerator = new SampleGenerator(function, seed, false);
			samples = sampleGenerator->randomLHS(testSampleSize);
			delete sampleGenerator;
		}
		trueVals = function->evaluateMany(samples);
		
	}
	vector<double> oneWeights(testSampleSize, 1);
	vector<double> modelVals = model->multipleSurfaceValues(samples);
	double minVal = model->unscaleObservation(*min_element(model->sampledPointsValues_.begin(), model->sampledPointsValues_.end()));
	double maxVal = model->unscaleObservation(*max_element(model->sampledPointsValues_.begin(), model->sampledPointsValues_.end()));
	double performanceError = relativeRootMeanSquaredError(trueVals, modelVals);
	double performanceCorrelation = weightedCorrelationCoefficient(trueVals, modelVals, oneWeights, false);

	if(printInfo){printf("Completed with model error %.4f and correlation %.4f\n", performanceError, performanceCorrelation);}
	if(printInfo){printf("Elapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000) / 1000);}

	// Store info
	ofstream outputFile;
	outputFile.open(outputFilename, ios_base::app);
	outputFile << problemType << 
					" " << functionName <<
					" " << to_string(highFiBudget) << 
					" " << to_string(lowFiBudget) << 
					" " << to_string(seed) << 
					" " << method <<
					" " << to_string(performanceError).substr(0,10) <<
					" " << to_string(performanceCorrelation).substr(0,10) <<
					" " << to_string(minVal).substr(0,10) << 
					" " << to_string(maxVal).substr(0,10) <<
					" " << (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000) << 
					"\n";

	outputFile.close();
	delete function;
}









#endif