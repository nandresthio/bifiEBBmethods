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
				printf("What??\n");
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
		samplePlanFilename = "../data/samplePlans/morrisMitchellLHS-dim" + to_string(function->d_) + "-n" + to_string(lowFiBudget) + "-s" + to_string(seed) + ".txt";
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
		samplePlanFilename = "../data/samplePlans/morrisMitchellSubset-dim" + to_string(function->d_) + "-nL" + to_string(lowFiBudget) + "-nH" + to_string(highFiBudget) + "-s" + to_string(seed) + ".txt";
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
		samplePlanFilename = "../data/samplePlans/morrisMitchellLHS-dim" + to_string(function->d_) + "-n" + to_string(highFiBudget) + "-s" + to_string(seed) + ".txt";
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
































// vector<double> executeExperiment(string filename, string functionName, string technique, int problemType, int budget, double lowFiBudgetORcostRatio, int seed,
// 									bool printInfo, bool printSolverInfo, int testSampleSize, int auxSolverIterations){
	
// 	int lowFiBudget = (int)lowFiBudgetORcostRatio;
// 	if(technique.compare("generateSample") == 0){
// 		generateAndSaveSample(functionName, budget, lowFiBudget, seed, printInfo);
// 		return {};
	
// 	}else if(technique.compare("functionPlot") == 0){
// 		plotAndAnalysieFunction(functionName, budget, lowFiBudget, seed, printInfo);
// 		return {};

// 	}else if(technique.compare(0, 7, "kriging") == 0 ||
// 				technique.compare(0, 9, "cokriging") == 0 ||
// 				technique.compare(0, 17, "adaptivecokriging") == 0){
// 				//  ||
				
// 				// technique.compare("krigingARSbasic") == 0 ||
// 				// technique.compare("krigingARSlocal") == 0 ||
// 				// technique.compare("krigingARSlocalmaxmin") == 0 ||
// 				// technique.compare("krigingRHTbasic") == 0 ||
// 				// technique.compare("krigingRHTlocal") == 0 ||
// 				// technique.compare("krigingRHTlocalmaxmin") == 0 ||
// 				// technique.compare("krigingHTbasic") == 0 ||
// 				// technique.compare("krigingHTlocal") == 0 ||
// 				// technique.compare("krigingHTlocalmaxmin") == 0 ||
// 				// technique.compare("krigingCMAESbasic") == 0 ||
// 				// technique.compare("krigingCMAESlocal") == 0 ||
// 				// technique.compare("krigingCMAESlocalmaxmin") == 0 ||
// 				// technique.compare("cokriging") == 0 ||
// 				// technique.compare("cokrigingcheap-0") == 0 ||
// 				// technique.compare("cokrigingcheap-s") == 0 ||
// 				// technique.compare("predictionKfold") == 0){

// 		clock_t tstart;
// 		tstart = clock();
// 		// First initialise function generator and solver
// 		BiFidelityFunction* function = processFunctionName(functionName);
// 		// Future research should define what kind of default function we want to use
// 		AuxSolver* auxSolver = new RandomHeavyTuneSolver(function, 10000, 10, true, seed, printSolverInfo, false, "maxmin", 2, 1000);
// 		// AuxSolver* auxSolver = new RandomHeavyTuneSolver(function, 10000, 10, true, seed, printSolverInfo, false, "maxmin", 0, 1000);
// 		// Get the sample
// 		int highFiBudget;
// 		if(problemType == 1){
// 			highFiBudget = budget;
// 		}else{
// 			highFiBudget = round(0.2*budget);
// 			lowFiBudget = 5 * round(0.2*budget);
// 		}

// 		printf("high low budget = %d-%d\n", highFiBudget, lowFiBudget);

// 		pair<vector<VectorXd>, vector<VectorXd> > points = readInOrGenerateInitialSample(function, highFiBudget, lowFiBudget, seed, printInfo);
// 		vector<VectorXd> sampledPoints = points.first;
// 		vector<VectorXd> sampledPointsLow = points.second;
// 		vector<double> sampledPointsValues = function->evaluateMany(sampledPoints);
// 		vector<double> sampledPointsValuesLow = function->evaluateManyLow(sampledPointsLow);

// 		// Generate analysis sample
// 		vector<VectorXd> samples;
// 		vector<double> trueVals;
// 		if(functionName.compare(0, 5, "SOLAR") == 0){
// 			// If we are here, read in the known values; taking a large sample will take forever
// 			string fileLocation = "../solar/points/evaluatedPointsReplaced.txt";
// 			samples = readPointsFromFile(fileLocation, 5000, 5);
// 			fileLocation = "../solar/points/evaluatedPointsReplacedValues.txt";
// 			vector<VectorXd> trueVals2 = readPointsFromFile(fileLocation, 5000, 1);
// 			for(int i = 0; i < (int)trueVals2.size(); i++){trueVals.push_back(trueVals2[i](0));}
// 		}else{
// 			SampleGenerator* sampleGenerator = new SampleGenerator(function, seed, false);
// 			samples = sampleGenerator->randomLHS(testSampleSize);
// 			trueVals = function->evaluateMany(samples);
// 			delete sampleGenerator;
// 		}
// 		vector<double> oneWeights(testSampleSize, 1);


// 		// Define performances
// 		double performanceError = -DBL_MAX;
// 		double performanceCorrelation = -DBL_MAX;
// 		// Create model, add sample, train and get predictions
// 		vector<double> modelVals;
// 		double minVal, maxVal;

// 		if(technique.compare(0, 7, "kriging") == 0){
// 			Kriging* krigingModel = new Kriging(function, auxSolver, highFiBudget, seed, printInfo);
// 			krigingModel->saveSample(sampledPoints, sampledPointsValues);
// 			krigingModel->trainModel();

// 			if(problemType == 2){
// 				int oldUsedBudget = -1;
// 				int usedBudget = (int)(krigingModel->sampledPoints_.size());
// 				while(usedBudget < budget){
// 					printf("Used budget: %d / %d.\n", usedBudget, budget);
// 					VectorXd point;
// 					if(abs(usedBudget - oldUsedBudget) < TOL || technique.compare(0, 13, "krigingMaxVar") == 0){
// 						point = krigingModel->chooseNextSampleSite("maximumVariance");
// 					}else if(technique.compare(0, 14, "krigingMaxImpr") == 0){
// 						point = krigingModel->chooseNextSampleSite("expectedImprovement");
// 					}else{
// 						printf("Don't understand kriging sampling technique! Stopping now...\n");
// 						exit(0);
// 					}
// 					oldUsedBudget = usedBudget;
// 					krigingModel->sampleAtLocation(point);
// 					usedBudget = (int)(krigingModel->sampledPoints_.size());
// 					krigingModel->trainModel();
// 				}
// 			}

// 			modelVals = krigingModel->multipleSurfaceValues(samples);
// 			minVal = krigingModel->unscaleObservation(*min_element(krigingModel->sampledPointsValues_.begin(), krigingModel->sampledPointsValues_.end()));
// 			maxVal = krigingModel->unscaleObservation(*max_element(krigingModel->sampledPointsValues_.begin(), krigingModel->sampledPointsValues_.end()));
// 			delete krigingModel;


// 		}else if(technique.compare(0, 9, "cokriging") == 0){

// 			CoKriging* cokrigingModel = new CoKriging(function, auxSolver, highFiBudget, lowFiBudget, seed, printInfo);
// 			cokrigingModel->saveSample(sampledPoints, sampledPointsValues, sampledPointsLow, sampledPointsValuesLow);
// 			cokrigingModel->trainModel();

// 			if(problemType == 2){
// 				double oldUsedBudget = -1;
// 				double usedBudget = (int)(cokrigingModel->sampledPoints_.size()) + lowFiBudgetORcostRatio * (int)(cokrigingModel->lowFiKriging_->sampledPoints_.size());
// 				while(usedBudget < budget){
// 					printf("Used budget: %.2f / %d.\n", usedBudget, budget);
// 					VectorXd point;
// 					if(abs(usedBudget - oldUsedBudget) < TOL || technique.compare(0, 15, "cokrigingMaxVar") == 0){
// 						point = cokrigingModel->chooseNextSampleSite("maximumVariance");
// 					}else if(technique.compare(0, 16, "cokrigingMaxImpr") == 0){
// 						point = cokrigingModel->chooseNextSampleSite("expectedImprovement");
// 					}else{
// 						printf("Don't understand kriging sampling technique! Stopping now...\n");
// 						exit(0);
// 					}
// 					oldUsedBudget = usedBudget;
// 					cokrigingModel->sampleAtLocation(point);
// 					usedBudget = (int)(cokrigingModel->sampledPoints_.size()) + lowFiBudgetORcostRatio * (int)(cokrigingModel->lowFiKriging_->sampledPoints_.size());
// 					cokrigingModel->trainModel();
// 				}
// 			}
// 			modelVals = cokrigingModel->multipleSurfaceValues(samples);
// 			minVal = cokrigingModel->unscaleObservation(*min_element(cokrigingModel->sampledPointsValues_.begin(), cokrigingModel->sampledPointsValues_.end()));
// 			maxVal = cokrigingModel->unscaleObservation(*max_element(cokrigingModel->sampledPointsValues_.begin(), cokrigingModel->sampledPointsValues_.end()));
// 			delete cokrigingModel;

// 		}else if(technique.compare(0, 17, "cokrigingAdaptive") == 0){
// 			// Need to read method, r and p values
// 			stringstream ss(technique);
// 			string line;
// 			// Don't care about the first line
// 			getline(ss, line, '-');
// 			getline(ss, line, '-');
// 			string strategy = line;
// 			getline(ss, line, '-');
// 			double pFactor = stof(line);
// 			getline(ss, line, '-');
// 			double rFactor = stof(line);
// 			AdaptiveCoKriging* adaptiveCokrigingModel = new AdaptiveCoKriging(function, strategy, pFactor, rFactor, auxSolver, highFiBudget, lowFiBudget, seed, printInfo);
// 			adaptiveCokrigingModel->saveSample(sampledPoints, sampledPointsValues, sampledPointsLow, sampledPointsValuesLow);
// 			adaptiveCokrigingModel->trainModel();
// 			modelVals = adaptiveCokrigingModel->multipleSurfaceValues(samples);
// 			minVal = adaptiveCokrigingModel->unscaleObservation(*min_element(adaptiveCokrigingModel->sampledPointsValues_.begin(), adaptiveCokrigingModel->sampledPointsValues_.end()));
// 			maxVal = adaptiveCokrigingModel->unscaleObservation(*max_element(adaptiveCokrigingModel->sampledPointsValues_.begin(), adaptiveCokrigingModel->sampledPointsValues_.end()));
// 			delete adaptiveCokrigingModel;
// 		}

// 		performanceError = relativeRootMeanSquaredError(trueVals, modelVals);
// 		performanceCorrelation = weightedCorrelationCoefficient(trueVals, modelVals, oneWeights, false);
		
// 		if(printInfo){printf("Elapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000) / 1000);}
		

// 		// Store info
// 		ofstream outputFile;
// 		outputFile.open(filename, ios_base::app);
// 		outputFile << functionName << 
// 						" " << technique << 
// 						" " << to_string(highFiBudget) << 
// 						" " << to_string(lowFiBudget) << 
// 						" " << to_string(seed) << 
// 						" " << to_string(performanceError).substr(0,10) <<
// 						" " << to_string(performanceCorrelation).substr(0,10) <<
// 						" " << to_string(minVal).substr(0,10) << 
// 						" " << to_string(maxVal).substr(0,10) <<
// 						" " << (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000) << 
// 						"\n";

// 		outputFile.close();
// 		delete function;
// 		return {performanceError, performanceCorrelation, minVal, maxVal};

// 	}else{
// 		printf("Unkown technique %s! Ending here...\n", technique.c_str()); return {};
// 	}



// 	// 	// THESE THINGS ARE ADDITIONAL, USED FOR THE CURRENT WORK
// 	// 	// double det = DBL_MAX;
// 	// 	// double likelihood = DBL_MAX;
// 	// 	// double detLow = DBL_MAX;
// 	// 	// double likelihoodLow = DBL_MAX;
// 	// 	// double detMid = DBL_MAX;
// 	// 	// double likelihoodMid = DBL_MAX;
// 	// 	// double relativeMin = DBL_MAX;
// 	// 	// double relativeMax = DBL_MAX;
// 	// 	// double rho = DBL_MAX;
// 	// 	// string hyperparamters = "(";
// 	// 	// double minDist = DBL_MAX;

		
// 	// 	if(technique.compare("kriging") == 0 ||
// 	// 			technique.compare("krigingARSbasic") == 0 ||
// 	// 			technique.compare("krigingARSlocal") == 0 ||
// 	// 			technique.compare("krigingARSlocalmaxmin") == 0 ||
// 	// 			technique.compare("krigingRHTbasic") == 0 ||
// 	// 			technique.compare("krigingRHTlocal") == 0 ||
// 	// 			technique.compare("krigingRHTlocalmaxmin") == 0 ||
// 	// 			technique.compare("krigingHTbasic") == 0 ||
// 	// 			technique.compare("krigingHTlocal") == 0 ||
// 	// 			technique.compare("krigingHTlocalmaxmin") == 0 ||
// 	// 			technique.compare("krigingCMAESbasic") == 0 ||
// 	// 			technique.compare("krigingCMAESlocal") == 0 ||
// 	// 			technique.compare("krigingCMAESlocalmaxmin") == 0){

// 	// 		bool printAll = false;
// 	// 		if(technique.compare("krigingARSbasic") == 0){
// 	// 			delete auxSolver;
// 	// 			auxSolver = new ARSsolver(function, 10, 10000, true, seed, printAll);

// 	// 		}else if(technique.compare("krigingARSlocal") == 0){
// 	// 			delete auxSolver;
// 	// 			auxSolver = new ARSsolver(function, 10, 10000, true, seed, printAll, false, "maxmin", 1, 1000);

// 	// 		}else if(technique.compare("krigingARSlocalmaxmin") == 0){
// 	// 			delete auxSolver;
// 	// 			auxSolver = new ARSsolver(function, 10, 10000, true, seed, printAll, false, "maxmin", 1, 1000);

// 	// 		}else if(technique.compare("krigingRHTbasic") == 0){
// 	// 			delete auxSolver;
// 	// 			auxSolver = new RandomHeavyTuneSolver(function, 10000, 10, true, seed, printAll);

// 	// 		}else if(technique.compare("krigingRHTlocal") == 0){
// 	// 			delete auxSolver;
// 	// 			auxSolver = new RandomHeavyTuneSolver(function, 10000, 10, true, seed, printAll, false, "maxmin", 1, 1000);

// 	// 		}else if(technique.compare("krigingRHTlocalmaxmin") == 0){
// 	// 			delete auxSolver;
// 	// 			auxSolver = new RandomHeavyTuneSolver(function, 10000, 10, true, seed, printAll, false, "maxmin", 2, 1000);

// 	// 		}else if(technique.compare("krigingHTbasic") == 0){
// 	// 			delete auxSolver;
// 	// 			auxSolver = new HeavyTuneSolver(function, 10000, 50, true, seed, printAll);

// 	// 		}else if(technique.compare("krigingHTlocal") == 0){
// 	// 			delete auxSolver;
// 	// 			auxSolver = new HeavyTuneSolver(function, 10000, 50, true, seed, printAll, false, "maxmin", 1, 1000);

// 	// 		}else if(technique.compare("krigingHTlocalmaxmin") == 0){
// 	// 			delete auxSolver;
// 	// 			auxSolver = new HeavyTuneSolver(function, 10000, 50, true, seed, printAll, false, "maxmin", 2, 1000);

// 	// 		}else if(technique.compare("krigingCMAESbasic") == 0){
// 	// 			delete auxSolver;
// 	// 			auxSolver = new CMAESsolver(function, 10000, true, seed, printAll);

// 	// 		}else if(technique.compare("krigingCMAESlocal") == 0){
// 	// 			delete auxSolver;
// 	// 			auxSolver = new CMAESsolver(function, 10000, true, seed, printAll, false, "maxmin", 1, 1000);

// 	// 		}else if(technique.compare("krigingCMAESlocalmaxmin") == 0){
// 	// 			delete auxSolver;
// 	// 			auxSolver = new CMAESsolver(function, 10000, true, seed, printAll, false, "maxmin", 2, 1000);

// 	// 		}else{
// 	// 			printf("This is unexpected!!\n");
// 	// 		}

// 	// 		Kriging* krigingModel = new Kriging(function, auxSolver, highFiBudget, seed, printInfo);
// 	// 		krigingModel->saveSample(sampledPoints, sampledPointsValues);
// 	// 		krigingModel->trainModel();
// 	// 		modelVals = krigingModel->multipleSurfaceValues(samples);
// 	// 		minVal = krigingModel->unscaleObservation(*min_element(krigingModel->sampledPointsValues_.begin(), krigingModel->sampledPointsValues_.end()));
// 	// 		maxVal = krigingModel->unscaleObservation(*max_element(krigingModel->sampledPointsValues_.begin(), krigingModel->sampledPointsValues_.end()));

// 	// 		likelihood = krigingModel->concentratedLikelihoodFunction();
// 	// 		det = krigingModel->rMatrix_.determinant();
// 	// 		for(int i = 0; i < (int)krigingModel->theta_.size(); i++){
// 	// 			hyperparamters = hyperparamters + to_string(krigingModel->theta_[i]) + ",";
// 	// 		}
// 	// 		for(int i = 0; i < (int)krigingModel->pVector_.size(); i++){
// 	// 			if(i == (int)(krigingModel->pVector_.size() - 1)){hyperparamters = hyperparamters + to_string(krigingModel->pVector_[i]) + ")";}
// 	// 			else{hyperparamters = hyperparamters + to_string(krigingModel->pVector_[i]) + ",";}
// 	// 		}
// 	// 		for(int i = 0; i < (int)krigingModel->theta_.size(); i++){
// 	// 			double val = pow(krigingModel->theta_[i], -1.0 / krigingModel->pVector_[i]);
// 	// 			if(val < minDist){minDist = val;}
// 	// 		}
			

// 	// 		delete krigingModel;


// 	// 	}else if(technique.compare("cokriging") == 0){

// 	// 		CoKriging* cokrigingModel = new CoKriging(function, auxSolver, highFiBudget, lowFiBudget, seed, printInfo);
// 	// 		cokrigingModel->saveSample(sampledPoints, sampledPointsValues, sampledPointsLow, sampledPointsValuesLow);
// 	// 		cokrigingModel->trainModel();
// 	// 		modelVals = cokrigingModel->multipleSurfaceValues(samples);
// 	// 		minVal = cokrigingModel->unscaleObservation(*min_element(cokrigingModel->sampledPointsValues_.begin(), cokrigingModel->sampledPointsValues_.end()));
// 	// 		maxVal = cokrigingModel->unscaleObservation(*max_element(cokrigingModel->sampledPointsValues_.begin(), cokrigingModel->sampledPointsValues_.end()));

// 	// 		// detLow = cokrigingModel->lowFiKriging_->rMatrix_.determinant();
// 	// 		// likelihoodLow = cokrigingModel->lowFiKriging_->concentratedLikelihoodFunction();
// 	// 		// auto data = cokrigingModel->intermediateMuSigmaCalculator();
// 	// 		// detMid = (get<3>(data)).matrixL().determinant() * (get<3>(data)).vectorD().prod();
// 	// 		// likelihoodMid = cokrigingModel->intermediateConcentratedLikelihoodFunction();
// 	// 		// det = cokrigingModel->cMatrix_.determinant();
// 	// 		// relativeMin = cokrigingModel->scaleObservation(*min_element(modelVals.begin(), modelVals.end()));
// 	// 		// relativeMax = cokrigingModel->scaleObservation(*max_element(modelVals.begin(), modelVals.end()));
// 	// 		// rho = cokrigingModel->rho_;

// 	// 		delete cokrigingModel;
		
// 	// 	// }else if(technique.compare("predictionKfold") == 0){
// 	// 	// 	Predictor* predictor = new Predictor(function, sampledPoints, sampledPointsLow, sampledPointsValues, sampledPointsLowValues, 10, seed, printInfo);
// 	// 	// 	predictor->chooseModel();
// 	// 	// 	predictor->trainModel();
// 	// 	// 	modelVals = predictor->multipleSurfaceValues(samples);
// 	// 	// 	minVal = *min_element(predictor->cokrigingModel_->sampledPointsValues_.begin(), predictor->cokrigingModel_->sampledPointsValues_.end());
// 	// 	// 	maxVal = *max_element(predictor->cokrigingModel_->sampledPointsValues_.begin(), predictor->cokrigingModel_->sampledPointsValues_.end());

// 	// 	// 	delete predictor;

// 	// 	// }else if(technique.compare("cokrigingcheap-0") == 0){

// 	// 	// 	CoKrigingCheap* cokrigingModel = new CoKrigingCheap(function, auxSolver, highFiBudget, '0', seed, printInfo);
// 	// 	// 	cokrigingModel->sampledPoints_ = sampledPoints;
// 	// 	// 	cokrigingModel->sampledPointsValues_ = sampledPointsValues;
// 	// 	// 	cokrigingModel->scaleObservations(cokrigingModel->sampledPointsValues_);
// 	// 	// 	cokrigingModel->trainModel();
// 	// 	// 	modelVals = cokrigingModel->multipleSurfaceValues(samples);
// 	// 	// 	minVal = *min_element(cokrigingModel->sampledPointsValues_.begin(), cokrigingModel->sampledPointsValues_.end());
// 	// 	// 	maxVal = *max_element(cokrigingModel->sampledPointsValues_.begin(), cokrigingModel->sampledPointsValues_.end());

// 	// 	// 	delete cokrigingModel;

// 	// 	// }else if(technique.compare("cokrigingcheap-s") == 0){

// 	// 	// 	CoKrigingCheap* cokrigingModel = new CoKrigingCheap(function, auxSolver, highFiBudget, 's', seed, printInfo);
// 	// 	// 	cokrigingModel->sampledPoints_ = sampledPoints;
// 	// 	// 	cokrigingModel->sampledPointsValues_ = sampledPointsValues;
// 	// 	// 	cokrigingModel->scaleObservations(cokrigingModel->sampledPointsValues_);
// 	// 	// 	cokrigingModel->trainModel();
// 	// 	// 	modelVals = cokrigingModel->multipleSurfaceValues(samples);
// 	// 	// 	minVal = *min_element(cokrigingModel->sampledPointsValues_.begin(), cokrigingModel->sampledPointsValues_.end());
// 	// 	// 	maxVal = *max_element(cokrigingModel->sampledPointsValues_.begin(), cokrigingModel->sampledPointsValues_.end());

// 	// 	// 	delete cokrigingModel;

		
		
// 	// 	}
		
// 	// 	delete auxSolver;

// 	// 	performanceError = relativeRootMeanSquaredError(trueVals, modelVals);
// 	// 	performanceCorrelation = weightedCorrelationCoefficient(trueVals, modelVals, oneWeights, false);
		
// 	// 	if(printInfo){printf("Elapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000) / 1000);}
		

// 	// 	// Store info
// 	// 	ofstream outputFile;
// 	// 	outputFile.open(filename, ios_base::app);
// 	// 	outputFile << functionName << 
// 	// 					" " << technique << 
// 	// 					" " << to_string(highFiBudget) << 
// 	// 					" " << to_string(lowFiBudget) << 
// 	// 					" " << to_string(seed) << 
// 	// 					" " << to_string(performanceError).substr(0,10) <<
// 	// 					" " << to_string(performanceCorrelation).substr(0,10) <<
// 	// 					" " << to_string(minVal).substr(0,10) << 
// 	// 					" " << to_string(maxVal).substr(0,10) <<
// 	// 					" " << (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000) << 
// 	// 					// " " << to_string(detLow) << 
// 	// 					// " " << to_string(likelihoodLow) << 
// 	// 					// " " << to_string(detMid) << 
// 	// 					// " " << to_string(likelihoodMid) << 
// 	// 					// " " << to_string(det) << 
// 	// 					// " " << to_string(relativeMin) << 
// 	// 					// " " << to_string(relativeMax) << 
// 	// 					// " " << to_string(rho) <<
// 	// 					" " << likelihood << 
// 	// 					" " << minDist << 
// 	// 					" " << det << 
// 	// 					" " << hyperparamters << 
// 	// 					"\n";

// 	// 	outputFile.close();

// 	// 	delete function;
// 	// 	return {performanceError, performanceCorrelation, minVal, maxVal};

// 	// }else{
// 	// 	printf("Unkown technique %s! Ending here...\n", technique.c_str()); return {};
// 	// }


// 	// }else if(technique.compare("specialPlotting") == 0){

//  	// 	function = new LiuBraninFunction();
//  	// 	VectorXd point(2);

//  	// 	vector<VectorXd> points;
//  	// 	vector<double> pointsVals;
//  	// 	vector<double> theta;
//  	// 	vector<double> pVector;


//  	// 	point << 0, 5;
//  	// 	points.push_back(point);
//  	// 	pointsVals.push_back(function->evaluate(point));
//  	// 	point << 0, 10;
//  	// 	points.push_back(point);
//  	// 	pointsVals.push_back(function->evaluate(point));
//  	// 	point << 5, 5;
//  	// 	points.push_back(point);
//  	// 	pointsVals.push_back(function->evaluate(point));
//  	// 	point << 5, 10;
//  	// 	points.push_back(point);
//  	// 	pointsVals.push_back(function->evaluate(point));
 		

//  	// 	theta.push_back(0.5);
//  	// 	theta.push_back(0.5);
//  	// 	pVector.push_back(2);
//  	// 	pVector.push_back(2);
 		
// 	// 	Kriging* krigingModel = new Kriging(function, auxSolver, 0, 0, false, false);
// 	// 	krigingModel->sampledPoints_ = points;
// 	// 	krigingModel->sampledPointsValues_ = pointsVals;
// 	// 	krigingModel->theta_ = theta;
// 	// 	krigingModel->pVector_ = pVector;
// 	// 	krigingModel->saveMuSigma();

// 	// 	ofstream outputFile;
// 	// 	printf("Outputting points to file...\n");
// 	// 	outputFile.open("../data/plots/contourVals.csv");
// 	// 	outputFile << "X1,X2,kriging\n";
// 	// 	int numPoints = 100;
// 	// 	for(int i = 0; i < numPoints; i++){
// 	// 		double val1 = function->lowerBound_[0] + i*(function->upperBound_[0] - function->lowerBound_[0])/(numPoints-1);
// 	// 		for(int j = 0; j < numPoints; j++){
// 	// 			double val2 = function->lowerBound_[1] + j*(function->upperBound_[1] - function->lowerBound_[1])/(numPoints-1);
// 	// 			VectorXd point(2);
// 	// 			point(0) = val1;
// 	// 			point(1) = val2;
// 	// 			outputFile << to_string(val1) << ',' << 
// 	// 						to_string(val2) << ',' << 
// 	// 						to_string(krigingModel->surfaceValue(point)) << "\n";
				
// 	// 		}
// 	// 	}
// 	// 	outputFile.close();
// 	// 	outputFile.open("../data/plots/highFiPointsContour.csv");
// 	// 	outputFile << "X1,X2,Y\n";
// 	// 	for(int i = 0; i < (int)points.size(); i++){
// 	// 		outputFile << to_string(points[i](0)) << "," <<
// 	// 						to_string(points[i](1)) << "," <<
// 	// 						to_string(pointsVals[i]) << "\n";
// 	// 	}
// 	// 	outputFile.close();

// 	// 	return {};
// }


// void generateAndSaveSample(string functionName, int highFiBudget, int lowFiBudget, int seed, bool printInfo){
// 	clock_t tstart;
// 	tstart = clock();
// 	// First initialise function, sample generator, solver and model
// 	BiFidelityFunction* function = processFunctionName(functionName);
// 	SampleGenerator* sampleGenerator = new SampleGenerator(function, seed, printInfo);
// 	// Note will not really use this, just needed when initialising a model
// 	AuxSolver* auxSolver = new ARSsolver(function);
// 	Kriging* noModel = new Kriging(function, auxSolver, 0, seed, printInfo);
// 	vector<VectorXd> sampledPointsLow;
// 	vector<VectorXd> sampledPoints;

// 	string samplePlanFilename;

// 	if(lowFiBudget > 0){
// 		if(printInfo){
// 			printf("Generate low fi data sample: ");
// 			sampleGenerator->prePrint_ = "Generate low fi data sample: ";
// 		}
// 		// If here, specifically want to spend the time to find "good" sampling plans, saving them, and then stopping.
// 		// Saved the points scaled from 0 to 1!
// 		// First check whether this exists or not.
// 		samplePlanFilename = "../data/samplePlans/morrisMitchellLHS-dim" + to_string(function->d_) + "-n" + to_string(lowFiBudget) + "-s" + to_string(seed) + ".txt";
// 		sampledPointsLow = readPointsFromFile(samplePlanFilename, lowFiBudget, function->d_);
// 		// If this is empty, the points did not exist and need to create them
// 		if(sampledPointsLow.empty()){
// 			sampledPointsLow = sampleGenerator->randomLHS(lowFiBudget);
// 			noModel->scalePoints(sampledPointsLow);
// 			sampleGenerator->morrisMitchellLocalOptimal(sampledPointsLow);
// 			writePointsToFile(samplePlanFilename, sampledPointsLow, lowFiBudget, function->d_);
// 		}
// 		if(printInfo){
// 			printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/ 1000);
// 			printf("Generate high fi data sample: ");
// 			sampleGenerator->prePrint_ = "Generate high fi data sample: ";
// 		}
// 		// Choose a subset
// 		samplePlanFilename = "../data/samplePlans/morrisMitchellSubset-dim" + to_string(function->d_) + "-nL" + to_string(lowFiBudget) + "-nH" + to_string(highFiBudget) + "-s" + to_string(seed) + ".txt";
// 		sampledPoints = readPointsFromFile(samplePlanFilename, highFiBudget, function->d_);
// 		// Again, if empty need to find them and store them
// 		if(sampledPoints.empty()){
// 			sampledPoints = sampleGenerator->morrisMitchellSubset(sampledPointsLow, highFiBudget);
// 			writePointsToFile(samplePlanFilename, sampledPoints, highFiBudget, function->d_);
// 		}
// 	}else{
// 		if(printInfo){
// 			printf("Generate high fi data sample: ");
// 			sampleGenerator->prePrint_ = "Generate high fi data sample: ";
// 		}
// 		samplePlanFilename = "../data/samplePlans/morrisMitchellLHS-dim" + to_string(function->d_) + "-n" + to_string(highFiBudget) + "-s" + to_string(seed) + ".txt";
// 		sampledPoints = readPointsFromFile(samplePlanFilename, highFiBudget, function->d_);
// 		if(sampledPoints.empty()){
// 			sampledPoints = sampleGenerator->randomLHS(highFiBudget);
// 			noModel->scalePoints(sampledPoints);
// 		}
// 		// Gonna add something here to run a portion of the optimisation, then save it, then repeat until no improvement is found
// 		// First need to know the current "score" i.e. min distance between a pair of points.
// 		tuple<double, int, int> info = sampleGenerator->morrisMitchellCriterion(sampledPoints);
// 		double distance = get<0>(info);
// 		double nextDistance = distance;
// 		// No a do while loop
// 		while(true){
// 			sampleGenerator->morrisMitchellLocalOptimal(sampledPoints, 100);
// 			info = sampleGenerator->morrisMitchellCriterion(sampledPoints);
// 			nextDistance = get<0>(info);
// 			writePointsToFile(samplePlanFilename, sampledPoints, highFiBudget, function->d_);
// 			if(abs(distance - nextDistance) > TOL){
// 				if(printInfo){printf("\nPausing sample optimisation, improved from %.4f to %.4f. Elapsed time since start: %.2f seconds. Saving to file and continuing...\n", distance, nextDistance, (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);}
// 				// Improved, want to rewrite file to save the current points, then continue
// 				writePointsToFile(samplePlanFilename, sampledPoints, highFiBudget, function->d_);
// 				distance = nextDistance;
// 			}else{
// 				// Otherwise optimisation is complete, get out of here!
// 				if(printInfo){printf("\nSample plan is optimised, final min point distance is %.4f.", nextDistance);}
// 				break;
// 			}
// 		}
			
		
// 	}
// 	if(printInfo){printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);}
	
// 	delete sampleGenerator;
// 	delete noModel;
// 	delete auxSolver;
// 	delete function;
// }




// pair<vector<VectorXd>, vector<VectorXd> > readInOrGenerateInitialSample(Function* function, int highFiBudget, int lowFiBudget, int seed, bool printInfo){
// 	clock_t tstart;
// 	tstart = clock();
// 	string samplePlanFilename;
// 	SampleGenerator* sampleGenerator = new SampleGenerator(function, seed, printInfo);
// 	// Note will not really use this, just needed when initialising a model
// 	AuxSolver* auxSolver = new ARSsolver(function);
// 	Kriging* noModel = new Kriging(function, auxSolver, 0, seed, printInfo);

// 	vector<VectorXd> sampledPoints = {};
// 	vector<VectorXd> sampledPointsLow = {};
// 	pair<vector<VectorXd>, vector<VectorXd> > points;
	
// 	// First generate the sample, read it if the file exists
// 	if(lowFiBudget > 0){
// 		// First try to read both sets
// 		samplePlanFilename = "../data/samplePlans/morrisMitchellLHS-dim" + to_string(function->d_) + "-n" + to_string(lowFiBudget) + "-s" + to_string(seed) + ".txt";
// 		sampledPointsLow = readPointsFromFile(samplePlanFilename, lowFiBudget, function->d_);
// 		samplePlanFilename = "../data/samplePlans/morrisMitchellSubset-dim" + to_string(function->d_) + "-nL" + to_string(lowFiBudget) + "-nH" + to_string(highFiBudget) + "-s" + to_string(seed) + ".txt";
// 		sampledPoints = readPointsFromFile(samplePlanFilename, highFiBudget, function->d_);
		
// 		// If there is no low fidelity sample, need to choose both samples
// 		if(sampledPointsLow.empty()){
// 			if(printInfo){
// 				printf("Generating low fi data sample, choosing random LHS sample, NOT Morris-Mitchell locally optimal.\n");
// 				printf("Generate low fi data sample: ");
// 				sampleGenerator->prePrint_ = "Generate low fi data sample: ";
// 			}
// 			sampledPointsLow = sampleGenerator->randomLHS(lowFiBudget);
// 			if(printInfo){printf("\n");}
// 			noModel->scalePoints(sampledPointsLow);
// 			if(printInfo){
// 				printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);
// 				printf("Generate high fi data sample: ");
// 				sampleGenerator->prePrint_ = "Generate high fi data sample: ";
// 			}
// 			sampledPoints = sampleGenerator->morrisMitchellSubset(sampledPointsLow, highFiBudget);
// 			if(printInfo){printf("\n");}
// 		}else if(sampledPoints.empty()){
// 			if(printInfo){
// 				printf("Read in sample locations for low fidelity sample only.");
// 				printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);
// 				printf("Generate high fi data sample: ");
// 				sampleGenerator->prePrint_ = "Generate high fi data sample: ";
// 			}
// 			sampledPoints = sampleGenerator->morrisMitchellSubset(sampledPointsLow, highFiBudget);
// 			if(printInfo){printf("\n");}

// 			// for(int i = 0; i < (int)sampledPoints.size(); i++){
// 			// 	printPoint(sampledPoints[i]);
// 			// 	printf("\n");
// 			// }
// 		}else{
// 			if(printInfo){
// 				printf("Read in sample locations for low and high fidelity samples.");
// 				printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);
// 			}
// 		}

// 		// Unscale both samples and return!
// 		noModel->unscalePoints(sampledPointsLow);
// 		noModel->unscalePoints(sampledPoints);
		
	
// 	}else{
// 		samplePlanFilename = "../data/samplePlans/morrisMitchellLHS-dim" + to_string(function->d_) + "-n" + to_string(highFiBudget) + "-s" + to_string(seed) + ".txt";
// 		sampledPoints = readPointsFromFile(samplePlanFilename, highFiBudget, function->d_);
// 		if(sampledPoints.empty()){
// 			if(printInfo){
// 				printf("Generating high fi data sample, choosing random LHS sample, NOT Morris-Mitchell locally optimal.\n");
// 				printf("Generate high fi data sample: ");
// 				sampleGenerator->prePrint_ = "Generate high fi data sample: ";
// 			}
// 			sampledPoints = sampleGenerator->randomLHS(highFiBudget);
// 			if(printInfo){
// 				printf("Read in sample locations for low and high fidelity samples.");
// 				printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);
// 			}
// 		}else{
// 			noModel->unscalePoints(sampledPoints);
// 			if(printInfo){
// 				printf("Read in sample locations for high fidelity samples.");
// 				printf("\nElapsed time since start: %.2f seconds\n", (clock() - tstart) / (double)(CLOCKS_PER_SEC / 1000)/1000);
// 			}
// 		}
// 	}

// 	delete sampleGenerator;
// 	delete noModel;
// 	delete auxSolver;

// 	points = make_pair(sampledPoints, sampledPointsLow);

// 	return points;
// }



// void plotAndAnalysieFunction(string functionName, int highFiBudget, int lowFiBudget, int seed, bool printInfo, bool scaling){
// 		BiFidelityFunction* function = processFunctionName(functionName);
// 		// Generate analysis sample
// 		SampleGenerator* sampleGenerator = new SampleGenerator(function, seed, false);
// 		vector<VectorXd> samples = sampleGenerator->randomLHS(5000);
// 		vector<double> trueVals = function->evaluateMany(samples);
// 		vector<double> oneWeights(5000, 1);

// 		// Get the sample
// 		pair<vector<VectorXd>, vector<VectorXd> > allPoints = readInOrGenerateInitialSample(function, highFiBudget, lowFiBudget, seed, printInfo);
// 		vector<VectorXd> sampledPoints = allPoints.first;
// 		vector<VectorXd> sampledPointsLow = allPoints.second;
// 		vector<double> sampledPointsValues = function->evaluateMany(sampledPoints);
// 		vector<double> sampledPointsValuesLow = function->evaluateManyLow(sampledPointsLow);

// 		// Save it into the models
// 		// AuxSolver* auxSolver = new ARSsolver(function, 10, 1000, true, seed, printInfo);
// 		// AuxSolver* auxSolverTieBreaker = new ARSsolver(function, 10, 1000, true, seed, printInfo);
// 		bool printSolverInfo = false;
// 		AuxSolver* auxSolver = new RandomHeavyTuneSolver(function, 10000, 10, true, seed, printSolverInfo, false, "maxmin", 2, 1000);
		
// 		printf("Train Kriging model");
// 		Kriging* krigingModel = new Kriging(function, auxSolver, highFiBudget, seed, printInfo, scaling);
// 		krigingModel->saveSample(sampledPoints, sampledPointsValues);
// 		krigingModel->trainModel();
// 		printf(" - done\n");

// 		// Kriging* krigingModelTieBreaker = new Kriging(function, auxSolverTieBreaker, highFiBudget, seed, printInfo, scaling);
// 		// krigingModelTieBreaker->saveSample(sampledPoints, sampledPointsValues);
// 		// krigingModelTieBreaker->trainModel();
// 		// vector<double> krigModelTieBreakerVals = krigingModelTieBreaker->multipleSurfaceValues(samples);
// 		// double performanceErrorKrigTieBreaker = relativeRootMeanSquaredError(trueVals, krigModelTieBreakerVals);
// 		// double performanceCorrelationKrigTieBreaker = weightedCorrelationCoefficient(trueVals, krigModelTieBreakerVals, oneWeights, false);
// 		printf("Train CoKriging model");
// 		CoKriging* cokrigingModel = new CoKriging(function, auxSolver, highFiBudget, lowFiBudget, seed, printInfo, scaling);
// 		cokrigingModel->saveSample(sampledPoints, sampledPointsValues, sampledPointsLow, sampledPointsValuesLow);
// 		cokrigingModel->trainModel();
// 		printf(" - done\n");

// 		printf("Train adaptive CoKriging model");
// 		AdaptiveCoKriging* adaptiveCokrigingModelHighLoc = new AdaptiveCoKriging(function, "lowLoc", 0.5, 0.1, auxSolver, highFiBudget, lowFiBudget, seed, printInfo, scaling);
// 		adaptiveCokrigingModelHighLoc->saveSample(sampledPoints, sampledPointsValues, sampledPointsLow, sampledPointsValuesLow);
// 		adaptiveCokrigingModelHighLoc->trainModel();
// 		printf(" - done\n");

// 		VectorXd point;
// 		for(int i = 0; i < 1; i++){
// 			vector<double> krigModelVals = krigingModel->multipleSurfaceValues(samples);
// 			vector<double> cokrigModelVals = cokrigingModel->multipleSurfaceValues(samples);
// 			vector<double> cokrigModelAdaptiveVals = adaptiveCokrigingModelHighLoc->multipleSurfaceValues(samples);
// 			double performanceErrorKrig = relativeRootMeanSquaredError(trueVals, krigModelVals);
// 			double performanceCorrelationKrig = weightedCorrelationCoefficient(trueVals, krigModelVals, oneWeights, false);
// 			double performanceErrorCoKrig = relativeRootMeanSquaredError(trueVals, cokrigModelVals);
// 			double performanceCorrelationCoKrig = weightedCorrelationCoefficient(trueVals, cokrigModelVals, oneWeights, false);
// 			double performanceErrorCoKrigAdaptiveHigh = relativeRootMeanSquaredError(trueVals, cokrigModelAdaptiveVals);
// 			double performanceCorrelationCoKrigAdaptiveHigh = weightedCorrelationCoefficient(trueVals, cokrigModelAdaptiveVals, oneWeights, false);
// 			printf("Performance Kriging: Error %.4f, Correlation %.4f\nCoKriging no tie breaker: Error %.4f, Correlation %.4f\n", performanceErrorKrig,
// 																															performanceCorrelationKrig,
// 																															performanceErrorCoKrig,
// 																															performanceCorrelationCoKrig);
// 			printf("Performance Adaptive CoKriging: Error %.4f, Correlation %.4f\n", 											performanceErrorCoKrigAdaptiveHigh,
// 																															performanceCorrelationCoKrigAdaptiveHigh);
		
// 			point = krigingModel->chooseNextSampleSite("maximumVariance");
// 			krigingModel->sampleAtLocation(point);
// 			point = cokrigingModel->chooseNextSampleSite("maximumVariance");
// 			cokrigingModel->sampleAtLocation(point);
// 			point = adaptiveCokrigingModelHighLoc->chooseNextSampleSite("maximumVariance");
// 			adaptiveCokrigingModelHighLoc->sampleAtLocation(point);
// 			krigingModel->trainModel();
// 			cokrigingModel->trainModel();
// 			adaptiveCokrigingModelHighLoc->trainModel();
			
			
			
// 			}


// 		// printf("Choose next sample site for each model.\n");
// 		// printf("Kriging");
// 		// VectorXd point = krigingModel->chooseNextSampleSite("expectedImprovement");
// 		// printf("Expected Improvement, point ");
// 		// printPoint(point);
// 		// printf("\n");
// 		// point = krigingModel->chooseNextSampleSite("maximumVariance");
// 		// printf("MaximumVar Vance, point ");
// 		// printPoint(point);
// 		// printf("\n");

// 		// printf("CoKriging");
// 		// point = cokrigingModel->chooseNextSampleSite("expectedImprovement");
// 		// printf("Expected Improvement, point ");
// 		// printPoint(point);
// 		// printf("\n");
// 		// point = cokrigingModel->chooseNextSampleSite("maximumVariance");
// 		// printf("MaximumVar Vance, point ");
// 		// printPoint(point);
// 		// printf("\n");

// 		// printf("Adaptive CoKriging");
// 		// point = adaptiveCokrigingModelHighLoc->chooseNextSampleSite("expectedImprovement");
// 		// printf("Expected Improvement, point ");
// 		// printPoint(point);
// 		// printf("\n");
// 		// point = adaptiveCokrigingModelHighLoc->chooseNextSampleSite("maximumVariance");
// 		// printf("MaximumVar Vance, point ");
// 		// printPoint(point);
// 		// printf("\n");
		
// 		vector<double> krigModelVals = krigingModel->multipleSurfaceValues(samples);
// 		vector<double> cokrigModelVals = cokrigingModel->multipleSurfaceValues(samples);
// 		vector<double> cokrigModelAdaptiveVals = adaptiveCokrigingModelHighLoc->multipleSurfaceValues(samples);
		
// 		double performanceErrorKrig = relativeRootMeanSquaredError(trueVals, krigModelVals);
// 		double performanceCorrelationKrig = weightedCorrelationCoefficient(trueVals, krigModelVals, oneWeights, false);
// 		double performanceErrorCoKrig = relativeRootMeanSquaredError(trueVals, cokrigModelVals);
// 		double performanceCorrelationCoKrig = weightedCorrelationCoefficient(trueVals, cokrigModelVals, oneWeights, false);
// 		double performanceErrorCoKrigAdaptiveHigh = relativeRootMeanSquaredError(trueVals, cokrigModelAdaptiveVals);
// 		double performanceCorrelationCoKrigAdaptiveHigh = weightedCorrelationCoefficient(trueVals, cokrigModelAdaptiveVals, oneWeights, false);
		
// 		// printf("Performance Kriging: Error %.4f, Correlation %.4f\nPerformance CoKriging: Error %.4f, Correlation %.4f\n", performanceErrorKrig,
// 		// 																													performanceCorrelationKrig,
// 		// 																													performanceErrorCoKrig,
// 		// 																													performanceCorrelationCoKrig);

// 		printf("Performance Kriging: Error %.4f, Correlation %.4f\nCoKriging no tie breaker: Error %.4f, Correlation %.4f\n", performanceErrorKrig,
// 																															performanceCorrelationKrig,
// 																															performanceErrorCoKrig,
// 																															performanceCorrelationCoKrig);
// 		printf("Performance Adaptive CoKriging: Error %.4f, Correlation %.4f\n", 											performanceErrorCoKrigAdaptiveHigh,
// 																															performanceCorrelationCoKrigAdaptiveHigh);
		
		
// 		ofstream outputFile;
// 		printf("Outputting points to file...\n");
// 		outputFile.open("../data/plots/plottingVals.csv");
// 		outputFile << "X,chosenIndex,fHigh,fLow,kriging,cokrigingHigh,cokrigingLow,adaptiveCokrigingHigh,adaptiveCokrigingLow,krigingVar,krigingImprov,cokrigingVar,cokrigingImprov,adaptiveCokrigingVar,adaptiveCokrigingImprob,cokriglowvar,adaptivecokriglowvar\n";
// 		// outputFile << "X,chosenIndex,fHigh,kriging,krigingTieBreaker\n";
// 		outputFile.close();
// 		outputFile.open("../data/plots/highFiPoints.csv");
// 		outputFile << "X,Y,chosenIndex\n";
// 		outputFile.close();
// 		outputFile.open("../data/plots/lowFiPoints.csv");
// 		outputFile << "X,Y,chosenIndex\n";
// 		outputFile.close();
// 		outputFile.open("../data/plots/lowFiPointsAdaptive.csv");
// 		outputFile << "X,Y,chosenIndex\n";
// 		outputFile.close();
		
// 		int points = 1000;
// 		for(int chosenIndex = 0; chosenIndex < function->d_; chosenIndex++){
// 			printf("%d chosen index\n", chosenIndex);
// 			for(int i = 0; i < points; i++){
// 				double val = function->lowerBound_[chosenIndex] + i*(function->upperBound_[chosenIndex] - function->lowerBound_[chosenIndex])/(points-1);
// 				if(i == 0){val += 0.0001;}
// 				VectorXd point(function->d_);
// 				point(chosenIndex) = val;
// 				for(int j = 0; j < function->d_; j++){
// 					if(j == chosenIndex){continue;}
// 					point(j) = function->lowerBound_[j] + (function->upperBound_[j] - function->lowerBound_[j])/2;
// 				}
// 				outputFile.open("../data/plots/plottingVals.csv", ios_base::app);
// 				outputFile << to_string(val) << "," <<
// 								to_string(chosenIndex) << "," <<
// 								to_string(function->evaluate(point)) << "," <<
// 								to_string(function->evaluateLow(point)) << "," <<
// 								to_string(krigingModel->surfaceValue(point)) << "," <<
// 								to_string(cokrigingModel->surfaceValue(point)) << "," <<
// 								to_string(cokrigingModel->lowFiKriging_->surfaceValue(point)) << "," <<
// 								to_string(adaptiveCokrigingModelHighLoc->surfaceValue(point)) << "," <<
// 								to_string(adaptiveCokrigingModelHighLoc->lowFiKriging_->surfaceValue(point)) << "," <<
// 								to_string(krigingModel->uncertainty(point)) << "," <<
// 								to_string(krigingModel->expectedImprovement(point)) << "," <<
// 								to_string(cokrigingModel->uncertainty(point)) << "," <<
// 								to_string(cokrigingModel->expectedImprovement(point)) << "," <<
// 								to_string(adaptiveCokrigingModelHighLoc->uncertainty(point)) << "," <<
// 								to_string(adaptiveCokrigingModelHighLoc->expectedImprovement(point)) << "," << 
// 								to_string(cokrigingModel->lowFiKriging_->uncertainty(point)) << "," <<
// 								to_string(adaptiveCokrigingModelHighLoc->lowFiKriging_->uncertainty(point)) << "\n";
// 				// outputFile << to_string(val) << "," <<
// 				// 				to_string(chosenIndex) << "," <<
// 				// 				to_string(function->evaluate(point)) << "," <<
// 				// 				to_string(krigingModel->surfaceValue(point)) << "," <<
// 				// 				to_string(krigingModelTieBreaker->surfaceValue(point)) << "\n";
// 				outputFile.close();
// 			}
			
// 			outputFile.open("../data/plots/highFiPoints.csv", ios_base::app);
// 			for(int i = 0; i < (int)sampledPoints.size(); i++){
// 				outputFile << to_string(sampledPoints[i](chosenIndex)) << "," <<
// 								to_string(sampledPointsValues[i]) << "," <<
// 								to_string(chosenIndex) << "\n";
// 			}
// 			outputFile.close();
// 			outputFile.open("../data/plots/lowFiPoints.csv", ios_base::app);
// 			for(int i = 0; i < (int)sampledPointsLow.size(); i++){
// 				outputFile << to_string(sampledPointsLow[i](chosenIndex)) << "," << 
// 								to_string(sampledPointsValuesLow[i]) << "," <<
// 								to_string(chosenIndex) << "\n";
// 			}
// 			outputFile.close();
// 			adaptiveCokrigingModelHighLoc->unscaleObservations(adaptiveCokrigingModelHighLoc->usedSampledPointsValuesLow_);
// 			adaptiveCokrigingModelHighLoc->unscalePoints(adaptiveCokrigingModelHighLoc->usedSampledPointsLow_);
			
// 			outputFile.open("../data/plots/lowFiPointsAdaptive.csv", ios_base::app);
// 			for(int i = 0; i < (int)adaptiveCokrigingModelHighLoc->usedSampledPointsLow_.size(); i++){
// 				outputFile << to_string(adaptiveCokrigingModelHighLoc->usedSampledPointsLow_[i](chosenIndex)) << "," << 
// 								to_string(adaptiveCokrigingModelHighLoc->usedSampledPointsValuesLow_[i]) << "," <<
// 								to_string(chosenIndex) << "\n";
// 			}
// 			outputFile.close();
// 		}
// 		points = 100;
// 		// if(function->d_ == 2){
// 		// 	outputFile.open("../data/plots/contourVals.csv");
// 		// 	// outputFile << "X1,X2,fHigh,fLow,kriging,cokrigingHigh,cokrigingLow\n";
// 		// 	outputFile << "X1,X2,fHigh,kriging,krigingTieBreaker\n";
// 		// 	for(int i = 0; i < points; i++){
// 		// 		double val1 = function->lowerBound_[0] + i*(function->upperBound_[0] - function->lowerBound_[0])/(points-1);
// 		// 		for(int j = 0; j < points; j++){
// 		// 			double val2 = function->lowerBound_[1] + j*(function->upperBound_[1] - function->lowerBound_[1])/(points-1);
// 		// 			VectorXd point(2);
// 		// 			point(0) = val1;
// 		// 			point(1) = val2;
// 		// 			outputFile << to_string(val1) << ',' << 
// 		// 						to_string(val2) << ',' << 
// 		// 						// to_string(function->evaluate(point)) << "," <<
// 		// 						// to_string(function->evaluateLow(point)) << "," <<
// 		// 						// to_string(krigingModel->surfaceValue(point)) << "," <<
// 		// 						// to_string(cokrigingModel->surfaceValue(point)) << "," <<
// 		// 						// to_string(cokrigingModel->lowFiKriging_->surfaceValue(point)) << "\n";
// 		// 						to_string(function->evaluate(point)) << "," <<
// 		// 						to_string(krigingModel->surfaceValue(point)) << "," <<
// 		// 						to_string(krigingModelTieBreaker->surfaceValue(point)) << "\n";
// 		// 		}
// 		// 	}
// 		// 	outputFile.close();
// 		// 	outputFile.open("../data/plots/highFiPointsContour.csv");
// 		// 	outputFile << "X1,X2,Y\n";
// 		// 	for(int i = 0; i < (int)sampledPoints.size(); i++){
// 		// 		outputFile << to_string(sampledPoints[i](0)) << "," <<
// 		// 						to_string(sampledPoints[i](1)) << "," <<
// 		// 						to_string(sampledPointsValues[i]) << "\n";
// 		// 	}
// 		// 	outputFile.close();
// 		// 	// outputFile.open("../data/plots/lowFiPointsContour.csv");
// 		// 	// outputFile << "X1,X2,Y\n";
// 		// 	// for(int i = 0; i < (int)sampledPointsLow.size(); i++){
// 		// 	// 	outputFile << to_string(sampledPointsLow[i](0)) << "," << 
// 		// 	// 					to_string(sampledPointsLow[i](1)) << "," << 
// 		// 	// 					to_string(sampledPointsLowValues[i]) << "\n";
// 		// 	// }
// 		// 	// outputFile.close();

// 		// }


// 		delete function;
// 		delete cokrigingModel;
// 		delete krigingModel;
// 		delete adaptiveCokrigingModelHighLoc;
// 		// delete krigingModelTieBreaker;
// 		delete auxSolver;
// 		// delete auxSolverTieBreaker;
// 		delete sampleGenerator;
// }












#endif