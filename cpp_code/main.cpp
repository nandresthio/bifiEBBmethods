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

#ifndef MAIN_CPP
#define MAIN_CPP

#include "libraries.hpp"
#include "run_experiments.hpp"


int main(int argc, char *argv[]){


	// BiFidelityFunction* function = processFunctionName("ToalBranin0.00");
	// pair<vector<VectorXd>, vector<VectorXd> > points = readInOrGenerateInitialSample(function, 5, 10, 1, true);
	// vector<VectorXd> sampledPoints = points.first;
	// vector<VectorXd> sampledPointsLow = points.second;
	// vector<double> sampledPointsValues = function->evaluateMany(sampledPoints);
	// vector<double> sampledPointsValuesLow = function->evaluateManyLow(sampledPointsLow);
	// SurrogateModel* model = processModelName("cokriging", function, 1, true, true);

	// model->saveSample(sampledPoints, sampledPointsLow, sampledPointsValues, sampledPointsValuesLow);
	// model->trainModel();
	// VectorXd point;

	// // model->setAquisitionFunction("surface");
	// // point = model->findNextSampleSite();
	// // model->setAquisitionFunction("globalVariance");
	// // point = model->findNextSampleSite();
	// // model->setAquisitionFunction("surface");
	// // point = model->findNextSampleSite();

	// printf("\nDone testing\n");


	// printf("Points\n");
	// for(int i = 0; i < (int)model->sampledPoints_.size(); i++){
	// 	printPoint(model->sampledPoints_[i]);
	// 	printf(" %.2f\n", model->sampledPointsValues_[i]);
	// }
	// printf("Points low\n");
	// for(int i = 0; i < (int)model->sampledPointsLow_.size(); i++){
	// 	printPoint(model->sampledPointsLow_[i]);
	// 	printf(" %.2f\n", model->sampledPointsValuesLow_[i]);
	// }

	// printf("\nGlobal variance based sample\n");
	// model->setAquisitionFunction("globalVariance");
	// point = model->findNextSampleSite();
	// model->addSample(point, true, true);
	// printf("Points\n");
	// for(int i = 0; i < (int)model->sampledPoints_.size(); i++){
	// 	printPoint(model->sampledPoints_[i]);
	// 	printf(" %.2f\n", model->sampledPointsValues_[i]);
	// }
	// printf("Points low\n");
	// for(int i = 0; i < (int)model->sampledPointsLow_.size(); i++){
	// 	printPoint(model->sampledPointsLow_[i]);
	// 	printf(" %.2f\n", model->sampledPointsValuesLow_[i]);
	// }
	// model->trainModel();
	// point = model->findNextSampleSite();
	// model->addSample(point, true, true);
	// model->trainModel();
	// printf("Points\n");
	// for(int i = 0; i < (int)model->sampledPoints_.size(); i++){
	// 	printPoint(model->sampledPoints_[i]);
	// 	printf(" %.2f\n", model->sampledPointsValues_[i]);
	// }
	// printf("Points low\n");
	// for(int i = 0; i < (int)model->sampledPointsLow_.size(); i++){
	// 	printPoint(model->sampledPointsLow_[i]);
	// 	printf(" %.2f\n", model->sampledPointsValuesLow_[i]);
	// }
	// printf("\nDone with global variance based\n");

	// printf("\nSurface based sample\n");
	// model->setAquisitionFunction("surface");
	// point = model->findNextSampleSite();
	// model->addSample(point, true, true);
	// model->trainModel();
	// point = model->findNextSampleSite();
	// model->addSample(point, true, true);
	// model->trainModel();
	// printf("\nDone with surface based\n");

	// printf("\nVariance based sample\n");
	// model->setAquisitionFunction("variance");
	// point = model->findNextSampleSite();
	// model->addSample(point, true, true);
	// model->trainModel();
	// point = model->findNextSampleSite();
	// model->addSample(point, true, true);
	// model->trainModel();
	// printf("\nDone with variance based\n");

	// printf("\nExpectedImprovementbased based sample\n");
	// model->setAquisitionFunction("expectedImprovement");
	// point = model->findNextSampleSite();
	// model->addSample(point, true, true);
	// model->trainModel();
	// point = model->findNextSampleSite();
	// model->addSample(point, true, true);
	// model->trainModel();
	// printf("\nDone with Expected improvement based\n\n\n");

	if(argc == 1){
		// If no inputs are given, show some examples on how to run this code
		printf("Note: This code is intended to run based on experiment information.\n");
		printf("Example run files are stored in data/runScripts. For formatting details,\n");
		printf("consult cpp_code/run_experiments.hpp In order to run an experiment please\n");
		printf("specify file containing information run, followed by ROWSTART, NUMROWS, ROWADD, with:\n\n");
		printf("\t-ROWSTART: Row of the first experiment to run in the run file.\n");
		printf("\t-NUMROWS: (Optional) Number of experiments to run, default is 1.\n");
		printf("\t-ROWADD: (Optional) Addition to ROWSTART, needed when running experiments in the SPARTAN cluster.\n\n");
		printf("The experiments in the specified which will run will go from line ROWADD + (ROWSTART - 1) * NUMROWS + 1 to line ROWADD + ROWSTART * NUMROWS\n\n");
		printf("Example run: main sampleCreationRun 1 10\n\n");


		string functionName = "ToalBranin0.10";
		string modelName = "cokriging";
		int seed = 1;
		int highFiBudget = 10;
		int lowFiBudget = 20;
		bool printInfo = true;

		printf("Next an example on how to use this code is given. The function %s is modelledby the surrogate model %s\n", functionName.c_str(),modelName.c_str()); 
		printf("using %d high-fidelity samples, %d low-fidelity samples, and the random seed %d\n", highFiBudget, lowFiBudget, seed);
		
		// This line creates a bifidelity function class which can be queried for both high and low fidelity function values
		BiFidelityFunction* function = processFunctionName(functionName);
		
		// This line creates a surrogate model based on the function, currently the strings "kriging" and "cokriging" are supported
		// Setting a non-zero seed allows reproducible results.
		SurrogateModel* model = processModelName(modelName, function, seed);

		// The following code reads in an already generated sample, or if a file with these details does not exist,
		// one is generated
		printf("Obtain sample.\n");
		pair<vector<VectorXd>, vector<VectorXd> > points = readInOrGenerateInitialSample(function, highFiBudget, lowFiBudget, seed, printInfo);
		// The next few lines extract the sample and obtain the high- and low-fidelity objective function values
		vector<VectorXd> sampledPoints = points.first;
		vector<VectorXd> sampledPointsLow = points.second;
		vector<double> sampledPointsValues = function->evaluateMany(sampledPoints);
		vector<double> sampledPointsValuesLow = function->evaluateManyLow(sampledPointsLow);

		// The next two lines save the sample to the model, and train the hyperparameters
		model->saveSample(sampledPoints, sampledPointsLow, sampledPointsValues, sampledPointsValuesLow);
		printf("\nTrain model.\n");
		model->trainModel();


		// Generate a large sample to assess model accuracy
		printf("\nObtain large sample to measure model accuracy.\n");
		SampleGenerator* sampleGenerator = new SampleGenerator(function, seed, printInfo);
		vector<VectorXd> testPoints = sampleGenerator->randomLHS(1000);
		vector<double> testPointsValues = function->evaluateMany(testPoints);

		vector<double> oneWeights(1000, 1);
		vector<double> modelVals = model->multipleSurfaceValues(testPoints);
		double performanceError = relativeRootMeanSquaredError(testPointsValues, modelVals);
		double performanceCorrelation = weightedCorrelationCoefficient(testPointsValues, modelVals, oneWeights, false);

		printf("\nTrained model has an error of %.2f and a correlation of %.4f with the true function.\n\n", performanceError, performanceCorrelation);
		return 0;
	}


	if(argc < 3 || argc > 5){
		printf("Error: Please specify file containing information run, followed by ROWSTART, NUMROWS, ROWADD, with:\n\n");
		printf("\t-ROWSTART: Row of the first experiment to run in the run file.\n");
		printf("\t-NUMROWS: (Optional) Number of experiments to run, default is 1.\n");
		printf("\t-ROWADD: (Optional) Addition to ROWSTART, needed when running experiments in the SPARTAN cluster.\n\n");
		printf("The experiments in the specified which will run will go from line ROWADD + (ROWSTART - 1) * NUMROWS + 1 to line ROWADD + ROWSTART * NUMROWS\n\n");
		printf("Example run: main sampleCreationRun 1 10\n\n");
		return 0;
	}


	string filename = argv[1];
	int indexNum = stoi(argv[2]);
	int indexMult, indexAdd;
	if(argc >= 4){indexMult = stoi(argv[3]);}
	else{indexMult = 1;}
	if(argc >= 5){indexAdd = stoi(argv[4]);}
	else{indexAdd = 0;}

	int arrayNumStart = indexAdd + (indexNum - 1) * indexMult + 1;
	int arrayNumEnd = indexAdd + indexNum * indexMult;
	string inputFilename = "data/runScripts/" + filename + ".txt";
	ifstream inputFile(inputFilename, ios::in);
	if(!inputFile.is_open()){printf("Error: Could not open input file %s to read experiment specifications! Stopping now...\n", inputFilename.c_str()); return 0;}
	printf("Reading from file %s, will run jobs from lines %d to %d (inclusive).\n", inputFilename.c_str(), arrayNumStart, arrayNumEnd);

	string outputFilename;
	if(indexMult == 1){outputFilename = "data/clusterResults/" + filename + "_arrayJob" + to_string(arrayNumStart) + ".txt";}
	else{outputFilename = "data/clusterResults/" + filename + "_arrayJob" + to_string(arrayNumStart) + "-" + to_string(arrayNumEnd) + ".txt";}

	string line;
	string header;
	int lineNum = 0;
	while(getline(inputFile, line)){
		if(lineNum >= arrayNumStart && lineNum <= arrayNumEnd){
			// Found the right spot. If this is the first entry, figure out the type of problem instance and save the header in the output file
			// It is assumed the problem type is always the same
			if(lineNum == arrayNumStart){
				stringstream ss(line);
				string token;
				getline(ss, token, ' ');
				string problemType = token;
				if(problemType.compare("sampleCreation") == 0){
					header = "problemType functionName nh nl seed method time\n";
				}else if(problemType.compare("surrogateModelWithFixedSample") == 0){
					header = "problemType functionName nh nl seed method modelError modelCorrelation minVal maxVal time\n";
				}else if(problemType.compare("surrogateModelWithBudget") == 0){
					header = "problemType functionName budget costRatio seed method usedBudget modelError modelCorrelation minVal maxVal time\n";
				}else if(problemType.compare("optimisationWithBudget") == 0){
					header = "problemType functionName budget costRatio seed method modelError modelCorrelation minVal maxVal time\n";
				}else{
					printf("Problem, do not recognise problem type in run file! Must be one of sampleCreation, surrogateModelWithFixedSample, surrogateModelWithBudget or optimisationWithBudget. Stopping now...\n");
					return 0;
				}
				ofstream outputFile(outputFilename);
				if(!outputFile.is_open()){printf("Error: Could not open output file %s to output experiment results! Stopping now...\n", outputFilename.c_str()); return 0;}
				outputFile << header;
				outputFile.close();
			}
			// Process each of the experimental lines
			processExperiment(outputFilename, line);
		}
		lineNum++;
	}
	inputFile.close();
	// Print warning if did not get all the expected lines
	if(lineNum <= arrayNumStart){printf("Could not find expected array jobs, wanted lines %d-%d but file only had %d lines!\n", arrayNumStart, arrayNumEnd, lineNum - 1);}
	else if(lineNum <= arrayNumEnd){printf("Looks like did not get as many array jobs as expected, wanted lines %d-%d but file only had %d lines!\n", arrayNumStart, arrayNumEnd, lineNum - 1);}
}


#endif