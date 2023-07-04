#ifndef MAIN_CPP
#define MAIN_CPP

#include "libraries.hpp"
#include "run_experiments.hpp"



// WHAT I WANT THE MAIN FUNCTION TO END UP LOOKING LIKE!
int main(int argc, char *argv[]){

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
					header = "problemType functionName budget costRatio seed method modelError modelCorrelation minVal maxVal time\n";
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