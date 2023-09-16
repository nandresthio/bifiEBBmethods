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

#ifndef HYPERPARAMETER_SOLVERS_CPP
#define HYPERPARAMETER_SOLVERS_CPP

#include "hyperparameter_solvers.hpp"





AmoebaSolver::AmoebaSolver(Function* function, bool min, int randomSeed, bool printInfo) :
	AuxSolver(function, min, randomSeed, printInfo),
	maxEval_(50 * function->d_){}

AmoebaSolver::AmoebaSolver(Function* function, int maxEval, bool min, int randomSeed, bool printInfo):
	AuxSolver(function, min, randomSeed, printInfo),
	maxEval_(maxEval){}

AmoebaSolver::~AmoebaSolver(){
}

VectorXd AmoebaSolver::optimise(){
	int d = function_->d_;
	// Start with a random set of d+1 points
	vector<VectorXd> points = sampleGenerator_->randomSample(d + 1);
	vector<double> pointVals;
	for(int i = 0; i < d + 1; i++){
		pointVals.push_back(function_->evaluate(points[i]));
	}
	return optimise(points, pointVals);
}


VectorXd AmoebaSolver::optimise(vector<VectorXd> simplexPoints, vector<double> simplexValues){
	double alpha = 1;
	double gamma = 2;
	double rho = 0.5;
	double sigma = 0.5;
	int d = function_->d_;
	VectorXd centroid = VectorXd(d);
	VectorXd reflectionPoint, expansionPoint, contractionPoint;
	double reflectionVal, expansionVal, contractionVal;
	int indexWorst, indexSecondWorst, indexBest;
	bool reflectionWithinBounds, expansionWithinBounds, contractionWithinBounds;

	vector<VectorXd> points = simplexPoints;
	vector<double> pointVals = simplexValues;
	
	int evals = d + 1;
	int maxEval = maxEval_;
	do{
		// Need to find the highest and second highest point; assuming we are minimising
		// i.e. find the worst and second worst points
		indexWorst = 0;
		indexBest = 0;
		for(int i = 1; i < d + 1; i++){
			// If the point at index worst is strictly better than the index being checked, reassign
			if(betterPoint(points[i], pointVals[i], points[indexWorst], pointVals[indexWorst]) == 1){indexWorst = i;}
			if(betterPoint(points[i], pointVals[i], points[indexBest], pointVals[indexBest]) == -1){indexBest = i;}
			
		}
		if(printInfo_){
			printf("\r%sRunning Amoeba, evaluations %d/%d, best point has value %.20f at point (", prePrint_.c_str(), evals, maxEval, pointVals[indexBest]);
			for(int i = 0; i < d; i++){
				printf("%.4f", points[indexBest](i));
				if(i < d-1){printf(", ");}
				else{printf(")");}
			}
		}
		// If the best and worst indexes are the same, it is time to stop
		if(indexWorst == indexBest){
			if(printInfo_){printf(" - completed with all points with the same objective value.\n");}
			break;
		}else if(evals >= maxEval){
			if(printInfo_){printf(" - stopping as reached max number of evaluations.\n");}
			break;
		}else if(printInfo_){printf("                         ");}

		if(indexWorst == 0){indexSecondWorst = 1;}
		else{indexSecondWorst = 0;}
		for(int i = indexSecondWorst + 1; i < d + 1; i++){
			if(indexWorst == i){continue;}
			// If the point at index worst is strictly better than the index being checked, reassign
			if(betterPoint(points[i], pointVals[i], points[indexSecondWorst], pointVals[indexSecondWorst]) == 1){indexSecondWorst = i;}
		}

		// Calculate the centroid of all points except the largest
		centroid = VectorXd(d);
		for(int i = 0; i < d; i++){
			centroid(i) = 0;
			for(int j = 0; j < d + 1; j++){
				if(j == indexWorst){continue;}
				centroid(i) += points[j](i);
			}
			centroid(i) = centroid(i) / d;
		}
		// Reflection
		reflectionPoint = centroid + alpha * (centroid - points[indexWorst]);
		reflectionWithinBounds = function_->pointWithinBounds(reflectionPoint);
		if(reflectionWithinBounds){
			reflectionVal = function_->evaluate(reflectionPoint);
			evals++;
			if(evals >= maxEval){continue;}
		}
		// If the reflected point is strictly better than the second worst, but not strictly better than the best,
		// replace the worst point with the reflection point and go back to the top
		if(reflectionWithinBounds &&
			betterPoint(reflectionPoint, reflectionVal, points[indexSecondWorst], pointVals[indexSecondWorst]) == -1 &&
			betterPoint(reflectionPoint, reflectionVal, points[indexBest], pointVals[indexBest]) != -1){

			points[indexWorst] = reflectionPoint;
			pointVals[indexWorst] = reflectionVal;
			continue;			
		}

		// Expansion
		// If the reflected point is the best point so far
		if(reflectionWithinBounds && betterPoint(reflectionPoint, reflectionVal, points[indexBest], pointVals[indexBest]) == -1){
			expansionPoint = centroid + gamma * (reflectionPoint - centroid);
			expansionWithinBounds = function_->pointWithinBounds(expansionPoint);
			if(expansionWithinBounds){
				expansionVal = function_->evaluate(expansionPoint);
				evals++;
				if(evals >= maxEval){continue;}
			}
			// If the expanded point is better than the reflected point, use it to replace the worst point.
			// Otherwise, use the reflected point to replave the worts point. Go back to top.
			if(expansionWithinBounds && betterPoint(reflectionPoint, reflectionVal, expansionPoint, expansionVal) == 1){
				points[indexWorst] = expansionPoint;
				pointVals[indexWorst] = expansionVal;
			}else{
				points[indexWorst] = reflectionPoint;
				pointVals[indexWorst] = reflectionVal;
			}
			continue;
		}
		// Contraction
		// Note if we are here, the reflection point is at least as bad as the second worst point
		if(reflectionWithinBounds && betterPoint(reflectionPoint, reflectionVal, points[indexWorst], pointVals[indexWorst]) == -1){
			// Contracted point outside
			contractionPoint = centroid + rho * (reflectionPoint - centroid);
			contractionWithinBounds = function_->pointWithinBounds(contractionPoint);
			if(contractionWithinBounds){
				contractionVal = function_->evaluate(contractionPoint);
				evals++;
			}
			// If the contracted point is better than the reflected point
			if(contractionWithinBounds && betterPoint(reflectionPoint, reflectionVal, contractionPoint, contractionVal) == 1){
				// Replace the worst point with the contraction point, go back
				points[indexWorst] = contractionPoint;
				pointVals[indexWorst] = contractionVal;
				continue;
			}
		}else{
			// Contracted point inside
			contractionPoint = centroid + rho * (points[indexWorst] - centroid);
			contractionWithinBounds = function_->pointWithinBounds(contractionPoint);
			if(contractionWithinBounds){
				contractionVal = function_->evaluate(contractionPoint);
				evals++;
				if(evals >= maxEval){continue;}
			}
			// If the contracted point is better than the worst point
			if(contractionWithinBounds && betterPoint(contractionPoint, contractionVal, points[indexWorst], pointVals[indexWorst]) == -1){
				// Replace the worst point with the contraction point, go back
				points[indexWorst] = contractionPoint;
				pointVals[indexWorst] = contractionVal;
				continue;
			}
		}
		// Shrink
		for(int i = 0; i < d + 1; i++){
			if(i == indexBest){continue;}
			points[i] = points[indexBest] + sigma * (points[i] - points[indexBest]);
			pointVals[i] = function_->evaluate(points[i]);
			evals++;
			if(evals >= maxEval){break;}
		}

	}while(true);

	// Storing how many evals were used so it can be accessed by heavy tune and random heavy tune
	evalsUsed_ = evals;

	VectorXd bestPoint = points[indexBest];
	return bestPoint;
}



GASolver::GASolver(Function* function, bool min, int randomSeed, bool printInfo) :
	AuxSolver(function, min, randomSeed, printInfo),
	maxEval_(10000),
	populationSize_(50){}

GASolver::GASolver(Function* function, int maxEval, int populationSize, bool min, int randomSeed, bool printInfo):
	AuxSolver(function, min, randomSeed, printInfo),
	maxEval_(maxEval),
	populationSize_(populationSize){}

GASolver::~GASolver(){
}

VectorXd GASolver::optimise(){
	// This is the main call, but as there is a function which takes in the initial population, 
	// this will simply create the initial population and call the next function
	vector<VectorXd> population = sampleGenerator_->randomSample(populationSize_);
	vector<double> populationVals;
	for(int i = 0; i < populationSize_; i++){
		populationVals.push_back(function_->evaluate(population[i]));
	}
	return optimise(population, populationVals);
}

VectorXd GASolver::optimise(vector<VectorXd> population, vector<double> populationVals){
	double recombination = 0.7;
	// Mutation is chosen to be between 0.5 and 1
	uniform_real_distribution<> mutationDist(0.5, 1);
	uniform_real_distribution<> selectionDist(0, 1);
	// Need a uniform integer distribution to choose members of the population
	uniform_int_distribution<> distrib(0, populationSize_ - 1);
	// Using the best1bin strategy
	int indexBest = 0;
	for(int i = 1; i < populationSize_; i++){
		// If the point at index worst is strictly better than the index being checked, reassign
		if(betterPoint(population[i], populationVals[i], population[indexBest], populationVals[indexBest]) == -1){indexBest = i;}
	}
	int evals = populationSize_;
	int maxEval = maxEval_;
	do{
		if(printInfo_){
			printf("\r%sRunning GA, evaluations %d/%d, best point has value %.20f at point (", prePrint_.c_str(), evals, maxEval, populationVals[indexBest]);
			for(int i = 0; i < function_->d_; i++){
				printf("%.4f", population[indexBest](i));
				if(i < function_->d_-1){printf(", ");}
				else{printf(")                         ");}
			}
		}
		if(evals >= maxEval){if(printInfo_){printf("                                            \n");}break;}
		// Mutate each candidate in the population
		// Mutation constant changes each distribution
		double mutation = mutationDist(randomGenerator_);
		for(int i = 0; i < populationSize_; i++){
			VectorXd mutated;
			// Want to get a mutated point that lies in the search space, repeat until one is found
			do{
				// Randomly choose two members in the population
				int firstMemberIndex = distrib(randomGenerator_);
				int secondMemberIndex;
				do{
					secondMemberIndex = distrib(randomGenerator_);
				}while(firstMemberIndex == secondMemberIndex);
				mutated = population[indexBest] + mutation * (population[firstMemberIndex] - population[secondMemberIndex]);
			}while(!function_->pointWithinBounds(mutated));
			
			VectorXd trialMember = population[i];
			for(int j = 0; j < function_->d_ - 1; j++){
				if(selectionDist(randomGenerator_) < recombination){trialMember(j) = mutated(j);}
			}
			// Always mutate the last entry
			trialMember(function_->d_ - 1) = mutated(function_->d_ - 1);
			// Check if it is an improvement!
			double trialVal = function_->evaluate(trialMember);
			evals++;
			if(betterPoint(population[i], populationVals[i], trialMember, trialVal) == 1){
				population[i] = trialMember;
				populationVals[i] = trialVal;
			}
			if(betterPoint(population[indexBest], populationVals[indexBest], trialMember, trialVal) == 1){
				indexBest = i;
			}
			if(evals >= maxEval){break;}

		}
	}while(true);

	VectorXd bestPoint = population[indexBest];

	return bestPoint;
}







DynamicHillClimberSolver::DynamicHillClimberSolver(Function* function, bool min, int randomSeed, bool printInfo) :
	AuxSolver(function, min, randomSeed, printInfo),
	maxEval_(5000){}

DynamicHillClimberSolver::DynamicHillClimberSolver(Function* function, int maxEval, bool min, int randomSeed, bool printInfo):
	AuxSolver(function, min, randomSeed, printInfo),
	maxEval_(maxEval){}

DynamicHillClimberSolver::~DynamicHillClimberSolver(){
}

VectorXd DynamicHillClimberSolver::optimise(){
	// This is the main call, but as there is a function which takes in the starting point, 
	// this will simply choose a random initial point and call the main function
	VectorXd startingPoint = (sampleGenerator_->randomSample(1))[0];
	double startingPointVal = function_->evaluate(startingPoint);
	return optimise(startingPoint, startingPointVal);
}

VectorXd DynamicHillClimberSolver::optimise(VectorXd startingPoint, double startingPointVal){
	int d = function_->d_;
	// Initialise all of the vectors as the potential directions
	vector<VectorXd> potentialDirections = initialiseDirections();
	VectorXd currentPoint = startingPoint;
	double currentPointValue = startingPointVal;
	VectorXd bestPoint = startingPoint;
	double bestVal = startingPointVal;
	int evals = 1;
	int maxEval = maxEval_;
	vector<VectorXd> visitedPoints;
	visitedPoints.reserve(maxEval);
	visitedPoints.push_back(currentPoint);
	int previouslyChosenIndex = -1;
	int stepsCurrentClimb = 0;
	int numVisited = 1;
	do{
		if(printInfo_){
			printf("\r%sRunning DHC, evlauted %d/%d, visited %d points, best point has value %.20f at point (", prePrint_.c_str(), evals, maxEval, numVisited, bestVal);
			for(int i = 0; i < function_->d_; i++){
				printf("%.4f", bestPoint(i));
				if(i < function_->d_-1){printf(", ");}
				else{printf("), ");}
			}
			printf("current point has value %.20f at point (", currentPointValue);
			for(int i = 0; i < function_->d_; i++){
				printf("%.4f", currentPoint(i));
				if(i < function_->d_-1){printf(", ");}
				else{printf(")                         ");}
			}
		}
		if(evals >= maxEval){if(printInfo_){printf("                                            \n");}break;}
		// First thing is to find the direction with the largest norm
		int indexChosenDirection = 0;
		for(int i = 1; i < 2*(d + 1); i++){
			if(potentialDirections[i].norm() > potentialDirections[indexChosenDirection].norm()){
				indexChosenDirection = i;
			}
		}
		// If the largest direction has a norm close enough to 0, start at a new point.
		// Create a large set of candidates, and choose the one that is farthest from all the point evaluated so far
		if(potentialDirections[indexChosenDirection].norm() < TOL){
			stepsCurrentClimb = 0;
			potentialDirections = initialiseDirections();
			vector<VectorXd> newCandidatePoints = sampleGenerator_->randomLHS(d * 100);
			VectorXd chosenNewCandidate;
			double maxMinDistance = -1;
			// printf("Looking for a restart point\n");
			for(int i = 0; i < (int)newCandidatePoints.size(); i++){
				double minDistance = DBL_MAX;
				for(int j = 0; j < (int)visitedPoints.size(); j++){
					double dist = (newCandidatePoints[i] - visitedPoints[j]).norm();
					if(dist < minDistance){minDistance = dist;}
				}
				if(minDistance > maxMinDistance){
					// printf("Found better point with dist %.5f\n", maxMinDistance);
					maxMinDistance = minDistance;
					chosenNewCandidate = newCandidatePoints[i];
				}
			}
			currentPoint = chosenNewCandidate;
			currentPointValue = function_->evaluate(currentPoint);
			evals++;
			visitedPoints.push_back(currentPoint);
			if(betterPoint(currentPoint, currentPointValue, bestPoint, bestVal) == -1){
				bestPoint = currentPoint;
				bestVal = currentPointValue;
			}
			continue;

		}
		// Use this direction to create a candidate point
		VectorXd candidatePoint = currentPoint + potentialDirections[indexChosenDirection];
		// If out of bounds, count as if the point is worse, and halve the direction
		if(!function_->pointWithinBounds(candidatePoint)){
			potentialDirections[indexChosenDirection] = potentialDirections[indexChosenDirection] / 2;
			continue;
		}
		// Get value and compare with point
		double candidatePointValue = function_->evaluate(candidatePoint);
		evals++;
		// See if this is an improvement. If not, halve the vector and go back
		if(betterPoint(candidatePoint, candidatePointValue, currentPoint, currentPointValue) != -1){
			potentialDirections[indexChosenDirection] = potentialDirections[indexChosenDirection] / 2;
			continue;
		}
		// If here, the candidate point is better than the current point, move there
		stepsCurrentClimb++;
		currentPoint = candidatePoint;
		currentPointValue = candidatePointValue;
		if(betterPoint(currentPoint, currentPointValue, bestPoint, bestVal) == -1){
			bestPoint = currentPoint;
			bestVal = currentPointValue;
		}
		// If the same direction was chosen twice in a row, double the direction
		if(previouslyChosenIndex == indexChosenDirection){
			potentialDirections[indexChosenDirection] = potentialDirections[indexChosenDirection] * 2;
		}else{
			previouslyChosenIndex = indexChosenDirection;
		}
		// Store the visited point, and if this is at least the third visited point, store the two extra directon vectors
		visitedPoints.push_back(currentPoint);
		numVisited = visitedPoints.size();
		if(stepsCurrentClimb >= 2){
			potentialDirections[2*d] = visitedPoints[numVisited - 1] - visitedPoints[numVisited - 3];
			potentialDirections[2*d] = -(visitedPoints[numVisited - 1] - visitedPoints[numVisited - 3]);
		}

	}while(true);

	return bestPoint;
}

vector<VectorXd> DynamicHillClimberSolver::initialiseDirections(){
	int d = function_->d_;
	vector<VectorXd> potentialDirections(2*(d + 1));
	for(int i = 0; i < d; i++){
		VectorXd positiveVector = VectorXd::Zero(d);
		VectorXd negativeVector = VectorXd::Zero(d);
		positiveVector(i) = (function_->upperBound_[i] - function_->lowerBound_[i]) / 2;
		negativeVector(i) = -(function_->upperBound_[i] - function_->lowerBound_[i]) / 2;
		potentialDirections[i] = positiveVector;
		potentialDirections[d + i] = negativeVector;
	}
	potentialDirections[2*d] = VectorXd::Zero(d);
	potentialDirections[2*d + 1] = VectorXd::Zero(d);
	return potentialDirections;
}









HeavyTuneSolver::HeavyTuneSolver(Function* function, bool min, int randomSeed, bool printInfo) :
	AuxSolver(function, min, randomSeed, printInfo),
	maxEval_(10000),
	popSizeGA_(50){}

HeavyTuneSolver::HeavyTuneSolver(Function* function, int maxEval, int popSizeGA, bool min, int randomSeed, bool printInfo):
	AuxSolver(function, min, randomSeed, printInfo),
	maxEval_(maxEval),
	popSizeGA_(popSizeGA){}

HeavyTuneSolver::~HeavyTuneSolver(){
}

VectorXd HeavyTuneSolver::optimise(){
	int d = function_->d_;
	int usedBudget = 0;
	// First run the simplex 
	AmoebaSolver* amoeba = new AmoebaSolver(function_, 50*d, min_, randomSeed_, printInfo_);
	amoeba->prePrint_ = "Running Heavy Tune - amoeaba first pass:  ";
	VectorXd currentBestPoint = amoeba->optimise();
	usedBudget += amoeba->evalsUsed_;
	// Now perturb the best point by up to 10% of the range of each dimension to create a new amoeba run
	vector<VectorXd> secondPassAmoebaPoints;
	vector<double> secondPassAmoebaVals;
	uniform_real_distribution<> perturbationDist(-0.1, 0.1);
	for(int i = 0; i < d + 1; i++){
		VectorXd perturbedPoint;
		do{
			perturbedPoint = currentBestPoint;
			for(int j = 0; j < d; j++){
				perturbedPoint(j) = perturbedPoint(j) + perturbationDist(randomGenerator_) * (function_->upperBound_[j] - function_->lowerBound_[j]);
			}
		}while(!function_->pointWithinBounds(perturbedPoint));
		secondPassAmoebaPoints.push_back(perturbedPoint);
		secondPassAmoebaVals.push_back(function_->evaluate(perturbedPoint));
	}
	amoeba->prePrint_ = "Running Heavy Tune - amoeaba second pass: ";
	currentBestPoint = amoeba->optimise(secondPassAmoebaPoints, secondPassAmoebaVals);
	usedBudget += amoeba->evalsUsed_;
	delete amoeba;

	// Work out how to divide the rest of the sample budget
	int gaBudget, dhcBudget;
	// If will do a local search with the tie breaker, must allocate less budget to dinamic hill climber, otherwise simply divide by two
	gaBudget = floor((maxEval_ - usedBudget)/2.0);
	dhcBudget = ceil((maxEval_ - usedBudget)/2.0);

	// Use this current best point as part of the population in the genetic algorithm
	GASolver* gaSolver = new GASolver(function_, gaBudget, popSizeGA_, min_, randomSeed_, printInfo_);
	gaSolver->prePrint_ = "Running Heavy Tune - GA pass:             ";
	vector<VectorXd> population;
	population = sampleGenerator_->randomSample(popSizeGA_);
	
	vector<double> populationVals;
	for(int i = 0; i < popSizeGA_; i++){
		populationVals.push_back(function_->evaluate(population[i]));
	}
	currentBestPoint = gaSolver->optimise(population, populationVals);
	delete gaSolver;

	// Finally, use this point as the starting point for the hill climber algorithm
	DynamicHillClimberSolver* dhsSolver = new DynamicHillClimberSolver(function_, dhcBudget, min_, randomSeed_, printInfo_);
	dhsSolver->prePrint_ = "Running Heavy Tune - DHS pass:            ";
	currentBestPoint = dhsSolver->optimise(currentBestPoint, function_->evaluate(currentBestPoint));
	delete dhsSolver;

	
	return currentBestPoint;
}






RandomHeavyTuneSolver::RandomHeavyTuneSolver(Function* function, bool min, int randomSeed, bool printInfo) :
	AuxSolver(function, min, randomSeed, printInfo),
	maxEval_(10000),
	numSphereARS_(10){}

RandomHeavyTuneSolver::RandomHeavyTuneSolver(Function* function, int maxEval, int numSphereARS, bool min, int randomSeed, bool printInfo) :
	AuxSolver(function, min, randomSeed, printInfo),
	maxEval_(maxEval),
	numSphereARS_(numSphereARS){}

RandomHeavyTuneSolver::~RandomHeavyTuneSolver(){
}

VectorXd RandomHeavyTuneSolver::optimise(){
	int d = function_->d_;
	int usedBudget = 0;
	// First run the simplex 
	AmoebaSolver* amoeba = new AmoebaSolver(function_, 50*d, min_, randomSeed_, printInfo_);
	amoeba->prePrint_ = "Running Random Heavy Tune - amoeaba first pass:  ";
	VectorXd currentBestPoint = amoeba->optimise();
	usedBudget += amoeba->evalsUsed_;
	// Now perturb the best point by up to 10% of the range of each dimension to create a new amoeba run
	vector<VectorXd> secondPassAmoebaPoints;
	vector<double> secondPassAmoebaVals;
	uniform_real_distribution<> perturbationDist(-0.1, 0.1);
	for(int i = 0; i < d + 1; i++){
		VectorXd perturbedPoint;
		do{
			perturbedPoint = currentBestPoint;
			for(int j = 0; j < d; j++){
				perturbedPoint(j) = perturbedPoint(j) + perturbationDist(randomGenerator_) * (function_->upperBound_[j] - function_->lowerBound_[j]);
			}
		}while(!function_->pointWithinBounds(perturbedPoint));
		secondPassAmoebaPoints.push_back(perturbedPoint);
		secondPassAmoebaVals.push_back(function_->evaluate(perturbedPoint));
	}
	amoeba->prePrint_ = "Running Random Heavy Tune - amoeaba second pass: ";
	currentBestPoint = amoeba->optimise(secondPassAmoebaPoints, secondPassAmoebaVals);
	usedBudget += amoeba->evalsUsed_;
	delete amoeba;

	// Work out how to divide the rest of the sample budget
	int arsBudget, dhcBudget;
	// If will do a local search with the tie breaker, must allocate less budget to dinamic hill climber, otherwise simply divide by two
	arsBudget = floor((maxEval_ - usedBudget)*0.8);
	dhcBudget = ceil((maxEval_ - usedBudget)*0.2);

	// Use this current best point as part of the population in the genetic algorithm
	ARSsolver* arsSolver = new ARSsolver(function_, numSphereARS_, arsBudget, min_, randomSeed_, printInfo_);
	arsSolver->prePrint_ = "Running Random Heavy Tune - ARS pass:            ";

	// Get the centers, and replace one with the current best point
	vector<VectorXd> centers = sampleGenerator_->randomLHS(numSphereARS_);
	centers[0] = currentBestPoint;
	currentBestPoint = arsSolver->optimise(centers);
	delete arsSolver;

	// Finally, use this point as the starting point for the hill climber algorithm
	DynamicHillClimberSolver* dhsSolver = new DynamicHillClimberSolver(function_, dhcBudget, min_, randomSeed_, printInfo_);
	dhsSolver->prePrint_ = "Running Random Heavy Tune - DHS pass:            ";
	currentBestPoint = dhsSolver->optimise(currentBestPoint, function_->evaluate(currentBestPoint));
	delete dhsSolver;
	return currentBestPoint;
}














// CMAESsolver::CMAESsolver(Function* function, bool min, int randomSeed, bool printInfo, bool useTieBreaker, string tieBreakerFlag, int finishWithLocalSearch, int additionalEvaluationsInLocalSearch) :
// 	AuxSolver(function, min, randomSeed, printInfo, useTieBreaker, tieBreakerFlag, finishWithLocalSearch, additionalEvaluationsInLocalSearch),
// 	maxEval_(10000){
// }

// CMAESsolver::CMAESsolver(Function* function, int maxEval, bool min, int randomSeed, bool printInfo, bool useTieBreaker, string tieBreakerFlag, int finishWithLocalSearch, int additionalEvaluationsInLocalSearch) :
// 	AuxSolver(function, min, randomSeed, printInfo, useTieBreaker, tieBreakerFlag, finishWithLocalSearch, additionalEvaluationsInLocalSearch),
// 	maxEval_(maxEval){
// }



// CMAESsolver::~CMAESsolver(){
// }

// VectorXd CMAESsolver::optimise(){
// 	cmaes_t evo; /* an CMA-ES type struct or "object" */
// 	cmaes_boundary_transformation_t boundaries;
// 	double *arFunvals, *x_in_bounds, *const*pop;
// 	int nb_bounds = function_->d_;
// 	double lowerBounds[nb_bounds];
// 	double upperBounds[nb_bounds];
// 	copy(function_->lowerBound_.begin(), function_->lowerBound_.end(), lowerBounds);
// 	copy(function_->upperBound_.begin(), function_->upperBound_.end(), upperBounds);
// 	unsigned long dimension;
// 	int i; 
// 	/* initialize boundaries, be sure that initialSigma is smaller than upper minus lower bound */
// 	cmaes_boundary_transformation_init(&boundaries, lowerBounds, upperBounds, nb_bounds);
	
// 	// Initialize starting point to the middle of the space and standard deviation to 1/4 of the range of that dimension
// 	double xStart[nb_bounds];
// 	double stds[nb_bounds];
// 	for(int i = 0; i < nb_bounds; i++){
// 		xStart[i] = function_->lowerBound_[i] + (function_->upperBound_[i] - function_->lowerBound_[i])/2.0;
// 		stds[i] = (function_->upperBound_[i] - function_->lowerBound_[i])/4.0;
// 	}

// 	// Keep track of of the points visited, and how many iterations we have gone through
// 	// Do this in case of need to reinitiallisation
// 	vector<VectorXd> visitedPoints;
// 	VectorXd bestPoint;
// 	// Will always work as if minimising; don't care about the actual best value so much as its location
// 	double bestPointValue;

// 	int evals = 1;
// 	int maxEval = maxEval_;
// 	if(finishWithLocalSearch_ == 1 || finishWithLocalSearch_ == 2){
// 		maxEval = maxEval - additionalEvaluationsInLocalSearch_;
// 	}

// 	do{	
// 		arFunvals = cmaes_init(&evo, function_->d_, xStart, stds, randomSeed_, 0, "");
// 		dimension = (unsigned long)cmaes_Get(&evo, "dimension");
// 		if(printInfo_){printf("%s\n", cmaes_SayHello(&evo));}
// 		x_in_bounds = cmaes_NewDouble(dimension); /* calloc another vector */
// 		evo.sp.stopMaxFunEvals = maxEval - evals;
// 		VectorXd currPoint;
// 		while(!cmaes_TestForTermination(&evo)){ 
// 			// printf("Iteration %d\n", completedIterations);
// 			pop = cmaes_SamplePopulation(&evo);
// 			for (i = 0; i < cmaes_Get(&evo, "lambda"); ++i){
// 				/* transform into bounds and evaluate the new search points */
// 				cmaes_boundary_transformation(&boundaries, pop[i], x_in_bounds, dimension);
// 				// Working with x_in_bounds, map it to a VectorXd
// 	        	currPoint = Map<VectorXd>(x_in_bounds, function_->d_);
// 	        	/* this loop can be omitted if is_feasible is invariably true */
// 		        while(!function_->pointWithinBounds(currPoint)) { /* is_feasible needs to be user-defined, in case, and can change/repair x */
// 		            cmaes_ReSampleSingle(&evo, i); 
// 		            cmaes_boundary_transformation(&boundaries, pop[i], x_in_bounds, dimension);
// 		            currPoint = Map<VectorXd>(x_in_bounds, function_->d_);
// 		        }
// 		        arFunvals[i] = round(arFunvals[i] / TOL) * TOL;
// 		        double val = function_->evaluate(currPoint);
// 		        // Looks like floating point error might be a thing here, round objective function value to defined TOL
// 		        val = round(val / TOL) * TOL;
// 		        // printf("(");
// 				// for(int i = 0; i < currPoint.size(); i++){
// 				// 	printf("%.10f",currPoint(i));
// 				// 	if(i < currPoint.size()-1){printf(", ");}
// 				// 	else{printf(")");}
// 				// }
// 		        // printf(": %.10f\n", val);

// 		        if(min_){arFunvals[i] = val;}
// 	        	else{arFunvals[i] = -1 * val;}

// 	        	// Keep track of what has happened
// 	        	visitedPoints.push_back(currPoint);
// 	        	evals++;
// 	        	// If found a better point, remember it
// 	        	if(betterPoint(currPoint, val, bestPoint, bestPointValue) == -1){
// 	        		bestPoint = currPoint;
// 	        		bestPointValue = val;
// 	        	}
// 	        }
// 	        // printf("STDs: ");
// 	        // for(int i = 0; i < nb_bounds; i++){
// 	        // 	printf("%.5f ", stds[i]);
// 			// }
// 			// printf("\n");
// 			// printf("Best so far:\n");
// 			// printf("(");
// 			// for(int i = 0; i < bestPoint.size(); i++){
// 			// 	printf("%.10f",bestPoint(i));
// 			// 	if(i < bestPoint.size()-1){printf(", ");}
// 			// 	else{printf(")");}
// 			// }
// 	        // printf(": %.10f\n", bestPointValue);

// 	        cmaes_UpdateDistribution(&evo, arFunvals);
// 	        /* read instructions for printing output or changing termination conditions */ 
// 		    if(printInfo_){cmaes_ReadSignals(&evo, "cmaes/cmaes_signals.par");}
// 		    if(printInfo_){fflush(stdout); /* useful in MinGW */}
// 		}
// 		if(printInfo_){printf("Stop: %s",  cmaes_TestForTermination(&evo)); /* print termination reason */}
// 		cmaes_WriteToFile(&evo, "all", "cmaes/allcmaes.dat");         /* write final results */		

// 		// Release memory
// 		cmaes_exit(&evo); /* release memory */
// 		free(x_in_bounds);

// 		if(evals < maxEval){
// 			if(printInfo_){
// 				printf("Best point found so far has value %.4f at point ", bestPointValue);
// 				printPoint(bestPoint);
// 				printf(", still have %d/%d samples left so continuing!\n", maxEval - evals, maxEval);
// 			}
// 			// Need to choose the new starting point!
// 			vector<VectorXd> newCandidatePoints = sampleGenerator_->randomLHS(function_->d_ * 100);
// 			VectorXd chosenNewCandidate;
// 			double maxMinDistance = -1;
// 			for(int i = 0; i < (int)newCandidatePoints.size(); i++){
// 				double minDistance = DBL_MAX;
// 				for(int j = 0; j < (int)visitedPoints.size(); j++){
// 					double dist = (newCandidatePoints[i] - visitedPoints[j]).norm();
// 					if(dist < minDistance){minDistance = dist;}
// 				}
// 				if(minDistance > maxMinDistance){
// 					maxMinDistance = minDistance;
// 					chosenNewCandidate = newCandidatePoints[i];
// 				}
// 			}
// 			for(int i = 0; i < function_->d_; i++){
// 				xStart[i] = chosenNewCandidate(i);
// 			}

// 		}else{
// 			if(printInfo_){
// 				printf("Best point found has value %.4f at point ", bestPointValue);
// 				printPoint(bestPoint);
// 				printf("\n");
// 			}
// 			break;
// 		}
// 	}while(true);

// 	cmaes_boundary_transformation_exit(&boundaries); /* release memory */
	
// 	if(finishWithLocalSearch_ > 0){
// 		DynamicHillClimberSolver* localSolver;
// 		if(finishWithLocalSearch_ == 1 || finishWithLocalSearch_ == 3){
// 			localSolver = new DynamicHillClimberSolver(function_, additionalEvaluationsInLocalSearch_, min_, randomSeed_, printInfo_, false, tieBreakerFlag_);
// 			localSolver->prePrint_ = "Running CMAES - local search without tie breaker: ";
// 		}else{
// 			localSolver = new DynamicHillClimberSolver(function_, additionalEvaluationsInLocalSearch_, min_, randomSeed_, printInfo_, true, tieBreakerFlag_);
// 			localSolver->prePrint_ = "Running CMAES - local search with tie breaker: ";
// 		}
// 		bestPoint = localSolver->optimise(bestPoint, bestPointValue);
// 		delete localSolver;
// 	}

// 	return bestPoint;
// }











#endif