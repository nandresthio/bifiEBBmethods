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

#ifndef SURROGATE_MODELS_CPP
#define SURROGATE_MODELS_CPP


#include "surrogate_models.hpp"



SurrogateModel::SurrogateModel(BiFidelityFunction* biFunction, AuxSolver* auxSolver, int randomSeed, bool printInfo, bool functionScaling):
	biFunction_(biFunction),
	auxSolver_(auxSolver),
	randomSeed_(randomSeed),
	printInfo_(printInfo),
	functionScaling_(functionScaling){
	if(randomSeed != 0){
		mt19937 gen(randomSeed);
		randomGenerator_ = gen;
	}else{
		random_device rd;
		mt19937 gen(rd());
		randomGenerator_ = gen;
	}
	auxSolver_->reseedRandom(randomSeed_);
	sampleGenerator_ = new SampleGenerator(biFunction_, randomSeed_, printInfo_);
	trainedModel_ = false;

	maxObservation_ = DBL_MAX;
	minObservation_ = DBL_MAX;

	chosenAcquisiton_ = "";

	costRatio_ = -1;

	startIntervalTraining_ = 100;
	intervalForTraining_ = 10;
	startLargeIntervalTraining_ = 200;
	intervalForLargeTraining_ = 50;
	startVeryLargeIntervalTraining_ = 500;
	intervalForVeryLargeTraining_ = 100;
	sampleSizeAtLastTrain_ = 0;
}

SurrogateModel::~SurrogateModel(){
	if(sampleGenerator_ != NULL){delete sampleGenerator_;}
}

void SurrogateModel::setCostRatio(double costRatio){
	if(abs(costRatio) > TOL &&
		abs(costRatio - 1) > TOL &&
		(costRatio < 0 || costRatio_ > 1)){
		printf("Error: Trying to set a cost ratio which is outside of the range [0,1], not saving this value.");
		return;
	}
	costRatio_ = costRatio;
}


void SurrogateModel::scalePoint(VectorXd &point){
	if(functionScaling_){
		for(int i = 0; i < (int)point.size(); i++){
			point(i) = 1*(point(i) - biFunction_->lowerBound_[i]) / (biFunction_->upperBound_[i] - biFunction_->lowerBound_[i]);
		}
	}
}

void SurrogateModel::scalePoints(vector<VectorXd> &points){
	for(int i = 0; i < (int)points.size(); i++){
		scalePoint(points[i]);
	}
}

void SurrogateModel::unscalePoint(VectorXd &point){
	if(functionScaling_){
		for(int i = 0; i < (int)point.size(); i++){
			point(i) = 1*point(i) * (biFunction_->upperBound_[i] - biFunction_->lowerBound_[i]) + biFunction_->lowerBound_[i];
		}
	}
}

void SurrogateModel::unscalePoints(vector<VectorXd> &points){
	for(int i = 0; i < (int)points.size(); i++){
		unscalePoint(points[i]);
	}
}

double SurrogateModel::scaleObservation(double value){
	if(!functionScaling_){return value;}
	if(maxObservation_ == DBL_MAX || minObservation_ == DBL_MAX){
		printf("Asking to scale observations without having had the model find the fmin and fmax first! Stopping now...\n");
		exit(0);
	}
	if(abs(maxObservation_ - minObservation_) < TOL){return value - minObservation_;}
	else{return 1*(value - minObservation_) / (maxObservation_ - minObservation_);}
}

void SurrogateModel::scaleObservations(vector<double> &observations){
	for(int i = 0; functionScaling_ && i < (int)observations.size(); i++){
		observations[i] = scaleObservation(observations[i]);
	}
}

double SurrogateModel::unscaleObservation(double value){
	if(!functionScaling_){return value;}
	if(maxObservation_ == DBL_MAX || minObservation_ == DBL_MAX){
		printf("Asking to unscale observations without having had the model scale them first! Stopping now...\n");
		exit(0);
	}
	if(abs(maxObservation_ - minObservation_) < TOL){return value + minObservation_;}
	else{return 1 * value * (maxObservation_ - minObservation_) + minObservation_;}
}

void SurrogateModel::unscaleObservations(vector<double> &observations){
	for(int i = 0; functionScaling_ && i < (int)observations.size(); i++){
		observations[i] = unscaleObservation(observations[i]);
	}
}

void SurrogateModel::scaleTwoObservations(vector<double> &observations, vector<double> &observationsLow){
	maxObservation_ = *max_element(observations.begin(), observations.end());
	minObservation_ = *min_element(observations.begin(), observations.end());
	scaleObservations(observations);
	scaleObservations(observationsLow);
}

void SurrogateModel::saveSample(vector<VectorXd> points, vector<VectorXd> pointsLow, vector<double> observations, vector<double> observationsLow){
	scalePoints(points);
	scalePoints(pointsLow);
	scaleTwoObservations(observations, observationsLow);
	sampledPoints_ = points;
	sampledPointsLow_ = pointsLow;
	sampledPointsValues_ = observations;
	sampledPointsValuesLow_ = observationsLow;
}

void SurrogateModel::addSample(VectorXd point, bool sampleHigh, bool sampleLow){
	// First evaluate each of the sources if relevant
	if(!sampleHigh && !sampleLow){
		printf("Note: calling this addSample method is doing nothing as both sampleHigh and sampleLow are set as false!\n");
	}
	// Just unscale everything, and save it again!
	unscalePoints(sampledPoints_);
	unscalePoints(sampledPointsLow_);
	unscaleObservations(sampledPointsValues_);
	unscaleObservations(sampledPointsValuesLow_);

	if(sampleHigh){
		double valHigh = biFunction_->evaluate(point);
		sampledPoints_.push_back(point);
		sampledPointsValues_.push_back(valHigh);
	}

	if(sampleLow){
		double valLow = biFunction_->evaluateLow(point);
		sampledPointsLow_.push_back(point);
		sampledPointsValuesLow_.push_back(valLow);
	}

	// And now save it again!
	saveSample(sampledPoints_, sampledPointsLow_, sampledPointsValues_, sampledPointsValuesLow_);
}



void SurrogateModel::trainModel(bool forcedOptimiseHyperparameters){
	printf("Should not get to this train model call! Exiting now...\n");
	exit(0);
}

vector<double> SurrogateModel::multipleSurfaceValues(vector<VectorXd> &points, bool pointIsScaled, bool unscaleOutput){
	vector<double> values((int)points.size(), 0.0);
	for(int i = 0; i < (int)points.size(); i++){
		values[i] = surfaceValue(points[i], pointIsScaled, unscaleOutput);
	}
	return values;
}

double SurrogateModel::surfaceValue(VectorXd &x, bool pointIsScaled, bool unscaleOutput){
	printf("Should not get to this surface value call! Exiting now...\n");
	exit(0);
}


double SurrogateModel::evaluateAcquisitionFunction(VectorXd x){
	if(chosenAcquisiton_.compare("surface") == 0){
		return surfaceValue(x);
	}else if(chosenAcquisiton_.compare("") == 0){
		printf("Error, asking to evaluate acquisition function without first definition what function to use! Stopping now...\n");
		exit(0);
	}else{
		printf("Error, asking to evaluate acquisition function with the unkown acquisition %s! Stopping now...\n", chosenAcquisiton_.c_str());
		exit(0);
	}
}

void SurrogateModel::setAcquisitionFunction(string chosenAcquisiton){
	if(chosenAcquisiton.compare("surface") == 0){
		chosenAcquisiton_ = chosenAcquisiton;
		acquisitionIsMin_ = true;
	}else{
		printf("Error, unkown acquisiton function %s specified! Stopping now...\n", chosenAcquisiton.c_str());
		exit(0);
	}
}

tuple<VectorXd, double, bool, bool> SurrogateModel::findNextSampleSite(){
	if(printInfo_){printf("Find next sample site with acquisition function %s.\n", chosenAcquisiton_.c_str());}
	// Ok, so idea is simply to initialise solver and run it!
	AcquisitionFunction* function = new AcquisitionFunction(this);
	auxSolver_->updateProblem(function, acquisitionIsMin_);
	VectorXd nextSampleSite = auxSolver_->optimise();
	double value = evaluateAcquisitionFunction(nextSampleSite);
	delete function;
	// Return false false indicating neither source should be sampled. Passing this on to add sample should
	// cause an error; these bools should be overriden by latter calls. 
	return make_tuple(nextSampleSite, value, false, false);
}


SurrogateModel::AcquisitionFunction::AcquisitionFunction(SurrogateModel* surrogateModel):
	Function(surrogateModel->biFunction_->d_, surrogateModel->biFunction_->lowerBound_, surrogateModel->biFunction_->upperBound_),
	surrogateModel_(surrogateModel){}

SurrogateModel::AcquisitionFunction::~AcquisitionFunction(){}

double SurrogateModel::AcquisitionFunction::evaluate(VectorXd &point){
	// Simply call the evaluate acquisition function, and return that
	return surrogateModel_->evaluateAcquisitionFunction(point);
}







Kriging::Kriging(BiFidelityFunction* biFunction, AuxSolver* auxSolver, int randomSeed, bool printInfo, bool functionScaling):
	SurrogateModel(biFunction, auxSolver, randomSeed, printInfo, functionScaling){

	theta_.reserve(biFunction_->d_);
	pVector_.reserve(biFunction_->d_);
	
}

Kriging::~Kriging(){}


void Kriging::trainModel(bool forcedOptimiseHyperparameters){
	int n = (int)sampledPoints_.size();
	if(printInfo_){printf("Training kriging model with %d samples - ", n);}
	if(forcedOptimiseHyperparameters){
		if(printInfo_){printf("optimising hyperparameters as user requested to regardless of other factors.\n");}
		trainHyperparameters();
		sampleSizeAtLastTrain_ = n;
	
	}else if(n == sampleSizeAtLastTrain_){
		if(printInfo_){printf("skipping  hyperparameter optimisation as sample size has not changed.\n");}
	
	}else if(startIntervalTraining_ >= n){
		if(printInfo_){printf("optimising hyperparameters as working with manageable sample size.\n");}
		trainHyperparameters();
		sampleSizeAtLastTrain_ = n;

	}else if((startLargeIntervalTraining_ > n && (n - sampleSizeAtLastTrain_) >= intervalForTraining_) ||
				(startVeryLargeIntervalTraining_ > n && (n - sampleSizeAtLastTrain_) >= intervalForLargeTraining_) ||
				(n - sampleSizeAtLastTrain_) >= intervalForVeryLargeTraining_){
		if(printInfo_){printf("optimising hyperparameters despite having large sample as it's been %d samples since last optimisation.\n", n - sampleSizeAtLastTrain_);}
		trainHyperparameters();
		sampleSizeAtLastTrain_ = n;
		
	}else{
		if(printInfo_){printf("skipping  hyperparameter optimisation.\n");}
	}
	
	saveMuSigma();	
	trainedModel_ = true;
	
	// It is possible that the predictions are now undefined,
	// if the hyperparameters where not optimised.
	// This can happen if there are lots of points.
	// Check that a random prediction can be given; if not call the function again
	// and force the training of hyperparameters
	// Get a random point
	VectorXd point = (sampleGenerator_->randomSample(1))[0];
	double val = surfaceValue(point);
	if(isnan(val)){
		if(printInfo_){printf("Skipped hyperparameter optimisation but now predictions give NaN. Forcing retuning of hyperparameters.\n");}
		trainModel(true);
	}
	return;
	
	return;
}

void Kriging::trainHyperparameters(){
	int d = biFunction_->d_;
	// First define the bounds for the function
	vector<double> lowerBound(2*d, 0.0);
	vector<double> upperBound(2*d, 0.0);
	// Bounds should be dependent on whether data has been scaled; if it has the bounds can be more restrictive
	for(int i = 0; i < d; i++){
		if(functionScaling_){
			lowerBound[i] = -3;
			// Increased bounds to 5 as it is possible to have so much data that for smaller values this is undefined
			upperBound[i] = 5;
		}else{
			lowerBound[i] = -10;
			upperBound[i] = 3;
		}
		lowerBound[d+i] = 0.1;
		upperBound[d+i] = 2.0;
	}
	maxLogLikelihood_ = -DBL_MAX;
	ConcentratedLikelihoodFunction* function;
	if(fixedP_){
		// If using fixed p hyperparameters, only need to optimise the thetas
		function = new ConcentratedLikelihoodFunction(d, lowerBound, upperBound, this);
	}else{
		function = new ConcentratedLikelihoodFunction(2*d, lowerBound, upperBound, this);
	}
	auxSolver_->updateProblem(function, false);
	VectorXd hyperparameters = auxSolver_->optimise();
	delete function;
	// Extract two vectors, and store!
	for(int i = 0; i < d; i++){
		theta_[i] = pow(10,hyperparameters(i));
		if(fixedP_){
			pVector_[i] = 2;
		}else{
			pVector_[i] = hyperparameters(d + i);
		}
	}
}

void Kriging::saveMuSigma(){
	auto data = muSigmaCalculator();
	mu_ = get<0>(data);
	sigma_ = get<1>(data);
	rMatrix_ = get<2>(data);
	rMatrixDecomposition_ = get<3>(data);
	modelWeights_ = get<4>(data);

	// rMatrixInverse_ = rMatrix_.inverse();
	int n = rMatrix_.cols();
	// I think that since I only find the inverse once here, I want to 
	// do a decomposition once again. Using LU decomposition from Eigen,
	// as it appears it is the most stable and accurate when finding an inverse.
	// Note this should work even when the matrix is non invertible, which is a huge plus,
	// even though this should not happen.
	rMatrixInverse_ = rMatrix_.fullPivLu().solve(MatrixXd::Identity(n,n));
	return;
}


tuple<double, double, MatrixXd, LDLT<MatrixXd>, MatrixXd> Kriging::muSigmaCalculator(){
	int n = (int)sampledPoints_.size();
	int d = (int)sampledPoints_[0].size();

	if(n < 1){
		printf("Trying to train Kriging model with less points than needed! Have %d but need %d. Exiting now...\n", n, 1);
		exit(0);
	}
	MatrixXd rMatrix(n,n);
	// Define matrix
	auto it1 = sampledPoints_.begin();
	for(int i = 0; i < n; i++, it1++){
		auto it2 = it1;
		for(int j = i; j < n; j++, it2++){
			double sum = 0;
			for(int k = 0; k < d; k++){
				sum += theta_[k] * pow(abs((*it1)(k) - (*it2)(k)), pVector_[k]);
			}
			rMatrix(i,j) = exp(-sum);
			rMatrix(j,i) = rMatrix(i,j);
		}
	}
	VectorXd one(n);
	VectorXd y(n);
	for(int i = 0; i < n; i++){
		one(i) = 1;
		y(i) = sampledPointsValues_[i];
	}
	// Save matrix decomposition
	LDLT<MatrixXd> rMatrixDecomposition = rMatrix.ldlt();
	MatrixXd mu_top = one.transpose() * rMatrixDecomposition.solve(y);
	MatrixXd mu_bottom = one.transpose() * rMatrixDecomposition.solve(one);
	double mu = mu_top(0,0)/ mu_bottom(0,0);
	MatrixXd modelWeights =  rMatrixDecomposition.solve(y - one*mu);
	MatrixXd sigma_m = (y - one*mu).transpose() * modelWeights;
	double sigma = sigma_m(0,0)/n;
	
	// Deal with special case when n = 1
	if(n == 1){
		mu = y(0);
		sigma = 0;
	}
	// return make_tuple(mu, sigma, rMatrix, rMatrixInverse, rMatrixDecomposition);
	return make_tuple(mu, sigma, rMatrix, rMatrixDecomposition, modelWeights);
}

double Kriging::surfaceValue(VectorXd &x, bool pointIsScaled, bool unscaleOutput){
	VectorXd xCopy = x;
	if(!pointIsScaled){scalePoint(xCopy);}
	auto data = meanVarianceCalculator(xCopy);
	double s_x = get<0>(data);
	if(unscaleOutput){s_x = unscaleObservation(s_x);}
	return s_x;
}

double Kriging::uncertainty(VectorXd &x, bool pointIsScaled, bool unscaleOutput){
	VectorXd xCopy = x;
	if(!pointIsScaled){scalePoint(xCopy);}
	auto data = meanVarianceCalculator(xCopy);
	double var = get<1>(data);
	// Not 100% sure, but for now just "unscale" it once
	// Unscaling it once is incorrect, as a variance of 0 will be displaced. Therefore only multiply by the original range
	// twice
	if(unscaleOutput){
		var = (maxObservation_ - minObservation_) * (maxObservation_ - minObservation_) * var;
	}
	return var;
}

double Kriging::expectedImprovement(VectorXd &x, bool pointIsScaled, bool unscaleOutput){
	VectorXd xCopy = x;
	if(!pointIsScaled){scalePoint(xCopy);}
	auto data = meanVarianceCalculator(xCopy);
	double s_x = get<0>(data);
	double var = get<1>(data);
	double value;
	if(unscaleOutput){
		s_x = unscaleObservation(s_x);
		var = (maxObservation_ - minObservation_) * (maxObservation_ - minObservation_) * var;
	}
	// If variance is close enough to 0, expected improvement is set to 0
	if(abs(var) < TOL){return 0;}
	// Need to decide whether we are working with a scaled model or not.
	// If yes, the min value is 0 (as it is scaled).
	// Otherwise need to use the true minimum
	if(unscaleOutput || !functionScaling_){
		value = (minObservation_ - s_x)/sqrt(var);
	}else{
		value = (0 - s_x)/sqrt(var);
	}
	double improvement = sqrt(var)*(value * normalCDF(value) + normalPDF(value));
	return improvement;
}

double Kriging::kernelProductIntegral(double x1, double x2, double theta1, double theta2){
	return ( 1 / ( 2 * sqrt(theta1 + theta2) ) ) * sqrt(M_PI) * exp( -1 * ( (x1 - x2) * (x1 - x2) * theta1 * theta2 ) / ( theta1 + theta2 ) ) * ( erf( ( theta1 - x1 * theta1 + theta2 - x2 * theta2 ) / ( sqrt(theta1 + theta2) ) ) + erf( ( x1 * theta1 + x2 * theta2 ) / sqrt(theta1 + theta2) ) );
}

void Kriging::setWmatrix(){
	int n = (int)sampledPoints_.size();
	int d = (int)sampledPoints_[0].size();
	MatrixXd wMatrix(n,n);
	// Define matrix
	auto it1 = sampledPoints_.begin();
	for(int i = 0; i < n; i++, it1++){
		auto it2 = it1;
		for(int j = i; j < n; j++, it2++){
			double product = 1;
			for(int k = 0; k < d; k++){
				product = product * 0.25 * sqrt(2 * M_PI * (1.0/theta_[k])) * exp(-1 * pow((*it1)(k) - (*it2)(k), 2) / (2 * (1.0/theta_[k]))) * (erf( (2 - ((*it1)(k) + (*it2)(k)) )/sqrt(2 * (1.0/theta_[k])) ) + erf( ((*it1)(k) + (*it2)(k))/sqrt(2 * (1.0/theta_[k])) ) ); 
			}
			wMatrix(i,j) = product;
			wMatrix(j,i) = product;
		}
	}
	wMatrix_ = wMatrix;
}

double Kriging::globalVarianceReduction(VectorXd &x, bool pointIsScaled, bool unscaleOutput){
	// Doing the approach detailed in "Replication or Exploration? Sequential Design for Stochastic Simulation Experiments"
	// Gives an O(n^2) value of reduction in global variance but is restricted to having the p hyperparameters equal to 2.
	VectorXd xCopy = x;
	if(!pointIsScaled){scalePoint(xCopy);}
	double value;
	int n = (int)sampledPoints_.size();
	int d = (int)sampledPoints_[0].size();
	// First want r vector
	VectorXd r(n);
	VectorXd w(n);
	auto it = sampledPoints_.begin();
	for(int i = 0; i < n; i++, it++){
		double sum = 0;
		double product = 1;
		for(int k = 0; k < d; k++){
			sum += theta_[k] * pow(abs((*it)(k) - xCopy(k)), pVector_[k]);
			product = product * kernelProductIntegral((*it)(k), xCopy(k), theta_[k], theta_[k]);
		}
		r(i) = exp(-sum);
		w(i) = product;
	}
	double wSelfProduct = 1;
	for(int k = 0; k < d; k++){
		wSelfProduct = wSelfProduct * kernelProductIntegral(xCopy(k), xCopy(k), theta_[k], theta_[k]); 
	}
	// This calculation is done in O(n^2)
	auto weights = rMatrixInverse_ * r;
	double val1 = (weights.transpose()  * wMatrix_ * weights);
	double val2 = 2 * (w.transpose() * weights)(0);
	
	value = sigma_ * (val1 - val2 + wSelfProduct) / (1 - r.transpose() * weights);
	if(unscaleOutput){
		value = (maxObservation_ - minObservation_) * (maxObservation_ - minObservation_) * value;
	}
	return value;
}



double Kriging::globalVarianceReductionApproximation(VectorXd &x, bool pointIsScaled, bool unscaleOutput){
	// Trying the approach given in "Gaussian Process Regression: Active Data Selection and Test Point Rejection"
	// First calculate the determinant
	VectorXd xCopy = x;
	if(!pointIsScaled){scalePoint(xCopy);}
	double denumerator;
	int n = (int)sampledPoints_.size();
	int d = (int)sampledPoints_[0].size();
	// First want r vector
	VectorXd r(n);
	auto it = sampledPoints_.begin();
	for(int i = 0; i < n; i++, it++){
		double sum = 0;
		for(int k = 0; k < d; k++){
			sum += theta_[k] * pow(abs((*it)(k) - xCopy(k)), pVector_[k]);
		}
		r(i) = exp(-sum);
	}
	auto weights = rMatrixInverse_ * r;
	// This will cause problems if the point is very close to an already sampled point. If so, do not use it
	denumerator = 1 - r.transpose() * weights;
	// if(denumerator < 0.001){return 0;}
	double sumNumerator = 0;
	for(int j = 0;  j < (int)varianceTestPoints_.size(); j++){
		VectorXd rOld(n);
		auto it = sampledPoints_.begin();
		for(int i = 0; i < n; i++, it++){
			double sum = 0;
			for(int k = 0; k < d; k++){
				sum += theta_[k] * pow(abs((*it)(k) - (varianceTestPoints_[j])(k)), pVector_[k]);
			}
			rOld(i) = exp(-sum);
		}
		double sum = 0;
		for(int k = 0; k < d; k++){
			sum += theta_[k] * pow(abs(xCopy(k) - (varianceTestPoints_[j])(k)), pVector_[k]);
		}

		double sumPointComp = exp(-sum);
		double result = sigma_ * (rOld.transpose() * weights - sumPointComp) * (rOld.transpose() * weights - sumPointComp);
		// It should not matter a whole lot but if want to unscale the output, need to multiply 
		if(unscaleOutput){
			result = (maxObservation_ - minObservation_) * (maxObservation_ - minObservation_) * result;
		}
		sumNumerator = sumNumerator + result;
	}
	return (sumNumerator / denumerator) / (int)varianceTestPoints_.size();
}


tuple<double, double> Kriging::meanVarianceCalculator(VectorXd &x){
	int n = (int)sampledPoints_.size();
	int d = (int)sampledPoints_[0].size();
	// First want r vector
	VectorXd r(n);
	auto it = sampledPoints_.begin();
	for(int i = 0; i < n; i++, it++){
		double sum = 0;
		for(int k = 0; k < d; k++){
			sum += theta_[k] * pow(abs((*it)(k) - x(k)), pVector_[k]);
		}
		r(i) = exp(-sum);
		
	}

	VectorXd rightHandSide = r.transpose() * modelWeights_;
	double s_x = mu_ + rightHandSide(0);
	if(n == 1){s_x = mu_;}
	MatrixXd var_mid = r.transpose() * rMatrixInverse_ * r;
	double variance = sigma_ * (1 - var_mid(0,0));
	
	return make_tuple(s_x, variance);
}

double Kriging::evaluateAcquisitionFunction(VectorXd x){
	if(chosenAcquisiton_.compare("variance") == 0){
		return uncertainty(x);
	}else if(chosenAcquisiton_.compare("globalVarianceApproximation") == 0){
		return globalVarianceReductionApproximation(x);
	}else if(chosenAcquisiton_.compare("globalVariance") == 0){
		return globalVarianceReduction(x);
	}else if(chosenAcquisiton_.compare("expectedImprovement") == 0){
		return expectedImprovement(x);
	}else{
		return SurrogateModel::evaluateAcquisitionFunction(x);
	}
}

void Kriging::setAcquisitionFunction(string chosenAcquisiton){
	// Adding this bool here in case the same model is assigned different acquisition functions
	// at different stages. Essentially it looks like if sampling is based on global variance,
	// need to restrict the p hyperparameters to be 2.
	if(chosenAcquisiton.compare("variance") == 0){
		chosenAcquisiton_ = chosenAcquisiton;
		acquisitionIsMin_ = false;
		fixedP_ = false;
	}else if(chosenAcquisiton.compare("globalVariance") == 0){
		chosenAcquisiton_ = chosenAcquisiton;
		acquisitionIsMin_ = false;
		fixedP_ = true;
	}else if(chosenAcquisiton.compare("globalVarianceApproximation") == 0){
		chosenAcquisiton_ = chosenAcquisiton;
		acquisitionIsMin_ = false;
		fixedP_ = false;	
	}else if(chosenAcquisiton.compare("expectedImprovement") == 0){
		chosenAcquisiton_ = chosenAcquisiton;
		acquisitionIsMin_ = false;
		fixedP_ = false;
	}else{
		fixedP_ = false;
		SurrogateModel::setAcquisitionFunction(chosenAcquisiton);
	}
}

tuple<VectorXd, double, bool, bool> Kriging::findNextSampleSite(){

	// If dealing with globalVarianceApproximation, choose a random set of points for monte-carlo integration
	if(chosenAcquisiton_.compare("globalVarianceApproximation") == 0){
		varianceTestPoints_ = sampleGenerator_->randomLHS(100);
		scalePoints(varianceTestPoints_);
	}
	// If dealing with globalVariance precalculate W matrix for speed
	if(chosenAcquisiton_.compare("globalVariance") == 0){
		setWmatrix();
	}
	tuple<VectorXd, double, bool, bool> result = SurrogateModel::findNextSampleSite();
	// Return the same, but specify that should only sample the high-fidelity source
	return make_tuple(get<0>(result), get<1>(result), true, false);
}


double Kriging::concentratedLikelihoodFunction(){
	// Get info needed, data contains mu, sigma, R and R^-1 in that order
	auto data = muSigmaCalculator();
	int n = (int)sampledPoints_.size();
	double logLikelihood = -n * logl(get<1>(data))/2 - logl((get<3>(data)).matrixL().determinant() * (get<3>(data)).vectorD().prod())/2;
	// If cannot calculate the second half, say it is a bad answer! Do so to not use hyperparamters which lead to errors in calculation
	if(isinf(logLikelihood) || isnan(logLikelihood)){
		return -DBL_MAX;
	}
	// Save this loglikelihood as the "best" one if it is the best one found
	if(maxLogLikelihood_ < logLikelihood){
		maxLogLikelihood_ = logLikelihood;
	}
	return logLikelihood;
}


Kriging::ConcentratedLikelihoodFunction::ConcentratedLikelihoodFunction(int d, vector<double> &lowerBound, vector<double> &upperBound, Kriging* krigingModel):
	Function(d, lowerBound, upperBound),
	krigingModel_(krigingModel){}

Kriging::ConcentratedLikelihoodFunction::~ConcentratedLikelihoodFunction(){}

double Kriging::ConcentratedLikelihoodFunction::evaluate(VectorXd &point){
	// First half of the entries of point contain the d theta variables, and the second half the p values
	// Should find out dimension and size, then extract theta and p values
	int d = (int)krigingModel_->sampledPoints_[0].size();
	vector<double> theta(d,0.0);
	vector<double> pVector(d,0.0);
	for(int i = 0; i < d; i++){
		theta[i] = pow(10,point(i));
		// theta[i] = point(i);
		if(krigingModel_->fixedP_){
			pVector[i] = 2;
		}else{
			pVector[i] = point(d + i);
		}
	}
	krigingModel_->theta_ = theta;
	krigingModel_->pVector_ = pVector;
	// Get info needed, data contains mu, sigma, R and R^-1 in that order
	return krigingModel_->concentratedLikelihoodFunction();
}


int Kriging::ConcentratedLikelihoodFunction::betterPoint(VectorXd &point1, double val1, VectorXd &point2, double val2, string flag){
	int d = (int)krigingModel_->sampledPoints_[0].size();
	double distance1, distance2;
	if(flag.compare("maxmin") == 0){
		distance1 = DBL_MAX;
		distance2 = DBL_MAX;
		for(int i = 0; i < d; i++){
			double val = pow(10, -point1(i) / point1(d + i));
			if(val < distance1){distance1 = val;}
			val = pow(10, -point2(i) / point2(d + i));
			if(val < distance2){distance2 = val;}
		}
	}else if(flag.compare("maxsum") == 0){
		distance1 = 0;
		distance2 = 0;
		for(int i = 0; i < d; i++){
			double val = pow(10, -point1(i) / point1(d + i));
			distance1 += sqrt(val);
			val = pow(10, -point2(i) / point2(d + i));
			distance2 += sqrt(val);
		}
	}else{
		printf("Problem, asking to compare loglikelihoods without specifying a flag with defined behaviour! Stopping now...\n");
		exit(0);
	}
	
	// If they are both close enough to the best likelihood found so far, compare the distances, prefer the larger one
	if(abs(val1 - krigingModel_->maxLogLikelihood_) < OBJTOL && abs(val2 - krigingModel_->maxLogLikelihood_) < OBJTOL){
		if(abs(distance1 - distance2) < TOL){return 0;}
		else if(distance1 > distance2){return -1;}
		else{return 1;}

	}else if(abs(val1 - krigingModel_->maxLogLikelihood_) < OBJTOL){
		// The first point is close enough to the max likelihood, but the second is not
		return -1;
	}else if(abs(val2 - krigingModel_->maxLogLikelihood_) < OBJTOL){
		// The second point is close enough to the max likelihood, but the first one is not
		return 1;
	}else{
		if(abs(val1 - val2) < TOL){
			return 0;
		}else if(val1 < val2){return 1;}
		else{return -1;}
	}
}















CoKriging::CoKriging(BiFidelityFunction* biFunction, AuxSolver* auxSolver, int randomSeed, bool printInfo, bool functionScaling):
	Kriging(biFunction, auxSolver, randomSeed, printInfo, functionScaling){

	lowFiKriging_ = new Kriging(biFunction_, auxSolver_, randomSeed_, printInfo_, functionScaling);
	
	thetaB_.reserve(biFunction_->d_);
	pBVector_.reserve(biFunction_->d_);

}

CoKriging::~CoKriging(){}


void CoKriging::saveSample(vector<VectorXd> points, vector<VectorXd> pointsLow, vector<double> observations, vector<double> observationsLow){
	// First to the usual saving of the data
	SurrogateModel::saveSample(points, pointsLow, observations, observationsLow);
	// And now save the relevant things to lowFiKriging
	lowFiKriging_->sampledPoints_ = sampledPointsLow_;
	lowFiKriging_->sampledPointsValues_ = sampledPointsValuesLow_;
	lowFiKriging_->maxObservation_ = maxObservation_;
	lowFiKriging_->minObservation_ = minObservation_;
	
}

void CoKriging::trainModel(bool forcedOptimiseHyperparameters){
	int n = (int)sampledPoints_.size();
	if(printInfo_){printf("Training cokriging model with %d high- and %d low-fidelity samples, start with low fi kriging.\n", n, (int)sampledPointsLow_.size());}
	lowFiKriging_->trainModel(forcedOptimiseHyperparameters);
	if(printInfo_){printf("Training cokriging model, train intermediate model - ");}

	if(forcedOptimiseHyperparameters){
		if(printInfo_){printf("optimising hyperparameters as user requested to regardless of other factors.\n");}
		trainHyperparameters();
		sampleSizeAtLastTrain_ = n;

	}else if(startIntervalTraining_ >= n){
		if(printInfo_){printf("optimising hyperparameters as working with manageable sample size.\n");}
		trainHyperparameters();
		sampleSizeAtLastTrain_ = n;

	}else if((startLargeIntervalTraining_ > n && (n - sampleSizeAtLastTrain_) >= intervalForTraining_) ||
				(startVeryLargeIntervalTraining_ > n && (n - sampleSizeAtLastTrain_) >= intervalForLargeTraining_) ||
				(n - sampleSizeAtLastTrain_) >= intervalForVeryLargeTraining_){
		if(printInfo_){printf("optimising hyperparameters despite having large sample as it's been %d samples since last optimisation.\n", n - sampleSizeAtLastTrain_);}
		trainHyperparameters();
		sampleSizeAtLastTrain_ = n;

	}else{
		if(printInfo_){printf("skipping  hyperparameter optimisation.\n");}
	}

	saveMuSigma();	

	// It is possible that the predictions are now undefined,
	// if the hyperparameters where not optimised.
	// This can happen if there are lots of points.
	// Check that a random prediction can be given; if not call the function again
	// and force the training of hyperparameters
	// Get a random point
	VectorXd point = (sampleGenerator_->randomSample(1))[0];
	double val = surfaceValue(point);
	if(isnan(val)){
		if(printInfo_){printf("Skipped hyperparameter optimisation but now predictions give NaN. Forcing retuning of hyperparameters.\n");}
		trainModel(true);
	}
	return;
}

void CoKriging::trainHyperparameters(){
	int d = biFunction_->d_;
	// First define the bounds for the function
	vector<double> lowerBound(2*d + 1, 0.0);
	vector<double> upperBound(2*d + 1, 0.0);
	// Bounds should be dependent on whether data has been scaled; if it has the bounds can be more restrictive
	for(int i = 0; i < d; i++){
		// Bounds for theta_b
		if(functionScaling_){
			lowerBound[i] = -3;
			upperBound[i] = 3;
		}else{
			lowerBound[i] = -10;
			upperBound[i] = 3;
		}
		// Bounds for p_b
		lowerBound[d+i] = 0.1;
		upperBound[d+i] = 2.0;
	}
	// Bounds for rho
	lowerBound[2*d] = -1000;
	upperBound[2*d] = 1000;

	// Store low fi surface values at high fi locations, clear vector if it was populated
	lowFiSurfaceValues_.clear();
	lowFiSurfaceValues_.reserve((int)sampledPoints_.size());
	for(int i = 0; i < (int)sampledPoints_.size(); i++){
		// Don't want to return scaled predictions from the low fi model, as internally everything is scaled (if it was requested by the user)
		lowFiSurfaceValues_.push_back(lowFiKriging_->surfaceValue(sampledPoints_[i], true, false));
	}
	maxLogLikelihood_ = -DBL_MAX;
	IntermediateConcentratedLikelihoodFunction* function;
	if(fixedP_){
		lowerBound[d] = -1000;
		upperBound[d] = 1000;
		function = new IntermediateConcentratedLikelihoodFunction(d + 1, lowerBound, upperBound, this);
	}else{
		function = new IntermediateConcentratedLikelihoodFunction(2*d + 1, lowerBound, upperBound, this);
	}
	auxSolver_->updateProblem(function, false);
	VectorXd hyperparameters = auxSolver_->optimise();
	delete function;
	// Extract two vectors, and store!
	for(int i = 0; i < d; i++){
		thetaB_[i] = pow(10,hyperparameters(i));

		if(fixedP_){
			pBVector_[i] = 2;
		}else{
			pBVector_[i] = hyperparameters(d + i);
		}
	}
	// Last value is rho
	if(fixedP_){
		rho_ = hyperparameters(d);
	}else{
		rho_ = hyperparameters(2*d);
	}
}

void CoKriging::saveMuSigma(){
	auto data = intermediateMuSigmaCalculator();
	sigmaL_ = lowFiKriging_->sigma_;
	sigmaB_ = get<1>(data);
	auto data1 = combinedMuCmatrixCalculator();
	muCombined_ = get<0>(data1);
	cMatrix_ = get<1>(data1);
	cMatrixDecomposition_ = get<2>(data1);
	
	// I think that since I only find the inverse once here, I want to 
	// do a decomposition once again. Using LU decomposition from Eigen,
	// as it appears it is the most stable and accurate when finding an inverse.
	// Note this should work even when the matrix is non invertible, which is a huge plus,
	// even though this should not happen.
	int n = cMatrix_.cols();
	cMatrixInverse_ = cMatrix_.fullPivLu().solve(MatrixXd::Identity(n,n));
	
	return;
}

tuple<double, double, MatrixXd, LDLT<MatrixXd>> CoKriging::intermediateMuSigmaCalculator(){	
	int n = (int)sampledPoints_.size();
	int d = (int)sampledPoints_[0].size();
	if(n < 1){
		printf("Trying to train CoKriging model with less points than needed! Have %d but need %d. Exiting now...\n", n, 1);
		exit(0);
	}

	MatrixXd bExpensiveMatrix(n,n);
	// First of all need to define matrix
	auto it1 = sampledPoints_.begin();
	for(int i = 0; i < n; i++, it1++){
		auto it2 = it1;
		for(int j = i; j < n; j++, it2++){
			double sum = 0;
			for(int k = 0; k < d; k++){
				sum += thetaB_[k] * pow(abs((*it1)(k) - (*it2)(k)), pBVector_[k]);
			}
			bExpensiveMatrix(i,j) = exp(-sum);
			bExpensiveMatrix(j,i) = bExpensiveMatrix(i,j);	
		}
	}
	VectorXd one(n);
	VectorXd b(n);
	for(int i = 0; i < n; i++){
		one(i) = 1;
		b(i) = sampledPointsValues_[i] - rho_ * lowFiSurfaceValues_[i];
	}

	LDLT<MatrixXd> bExpensiveMatrixDecomposition = bExpensiveMatrix.ldlt();
	MatrixXd mu_top = one.transpose() * bExpensiveMatrixDecomposition.solve(b);
	MatrixXd mu_bottom = one.transpose() * bExpensiveMatrixDecomposition.solve(one);
	double mu = mu_top(0,0)/ mu_bottom(0,0);

	MatrixXd sigma_m = (b - one*mu).transpose() * bExpensiveMatrixDecomposition.solve(b - one*mu);
	double sigma = sigma_m(0,0)/n;

	return make_tuple(mu, sigma, bExpensiveMatrix, bExpensiveMatrixDecomposition);
}

tuple<double, MatrixXd, LDLT<MatrixXd>> CoKriging::combinedMuCmatrixCalculator(){

	int nL = (int)lowFiKriging_->sampledPoints_.size();
	int nH = (int)sampledPoints_.size();
	int d = (int)sampledPoints_[0].size();

	MatrixXd cMatrix(nH + nL, nH + nL);
	// Start with Rl(Xl, Xl)
	auto it1 = lowFiKriging_->sampledPoints_.begin();
	for(int i = 0; i < nL; i++, it1++){
		auto it2 = it1;
		for(int j = i; j < nL; j++, it2++){
			double sum = 0;
			for(int k = 0; k < d; k++){
				sum += lowFiKriging_->theta_[k] * pow(abs((*it1)(k) - (*it2)(k)), lowFiKriging_->pVector_[k]);
			}
			sum = exp(-sum);
			sum = sum * sigmaL_;
			cMatrix(i,j) = sum;
			cMatrix(j,i) = cMatrix(i,j);	
		}
	}
	// Now go with Rl(Xl, Xh) and Rl(Xh, Xl)
	it1 = lowFiKriging_->sampledPoints_.begin();
	for(int i = 0; i < nL; i++, it1++){
		auto it2 = sampledPoints_.begin();
		for(int j = 0; j < nH; j++, it2++){
			double sum = 0;
			for(int k = 0; k < d; k++){
				sum += lowFiKriging_->theta_[k] * pow(abs((*it1)(k) - (*it2)(k)), lowFiKriging_->pVector_[k]);
			}
			sum = exp(-sum);
			sum = sum * rho_ * sigmaL_;
			cMatrix(i, nL + j) = sum;
			cMatrix(nL + j, i) = sum;
		}
	}

	// Now go with Rl(Xh, Xh) + Rb(Xh, Xh)
	it1 = sampledPoints_.begin();
	for(int i = 0; i < nH; i++, it1++){
		auto it2 = it1;
		for(int j = i; j < nH; j++, it2++){
			double sumL = 0;
			double sumB = 0;
			for(int k = 0; k < d; k++){
				sumL += lowFiKriging_->theta_[k] * pow(abs((*it1)(k) - (*it2)(k)), lowFiKriging_->pVector_[k]);
				sumB += thetaB_[k] * pow(abs((*it1)(k) - (*it2)(k)), pBVector_[k]);
			}
			sumL = exp(-sumL);
			sumB = exp(-sumB);
			double total = rho_ * rho_ * sigmaL_ * sumL + sigmaB_ * sumB;
			cMatrix(nL + i, nL + j) = total;
			cMatrix(nL + j, nL + i) = total;
		}
	}

	VectorXd one(nL + nH);
	VectorXd y(nL + nH);
	for(int i = 0; i < nL; i++){
		one(i) = 1;
		y(i) = lowFiKriging_->sampledPointsValues_[i];
	}
	for(int i = 0; i < nH; i++){
		one(nL + i) = 1;
		y(nL + i) = sampledPointsValues_[i];
	}

	LDLT<MatrixXd> cMatrixDecomposition = cMatrix.ldlt();
	// MatrixXd mu_top = one.transpose() * cMatrixDecomposition.solve(y);
	// MatrixXd mu_bottom = one.transpose() * cMatrixDecomposition.solve(one);

	MatrixXd mu_top = one.transpose() * cMatrix.fullPivLu().solve(y);
	MatrixXd mu_bottom = one.transpose() * cMatrix.fullPivLu().solve(one);

	double mu = mu_top(0,0)/ mu_bottom(0,0);
	return make_tuple(mu, cMatrix, cMatrixDecomposition);
}

double CoKriging::intermediateConcentratedLikelihoodFunction(){
	auto data = intermediateMuSigmaCalculator();
	int n = (int)sampledPoints_.size();
	double logLikelihood = -n * logl(get<1>(data))/2 - logl((get<3>(data)).matrixL().determinant()) - logl((get<3>(data)).vectorD().prod())/2;

	if(isinf(logLikelihood) || isnan(logLikelihood)){return -DBL_MAX;}
	if(maxLogLikelihood_ < logLikelihood){
		maxLogLikelihood_ = logLikelihood;
	}

	return logLikelihood;
}

tuple<double, double> CoKriging::meanVarianceCalculator(VectorXd &x){
	int nL = (int)lowFiKriging_->sampledPoints_.size();
	int nH = (int)sampledPoints_.size();
	// First want c vector
	VectorXd c = cVector(x);
	VectorXd one(nL + nH);
	VectorXd y(nL + nH);
	for(int i = 0; i < nL; i++){
		one(i) = 1;
		y(i) = lowFiKriging_->sampledPointsValues_[i];
	}
	for(int i = 0; i < nH; i++){
		one(nL + i) = 1;
		y(nL + i) = sampledPointsValues_[i];
	}

	// New way
	VectorXd rightHandSide = c.transpose() * cMatrixInverse_ * (y - one * muCombined_);
	
	double s_x = muCombined_ + rightHandSide(0);
	rightHandSide = c.transpose() * cMatrixInverse_ * c;
	double variance = rho_ * rho_ * sigmaL_ + sigmaB_ - rightHandSide(0);

	return make_tuple(s_x, variance);
}


VectorXd CoKriging::cVector(VectorXd &x){
	int nL = (int)lowFiKriging_->sampledPoints_.size();
	int nH = (int)sampledPoints_.size();
	int d = (int)sampledPoints_[0].size();
	VectorXd c(nL + nH);
	// First deal with Rl(Xl, x)
	auto it = lowFiKriging_->sampledPoints_.begin();
	for(int i = 0; i < nL; i++, it++){
		double sum = 0;
		for(int k = 0; k < d; k++){
			sum += lowFiKriging_->theta_[k] * pow(abs((*it)(k) - x(k)), lowFiKriging_->pVector_[k]);
		}
		c(i) = rho_ * sigmaL_ * exp(-sum);
	}
	// Deal with Rl(Xh, x) + Rb(Xh, x)
	it = sampledPoints_.begin();
	for(int i = 0; i < nH; i++, it++){
		double sumC = 0;
		double sumB = 0;
		for(int k = 0; k < d; k++){
			sumC += lowFiKriging_->theta_[k] * pow(abs((*it)(k) - x(k)), lowFiKriging_->pVector_[k]);
			sumB += thetaB_[k] * pow(abs((*it)(k) - x(k)), pBVector_[k]);
		}
		c(nL + i) = rho_ * rho_ * sigmaL_ * exp(-sumC) + sigmaB_ * exp(-sumB);
	}
	return c;
}

VectorXd CoKriging::cLowVector(VectorXd &x){
	int nL = (int)lowFiKriging_->sampledPoints_.size();
	int nH = (int)sampledPoints_.size();
	int d = (int)sampledPoints_[0].size();
	VectorXd c(nL + nH);
	// First deal with X_l
	auto it = lowFiKriging_->sampledPoints_.begin();
	for(int i = 0; i < nL; i++, it++){
		double sum = 0;
		for(int k = 0; k < d; k++){
			sum += lowFiKriging_->theta_[k] * pow(abs((*it)(k) - x(k)), lowFiKriging_->pVector_[k]);
		}
		c(i) = sigmaL_ * exp(-sum);
	}
	// Deal with X_h
	it = sampledPoints_.begin();
	for(int i = 0; i < nH; i++, it++){
		double sum = 0;
		for(int k = 0; k < d; k++){
			sum += lowFiKriging_->theta_[k] * pow(abs((*it)(k) - x(k)), lowFiKriging_->pVector_[k]);
		}
		c(nL + i) = rho_ * sigmaL_ * exp(-sum);
	}
	return c;
}

MatrixXd CoKriging::cBothMatrix(VectorXd &x){
	int nL = (int)lowFiKriging_->sampledPoints_.size();
	int nH = (int)sampledPoints_.size();
	VectorXd cLow = cLowVector(x);
	VectorXd c = cVector(x);
	MatrixXd cBoth(nL + nH, 2);
	for(int i = 0; i < nL + nH; i++){
		cBoth(i, 0) = cLow(i);
		cBoth(i, 1) = c(i);
	}
	return cBoth;

}

double CoKriging::w1product(VectorXd &point1, VectorXd &point2){
	double product = (rho_ * sigmaL_)*(rho_ * sigmaL_);
	for(int i = 0; i < point1.size(); i++){
		product = product * kernelProductIntegral(point1(i), point2(i), lowFiKriging_->theta_[i], lowFiKriging_->theta_[i]);
	}
	return product;
}

double CoKriging::w2product(VectorXd &point1, VectorXd &point2){	
	double productLow = rho_ * rho_ * rho_ * sigmaL_ * sigmaL_;
	double productDiff = rho_ * sigmaL_ * sigmaB_;
	for(int i = 0; i < point1.size(); i++){
		productLow = productLow * kernelProductIntegral(point1(i), point2(i), lowFiKriging_->theta_[i], lowFiKriging_->theta_[i]);
		productDiff = productDiff * kernelProductIntegral(point1(i), point2(i), lowFiKriging_->theta_[i], thetaB_[i]);	
	}
	return productLow + productDiff;
}

double CoKriging::w3product(VectorXd &point1, VectorXd &point2){
	double productLow = (rho_ * rho_ * sigmaL_) * (rho_ * rho_ * sigmaL_);
	double productDiffOne = rho_ * rho_ * sigmaL_ * sigmaB_;
	double productDiffTwo = rho_ * rho_ * sigmaL_ * sigmaB_;
	double productDiff = sigmaB_ * sigmaB_;
	for(int i = 0; i < point1.size(); i++){
		productLow = productLow * kernelProductIntegral(point1(i), point2(i), lowFiKriging_->theta_[i], lowFiKriging_->theta_[i]);
		productDiffOne = productDiffOne * kernelProductIntegral(point1(i), point2(i), lowFiKriging_->theta_[i], thetaB_[i]);
		productDiffTwo = productDiffTwo * kernelProductIntegral(point1(i), point2(i), thetaB_[i], lowFiKriging_->theta_[i]);
		productDiff = productDiff * kernelProductIntegral(point1(i), point2(i), thetaB_[i], thetaB_[i]);
	}
	return productLow + productDiffOne + productDiffTwo + productDiff;
}


void CoKriging::setWmatrix(){
	int nL = (int)lowFiKriging_->sampledPoints_.size();
	int nH = (int)sampledPoints_.size();

	MatrixXd wMatrix(nH + nL, nH + nL);
	// Start with W1
	auto it1 = lowFiKriging_->sampledPoints_.begin();
	for(int i = 0; i < nL; i++, it1++){
		auto it2 = it1;
		for(int j = i; j < nL; j++, it2++){
			wMatrix(i,j) = w1product((*it1), (*it2));
			wMatrix(j,i) = wMatrix(i,j);
		}
	}
	// Now go with W2
	it1 = lowFiKriging_->sampledPoints_.begin();
	for(int i = 0; i < nL; i++, it1++){
		auto it2 = sampledPoints_.begin();
		for(int j = 0; j < nH; j++, it2++){
			wMatrix(i,nL + j) = w2product((*it1), (*it2));
			wMatrix(nL + j,i) = w2product((*it1), (*it2));
		}
	}
	// Now go with W3
	it1 = sampledPoints_.begin();
	for(int i = 0; i < nH; i++, it1++){
		auto it2 = it1;
		for(int j = i; j < nH; j++, it2++){
			wMatrix(nL + i, nL + j) = w3product((*it1), (*it2));
			wMatrix(nL + j, nL + i) = wMatrix(nL + i,nL + j);
		}
	}
	wMatrix_ = wMatrix;
}

double CoKriging::globalVarianceReductionLowFiSample(VectorXd &x, bool pointIsScaled, bool unscaleOutput){
	// First calculate the denumerator
	int nL = (int)lowFiKriging_->sampledPoints_.size();
	int nH = (int)sampledPoints_.size();
	VectorXd xCopy = x;
	if(!pointIsScaled){scalePoint(xCopy);}
	double denumerator;
	// Will want cL vector at this point
	VectorXd cLow = cLowVector(xCopy);
	// Will need this vector times inverse of C matrix in a couple of places, work it out
	auto weights = cMatrixInverse_ * cLow;
	denumerator = sigmaL_ - cLow.transpose() * weights;
	// Populate wLow vector
	VectorXd wLow(nL + nH);
	auto it = lowFiKriging_->sampledPoints_.begin();
	for(int i = 0; i < nL; i++, it++){
		wLow(i) = w1product(xCopy, (*it));
	}
	it = sampledPoints_.begin();
	for(int i = 0; i < nH; i++, it++){
		wLow(nL + i) = w2product(xCopy, (*it));
	}

	double numerator = (weights.transpose() * wMatrix_ * weights)(0) - 2 * (wLow.transpose() * weights)(0) + w1product(xCopy, xCopy);
	double result =  numerator / denumerator;
	if(unscaleOutput){
		result = (maxObservation_ - minObservation_) * (maxObservation_ - minObservation_) * result;
	}

	return result;

}

double CoKriging::globalVarianceReductionHighFiSample(VectorXd &x, bool pointIsScaled, bool unscaleOutput){
	int nL = (int)lowFiKriging_->sampledPoints_.size();
	int nH = (int)sampledPoints_.size();
	VectorXd xCopy = x;
	if(!pointIsScaled){scalePoint(xCopy);}
	double denumerator;
	// Will want cL vector at this point
	VectorXd c = cVector(xCopy);
	// Will need this vector times inverse of C matrix in a couple of places, work it out
	auto weights = cMatrixInverse_ * c;
	denumerator = rho_ * rho_ * sigmaL_ + sigmaB_ - c.transpose() * weights;
	// Populate wLow vector
	VectorXd wHigh(nL + nH);
	auto it = lowFiKriging_->sampledPoints_.begin();
	for(int i = 0; i < nL; i++, it++){
		wHigh(i) = w2product((*it), xCopy);
	}
	it = sampledPoints_.begin();
	for(int i = 0; i < nH; i++, it++){
		wHigh(nL + i) = w3product((*it), xCopy);
	}
	double numerator = (weights.transpose() * wMatrix_ * weights)(0) - 2 * (wHigh.transpose() * weights)(0) + w3product(xCopy, xCopy);
	double result =  numerator / denumerator;
	if(unscaleOutput){
		result = (maxObservation_ - minObservation_) * (maxObservation_ - minObservation_) * result;
	}
	return result;
}

double CoKriging::globalVarianceReductionBothFiSample(VectorXd &x, bool pointIsScaled, bool unscaleOutput){

	// First check that there is some variance coming from both sources. If not from low-fi source, only 
	// reduction comes from sampling the high fi source. If not from high-fi source, only reduction comes from
	// sampling the low fi source. If both sources have no variance then adding a point will not change anything.
	if(sigmaL_ < TOL && sigmaB_ < TOL){return 0.0;}
	if(sigmaL_ < TOL){return globalVarianceReductionHighFiSample(x, pointIsScaled, unscaleOutput);}
	if(sigmaB_ < TOL){return globalVarianceReductionLowFiSample(x, pointIsScaled, unscaleOutput);}
	

	int nL = (int)lowFiKriging_->sampledPoints_.size();
	int nH = (int)sampledPoints_.size();
	VectorXd xCopy = x;
	if(!pointIsScaled){scalePoint(xCopy);}
	// Will want cL vector at this point
	MatrixXd cBothNext = cBothMatrix(xCopy);
	MatrixXd sigmaMatrix(2,2);
	sigmaMatrix(0,0) = sigmaL_;
	sigmaMatrix(1,0) = rho_ * sigmaL_;
	sigmaMatrix(0,1) = rho_ * sigmaL_;
	sigmaMatrix(1,1) = rho_ * rho_ * sigmaL_ + sigmaB_;
	
	MatrixXd wBmatrix(2,2);
	wBmatrix(0,0) = w1product(xCopy, xCopy);
	wBmatrix(1,0) = w2product(xCopy, xCopy);
	wBmatrix(0,1) = w2product(xCopy, xCopy);
	wBmatrix(1,1) = w3product(xCopy, xCopy);

	MatrixXd wBoth(nL + nH, 2);
	auto it = lowFiKriging_->sampledPoints_.begin();
	for(int i = 0; i < nL; i++, it++){
		wBoth(i, 0) = w1product(xCopy, (*it));
		wBoth(i, 1) = w2product((*it), xCopy);
	}
	it = sampledPoints_.begin();
	for(int i = 0; i < nH; i++, it++){
		wBoth(nL + i, 0) = w2product(xCopy, (*it));
		wBoth(nL + i, 1) = w3product((*it), xCopy);
	}

	
	auto weights = cMatrixInverse_ * cBothNext;
	MatrixXd alphaMatrix = (sigmaMatrix - cBothNext.transpose() * weights).inverse();
	MatrixXd topLeft = weights * alphaMatrix * weights.transpose();
	MatrixXd topRight = - weights * alphaMatrix;
	MatrixXd bottomLeft = - alphaMatrix * weights.transpose();
	
	double result = 0;
	for(int i = 0; i < topLeft.rows(); i++){
		for(int j = 0; j < topLeft.cols(); j++){
			result += topLeft(i,j) * wMatrix_(i,j);
		}
	}
	for(int i = 0; i < topRight.rows(); i++){
		for(int j = 0; j < topRight.cols(); j++){
			result += topRight(i,j) * wBoth(i,j);
		}
	}
	for(int i = 0; i < bottomLeft.rows(); i++){
		for(int j = 0; j < bottomLeft.cols(); j++){
			result += bottomLeft(i,j) * wBoth(j,i);
		}
	}
	for(int i = 0; i < alphaMatrix.rows(); i++){
		for(int j = 0; j < alphaMatrix.cols(); j++){
			result += alphaMatrix(i,j) * wBmatrix(i,j);
		}
	}
	
	if(unscaleOutput){
		result = (maxObservation_ - minObservation_) * (maxObservation_ - minObservation_) * result;
	}
	return result;


}


// double CoKriging::globalVarianceReductionApproximationLowFiSample(VectorXd &x, bool pointIsScaled, bool unscaleOutput){
// 	// Trying the approach given in "Gaussian Process Regression: Active Data Selection and Test Point Rejection"
// 	// First calculate the denumerator
// 	VectorXd xCopy = x;
// 	if(!pointIsScaled){scalePoint(xCopy);}
// 	double denumerator;
// 	int d = (int)sampledPoints_[0].size();
// 	// Will want cL vector at this point
// 	VectorXd cLow = cLowVector(xCopy);
// 	// Will need this vector times inverse of C matrix in a couple of places, work it out
// 	auto weights = cMatrixInverse_ * cLow;
// 	// auto weights = cMatrixDecomposition_.solve(cLow);
// 	denumerator = sigmaL_ - cLow.transpose() * weights;
// 	// if(sigmaL_ < TOL || (denumerator/sigmaL_) < TOL){return 0;}
// 	// Now for each trial point calculate the nominator
// 	double sumNumerator = 0;
// 	for(int j = 0;  j < (int)varianceTestPoints_.size(); j++){
// 		// Want c vector at this point
// 		VectorXd c = cVector(varianceTestPoints_[j]);
// 		// Here want Rl between the test point and the candidate point
// 		double sum = 0;
// 		for(int k = 0; k < d; k++){
// 			sum += lowFiKriging_->theta_[k] * pow(abs(xCopy(k) - (varianceTestPoints_[j](k))), lowFiKriging_->pVector_[k]);
// 		}
// 		double sumPointComp = exp(-sum);
// 		double result = (c.transpose() * weights - rho_ * sigmaL_ * sumPointComp) * (c.transpose() * weights - rho_ * sigmaL_ * sumPointComp);
// 		// It should not matter a whole lot but if want to unscale the output, need to multiply 
// 		if(unscaleOutput){
// 			result = (maxObservation_ - minObservation_) * (maxObservation_ - minObservation_) * result;
// 		}
// 		sumNumerator = sumNumerator + result;
// 	}
// 	// Note the actual change is the negative of this value, but will keep the absolute and maximise (so as to maximise the reduction)
// 	return (sumNumerator / denumerator) / (int)varianceTestPoints_.size();
// }

// double CoKriging::globalVarianceReductionApproximationHighFiSample(VectorXd &x, bool pointIsScaled, bool unscaleOutput){
// 	// Trying the approach given in "Gaussian Process Regression: Active Data Selection and Test Point Rejection"
// 	// First calculate the denumerator
// 	VectorXd xCopy = x;
// 	if(!pointIsScaled){scalePoint(xCopy);}
// 	double denumerator;
// 	int d = (int)sampledPoints_[0].size();
// 	// Will want cL vector at this point
// 	VectorXd cXnext = cVector(xCopy);
// 	// Will need this vector times inverse of C matrix in a couple of places, work it out
// 	auto weights = cMatrixInverse_ * cXnext;
// 	// auto weights = cMatrixDecomposition_.solve(cXnext);
// 	denumerator = rho_ * rho_ * sigmaL_ + sigmaB_ - cXnext.transpose() * weights;
// 	// if((rho_ * rho_ * sigmaL_ + sigmaB_) < TOL || (denumerator/(rho_ * rho_ * sigmaL_ + sigmaB_)) < TOL){return 0;}
// 	// Now for each trial point calculate the nominator
// 	double sumNumerator = 0;
// 	for(int j = 0;  j < (int)varianceTestPoints_.size(); j++){
// 		// Want c vector at this point
// 		VectorXd c = cVector(varianceTestPoints_[j]);
// 		// Here want Rl between the test point and the candidate point
// 		double sumRl = 0;
// 		double sumRdiff = 0;
// 		for(int k = 0; k < d; k++){
// 			sumRl += lowFiKriging_->theta_[k] * pow(abs(xCopy(k) - (varianceTestPoints_[j](k))), lowFiKriging_->pVector_[k]);
// 			sumRdiff += thetaB_[k] * pow(abs(xCopy(k) - (varianceTestPoints_[j](k))), pBVector_[k]);
			
// 		}
// 		double sumPointComp = rho_ * rho_ * sigmaL_ * exp(-sumRl) + sigmaB_ * exp(-sumRdiff);
// 		double result = ((c.transpose() * weights)(0) - sumPointComp) * ((c.transpose() * weights)(0) - sumPointComp);
// 		// It should not matter a whole lot but if want to unscale the output, need to multiply 
// 		if(unscaleOutput){
// 			result = (maxObservation_ - minObservation_) * (maxObservation_ - minObservation_) * result;
// 		}
// 		sumNumerator = sumNumerator + result;
		
// 	}
// 	// Note the actual change is the negative of this value, but will keep the absolute and maximise (so as to maximise the reduction)
// 	return (sumNumerator / denumerator) / (int)varianceTestPoints_.size();
// }

// double CoKriging::globalVarianceReductionApproximationBothFiSample(VectorXd &x, bool pointIsScaled, bool unscaleOutput){
// 	// First check that there is some variance coming from both sources. If not from low-fi source, only 
// 	// reduction comes from sampling the high fi source. If not from high-fi source, only reduction comes from
// 	// sampling the low fi source. If both sources have no variance then adding a point will not change anything.
// 	if(sigmaL_ < TOL && sigmaB_ < TOL){return 0.0;}
// 	if(sigmaL_ < TOL){return globalVarianceReductionApproximationHighFiSample(x, pointIsScaled, unscaleOutput);}
// 	if(sigmaB_ < TOL){return globalVarianceReductionApproximationLowFiSample(x, pointIsScaled, unscaleOutput);}

// 	// Trying the approach given in "Gaussian Process Regression: Active Data Selection and Test Point Rejection"
// 	// First calculate the denumerator
// 	VectorXd xCopy = x;
// 	if(!pointIsScaled){scalePoint(xCopy);}
// 	int d = (int)sampledPoints_[0].size();
// 	// Will want cL vector at this point
// 	MatrixXd cBothNext = cBothMatrix(xCopy);
// 	MatrixXd sigmaMatrix(2,2);
// 	sigmaMatrix(0,0) = sigmaL_;
// 	sigmaMatrix(1,0) = rho_ * sigmaL_;
// 	sigmaMatrix(0,1) = rho_ * sigmaL_;
// 	sigmaMatrix(1,1) = rho_ * rho_ * sigmaL_ + sigmaB_;
	
// 	auto weights = cMatrixInverse_ * cBothNext;

// 	MatrixXd alphaMatrix = (sigmaMatrix - cBothNext.transpose() * weights).inverse();
// 	// auto alphaMatrixDecomposition = (sigmaMatrix - cBothNext.transpose() * weights).ldlt();

// 	// Now for each trial point calculate the numerator
// 	double sumNumerator = 0;
// 	for(int j = 0;  j < (int)varianceTestPoints_.size(); j++){
// 		// Want c vector at this point
// 		VectorXd c = cVector(varianceTestPoints_[j]);
// 		// Here want Rl between the test point and the candidate point
// 		double sumRl = 0;
// 		double sumRdiff = 0;
// 		for(int k = 0; k < d; k++){
// 			sumRl += lowFiKriging_->theta_[k] * pow(abs(xCopy(k) - (varianceTestPoints_[j](k))), lowFiKriging_->pVector_[k]);
// 			sumRdiff += thetaB_[k] * pow(abs(xCopy(k) - (varianceTestPoints_[j](k))), pBVector_[k]);
			
// 		}
// 		VectorXd kernelVector(2);
// 		kernelVector(0) = rho_ * sigmaL_ * exp(-sumRl);
// 		kernelVector(1) = rho_ * rho_ * sigmaL_ * exp(-sumRl) + sigmaB_ * exp(-sumRdiff);

// 		auto matrixProduct = c.transpose() * weights;

// 		double result = (matrixProduct * alphaMatrix * matrixProduct.transpose())(0) - (kernelVector.transpose() * alphaMatrix * matrixProduct.transpose())(0) - (matrixProduct * alphaMatrix * kernelVector)(0) + (kernelVector.transpose() * alphaMatrix * kernelVector)(0);
		
// 		// It should not matter a whole lot but if want to unscale the output, need to multiply 
// 		if(unscaleOutput){
// 			result = (maxObservation_ - minObservation_) * (maxObservation_ - minObservation_) * result;
// 		}
// 		sumNumerator = sumNumerator + result;	
// 	}
// 	// Note the actual change is the negative of this value, but will keep the absolute and maximise (so as to maximise the reduction)
// 	return sumNumerator / (int)varianceTestPoints_.size();
// }

double CoKriging::evaluateAcquisitionFunction(VectorXd x){
	if(chosenAcquisiton_.compare("globalVarianceLow") == 0){
		return  globalVarianceReductionLowFiSample(x);
	}else if(chosenAcquisiton_.compare("globalVarianceHigh") == 0){
		return globalVarianceReductionHighFiSample(x);
	// }else if(chosenAcquisiton_.compare("globalVarianceApproximationLow") == 0){
	// 	return globalVarianceReductionApproximationLowFiSample(x);
	// }else if(chosenAcquisiton_.compare("globalVarianceApproximationHigh") == 0){
		// return globalVarianceReductionApproximationHighFiSample(x);
	}else if(chosenAcquisiton_.compare("globalVariance") == 0){
		return globalVarianceReductionBothFiSample(x);
	// }else if(chosenAcquisiton_.compare("globalVarianceApproximation") == 0){
	// 	return globalVarianceReductionApproximationBothFiSample(x);
	}else{
		return Kriging::evaluateAcquisitionFunction(x);
	}
}

void CoKriging::setAcquisitionFunction(string chosenAcquisiton){
	// Adding this bool here in case the same model is assigned different acquisition functions
	// at different stages. Essentially it looks like if sampling is based on global variance,
	// need to restrict the p hyperparameters to be 2.
	if(chosenAcquisiton.compare("globalVarianceWithChoice") == 0){
		chosenAcquisiton_ = chosenAcquisiton;
		acquisitionIsMin_ = false;
		fixedP_ = true;
	}else if(chosenAcquisiton.compare("globalVarianceApproximationWithChoice") == 0){
		chosenAcquisiton_ = chosenAcquisiton;
		acquisitionIsMin_ = false;
		fixedP_ = false;	
	}else{
		Kriging::setAcquisitionFunction(chosenAcquisiton);
	}

	// Need to apply the fixedP_ restriction to the low fi Kriging model as well!
	lowFiKriging_->fixedP_ = fixedP_;
}


tuple<VectorXd, double, bool, bool> CoKriging::findNextSampleSite(){
	// Kind of need to repeat the logic of Kriging and then skip the Kriging call, things run slightly different here
	if(chosenAcquisiton_.compare("globalVarianceApproximationWithChoice") == 0){
		varianceTestPoints_ = sampleGenerator_->randomLHS(100);
		scalePoints(varianceTestPoints_);
		// Change name to low and high, and find the best point for each
		chosenAcquisiton_ = "globalVarianceApproximationHigh";
		tuple<VectorXd, double, bool, bool> resultHigh = SurrogateModel::findNextSampleSite();
		VectorXd pointHigh = get<0>(resultHigh);
		double valueHigh = get<1>(resultHigh);

		// It is possible that we already have way too many low fi sample points.
		// If this is the case, do not bother looking for promising sample points.
		if((int)sampledPointsLow_.size() >= 50 * biFunction_->d_){
			if(printInfo_){
				printf("Low fi point is skipped as already have %d low-fi points ", (int)sampledPointsLow_.size());
				printf("\nHigh fi point has value %.4f at point ", valueHigh);
				printPoint(pointHigh);
				printf("\n");
			}
			chosenAcquisiton_ = "globalVarianceApproximationWithChoice";
			return make_tuple(pointHigh, valueHigh, true, false);
		}

		chosenAcquisiton_ = "globalVarianceApproximationLow";
		tuple<VectorXd, double, bool, bool> resultLow = SurrogateModel::findNextSampleSite();
		VectorXd pointLow = get<0>(resultLow);
		double valueLow = get<1>(resultLow);

		if(printInfo_){
			printf("Low fi point has value %.4f at point ", valueLow);
			printPoint(pointLow);
			printf("\nHigh fi point has value %.4f at point ", valueHigh);
			printPoint(pointHigh);
			printf("\n");
		}
		chosenAcquisiton_ = "globalVarianceApproximationWithChoice";

		// Here want to choose which source to sample, base it on cost ratio.
		// First check that this value has been set.
		// Default should be -1 (which is incorrect), so just check this
		if(costRatio_ < -0.5){
			printf("Cost ratio needs to be used but it has not yet been set! Stopping now...\n");
			exit(0);
		}
		// Now just divide low fi value by cost ratio for comparison. If the cost ratio is too close to 0, use 0.001 as a stand in
		double relativeValueLow;
		if(abs(costRatio_) < TOL){relativeValueLow = valueLow / 0.001;}
		else{relativeValueLow = valueLow / costRatio_;}

		if(valueHigh > relativeValueLow){
			return make_tuple(pointHigh, valueHigh, true, false);
		
		}else{
			return make_tuple(pointLow, relativeValueLow, false, true);
		}
		
	}
	// If dealing with globalVariance precalculate W matrix for speed
	if(chosenAcquisiton_.compare("globalVarianceWithChoice") == 0){
		setWmatrix();
		// Change name to low and high, and find the best point for each
		chosenAcquisiton_ = "globalVarianceHigh";
		tuple<VectorXd, double, bool, bool> resultHigh = SurrogateModel::findNextSampleSite();
		VectorXd pointHigh = get<0>(resultHigh);
		double valueHigh = get<1>(resultHigh);

		// It is possible that we already have way too many low fi sample points.
		// If this is the case, do not bother looking for promising sample points.
		if((int)sampledPointsLow_.size() >= 50 * biFunction_->d_){
			if(printInfo_){
				printf("Low fi point is skipped as already have %d low-fi points ", (int)sampledPointsLow_.size());
				printf("\nHigh fi point has value %.4f at point ", valueHigh);
				printPoint(pointHigh);
				printf("\n");
			}
			chosenAcquisiton_ = "globalVarianceWithChoice";
			return make_tuple(pointHigh, valueHigh, true, false);
		}

		chosenAcquisiton_ = "globalVarianceLow";
		tuple<VectorXd, double, bool, bool> resultLow = SurrogateModel::findNextSampleSite();
		VectorXd pointLow = get<0>(resultLow);
		double valueLow = get<1>(resultLow);

		if(printInfo_){
			printf("Low fi point has value %.4f at point ", valueLow);
			printPoint(pointLow);
			printf("\nHigh fi point has value %.4f at point ", valueHigh);
			printPoint(pointHigh);
			printf("\n");
		}
		chosenAcquisiton_ = "globalVarianceWithChoice";

		// Here want to choose which source to sample, base it on cost ratio.
		// First check that this value has been set.
		// Default should be -1 (which is incorrect), so just check this
		if(costRatio_ < -0.5){
			printf("Cost ratio needs to be used but it has not yet been set! Stopping now...\n");
			exit(0);
		}
		// Now just divide low fi value by cost ratio for comparison. If the cost ratio is too close to 0, use 0.001 as a stand in
		double relativeValueLow;
		if(abs(costRatio_) < TOL){relativeValueLow = valueLow / 0.001;}
		else{relativeValueLow = valueLow / costRatio_;}

		if(valueHigh > relativeValueLow){
			return make_tuple(pointHigh, valueHigh, true, false);
		
		}else{
			return make_tuple(pointLow, relativeValueLow, false, true);
		}

		return make_tuple(pointHigh, valueHigh, true, false);
	}
	tuple<VectorXd, double, bool, bool> result = Kriging::findNextSampleSite();
	// Return the same, but specify that should sample both sources
	return make_tuple(get<0>(result), get<1>(result), true, true);
}



CoKriging::IntermediateConcentratedLikelihoodFunction::IntermediateConcentratedLikelihoodFunction(int d, vector<double> &lowerBound, vector<double> &upperBound, CoKriging* cokrigingModel):
	Function(d, lowerBound, upperBound),
	cokrigingModel_(cokrigingModel){}

CoKriging::IntermediateConcentratedLikelihoodFunction::~IntermediateConcentratedLikelihoodFunction(){}

double CoKriging::IntermediateConcentratedLikelihoodFunction::evaluate(VectorXd &point){
	// Ok so here idea is that the first half of the entries of point contain the d theta variables, and the second half the p values
	// Should find out dimension and size, then extract theta and p values
	int d = (int)cokrigingModel_->sampledPoints_[0].size();
	vector<double> theta(d,0.0);
	vector<double> pVector(d,0.0);
	for(int i = 0; i < d; i++){
		theta[i] = pow(10,point(i));
		// theta[i] = point(i);
		if(cokrigingModel_->fixedP_){
			pVector[i] = 2;
		}else{
			pVector[i] = point(d + i);
		}
	}
	cokrigingModel_->thetaB_ = theta;
	cokrigingModel_->pBVector_ = pVector;
	if(cokrigingModel_->fixedP_){
		cokrigingModel_->rho_ = point(d);
	}else{
		cokrigingModel_->rho_ = point(2*d);
	}
	
	// Get info needed, data contains mu, sigma, R and R^-1 in that order
	return cokrigingModel_->intermediateConcentratedLikelihoodFunction();
}


int CoKriging::IntermediateConcentratedLikelihoodFunction::betterPoint(VectorXd &point1, double val1, VectorXd &point2, double val2, string flag){
	int d = (int)cokrigingModel_->sampledPoints_[0].size();
	double distance1, distance2;
	if(flag.compare("maxmin") == 0){
		distance1 = DBL_MAX;
		distance2 = DBL_MAX;
		for(int i = 0; i < d; i++){
			double val = pow(10, -point1(i) / point1(d + i));
			if(val < distance1){distance1 = val;}
			val = pow(10, -point2(i) / point2(d + i));
			if(val < distance2){distance2 = val;}
		}
	}else if(flag.compare("maxsum") == 0){
		distance1 = 0;
		distance2 = 0;
		for(int i = 0; i < d; i++){
			double val = pow(10, -point1(i) / point1(d + i));
			distance1 += sqrt(val);
			val = pow(10, -point2(i) / point2(d + i));
			distance2 += sqrt(val);
		}
	}else{
		printf("Problem, asking to compare loglikelihoods without specifying a flag with defined behaviour! Stopping now...\n");
		exit(0);
	}
	
	// If they are both close enough to the best likelihood found so far, compare the distances, prefer the larger one
	if(abs(val1 - cokrigingModel_->maxLogLikelihood_) < OBJTOL && abs(val2 - cokrigingModel_->maxLogLikelihood_) < OBJTOL){
		if(abs(distance1 - distance2) < TOL){return 0;}
		else if(distance1 > distance2){return -1;}
		else{return 1;}

	}else if(abs(val1 - cokrigingModel_->maxLogLikelihood_) < OBJTOL){
		// The first point is close enough to the max likelihood, but the second is not
		return -1;
	}else if(abs(val2 - cokrigingModel_->maxLogLikelihood_) < OBJTOL){
		// The second point is close enough to the max likelihood, but the first one is not
		return 1;
	}else{
		if(abs(val1 - val2) < TOL){
			return 0;
		}else if(val1 < val2){return 1;}
		else{return -1;}
	}
}









AdaptiveCoKriging::AdaptiveCoKriging(BiFidelityFunction* biFunction, AuxSolver* auxSolver, bool advanced, int randomSeed, bool printInfo, bool functionScaling):
	CoKriging(biFunction, auxSolver, randomSeed, printInfo, functionScaling){
	krigingModel_ = new Kriging(biFunction_, auxSolver_, randomSeed_, printInfo_, functionScaling);	
	cokrigingModel_ = new CoKriging(biFunction_, auxSolver_, randomSeed_, printInfo_, functionScaling);	
	lowFiDataIsUseful_ = false;
	advanced_ = advanced;
}

AdaptiveCoKriging::AdaptiveCoKriging(BiFidelityFunction* biFunction, AuxSolver* auxSolver, Kriging* krigingModel, CoKriging* cokrigingModel, bool advanced, int randomSeed, bool printInfo, bool functionScaling):
	CoKriging(biFunction, auxSolver, randomSeed, printInfo, functionScaling){
	krigingModel_ = krigingModel;	
	cokrigingModel_ = cokrigingModel;	
	lowFiDataIsUseful_ = false;
	advanced_ = advanced;
}

AdaptiveCoKriging::~AdaptiveCoKriging(){}

void AdaptiveCoKriging::assessLowFiFunction(){
	// Ok so first thing is to get the locations at which the high fidelity function is known,
	// the high fidelity function values, and the low fidelity function values (which are predicted 
	// if the true value is not known.)
	int n = (int)sampledPoints_.size();
	int d = biFunction_->d_;

	vector<VectorXd> locations = sampledPoints_;
	vector<double> highFiValues = sampledPointsValues_;
	vector<double> lowFiValues;
	for(int i = 0; i < n; i++){
		lowFiValues.push_back(cokrigingModel_->lowFiKriging_->surfaceValue(sampledPoints_[i], true, false));
	}
	vector<double> differences;
	for(int i = 0; i < n; i++){
		differences.push_back(highFiValues[i] - lowFiValues[i]);
	}

	if(printInfo_){printf("Assessing low fidelity function. Calculate Adjusted R^2 for a linear model\n");}
	leastSquaresFunction* function = new leastSquaresFunction(d + 1, locations, differences);
	auxSolver_->updateProblem(function, true);
	VectorXd coefficients = auxSolver_->optimise();
	// Now want the mean, sum of squares, and sum of residuals
	double mean = 0;
	double min = differences[0];
	double max = differences[0];
	for(int i = 0; i < n; i++){
		mean += differences[i];
		if(min > differences[i]){min = differences[i];}
		if(max < differences[i]){max = differences[i];}
	}
	mean = mean / (int)differences.size();
	double ssReg = 0;
	for(int i = 0; i < n; i++){
		ssReg += pow(differences[i] - mean, 2);
	}

	// for(int i = 0; i < n; i++){
	// 	// Get the prediction
	// 	double prediction = coefficients[0];
	// 	for(int j = 0; j < d; j++){
	// 		prediction += coefficients[j + 1] * locations[i](j);
	// 	}
	// 	printf("Prediction %.4f vs real %.4f\n", prediction, differences[i]);
	// 	// ssReg += pow(prediction - mean, 2);
	// }
	double ssRes = function->evaluate(coefficients);

	// printf("Got ssReg %.20f and ssRes %.20f\n", ssReg, ssRes);
	// double ssTot = ssRes + ssReg;
	// double adjustedR2 = 1 - (ssRes / ssTot) / (((double)n - d - 1) / ((double)n - 1));
	double adjustedR2;
	// If have very little data the model should be able to fit it exactly
	if(n <= d + 1){adjustedR2 = 1;}
	// If data is constant should also fit it perfectly (and will divide by 0 so need to set 1 by hand)
	else if(abs(min - max) < TOL){adjustedR2 = 1;}
	// Otherwise can use formula
	else{adjustedR2 = 1 - (((double)n - 1) / ((double)n - d - 1)) * (ssRes / ssReg);}
	
	

	// Ok so need sample to be unscaled when passing it to calculateLocalCorrelations, as the function scales internally
	unscalePoints(locations);
	unscaleObservations(highFiValues);
	unscaleObservations(lowFiValues);
	pair<vector<double>, vector<double> > lccs = calculateLocalCorrelations(biFunction_, 0.2, {0.4, 0.95}, locations, highFiValues, lowFiValues);
	
	double LCC_4_rel = lccs.second[0];
	double LCC_95_rel = lccs.second[1];
	double brh = locations.size() / (double)biFunction_->d_;
	double br = locations.size() / (double)sampledPointsLow_.size();

	// Ok I think I have it all! Print them and make a decision
	if(printInfo_){printf("Have feature values of Brh %.2f, Br %.2f, LCC_0.4 %.2f, LCC_0.95 %.2f and R2_L %.2f\n", brh, br, LCC_4_rel, LCC_95_rel, adjustedR2);}

	if(!advanced_){

		if(brh >= 18 || br >= 1 || LCC_4_rel <= 0.7){
			if(printInfo_){printf("Instance belongs to first group, should NOT use low fi function.\n");}
			lowFiDataIsUseful_ = false;
		}else if(LCC_95_rel >= 0.5 || adjustedR2 >= 0.4){
			if(printInfo_){printf("Instance belongs to second group, should use low fi function.\n");}
			lowFiDataIsUseful_ = true;
		}else if(brh <= 5){
			if(printInfo_){printf("Instance belongs to third group and have little high fi data, should use low fi function.\n");}
			lowFiDataIsUseful_ = true;
		}else{
			if(printInfo_){printf("Instance belongs to third group and have a fair amount of high fi data, should NOT use low fi function.\n");}
			lowFiDataIsUseful_ = false;
		}
	
	}else{
		// Ok so idea is to base decisions of brh also on previous decisions of low fi data usefulness
		if(br >= 1 || LCC_4_rel <= 0.7){
			if(printInfo_){printf("Instance belongs to first group (not based on brh), should NOT use low fi function.\n");}
			lowFiDataIsUseful_ = false;
		
		}else if(brh >= 18 && !lowFiDataIsUseful_){
			if(printInfo_){printf("Instance belongs to first group based on Brh and in previous iteration did not use low fi, should NOT use low fi function.\n");}
			lowFiDataIsUseful_ = false;
		
		}else if(LCC_95_rel >= 0.5 || adjustedR2 >= 0.4){
			if(printInfo_){printf("Instance belongs to second group, should use low fi function.\n");}
			lowFiDataIsUseful_ = true;
		}else{
			if(printInfo_){
				printf("Instance belongs to third group, will continue as before: ");
				if(lowFiDataIsUseful_){printf("Use low fi function.\n");}
				else{printf("Do NOT use low fi function.\n");}
			}
		}
	}
}

void AdaptiveCoKriging::setCostRatio(double costRatio){
	krigingModel_->setCostRatio(costRatio);
	cokrigingModel_->setCostRatio(costRatio);
}

void AdaptiveCoKriging::saveSample(vector<VectorXd> points, vector<VectorXd> pointsLow, vector<double> observations, vector<double> observationsLow){
	// Sounds like what I should do is simply call the save sample function on each of Kriging and CoKriging
	krigingModel_->saveSample(points, pointsLow, observations, observationsLow);
	cokrigingModel_->saveSample(points, pointsLow, observations, observationsLow);
	// Not going to call the parent function as don't want to save anything internally.
	SurrogateModel::saveSample(points, pointsLow, observations, observationsLow);
}


void AdaptiveCoKriging::trainModel(bool forcedOptimiseHyperparameters){
	// Ok so first need to train Co-Kriging model
	if(printInfo_){printf("Training cokriging model to estimate quality of low fi source.\n");}
	cokrigingModel_->trainModel(forcedOptimiseHyperparameters);
	// Now using the predicted low fi values at known high fi locations assess whether the low fi function looks usable
	assessLowFiFunction();
	// If the low fi function is usable, no need to train a Kriging model. If not usable, 
	// need to train a Kriging model. 
	// I think that whenever hyperparameters are optimised I should train the Kriging model just to keep it
	// "up to date".
	if(lowFiDataIsUseful_){
		if(printInfo_){printf("Not training kriging model as low-fi looks benefitial.\n");}
	}else{
		if(printInfo_){printf("Training kriging model as low-fi looks harmful.\n");}
		krigingModel_->trainModel(forcedOptimiseHyperparameters);
	}
	printf("Training complete!\n");
}


tuple<double, double> AdaptiveCoKriging::meanVarianceCalculator(VectorXd &x){
	if(lowFiDataIsUseful_){
		return cokrigingModel_->meanVarianceCalculator(x);
	}else{
		return krigingModel_->meanVarianceCalculator(x);
	}
}

void AdaptiveCoKriging::setAcquisitionFunction(string chosenAcquisiton){
	cokrigingModel_->setAcquisitionFunction(chosenAcquisiton);
	if(chosenAcquisiton.compare("globalVarianceWithChoice") == 0){
		krigingModel_->setAcquisitionFunction("globalVariance");
	}else if(chosenAcquisiton.compare("globalVarianceApproximationWithChoice") == 0){
		krigingModel_->setAcquisitionFunction("globalVarianceApproximation");
	}else{
		krigingModel_->setAcquisitionFunction(chosenAcquisiton);
	}
}

tuple<VectorXd, double, bool, bool> AdaptiveCoKriging::findNextSampleSite(){
	if(lowFiDataIsUseful_){
		return cokrigingModel_->findNextSampleSite();
	}else{
		return krigingModel_->findNextSampleSite();
	}
}


AdaptiveCoKriging::leastSquaresFunction::leastSquaresFunction(int d, vector<VectorXd> points, vector<double> values):
	Function(d, -100, 100),
	points_(points),
	values_(values){}

AdaptiveCoKriging::leastSquaresFunction::~leastSquaresFunction(){}

double AdaptiveCoKriging::leastSquaresFunction::evaluate(VectorXd &point){
	// Just want to find the prediction at each point, and get the error squared
	double sumError = 0;
	for(int i = 0; i < (int)points_.size(); i++){
		double prediction = point[0];
		for(int j = 0; j < (int)points_[0].size(); j++){
			prediction += point(j + 1) * points_[i](j);
		}
		sumError += pow(prediction - values_[i], 2);
	}
	return sumError;
}




#endif