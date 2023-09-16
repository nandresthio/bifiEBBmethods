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
}

SurrogateModel::~SurrogateModel(){
	if(sampleGenerator_ != NULL){delete sampleGenerator_;}
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

void SurrogateModel::saveSample(vector<VectorXd> &points, vector<VectorXd> &pointsLow, vector<double> &observations, vector<double> &observationsLow){
	scalePoints(points);
	scalePoints(pointsLow);
	scaleTwoObservations(observations, observationsLow);
	sampledPoints_ = points;
	sampledPointsLow_ = pointsLow;
	sampledPointsValues_ = observations;
	sampledPointsValuesLow_ = observationsLow;
	
	
}


void SurrogateModel::trainModel(){
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











Kriging::Kriging(BiFidelityFunction* biFunction, AuxSolver* auxSolver, int randomSeed, bool printInfo, bool functionScaling):
	SurrogateModel(biFunction, auxSolver, randomSeed, printInfo, functionScaling){

	theta_.reserve(biFunction_->d_);
	pVector_.reserve(biFunction_->d_);
	
}

Kriging::~Kriging(){}


void Kriging::trainModel(){
	if(printInfo_){printf("Training kriging model.\n");}
	trainHyperparameters();
	saveMuSigma();	
	trainedModel_ = true;
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
			upperBound[i] = 3;
		}else{
			lowerBound[i] = -10;
			upperBound[i] = 3;
		}
		lowerBound[d+i] = 0.1;
		upperBound[d+i] = 2.0;
	}
	maxLogLikelihood_ = -DBL_MAX;
	ConcentratedLikelihoodFunction* function = new ConcentratedLikelihoodFunction(2*d, lowerBound, upperBound, this);
	auxSolver_->updateProblem(function, false);
	VectorXd hyperparameters = auxSolver_->optimise();
	delete function;
	// Extract two vectors, and store!
	for(int i = 0; i < d; i++){
		theta_[i] = pow(10,hyperparameters(i));
		pVector_[i] = hyperparameters(d + i);
	}
}

void Kriging::saveMuSigma(){
	auto data = muSigmaCalculator();
	mu_ = get<0>(data);
	sigma_ = get<1>(data);
	rMatrix_ = get<2>(data);
	rMatrixDecomposition_ = get<3>(data);
	modelWeights_ = get<4>(data);

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
	// Need to decide whether we are working with a scaled model or not.
	// If yes, the min value is 0 (as it is scaled).
	// Otherwise need to use the true minimum
	if(unscaleOutput || !functionScaling_){
		value = (fMin_ - s_x)/sqrt(var);
	}else{
		value = (0 - s_x)/sqrt(var);
	}
	double improvement = sqrt(var)*(value * normalCDF(value) + normalPDF(value));
	return improvement;

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
	MatrixXd var_mid = r.transpose() * rMatrixDecomposition_.solve(r);
	double variance = sigma_ * (1 - var_mid(0,0));
	
	return make_tuple(s_x, variance);
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
		pVector[i] = point(d + i);
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



// VectorXd Kriging::chooseNextSampleSite(string technique){
// 	if(printInfo_){printf("Choosing next sample.");}
// 	bool min;
// 	// When this is called, first define any values that might be needed, then define next sample function, and optimise it!
// 	if(technique.compare("expectedImprovement") == 0){
// 		// Will need to know the minimum.
// 		// This should be known, as when adding a sample the values are automatically scaled
// 		min = false;
// 	}else if(technique.compare("maximumVariance") == 0){
// 		// Should not need any extra information
// 		min = false;
// 	}else{
// 		printf("Asking for next sample site with undefined technique! Stopping now...\n");
// 		exit(0);
// 	}

// 	// Define function, update aux solver, and find the point!
// 	nextSampleSiteFunction* function = new nextSampleSiteFunction(ebbFunction_->d_, ebbFunction_->lowerBound_, ebbFunction_->upperBound_, this, technique);
// 	auxSolver_->updateProblem(function, min);
// 	VectorXd nextSite = auxSolver_->optimise();
// 	if(printInfo_){
// 		printf(" Chosen point ");
// 		printPoint(nextSite);
// 		printf(" with value %.5f\n", function->evaluate(nextSite));
// 	}
// 	delete function;
// 	return nextSite;
// }




// Kriging::nextSampleSiteFunction::nextSampleSiteFunction(int d, vector<double> &lowerBound, vector<double> &upperBound, Kriging* model, string technique):
// 	Function(d, lowerBound, upperBound),
// 	model_(model),
// 	technique_(technique){}

// Kriging::nextSampleSiteFunction::~nextSampleSiteFunction(){}

// double Kriging::nextSampleSiteFunction::evaluate(VectorXd &point){
// 	// Simply query the right function from the model
// 	if(technique_.compare("expectedImprovement") == 0){
// 		return model_->expectedImprovement(point);
// 	}else if(technique_.compare("maximumVariance") == 0){
// 		return model_->uncertainty(point);
// 	}else{
// 		printf("Using next sample function with undefined sampling strategy! Stopping now...\n");
// 		exit(0);
// 	}
// 	return 0;
// }























CoKriging::CoKriging(BiFidelityFunction* biFunction, AuxSolver* auxSolver, int randomSeed, bool printInfo, bool functionScaling):
	Kriging(biFunction, auxSolver, randomSeed, printInfo, functionScaling){

	lowFiKriging_ = new Kriging(biFunction_, auxSolver_, randomSeed_, printInfo_, functionScaling);
	
	thetaB_.reserve(biFunction_->d_);
	pBVector_.reserve(biFunction_->d_);

}

CoKriging::~CoKriging(){}


void CoKriging::saveSample(vector<VectorXd> &points, vector<VectorXd> &pointsLow, vector<double> &observations, vector<double> &observationsLow){
	// First to the usual saving of the data
	SurrogateModel::saveSample(points, pointsLow, observations, observationsLow);
	// And now save the relevant things to lowFiKriging
	lowFiKriging_->sampledPoints_ = pointsLow;
	lowFiKriging_->sampledPointsValues_ = observationsLow;
	lowFiKriging_->maxObservation_ = maxObservation_;
	lowFiKriging_->minObservation_ = minObservation_;
	
}

void CoKriging::trainModel(){
	if(printInfo_){printf("Training cokriging model, start with low fi kriging.\n");}
	lowFiKriging_->trainModel();
	if(printInfo_){printf("Training cokriging model, train intermediate model.\n");}
	trainHyperparameters();
	saveMuSigma();	
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

	// Store low fi surface values at high fi locations
	lowFiSurfaceValues_.reserve((int)sampledPoints_.size());
	for(int i = 0; i < (int)sampledPoints_.size(); i++){
		// Don't want to return scaled predictions from the low fi model, as internally everything is scaled (if it was requested by the user)
		lowFiSurfaceValues_.push_back(lowFiKriging_->surfaceValue(sampledPoints_[i], true, false));
	}
	maxLogLikelihood_ = -DBL_MAX;
	IntermediateConcentratedLikelihoodFunction* function = new IntermediateConcentratedLikelihoodFunction(2*d + 1, lowerBound, upperBound, this);
	auxSolver_->updateProblem(function, false);
	VectorXd hyperparameters = auxSolver_->optimise();
	delete function;
	// Extract two vectors, and store!
	for(int i = 0; i < d; i++){
		thetaB_[i] = pow(10,hyperparameters(i));
		pBVector_[i] = hyperparameters(d + i);
	}
	// Last value is rho
	rho_ = hyperparameters(2*d);
}

void CoKriging::saveMuSigma(){
	auto data = intermediateMuSigmaCalculator();
	sigmaL_ = lowFiKriging_->sigma_;
	sigmaB_ = get<1>(data);
	auto data1 = combinedMuCmatrixCalculator();
	muCombined_ = get<0>(data1);
	cMatrix_ = get<1>(data1);
	cMatrixDecomposition_ = get<2>(data1);
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
		// b(i) = sampledPointsValues_[i] - rho_ * lowFiKriging_->surfaceValue(sampledPoints_[i]);
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
	MatrixXd mu_top = one.transpose() * cMatrixDecomposition.solve(y);
	MatrixXd mu_bottom = one.transpose() * cMatrixDecomposition.solve(one);
	double mu = mu_top(0,0)/ mu_bottom(0,0);
	return make_tuple(mu, cMatrix, cMatrixDecomposition);
}

double CoKriging::intermediateConcentratedLikelihoodFunction(){
	auto data = intermediateMuSigmaCalculator();
	int n = (int)sampledPoints_.size();
	double logLikelihood = -n * logl(get<1>(data))/2 - logl((get<3>(data)).matrixL().determinant()) - logl((get<3>(data)).vectorD().prod())/2;

	if(isinf(logLikelihood) || isnan(logLikelihood)){return -DBL_MAX;}
	// if(get<1>(data) < TOL){return -DBL_MAX;}
	if(maxLogLikelihood_ < logLikelihood){
		maxLogLikelihood_ = logLikelihood;
	}

	return logLikelihood;
}

tuple<double, double> CoKriging::meanVarianceCalculator(VectorXd &x){
	int nL = (int)lowFiKriging_->sampledPoints_.size();
	int nH = (int)sampledPoints_.size();
	int d = (int)sampledPoints_[0].size();
	// First want c vector
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
	VectorXd rightHandSide = c.transpose() * cMatrixDecomposition_.solve(y - one * muCombined_);
	double s_x = muCombined_ + rightHandSide(0);
	rightHandSide = c.transpose() * cMatrixDecomposition_.solve(c);
	double variance = rho_ * rho_ * sigmaL_ + sigmaB_ - rightHandSide(0);
	
	return make_tuple(s_x, variance);
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
		pVector[i] = point(d + i);
	}
	cokrigingModel_->thetaB_ = theta;
	cokrigingModel_->pBVector_ = pVector;
	cokrigingModel_->rho_ = point(2*d);
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










// CODE BELOW IS UNDER DEVELOPMENT!





// CoKrigingCheap::CoKrigingCheap(BiFidelityFunction* biFunction, AuxSolver* auxSolver, int highFiSampleBudget, char mode, int randomSeed, bool printInfo, bool functionScaling):
// 	Kriging(biFunction, auxSolver, 0, randomSeed, printInfo, functionScaling), 
// 	biFunction_(biFunction),
// 	highFiSampleBudget_(highFiSampleBudget),
// 	mode_(mode){

// 	if(mode_ != '0' && mode_ != 's'){
// 		printf("Initialised CoKriging cheap low fi model with undefined mode! Setting model to 's'...\n");
// 		mode_ = 's';
// 	}

// 	thetaB_.reserve(ebbFunction_->d_);
// 	pBVector_.reserve(ebbFunction_->d_);

// }

// CoKrigingCheap::~CoKrigingCheap(){}

// void CoKrigingCheap::trainHyperparameters(){
// 	int d = biFunction_->d_;
// 	// First define the bounds for the function
// 	vector<double> lowerBound(2*d + 1, 0.0);
// 	vector<double> upperBound(2*d + 1, 0.0);
// 	// Bounds should be dependent on whether data has been scaled; if it has the bounds can be more restrictive
// 	for(int i = 0; i < d; i++){
// 		// Bounds for theta_b
// 		if(functionScaling_){
// 			lowerBound[i] = -3;
// 			upperBound[i] = 3;
// 		}else{
// 			lowerBound[i] = -10;
// 			upperBound[i] = 3;
// 		}
// 		// Bounds for p_b
// 		lowerBound[d+i] = 0.1;
// 		upperBound[d+i] = 2.0;
// 	}
// 	// Bounds for rho
// 	lowerBound[2*d] = -1000;
// 	upperBound[2*d] = 1000;

// 	maxLogLikelihood_ = -DBL_MAX;
// 	IntermediateConcentratedLikelihoodFunction* function = new IntermediateConcentratedLikelihoodFunction(2*d + 1, lowerBound, upperBound, this);
// 	auxSolver_->updateProblem(function, false);
// 	VectorXd hyperparameters = auxSolver_->optimise();
// 	delete function;
// 	// Extract two vectors, and store!
// 	for(int i = 0; i < d; i++){
// 		thetaB_[i] = pow(10,hyperparameters(i));
// 		pBVector_[i] = hyperparameters(d + i);
// 	}
// 	// Last value is rho
// 	rho_ = hyperparameters(2*d);
// }

// void CoKrigingCheap::saveMuSigma(){
// 	auto data = intermediateMuSigmaCalculator();
// 	muIntermediate_ = get<0>(data);
// 	sigmaB_ = get<1>(data);
// 	bMatrix_ = get<2>(data);
// 	bMatrixDecomposition_ = get<3>(data);
// 	mu_ = muCalculator();
// 	return;
// }

// tuple<double, double, MatrixXd, LDLT<MatrixXd>> CoKrigingCheap::intermediateMuSigmaCalculator(){	
// 	int n = (int)sampledPoints_.size();
// 	int d = (int)sampledPoints_[0].size();
// 	if(n < 1){
// 		printf("Trying to train CoKriging model with less points than needed! Have %d but need %d. Exiting now...\n", n, 1);
// 		exit(0);
// 	}

// 	MatrixXd bExpensiveMatrix(n,n);
// 	// First of all need to define matrix
// 	auto it1 = sampledPoints_.begin();
// 	for(int i = 0; i < n; i++, it1++){
// 		auto it2 = it1;
// 		for(int j = i; j < n; j++, it2++){
// 			double sum = 0;
// 			for(int k = 0; k < d; k++){
// 				sum += thetaB_[k] * pow(abs((*it1)(k) - (*it2)(k)), pBVector_[k]);
// 			}
// 			bExpensiveMatrix(i,j) = exp(-sum);
// 			bExpensiveMatrix(j,i) = bExpensiveMatrix(i,j);	
// 		}
// 	}
// 	VectorXd one(n);
// 	VectorXd b(n);
// 	for(int i = 0; i < n; i++){
// 		one(i) = 1;
// 		// Next is key, low fi function is cheap so evaluate it
// 		double val = biFunction_->evaluateLow(sampledPoints_[i]);
// 		val = scaleObservation(val);
// 		b(i) = sampledPointsValues_[i] - rho_ * val; 
// 		// b(i) = sampledPointsValues_[i] - rho_ * biFunction_->evaluateLow(sampledPoints_[i]);
// 	}
// 	LDLT<MatrixXd> bExpensiveMatrixDecomposition = bExpensiveMatrix.ldlt();
// 	MatrixXd mu_top = one.transpose() * bExpensiveMatrixDecomposition.solve(b);
// 	MatrixXd mu_bottom = one.transpose() * bExpensiveMatrixDecomposition.solve(one);
// 	double mu = mu_top(0,0)/ mu_bottom(0,0);
// 	MatrixXd sigma_m = (b - one*mu).transpose() * bExpensiveMatrixDecomposition.solve(b - one*mu);
// 	double sigma = sigma_m(0,0)/n;

// 	return make_tuple(mu, sigma, bExpensiveMatrix, bExpensiveMatrixDecomposition);
// }

// double CoKrigingCheap::muCalculator(int nL){
// 	int n = (int)sampledPoints_.size();
// 	VectorXd one(n);
// 	VectorXd b(n);
// 	for(int i = 0; i < n; i++){
// 		one(i) = 1;
// 		double val = biFunction_->evaluateLow(sampledPoints_[i]);
// 		val = scaleObservation(val);
// 		b(i) = sampledPointsValues_[i] - rho_ * val; 
// 		// b(i) = sampledPointsValues_[i] - rho_ * biFunction_->evaluateLow(sampledPoints_[i]);
// 	}


// 	if(mode_ == '0'){
// 		MatrixXd mu_top = one.transpose() * bMatrixDecomposition_.solve(b);
// 		MatrixXd mu_bottom = one.transpose() * bMatrixDecomposition_.solve(one);
// 		double mu = mu_top(0,0) / ((1 - rho_) * mu_bottom(0,0));
		
// 		return mu;
	
// 	}else if(mode_ == 's'){
// 		vector<VectorXd> lowFiSample = sampleGenerator_->randomLHS(nL);
// 		vector<double> lowFiSampleVals = biFunction_->evaluateManyLow(lowFiSample);
// 		scaleObservations(lowFiSampleVals);

// 		double sumL = 0;
// 		for(int i = 0; i < nL; i++){
// 			sumL += lowFiSampleVals[i];
// 		}
// 		double muL = sumL / nL;
// 		double sigmaL = 0;
// 		for(int i = 0; i < nL; i++){
// 			sigmaL += (lowFiSampleVals[i] - muL)*(lowFiSampleVals[i] - muL);
// 		}
// 		sigmaL = sigmaL / nL;
// 		MatrixXd mu_top = one.transpose() * bMatrixDecomposition_.solve(b);
// 		MatrixXd mu_bottom = one.transpose() * bMatrixDecomposition_.solve(one);
// 		double mu;
// 		if(sigmaL < TOL){
// 			mu = ((1 - rho_) * mu_top(0,0) / sigmaB_)/((1 - rho_) * (1 - rho_) * mu_bottom(0,0) / sigmaB_);
// 		}else{
// 			mu = (sumL / sigmaL + (1 - rho_) * mu_top(0,0) / sigmaB_)/(nL / sigmaL + (1 - rho_) * (1 - rho_) * mu_bottom(0,0) / sigmaB_);
// 		}
		
// 		return mu;
	
// 	}else{
// 		printf("Weird, caluclating mu with undefined mode for cheap Co-Kriging. Stopping now...\n");
// 		exit(0);
// 	}

// 	return 0.0;
// }

// double CoKrigingCheap::intermediateConcentratedLikelihoodFunction(){
// 	auto data = intermediateMuSigmaCalculator();
// 	int n = (int)sampledPoints_.size();
// 	double logLikelihood = -n * logl(get<1>(data))/2 - logl((get<3>(data)).matrixL().determinant()) - logl((get<3>(data)).vectorD().prod())/2;

// 	if(isinf(logLikelihood) || isnan(logLikelihood)){return -DBL_MAX;}
// 	if(maxLogLikelihood_ < logLikelihood){
// 		maxLogLikelihood_ = logLikelihood;
// 	}
// 	// We are going to say that if the determinant of b is 0.999, it is assumed to be the identity, so nothing larger that 0.999 is allowed. 
// 	// This should lead to a log likelihood not impacted by the determinant of b, since 0 > log |b| > -0.000435 for |b| > 0.999
// 	// MatrixXd bMatrix = get<2>(data);
// 	// Same with R being the identity, if it is too close to the identity ignore it
// 	// This is done so that something very close to the identity can still be found, but anything "closer" to it
// 	// which lead to an almost identical model is ignored
// 	// bool isIdentity = true;
// 	// for(int i = 0; isIdentity && i < bMatrix.rows(); i++){
// 	// 	for(int j = i + 1; isIdentity && j < bMatrix.cols(); j++){
// 	// 		if(bMatrix(i,j) > 0.001){
// 	// 			isIdentity = false;
// 	// 		}
// 	// 	}
// 	// }
// 	// if(isIdentity){return -DBL_MAX;}
// 	// if(bMatrix.determinant() > 0.999){return -DBL_MAX;}
// 	return logLikelihood;
// }

// tuple<double, double> CoKrigingCheap::meanVarianceCalculator(VectorXd &x){
// 	int nH = (int)sampledPoints_.size();
// 	int d = (int)sampledPoints_[0].size();
// 	// First want c vector
// 	VectorXd c(nH);
// 	auto it = sampledPoints_.begin();
// 	for(int i = 0; i < nH; i++, it++){
// 		double sumB = 0;
// 		for(int k = 0; k < d; k++){
// 			sumB += thetaB_[k] * pow(abs((*it)(k) - x(k)), pBVector_[k]);
// 		}
// 		c(i) = exp(-sumB);
		
// 	}

// 	VectorXd one(nH);
// 	VectorXd b(nH);
// 	for(int i = 0; i < nH; i++){
// 		one(i) = 1;
// 		double val = biFunction_->evaluateLow(sampledPoints_[i]);
// 		val = scaleObservation(val);
// 		b(i) = sampledPointsValues_[i] - rho_ * val; 
// 		// b(i) = sampledPointsValues_[i] - rho_ * biFunction_->evaluateLow(sampledPoints_[i]);
// 	}

	
// 	VectorXd rightHandSide = c.transpose() * bMatrixDecomposition_.solve(b - one * (1 - rho_) * mu_);
// 	double s_x = (1 - rho_) * mu_ + rho_ * biFunction_->evaluateLow(x) + rightHandSide(0);
// 	rightHandSide = c.transpose() * bMatrixDecomposition_.solve(c);
// 	double variance = sigmaB_ * (1 - rightHandSide(0));
	
// 	return make_tuple(s_x, variance);
// }

// CoKrigingCheap::IntermediateConcentratedLikelihoodFunction::IntermediateConcentratedLikelihoodFunction(int d, vector<double> &lowerBound, vector<double> &upperBound, CoKrigingCheap* cokrigingCheapModel):
// 	Function(d, lowerBound, upperBound),
// 	cokrigingCheapModel_(cokrigingCheapModel){}

// CoKrigingCheap::IntermediateConcentratedLikelihoodFunction::~IntermediateConcentratedLikelihoodFunction(){}

// double CoKrigingCheap::IntermediateConcentratedLikelihoodFunction::evaluate(VectorXd &point){
// 	// Ok so here idea is that the first half of the entries of point contain the d theta variables, and the second half the p values
// 	// Should find out dimension and size, then extract theta and p values
// 	int d = (int)cokrigingCheapModel_->sampledPoints_[0].size();
// 	vector<double> theta(d,0.0);
// 	vector<double> pVector(d,0.0);
// 	for(int i = 0; i < d; i++){
// 		theta[i] = pow(10,point(i));
// 		// theta[i] = point(i);
// 		pVector[i] = point(d + i);
// 	}
// 	cokrigingCheapModel_->thetaB_ = theta;
// 	cokrigingCheapModel_->pBVector_ = pVector;
// 	cokrigingCheapModel_->rho_ = point(2*d);
// 	// Get info needed, data contains mu, sigma, R and R^-1 in that order
// 	return cokrigingCheapModel_->intermediateConcentratedLikelihoodFunction();
// }

// int CoKrigingCheap::IntermediateConcentratedLikelihoodFunction::betterPoint(VectorXd &point1, double val1, VectorXd &point2, double val2, string flag){
// 	int d = (int)cokrigingCheapModel_->sampledPoints_[0].size();
// 	double distance1, distance2;
// 	if(flag.compare("maxmin") == 0){
// 		distance1 = DBL_MAX;
// 		distance2 = DBL_MAX;
// 		for(int i = 0; i < d; i++){
// 			double val = pow(10, -point1(i) / point1(d + i));
// 			if(val < distance1){distance1 = val;}
// 			val = pow(10, -point2(i) / point2(d + i));
// 			if(val < distance2){distance2 = val;}
// 		}
// 	}else if(flag.compare("maxsum") == 0){
// 		distance1 = 0;
// 		distance2 = 0;
// 		for(int i = 0; i < d; i++){
// 			double val = pow(10, -point1(i) / point1(d + i));
// 			distance1 += sqrt(val);
// 			val = pow(10, -point2(i) / point2(d + i));
// 			distance2 += sqrt(val);
// 		}
// 	}else{
// 		printf("Problem, asking to compare loglikelihoods without specifying a flag with defined behaviour! Stopping now...\n");
// 		exit(0);
// 	}
	
// 	// If they are both close enough to the best likelihood found so far, compare the distances, prefer the larger one
// 	if(abs(val1 - cokrigingCheapModel_->maxLogLikelihood_) < OBJTOL && abs(val2 - cokrigingCheapModel_->maxLogLikelihood_) < OBJTOL){
// 		if(abs(distance1 - distance2) < TOL){return 0;}
// 		else if(distance1 > distance2){return -1;}
// 		else{return 1;}

// 	}else if(abs(val1 - cokrigingCheapModel_->maxLogLikelihood_) < OBJTOL){
// 		// The first point is close enough to the max likelihood, but the second is not
// 		return -1;
// 	}else if(abs(val2 - cokrigingCheapModel_->maxLogLikelihood_) < OBJTOL){
// 		// The second point is close enough to the max likelihood, but the first one is not
// 		return 1;
// 	}else{
// 		if(abs(val1 - val2) < TOL){
// 			return 0;
// 		}else if(val1 < val2){return 1;}
// 		else{return -1;}
// 	}
// }










// AdaptiveCoKriging::AdaptiveCoKriging(BiFidelityFunction* biFunction, string strategy, double pFactor, double rFactor, AuxSolver* auxSolver, int highFiSampleBudget, int lowFiSampleBudget, int randomSeed, bool printInfo, bool functionScaling):
// 	CoKriging(biFunction, auxSolver, highFiSampleBudget, lowFiSampleBudget, randomSeed, printInfo, functionScaling),
// 	strategy_(strategy),
// 	pFactor_(pFactor),
// 	rFactor_(rFactor){

// 	lowFiKriging_ = new Kriging(biFunction_, auxSolver_, lowFiSampleBudget_, randomSeed_, printInfo_, functionScaling);
// 	highFiKriging_ = new Kriging(biFunction_, auxSolver, highFiSampleBudget_, randomSeed_, printInfo_, functionScaling);
	
// 	thetaB_.reserve(ebbFunction_->d_);
// 	pBVector_.reserve(ebbFunction_->d_);

// }

// AdaptiveCoKriging::~AdaptiveCoKriging(){}


// void AdaptiveCoKriging::saveSample(vector<VectorXd> sample, vector<double> sampleVals, vector<VectorXd> sampleLow, vector<double> sampleValsLow){
// 	sampledPoints_ = sample;
// 	sampledPointsValues_ = sampleVals;
// 	sampledPointsLow_ = sampleLow;
// 	sampledPointsValuesLow_ = sampleValsLow;
// 	scalePoints(sampledPoints_);
// 	scalePoints(sampledPointsLow_);
// 	scaleTwoObservations(sampledPointsValues_, sampledPointsValuesLow_);
// 	highFiKriging_->maxObservation_ = maxObservation_;
// 	highFiKriging_->minObservation_ = minObservation_;
// }

// void AdaptiveCoKriging::sampleAtLocation(VectorXd point){
// 	// First of all scale the point
// 	VectorXd scaledCopy = point;
// 	scalePoint(scaledCopy);
// 	// First check that the point is not repeated
// 	bool isRepeated = checkIfRepeated(sampledPoints_, scaledCopy);
// 	bool isRepeatedLow = checkIfRepeated(sampledPointsLow_, scaledCopy);
// 	if(isRepeated){
// 		printf("Not saving high fi value of point as it is already used by the model!\n");
// 	}else{
// 		double valHigh = biFunction_->evaluate(point);
// 		sampledPoints_.push_back(scaledCopy);
// 		if(valHigh < minObservation_ || valHigh > maxObservation_){
// 			unscaleTwoObservations(sampledPointsValues_, lowFiKriging_->sampledPointsValues_);
// 			sampledPointsValues_.push_back(valHigh);
// 			scaleTwoObservations(sampledPointsValues_, lowFiKriging_->sampledPointsValues_);
// 			highFiKriging_->maxObservation_ = maxObservation_;
// 			highFiKriging_->minObservation_ = minObservation_;
// 		}else{
// 			valHigh = scaleObservation(valHigh);
// 			sampledPointsValues_.push_back(valHigh);
// 		}

// 	}
// 	if(isRepeatedLow){
// 		printf("Not saving high fi value of point as it is already used by the model!\n");
// 	}else{
// 		double valLow = biFunction_->evaluateLow(point);
// 		sampledPointsLow_.push_back(scaledCopy);
// 		valLow = scaleObservation(valLow);
// 		sampledPointsValuesLow_.push_back(valLow);
// 	}
// }



// void AdaptiveCoKriging::chooseReliableLowFidelityData(double p, double r){
// 	printf("Working with p %.2f and r %.2f\n", p, r);
// 	// NOTE: Will need to deal with the case where f_l = -f_h
	
// 	// Clear vectors containing chosen points, as doing the selection again.
// 	usedSampledPointsLow_.clear();
// 	usedSampledPointsValuesLow_.clear();
// 	usedSampledPointsLow_.reserve(sampledPointsLow_.size());
// 	usedSampledPointsValuesLow_.reserve(sampledPointsLow_.size());

// 	// First thing will want is to train the lowFiModel and the highFiModel with all of the available data
// 	lowFiKriging_->sampledPoints_ = sampledPointsLow_;
// 	lowFiKriging_->sampledPointsValues_ = sampledPointsValuesLow_;
// 	highFiKriging_->sampledPoints_ = sampledPoints_;
// 	highFiKriging_->sampledPointsValues_ = sampledPointsValues_;

// 	// The locations at which the high and low fidelity function values will be estimated
// 	vector<VectorXd> informationLocation;
// 	vector<double> informationValsHigh;
// 	vector<double> informationValsLow;
// 	// Seems like behaviour should always be to get LCC at the low fidelity sample sites.
// 	// Chose the locations which will be used to calculate LCC

// 	// Can be smart about what model to train, but for now to get it all working just train both always
// 	if(printInfo_){printf("Choosing reliable low fi data, train two kriging models on high and low fidelity data only.\n");}
// 	lowFiKriging_->trainModel();
// 	highFiKriging_->trainModel();
// 	if(strategy_.compare("highLoc") == 0){
// 		// if(printInfo_){
// 		// 	printf("Train low fi model to decide what low fidelity data is reliable:\n");
// 		// }
// 		informationLocation = sampledPoints_;
// 		unscalePoints(informationLocation);
// 	}else if(strategy_.compare("lowLoc") == 0){
// 		// if(printInfo_){
// 		// 	printf("Train high fi model to decide what low fidelity data is reliable:\n");
// 		// }
// 		informationLocation = sampledPointsLow_;
// 		unscalePoints(informationLocation);
	
// 	}else if(strategy_.compare("largeLoc") == 0){
// 		// if(printInfo_){
// 		// 	printf("Train both high and low fi models to decide what low fidelity data is reliable:\n");

// 		// }
// 		informationLocation = sampleGenerator_->randomLHS(100 * biFunction_->d_);
// 	}
// 	informationValsHigh = highFiKriging_->multipleSurfaceValues(informationLocation);
// 	informationValsLow = lowFiKriging_->multipleSurfaceValues(informationLocation);
// 	scalePoints(informationLocation);
// 	// Create necessary vectors to calculate LCC
// 	int locationSize = sampledPointsLow_.size();
// 	int informationSize = informationLocation.size();
// 	vector<double> localCorrValues(informationSize, 0.0);
// 	vector<double> weights;
// 	vector<double> localHighSample;
// 	vector<double> localLowSample;
// 	weights.reserve(informationSize);
// 	localHighSample.reserve(informationSize);
// 	localLowSample.reserve(informationSize);
// 	// Define neighbourhood distance
// 	double maxDist = 0.0;
// 	// If using fuction scaling, internally every point lies in [0,1]^d, so sqrt(d) is the max dist.
// 	// If no scaling, make it relative to to actual sample space
// 	if(functionScaling_){maxDist = sqrt(biFunction_->d_);}
// 	else{
// 		for(int i = 0; i < biFunction_->d_; i++){
// 			maxDist += pow(biFunction_->upperBound_[i] - biFunction_->lowerBound_[i], 2);
// 		}
// 	}
// 	maxDist = r *sqrt(maxDist);

// 	// For each location, calculate LCC and choose whether to use the point
// 	for(int i = 0; i < locationSize; i++){
// 		VectorXd point = sampledPointsLow_[i];
// 		int addedPoints = 0;
// 		// Cycle through all points, if closer than max dist, add to local vector
// 		printPoint(point);
// 		printf(" - real ");
// 		unscalePoint(point);
// 		printPoint(point);
// 		scalePoint(point);
// 		// printf("\n");
// 		for(int j = 0; j < informationSize; j++){
// 			double localDist = dist(point, informationLocation[j]);
// 			if(localDist > maxDist){continue;}
// 			addedPoints++;
// 			localHighSample.push_back(informationValsHigh[j]);
// 			localLowSample.push_back(informationValsLow[j]);
// 			weights.push_back(1.0 - localDist / maxDist);
// 			// unscalePoint(informationLocation[j]);
// 			// printPoint(informationLocation[j]);
// 			// scalePoint(informationLocation[j]);
// 			// printf(": %.2f %.2f %.2f %.2f\n", localDist, 1.0 - localDist / maxDist, informationValsHigh[j], informationValsLow[j]);
// 		}
// 		// Calculate local correlation if have at least 2 points, otherwise base case is use point (as there is no high fi data nearby, at least can rely on this)
// 		double localCorr;
// 		if(addedPoints < 2){
// 			localCorr = 1;
// 		}else{
// 			localCorr = weightedCorrelationCoefficient(localLowSample, localHighSample, weights, false);
// 		}
// 		printf(": Used %d points, local corr %.2f\n", addedPoints, localCorr);
// 		// If this is past a certain threshold, save it
// 		if(localCorr > p){
// 			usedSampledPointsLow_.push_back(sampledPointsLow_[i]);
// 			usedSampledPointsValuesLow_.push_back(sampledPointsValuesLow_[i]);
// 		}

// 		// Done! Clear vectors for next iteration
// 		weights.clear();
// 		localHighSample.clear();
// 		localLowSample.clear();
// 	}

// 	// I think we are done! Chosen points and values are saved in usedSampledPointsLow_ and usedSampledPointsValuesLow_, respectively.
// 	onlyUseHighFiData_ = ((int)usedSampledPointsLow_.size() == 0);
// }

// void AdaptiveCoKriging::trainModel(){
// 	if(printInfo_){printf("Training adaptive cokriging model.\n");}
// 	chooseReliableLowFidelityData(pFactor_, rFactor_);
// 	// If chosen low fi data is empty, do not train the low fi model
// 	if(onlyUseHighFiData_){
// 		// No need for further training, will only use the high fi model already trained
// 		if(printInfo_){printf("Chosen to ignore all low fi data, currently working exactly like Kriging.\n");}
// 	}else{
// 		lowFiKriging_->sampledPoints_ = usedSampledPointsLow_;
// 		lowFiKriging_->sampledPointsValues_ = usedSampledPointsValuesLow_;
// 		if(printInfo_){printf("Training adaptive cokriging model, train low fi kriging with chosen data.\n");}
// 		lowFiKriging_->trainModel();
// 		if(printInfo_){printf("Training adaptive cokriging model, train intermediate model.\n");}
// 		trainHyperparameters();
// 		saveMuSigma();
// 	}
// }


// tuple<double, double> AdaptiveCoKriging::meanVarianceCalculator(VectorXd &x){
// 	if(onlyUseHighFiData_){
// 		return highFiKriging_->meanVarianceCalculator(x);
// 	}else{
// 		return CoKriging::meanVarianceCalculator(x);
// 	}
// }










#endif