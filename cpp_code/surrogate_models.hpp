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

#ifndef SURROGATE_MODELS_HPP
#define SURROGATE_MODELS_HPP

#include "libraries.hpp"
#include "aux_solvers.hpp"
#include "functions.hpp"
#include "sample_generator.hpp"




// Parent class of a bi-fidelity surrogate model. This class should not be instantiated, instead one of the children classes
// should be used. It contains the functions which are shared with all surrogate models, which scale the data used to train the model
// to lie in the unit hypercube [0,1]^d, and the objective function values to lie within [0,1]. It also has 
// virtual functions for saving the training data within an instance of this class, training a model, and querying the model
// for a surface value and uncertainty prediction.
class SurrogateModel{
public:
	// The constructor requires the specification of the bi-fidelity source being modelled i.e. a pair f_h and f_l of high- and low- fidelity
	// functions. This should be a pointer to an instance of the class BiFidelityFunction. This instance also requires an black-box solver
	// passed as auxSolver. This solver is used for all auxiliary optimisation problems such as choosing optimal hyper-parameters
	// and choosing further sample locations. Note both BiFidelityFunction and AuxSolver are not defined in this repository,
	// rather they are defined in the repository https://github.com/nandresthio/bifiEBBmethods. In order to compile this code please
	// follow the instructions for linking the two repositories given in the README. For reproducibility purposes, the seed can also be
	// specified. Leaving the seed as 0 leads to new random behaviour at each run. Setting printInfo to true leads to all of the 
	// inner workings being printed to the terminal. Finally, it is recommended to scale the data being used to train the model for improved
	// accuracy and stability. Setting functionScaling to false turns this scalling off.
	SurrogateModel(BiFidelityFunction* biFunction, AuxSolver* auxSolver, int randomSeed = 0, bool printInfo = false, bool functionScaling = true);

	virtual ~SurrogateModel();

	virtual void setCostRatio(double costRatio);

	// Function which linearly scales a given point so that is lies in the unit hypercube [0,1]^d. It does so
	// based on the domain of the given biFunction.
	void scalePoint(VectorXd &point);

	// Function which linearly scales a set of points to lie in the unit hypercube [0,1]^d; calls the function scalePoint
	// on every point in the vector.
	void scalePoints(vector<VectorXd> &points);

	// Function which linearly scales a point in the unit hypercube [0,1]^d to the domain of the domain of the given function.
	void unscalePoint(VectorXd &point);

	// Function which linearly scales a set of points from the unity hypercube [0,1]^d to the domain of the biFunction
	void unscalePoints(vector<VectorXd> &points);

	// Function which linearly scales an objective function value from [minObservation_,maxObservation_] to [0,1]. Note
	// that this should lead to all f_h values to lie within [0,1], but f_l which lie outside of the range of f_h will also
	// lie outside of [0,1]. Note also that this function should be called only after calling scaleTwoObservations at least
	// once, as this finds the values for minObservation_ and maxObservation_ first.
	double scaleObservation(double value);

	// Function which linearly scales a set of objective function values by calling the scaleObservation function on each of the 
	// observations.
	void scaleObservations(vector<double> &observations);

	// Function which linearly scales an objective function value from [0,1] to [minObservation_,maxObservation_]. Note that
	// for scaled objective function values from f_l which lie outside of [0,1], scaling them will also take outside of
	// the range [minObservation_,maxObservation_].
	double unscaleObservation(double value);

	// Function which linearly scales a set of objective function values from [0,1] to [minObservation_,maxObservation_] 
	// by calling unscaleObservation on each of the observations.
	void unscaleObservations(vector<double> &observations);

	// Function which scales two sets of observations. The set "observations" should come from f_h, and the set "observationsLow"
	// should come from f_l. The function finds the min and max values from "observations", and uses these values
	// to linearly scale all the values to [0,1]. Note that some values in observationsLow might get mapped outside of [0,1]
	void scaleTwoObservations(vector<double> &observations, vector<double> &observationsLow);

	// Function which takes a gathered sample assumed to come from two sources f_h and f_l, and internally saves them inside the 
	// surrogate model instance. The sets "points" and "pointsLow" are the locations at which f_h and f_l have been sampled, respectively.
	// The sets "observations" and "observationsLow" are the objective function values of f_h and f_l, respectively. Note that the
	// default behaviour is to scale everything internally, that is the samples to lie in [0,1]^d, and the observations to lie
	// within [0,1]. In order to stop this scaling, set functionScaling to false when initialising the surrogate model.
	virtual void saveSample(vector<VectorXd> points, vector<VectorXd> pointsLow, vector<double> observations, vector<double> observationsLow);

	// Function which adds a point to the saved sample. It evaluates the high-fidelity source and/or the low-fidelity source, scales
	// the resulting values and adds them to the stored data.
	void addSample(VectorXd point, bool sampleHigh, bool sampleLow);

	// Parent function to train model. For this class this is empty and should not be called; each of the derived classes
	// implement this based on how each of the models is trained.
	// If forcedOptimiseHyperparameters is true, will re optimise the hyperparameters. If not, 
	// will optimise them "smartly" i.e. always for small samples, and only sometimes for large samples.
	virtual void trainModel(bool forcedOptimiseHyperparameters = false);

	// Returns the surrogate surface value prediction at a point. It is assumed the point is unscaled i.e. it lies in the domain
	// of the specified biFunction, and that the desired output should not be scaled i.e. it should lie near the range of the biFunction.
	// Internally these flags might be modified to account for internal scaling. Note this function should not
	// be called at this level, and instead should be implemented in each of the children classes.
	virtual double surfaceValue(VectorXd &x, bool pointIsScaled = false, bool unscaleOutput = true);

	// Function which returns a vector with the surrogate surgace value precitions for a set of points. For further details
	// look at the surfaceValue function above. 
	vector<double> multipleSurfaceValues(vector<VectorXd> &points, bool pointIsScaled = false, bool unscaleOutput = true);

	// Returns the evaluation of the acquisition function at a particular location. It is assumed the point is unscaled i.e. it lies in the domain
	// of the specified biFunction. Internally this flag might be modified to account for internal scaling. Note this function should not
	// be called at this level, and instead should be implemented in each of the children classes.
	virtual double evaluateAcquisitionFunction(VectorXd x);

	// Function which defines which acquisiton function should be used. To see what acquisition functions are
	// defined consult the evaluateAcquisitionFunction definition in derived classes.
	virtual void setAcquisitionFunction(string chosenAcquisiton);

	// Finds the location of the next sampling site, based on an already specified acquisition function
	// Returns the location and two booleans specifying whether the high / low sources should be sampled
	virtual tuple<VectorXd, double, bool, bool> findNextSampleSite();

	class AcquisitionFunction : public Function{
		public:
	
		AcquisitionFunction(SurrogateModel* surrogateModel);

		~AcquisitionFunction();

		virtual double evaluate(VectorXd &point) override;

		SurrogateModel* surrogateModel_; 	// Reference to surrogate model needed to evaluate the acquisition function
	};


	BiFidelityFunction* biFunction_;		// Two-source black box function for which a surrogate model is built.
	AuxSolver* auxSolver_;					// Auxiliary solver used in auxiliary optimisation problems i.e. hyperparameter optimisation.
	int randomSeed_;						// Random seed stored to allow for reproducible results.
	bool printInfo_;						// Indicates whether information should be printed as the model is being trained.
	mt19937 randomGenerator_;				// Random number generator used for various purposes.
	SampleGenerator* sampleGenerator_;		// Sample generator used when creating sample to train surrogate model.
	bool functionScaling_;					// Internal flag which determines whether the sample locations and objective function values are scaled
	
	vector<VectorXd> sampledPoints_;		// Points at which the function f_h was sampled, which will be used to train the model.
	vector<double> sampledPointsValues_;	// Function value of f_h at the sampled points.
	vector<VectorXd> sampledPointsLow_;		// Points at which the function f_l was sampled, which will be used to train the model.
	vector<double> sampledPointsValuesLow_;	// Function value of f_l at the sampled points.

	double maxObservation_;					// Highest f_h value seen so far. Used for internal scaling.
	double minObservation_;					// Lowest f_h value seen so far. Used for internal scaling.

	bool trainedModel_;						// bool which indicates whether model has been trained, used in order to assess whether model values can be calculated.

	string chosenAcquisiton_;				// String which defines what the acquisition function should be.
	bool acquisitionIsMin_;					// Bool which defines whether the acquisition function is being minimised or maximised.

	double costRatio_;						// The cost ratio of sampling the low-fidelity source; should be between 0 and 1.

	// Can't reoptimise hyperparameters at every iteration for large data samples as it takes forever.
	// Have a couple of checks to only reoptimise them every so often in this case
	int startIntervalTraining_;
	int intervalForTraining_;
	int startLargeIntervalTraining_;
	int intervalForLargeTraining_;
	int startVeryLargeIntervalTraining_;
	int intervalForVeryLargeTraining_;
	int sampleSizeAtLastTrain_;
};





// Implementation of a Kriging surrogate model of a single source function, which is normally an expensive black box.
// Once an instance of this class has been trained via a call to trainModel, it can be queried for 
// surface value and variance (i.e. uncertainty) in the sample space.
// This implementation follows the details given in the appendix folder.
class Kriging : public SurrogateModel{
	public:

	// For information on the inputs of the constructor, look at the parent constructor
	Kriging(BiFidelityFunction* biFunction, AuxSolver* auxSolver, int randomSeed = 0, bool printInfo = false, bool functionScaling = true);
	
	virtual ~Kriging();

	// Main method which trains the Kriging model. Undergoes an auxiliary optimisation to find the best hyperparameters,
	// then stores some useful values (i.e. mu and sigma) and matrices.
	virtual void trainModel(bool forcedOptimiseHyperparameters = false) override;

	// Method which uses ARS to find the best hyperparameter values by maximising the concentrated likelihood function.
	virtual void trainHyperparameters();

	// Method which saves mu, sigma and two useful matrices for further calculations once the optimal hyperparameters have been found.
	virtual void saveMuSigma();

	// Calculates mu, sigma, R and a decomposition of R using the storeds hyperparameters. Used when calculating 
	// the concentrated likelihood function, as well as when saving these values after the hyperparameters have been found.
	tuple<double, double, MatrixXd, LDLT<MatrixXd>, MatrixXd> muSigmaCalculator();

	// Returns the surrogate surface value at a point
	double surfaceValue(VectorXd &x, bool pointIsScaled = false, bool unscaleOutput = true) override;

	// Returns the surrogate surface variance, or uncertainty, at a point
	double uncertainty(VectorXd &x, bool pointIsScaled = false, bool unscaleOutput = true);

	// Returns the expected improvement at a point; assumption is we are minimising
	double expectedImprovement(VectorXd &x, bool pointIsScaled = false, bool unscaleOutput = true);

	double kernelProductIntegral(double x1, double x2, double theta1, double theta2);

	virtual void setWmatrix();

	double globalVarianceReduction(VectorXd &x, bool pointIsScaled = false, bool unscaleOutput = true);

	double globalVarianceReductionApproximation(VectorXd &x, bool pointIsScaled = false, bool unscaleOutput = true);

	// Calculates the surrogate surface value and variance at a particular point
	virtual tuple<double, double> meanVarianceCalculator(VectorXd &x);

	// Returns the evaluation of the acquisition function at a particular location. It is assumed the point is unscaled i.e. it lies in the domain
	// of the specified biFunction. Internally this flag might be modified to account for internal scaling. Note this function should not
	// be called at this level, and instead should be implemented in each of the children classes.
	virtual double evaluateAcquisitionFunction(VectorXd x) override;

	// Function which defines which acquisiton function should be used. To see what acquisition functions are
	// defined consult the evaluateAcquisitionFunction definition in derived classes.
	virtual void setAcquisitionFunction(string chosenAcquisiton) override;

	// Finds the location of the next sampling site, based on an already specified acquisition function
	virtual tuple<VectorXd, double, bool, bool> findNextSampleSite() override;

	// Calculates the concentrated likelihood function using the stored hyperparameters
	double concentratedLikelihoodFunction();

	// Function defined to be used for the auxiliary optimisation problem of optimising the hyperparameters.
	// A point passed to the evaluate function is a set of hyperparameters. This point is evaluated by storing 
	// these hyperparameters in the Kriging model and calculating the concentrated likelihood function. 
	class ConcentratedLikelihoodFunction : public Function{
		public:
	
		ConcentratedLikelihoodFunction(int d, vector<double> &lowerBound, vector<double> &upperBound, Kriging* krigingModel);

		~ConcentratedLikelihoodFunction();

		virtual double evaluate(VectorXd &point) override;

		virtual int betterPoint(VectorXd &point1, double val1, VectorXd &point2, double val2, string flag = "") override; 

		Kriging* krigingModel_; 			// Reference to Kriging model needed to calculate loglikelihood for a given set of hyperparameters
	};


	vector<double> theta_;					// First set of d hyperparameters, where d is the problem dimension.
	vector<double> pVector_;				// Second set of d hyperparameters, where d is the problem dimension.

	double mu_;								// mu value of the surrogate model, stored for speed when calculating surface values.
	double sigma_;							// sigma value of the surrogate model, stored for speed when calculating surface values.
	MatrixXd rMatrix_;						// R matrix, stored for speed
	LDLT<MatrixXd> rMatrixDecomposition_;	// R matrix decomposition, stored for speed instead of the inverse as this is faster.
	MatrixXd modelWeights_;					// To save time when needing predictions, save the multiplication of the system solved when training the model.
	MatrixXd rMatrixInverse_;				// Again, to save time this matrix is stored once the hyperparameters have been chosen. It is calculated once and it is faster than using solve with the decomposition (if the inverse is only calculated once).

	double maxLogLikelihood_;				// Maximum loglikelihood observed so far, stored for tie breaker purposes for certain training approaches.

	vector<VectorXd> varianceTestPoints_;	// Locations used to estimate the integral of the variance, used for acquisition function based on overall variance reduction

	bool fixedP_;							// Bool which indicates whether the hyperparameters p_i are all fixed to 2.
	MatrixXd wMatrix_;						// Matrix for variance integration
};








// Implementation of a CoKriging surrogate model of a two-source function, which is normally a high fidelity expensive black box function,
// and a cheaper, low fidelity but also expensive black box function. Once an instance of this class has been trained via a call to 
// trainModel, it can be queried for surface value and variance (i.e. uncertainty) in the sample space.
// This implementation follows the details given in the appendix folder.
class CoKriging: public Kriging{
	public:
	// For details on this constructed, look at the constructor of the parent class SurrogateModel
	CoKriging(BiFidelityFunction* biFunction, AuxSolver* auxSolver, int randomSeed = 0, bool printInfo = false, bool functionScaling = true);

	~CoKriging() override;

	// Override function which also saves data to the low fi krig model
	void saveSample(vector<VectorXd> points, vector<VectorXd> pointsLow, vector<double> observations, vector<double> observationsLow) override;

	// Main method which trains the CoKriging model. First trains a Kriging model on the low fidelity data,
	// then undergoes an auxiliary optimisation to find the best hyperparameters of the difference model including rho.
	// Finally it stores some useful values (i.e. mu and sigma) and matrices.
	virtual void trainModel(bool forcedOptimiseHyperparameters = false) override;

	// Method which finds the optimal hyperparameters for the difference model (i.e. theta_b's and p_b's) as well as rho
	// using the ARS by maximising the intermediate concentrated likelihood function.
	void trainHyperparameters() override;

	// Saves important mu, sigma and matrix values for speed once hyperparamters have been calculated, inculding sigma_b, sigma_l, mu,
	// C matrix and its decomposition.
	void saveMuSigma() override;

	// Calculates mu_b, sigma_b, R_b and a decomposition of R_b using the storeds hyperparameters. Used when calculating 
	// the intermediate concentrated likelihood function, as well as when saving sigma_b after the hyperparameters have been found.
	tuple<double, double, MatrixXd, LDLT<MatrixXd>> intermediateMuSigmaCalculator();

	// Calculates the combined mu value as well as the C matrix and its decomposition. Used to store these values for speed.
	// when calculating surrogate model values.
	tuple<double, MatrixXd, LDLT<MatrixXd>> combinedMuCmatrixCalculator();

	// Calculates the surrogate surface value and variance at a particular point
	virtual tuple<double, double> meanVarianceCalculator(VectorXd &x) override;

	// Function which populates the vector c(x) as defined in the appendices for x = point
	VectorXd cVector(VectorXd &point);

	// Function which populates the vector c_l(x) as defined in the appendices for x = point
	VectorXd cLowVector(VectorXd &point);

	// Function which populates the vector c_b(x) as defined in the appendices for x = point
	MatrixXd cBothMatrix(VectorXd &x);

	// Function which returns the value W^1(x_1, x_2) as defined in Chapter 5.1.3.2 given in the
	// repo for x_1 = point1, x_2 = point2
	double w1product(VectorXd &point1, VectorXd &point2);

	// Function which returns the value W^2(x_1, x_2) as defined in Chapter 5.1.3.2 given in the
	// repo for x_1 = point1, x_2 = point2
	double w2product(VectorXd &point1, VectorXd &point2);

	// Function which returns the value W^3(x_1, x_2) as defined in Chapter 5.1.3.2 given in the
	// repo for x_1 = point1, x_2 = point2
	double w3product(VectorXd &point1, VectorXd &point2);

	// Function which populates the matrix W^cokrig as defined in Chapter 5.1.3.2 given the saved
	// sample
	virtual void setWmatrix() override;

	// Function which evaluates Delta^l_{x} IMSPE as defined in Chapter 5.1.3.2
	double globalVarianceReductionLowFiSample(VectorXd &x, bool pointIsScaled = false, bool unscaleOutput = true);

	// Function which evaluates Delta^h_{x} IMSPE as defined in Chapter 5.1.3.2
	double globalVarianceReductionHighFiSample(VectorXd &x, bool pointIsScaled = false, bool unscaleOutput = true);

	// Function which evaluates Delta^b_{x} IMSPE as defined in Chapter 5.1.3.2
	double globalVarianceReductionBothFiSample(VectorXd &x, bool pointIsScaled = false, bool unscaleOutput = true);

	// double globalVarianceReductionApproximationLowFiSample(VectorXd &x, bool pointIsScaled = false, bool unscaleOutput = true);

	// double globalVarianceReductionApproximationHighFiSample(VectorXd &x, bool pointIsScaled = false, bool unscaleOutput = true);

	// double globalVarianceReductionApproximationBothFiSample(VectorXd &x, bool pointIsScaled = false, bool unscaleOutput = true);

	// Returns the evaluation of the acquisition function at a particular location. It is assumed the point is unscaled i.e. it lies in the domain
	// of the specified biFunction. Internally this flag might be modified to account for internal scaling.
	virtual double evaluateAcquisitionFunction(VectorXd x) override;

	// Function which saves a string indicating to the model which acquisition function should be used.
	virtual void setAcquisitionFunction(string chosenAcquisiton) override;

	// Finds the location of the next sampling site, based on an already specified acquisition function.
	virtual tuple<VectorXd, double, bool, bool> findNextSampleSite() override;

	// Calculates the log likelihood function of the intermediate (i.e. error) model using the stored hyperparameters.
	double intermediateConcentratedLikelihoodFunction();
	
	// Function defined to be used for the auxiliary optimisation problem of optimising the hyperparameters of the intermediate (i.e. difference).
	// A point passed to the evaluate function is a set of hyperparameters. This point is evaluated by storing 
	// these hyperparameters in the CoKriging model and calculating the concentrated likelihood function. 
	class IntermediateConcentratedLikelihoodFunction : public Function{
	public:

		IntermediateConcentratedLikelihoodFunction(int d, vector<double> &lowerBound, vector<double> &upperBound, CoKriging* cokrigingModel);

		~IntermediateConcentratedLikelihoodFunction();

		virtual double evaluate(VectorXd &point) override;

		virtual int betterPoint(VectorXd &point1, double val1, VectorXd &point2, double val2, string flag = "") override; 

		CoKriging* cokrigingModel_;
	};

	Kriging* lowFiKriging_;					// Kriging model trained on the low fidelity data. 
	vector<double> thetaB_;					// Intermediate set of first hyperparameters.
	vector<double> pBVector_;				// Intermediate set of second hyperparameters.
	double rho_;							// Multiplier of low fidelity model, also a hyperparameter.
	vector<double> lowFiSurfaceValues_;		// Surface value prediction of low fidelity kriging at high fidelity locations, stored for speed.

	double sigmaL_;							// Sigma of the low fidelity model.
	double sigmaB_;							// Sigma of the intermediate model.
	double muCombined_;						// Mu of the combined model.
	MatrixXd cMatrix_;						// C matrix, stored for speed of future calculations.
	LDLT<MatrixXd> cMatrixDecomposition_;	// C matrix decomposition, stored for speed of future calculations.
	MatrixXd cMatrixInverse_;				// Again, to save time this matrix is stored once the hyperparameters have been chosen. It is calculated once and it is faster than using solve with the decomposition (if the inverse is only calculated once).

	MatrixXd pMatrix_;
	MatrixXd pMatrixTranspose_;
	MatrixXd lMatrixInverse_;
	MatrixXd lMatrixStarInverse_;
	MatrixXd dMatrixInverse_;

};




// Adaptive surrogate model which, when being trained, chooses whether the low-fidelity source should be used.
// If not, a Kriging model is used (trained only on high-fidelity data).
// If yet, a Co-Kriging model is used (trained on low and high fidelity data).
// It is a derived class of Co-Kriging since you can call exptected improvement, uncertainty, etc...
class AdaptiveCoKriging: public CoKriging{
	public:
	// For details on this constructed, look at the constructor of the parent class SurrogateModel
	AdaptiveCoKriging(BiFidelityFunction* biFunction, AuxSolver* auxSolver, bool advanced, int randomSeed = 0, bool printInfo = false, bool functionScaling = true);

	// Additional constructed where the kriging and cokriging models are isntantiated separately and passed to the constructor
	AdaptiveCoKriging(BiFidelityFunction* biFunction, AuxSolver* auxSolver, Kriging* krigingModel, CoKriging* cokrigingModel, bool advanced, int randomSeed = 0, bool printInfo = false, bool functionScaling = true);


	~AdaptiveCoKriging() override;

	// Stores the cost ratio in the surrogate model; used for comparison of adaptive function.
	void setCostRatio(double costRatio) override;

	// Function which assesses the low fi function i.e. checks its quality to decide whether the 
	// model should rely on it. The rules followed are those derived in the provided thesis in the repo
	// in Chapter 5.1.2
	void assessLowFiFunction();

	// Override function which also saves data to the low fi krig model
	void saveSample(vector<VectorXd> points, vector<VectorXd> pointsLow, vector<double> observations, vector<double> observationsLow) override;

	// Main method which trains the CoKriging model. First trains a Kriging model on the low fidelity data,
	// then undergoes an auxiliary optimisation to find the best hyperparameters of the difference model including rho.
	// Finally it stores some useful values (i.e. mu and sigma) and matrices.
	virtual void trainModel(bool forcedOptimiseHyperparameters = false) override;

	// Calculates the surrogate surface value and variance at a particular point
	virtual tuple<double, double> meanVarianceCalculator(VectorXd &x) override;

	// Function which saves a string indicating to the model which acquisition function should be used.
	virtual void setAcquisitionFunction(string chosenAcquisiton) override;

	// Finds the location of the next sampling site, based on an already specified acquisition function
	virtual tuple<VectorXd, double, bool, bool> findNextSampleSite() override;

	// Function used to train a simple regression model by least squares regression.
	class leastSquaresFunction : public Function{
	public:

		leastSquaresFunction(int d, vector<VectorXd> points, vector<double> values);

		~leastSquaresFunction();

		virtual double evaluate(VectorXd &point) override;

		vector<VectorXd> points_;
		vector<double> values_;
	};


	Kriging* krigingModel_;					// Kriging model trained on high fidelity data
	CoKriging* cokrigingModel_;				// Kriging model trained on high and low fidelity data
	bool lowFiDataIsUseful_;				// Bool which stated whether the low fi data should be used

	bool advanced_;							// Aditional flag which only does switch usefulness
											// based on sample size

};



#endif