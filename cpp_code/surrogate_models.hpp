#ifndef SURROGATE_MODELS_HPP
#define SURROGATE_MODELS_HPP

#include "libraries.hpp"
#include "aux_solvers.hpp"
#include "functions.hpp"
#include "sample_generator.hpp"



// Implementation of a Kriging surrogate model of a single source function, which is normally an expensive black box.
// Includes methods to place down a sample using Latin Hypercube Sampling and the training
// of its hyperparameters. Uses Accelerated Random Search (ARS) for all auxiliary optimisation problems. 
// Once an instance of this class has been trained via a call to createSurrogateModel, it can be queried for 
// surface value and variance (i.e. uncertainty) in the sample space.
// This implementation follows the details given in the Appendix of 
// "Bi-fidelity Surrogate Modelling: Showcasing the need for new test instances" by Nicolau Andres-Thio, Mario Andres Munoz and Kate Smith-Miles.
class Kriging{
	public:

	Kriging(Function* ebbFunction, AuxSolver* auxSolver, int sampleBudget, int randomSeed = 0, bool printInfo = false, bool functionScaling = true);
	
	virtual ~Kriging();


	void scalePoint(VectorXd &point);

	void scalePoints(vector<VectorXd> &points);

	void unscalePoint(VectorXd &point);

	void unscalePoints(vector<VectorXd> &points);

	double scaleObservation(double value);

	void scaleObservations(vector<double> &observations);

	double unscaleObservation(double value);

	void unscaleObservations(vector<double> &observations);

	void saveSample(vector<VectorXd> sample, vector<double> sampleVals);

	virtual void sampleAtLocation(VectorXd point);

	// Main method which queries the single source function and then trains the Kriging model.
	void createSurrogateModel();

	// Method which uses the sample generator to query the source function, using a Latin Hypercube Sampling plan
	virtual void generateSample();

	// Main method which trains the Kriging model. Undergoes an auxiliary optimisation to find the best hyperparameters,
	// then stores some useful values (i.e. mu and sigma) and matrices.
	virtual void trainModel();

	// Method which uses ARS to find the best hyperparameter values by maximising the concentrated likelihood function.
	virtual void trainHyperparameters();

	// Method which saves mu, sigma and two useful matrices for further calculations once the optimal hyperparameters have been found.
	virtual void saveMuSigma();

	// Calculates mu, sigma, R and a decomposition of R using the storeds hyperparameters. Used when calculating 
	// the concentrated likelihood function, as well as when saving these values after the hyperparameters have been found.
	tuple<double, double, MatrixXd, LDLT<MatrixXd>, MatrixXd> muSigmaCalculator();

	// Returns the surrogate surface value at a point
	double surfaceValue(VectorXd &x, bool pointIsScaled = false, bool scaleOutput = true);

	// Returns the surrogate surface variance, or uncertainty, at a point
	double uncertainty(VectorXd &x, bool pointIsScaled = false, bool scaleOutput = true);

	// Returns the expected improvement at a point; assumption is we are minimising
	double expectedImprovement(VectorXd &x, bool pointIsScaled = false, bool scaleOutput = true);

	// Returns the surrogate surface value at multiple points
	vector<double> multipleSurfaceValues(vector<VectorXd> &points, bool pointIsScaled = false, bool scaleOutput = true);

	// Calculates the surrogate surface value and variance at a particular point
	virtual tuple<double, double> meanVarianceCalculator(VectorXd &x);

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

		Kriging* krigingModel_;
	};


	VectorXd chooseNextSampleSite(string technique);


	class nextSampleSiteFunction : public Function{
		public:
	
		nextSampleSiteFunction(int d, vector<double> &lowerBound, vector<double> &upperBound, Kriging* model, string technique);

		~nextSampleSiteFunction();

		virtual double evaluate(VectorXd &point) override;

		Kriging* model_;
		string technique_;

	};




	Function* ebbFunction_;					// Function for which a surrogate model will be built.
	AuxSolver* auxSolver_;					// Auxiliary solver used in auxiliary optimisation problems i.e. hyperparameter optimisation.
	int sampleBudget_;						// Sample size to be used to train the surrogate model.
	int randomSeed_;						// Random seed stored to allow for reproducible results.
	bool printInfo_;						// Indicates whether information should be printed as the model is being trained.
	mt19937 randomGenerator_;				// Random number generator used for various purposes.
	SampleGenerator* sampleGenerator_;		// Sample generator used when creating sample to train surrogate model.
	vector<double> theta_;					// First set of d hyperparameters, where d is the problem dimension.
	vector<double> pVector_;				// Second set of d hyperparameters, where d is the problem dimension.

	vector<VectorXd> sampledPoints_;		// Points at which the function was sampled, which will be used to train the model.
	vector<double> sampledPointsValues_;	// Function value at the sampled points.
	double mu_;								// mu value of the surrogate model, stored for speed when calculating surface values.
	double sigma_;							// sigma value of the surrogate model, stored for speed when calculating surface values.
	MatrixXd rMatrix_;						// R matrix, stored for speed
	LDLT<MatrixXd> rMatrixDecomposition_;	// R matrix decomposition, stored for speed instead of the inverse as this is faster.
	MatrixXd modelWeights_;					// To save time when needing predictions, save the multiplication of the system solved when training the model.

	bool trainedModel_;						// bool which indicates whether model has been trained, used in order to assess whether model values can be calculated.
	bool functionScaling_;

	double maxObservation_;
	double minObservation_;

	double maxLogLikelihood_;

	double fMin_;
};











// Implementation of a CoKriging surrogate model of a two-source function, which is normally a high fidelity expensive black box function,
// and a cheaper, low fidelity but also expensive black box function.
// Includes methods to place down a sample of both the high and low fidelity functions using Latin Hypercube Sampling and 
// the Morris-Mitchell criterion, followed by the training of its hyperparameters. Uses Accelerated Random Search (ARS) 
// for all auxiliary optimisation problems. Once an instance of this class has been trained via a call to createSurrogateModel, it can be queried for 
// surface value and variance (i.e. uncertainty) in the sample space.
// This implementation follows the details given in the Appendix of 
// "Bi-fidelity Surrogate Modelling: Showcasing the need for new test instances" by Nicolau Andres-Thio, Mario Andres Munoz and Kate Smith-Miles.
class CoKriging: public Kriging{
	public:
	CoKriging(BiFidelityFunction* biFunction, AuxSolver* auxSolver, int highFiSampleBudget, int lowFiSampleBudget, int randomSeed = 0, bool printInfo = false, bool functionScaling = true);

	~CoKriging() override;

	void scaleTwoObservations(vector<double> &observations, vector<double> &observationsLow);

	void unscaleTwoObservations(vector<double> &observations, vector<double> &observationsLow);

	virtual void saveSample(vector<VectorXd> sample, vector<double> sampleVals, vector<VectorXd> sampleLow, vector<double> sampleValsLow);

	virtual void sampleAtLocation(VectorXd point) override;

	// Method which uses the sample generator to query the source function. First puts down a Latin Hypercube Sampling plan,
	// followed by choosing a subset which is locally maximal in terms of the Morris-Mitchell criterion.
	void generateSample() override;

	// Main method which trains the CoKriging model. First trains a Kriging model on the low fidelity data,
	// then undergoes an auxiliary optimisation to find the best hyperparameters of the difference model including rho.
	// Finally it stores some useful values (i.e. mu and sigma) and matrices.
	virtual void trainModel() override;

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

	BiFidelityFunction* biFunction_;		// Two-source black box function for which a surrogate model is built.
	int highFiSampleBudget_;				// Number of samples of the high fidelity source used when training the model.
	int lowFiSampleBudget_;					// Number of samples of the low fidelity source used when training the model.
	Kriging* lowFiKriging_;					// Kriging model trained on the low fidelity data. 
	vector<double> thetaB_;					// Intermediate set of first hyperparameters.
	vector<double> pBVector_;				// Intermediate set of second hyperparameters.
	double rho_;							// Multiplier of low fidelity model, also a hyperparameter.
	vector<double> lowFiSurfaceValues_;		// Surface values of low fidelity kriging at high fidelity locations, stored for speed.

	double sigmaL_;							// Sigma of the low fidelity model.
	double sigmaB_;							// Sigma of the intermediate model.
	double muCombined_;						// Mu of the combined model.
	MatrixXd cMatrix_;						// C matrix, stored for speed of future calculations.
	LDLT<MatrixXd> cMatrixDecomposition_;	// C matrix decomposition, stored for speed of future calculations.

};





// Implementation of CoKriging when the low fidelity function is "cheap", i.e. there is no cost to sample it and is therefore sampled freely.
class CoKrigingCheap: public Kriging{
	public:

	CoKrigingCheap(BiFidelityFunction* biFunction, AuxSolver* auxSolver, int highFiSampleBudget, char mode = 's', int randomSeed = 0, bool printInfo = false, bool functionScaling = true);

	~CoKrigingCheap() override;

	void trainHyperparameters() override;

	void saveMuSigma() override;

	tuple<double, double, MatrixXd, LDLT<MatrixXd>> intermediateMuSigmaCalculator();

	double muCalculator(int nL = 10000);

	tuple<double, double> meanVarianceCalculator(VectorXd &x) override;

	double intermediateConcentratedLikelihoodFunction();
	
	class IntermediateConcentratedLikelihoodFunction : public Function{
	public:

		IntermediateConcentratedLikelihoodFunction(int d, vector<double> &lowerBound, vector<double> &upperBound, CoKrigingCheap* cokrigingCheapModel);

		~IntermediateConcentratedLikelihoodFunction();

		virtual double evaluate(VectorXd &point) override;

		virtual int betterPoint(VectorXd &point1, double val1, VectorXd &point2, double val2, string flag = "") override; 

		CoKrigingCheap* cokrigingCheapModel_;
	};

	BiFidelityFunction* biFunction_;		// Two-source black box function for which a surrogate model is built.
	int highFiSampleBudget_;				// Number of samples of the high fidelity source used when training the model.
	vector<double> thetaB_;					// Intermediate set of first hyperparameters.
	vector<double> pBVector_;				// Intermediate set of second hyperparameters.
	double rho_;							// Multiplier of low fidelity model, also a hyperparameter.
	
	double sigmaB_;							// Sigma of the intermediate model.
	double muIntermediate_;					// Mu of the intermediate model.
	double mu_;
	MatrixXd bMatrix_;						
	LDLT<MatrixXd> bMatrixDecomposition_;	

	char mode_;

};



class AdaptiveCoKriging: public CoKriging{
	public:
	AdaptiveCoKriging(BiFidelityFunction* biFunction, string strategy, double pFactor, double rFactor, AuxSolver* auxSolver, int highFiSampleBudget, int lowFiSampleBudget, int randomSeed = 0, bool printInfo = false, bool functionScaling = true);

	~AdaptiveCoKriging() override;

	void saveSample(vector<VectorXd> sample, vector<double> sampleVals, vector<VectorXd> sampleLow, vector<double> sampleValsLow) override;

	virtual void sampleAtLocation(VectorXd point) override;

	void chooseReliableLowFidelityData(double p = 0.5, double r = 0.2);

	void trainModel() override;

	tuple<double, double> meanVarianceCalculator(VectorXd &x) override;

		
	string strategy_;
	double pFactor_;
	double rFactor_;

	Kriging* highFiKriging_;

	vector<VectorXd> sampledPointsLow_;				// Points at which the function was sampled, which will be used to train the model.
	vector<double> sampledPointsValuesLow_;			// Function value at the sampled points.

	vector<VectorXd> usedSampledPointsLow_;			// Points which have been chosen as "trustworthy" by the selection process to train a low kriging model.
	vector<double> usedSampledPointsValuesLow_; 	// Low fi values at the chosen points.

	bool onlyUseHighFiData_;						// Bool on whether all low fi data has been discarded.
};


#endif