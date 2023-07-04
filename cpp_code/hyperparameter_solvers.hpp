#ifndef HYPERPARAMETER_SOLVERS_HPP
#define HYPERPARAMETER_SOLVERS_HPP

#include "libraries.hpp"
#include "aux_solvers.hpp"
#include "aux_functions.hpp"
#include "functions.hpp"
#include "sample_generator.hpp"
// extern "C" {
// 	#include "cmaes/src/boundary_transformation.h"
// 	#include "cmaes/src/cmaes.h"
// 	#include "cmaes/src/cmaes_interface.h"
// }





class AmoebaSolver : public AuxSolver{
	public:

	AmoebaSolver(Function* function, bool min = true, int randomSeed = 0, bool printInfo = false);

	AmoebaSolver(Function* function, int maxEval, bool min = true, int randomSeed = 0, bool printInfo = false);

	virtual ~AmoebaSolver() override;

	// Main method which finds best possible point given budget; returns best point found.
	virtual VectorXd optimise() override;

	// Extra function so that the initial simplex and values are given, with heavy tune in mind
	VectorXd optimise(vector<VectorXd> simplexPoints, vector<double> simplexValues);

	int maxEval_;							// Number of evaluations the search is allowed
	int evalsUsed_;
	
};

class GASolver : public AuxSolver{
	public:

	GASolver(Function* function, bool min = true, int randomSeed = 0, bool printInfo = false);

	GASolver(Function* function, int maxEval, int populationSize, bool min = true, int randomSeed = 0, bool printInfo = false);

	virtual ~GASolver() override;

	// Main method which finds best possible point given budget; returns best point found.
	virtual VectorXd optimise() override;

	// Extra function so that the initial population and values are given, with heavy tune in mind
	VectorXd optimise(vector<VectorXd> population, vector<double> populationVals);

	int maxEval_;							// Number of evaluations the search is allowed
	int populationSize_;					// Population at each iteration
	
};


class DynamicHillClimberSolver : public AuxSolver{
	public:

	DynamicHillClimberSolver(Function* function, bool min = true, int randomSeed = 0, bool printInfo = false);

	DynamicHillClimberSolver(Function* function, int maxEval, bool min = true, int randomSeed = 0, bool printInfo = false);

	virtual ~DynamicHillClimberSolver() override;

	// Main method which finds best possible point given budget; returns best point found.
	virtual VectorXd optimise() override;

	// Extra function so that the initial point is given, with heavy tune in mind
	VectorXd optimise(VectorXd startingPoint, double startingPointVal);

	vector<VectorXd> initialiseDirections();

	int maxEval_;							// Number of evaluations the search is allowed
	
};


class HeavyTuneSolver : public AuxSolver{
	public:

	HeavyTuneSolver(Function* function, bool min = true, int randomSeed = 0, bool printInfo = false);

	// HeavyTuneSolver(Function* function, int maxIterSimplex, int maxIterGA, int popSizeGA, int maxIterDHC, bool min = true, int randomSeed = 0, bool printInfo = false, bool useTieBreaker = false, string tieBreakerFlag = "", int finishWithLocalSearch = 0, int additionalEvaluationsInLocalSearch = 1000);
	HeavyTuneSolver(Function* function, int maxEval, int popSizeGA, bool min = true, int randomSeed = 0, bool printInfo = false);

	virtual ~HeavyTuneSolver() override;

	// Main method which finds best possible point given budget; returns best point found.
	virtual VectorXd optimise() override;

	int maxEval_;							// Number of evaluations the search is allowed

	// int maxIterSimplex_;					// Max number of iterations for the simplex solver
	// int maxIterGA_;							// Max number of iterations for the GA solver
	int popSizeGA_;							// Population size of GA solver
	// int maxIterDHC_;						// Max number of iterations for the DHC solver

};


class RandomHeavyTuneSolver : public AuxSolver{
	public:

	RandomHeavyTuneSolver(Function* function, bool min = true, int randomSeed = 0, bool printInfo = false);

	RandomHeavyTuneSolver(Function* function, int maxEval, int numSphereARS, bool min = true, int randomSeed = 0, bool printInfo = false);

	virtual ~RandomHeavyTuneSolver() override;

	// Main method which finds best possible point given budget; returns best point found.
	virtual VectorXd optimise() override;

	int maxEval_;							// Number of evaluations the search is allowed
	int numSphereARS_;						// Number of spheres when running ARS solver
	// int maxIterARS_;						// Max number of iterations for the ARS solver
	// int maxIterDHC_;						// Max number of iterations for the DHC solver
};


// class CMAESsolver : public AuxSolver{
// 	public:

// 	CMAESsolver(Function* function, bool min = true, int randomSeed = 0, bool printInfo = false, bool useTieBreaker = false, string tieBreakerFlag = "", int finishWithLocalSearch = 0, int additionalEvaluationsInLocalSearch = 1000);

// 	CMAESsolver(Function* function, int maxEval, bool min = true, int randomSeed = 0, bool printInfo = false, bool useTieBreaker = false, string tieBreakerFlag = "", int finishWithLocalSearch = 0, int additionalEvaluationsInLocalSearch = 1000);

// 	virtual ~CMAESsolver() override;

// 	// Main method which finds best possible point given budget; returns best point found.
// 	virtual VectorXd optimise() override;

// 	int maxEval_;							// Max number of iterations
	
// };



#endif