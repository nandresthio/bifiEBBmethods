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

// This file contains a variety of black-box solvers, used mainly for solving 
// the auxiliary optimisation problems of finding good hyper-parameters and 
// choosing the next sample point. They are child classes of the AuxSolver
// defined in the BiFiEBBbenchmarks repo.


#ifndef HYPERPARAMETER_SOLVERS_HPP
#define HYPERPARAMETER_SOLVERS_HPP

#include "libraries.hpp"
#include "aux_solvers.hpp"
#include "aux_functions.hpp"
#include "functions.hpp"
#include "sample_generator.hpp"


// Aomeba solver implementation (also known as Nelder–Mead method, see https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method)
// A local optimiser which stops after converging to a local minimum.
class AmoebaSolver : public AuxSolver{
	public:

	// Constructor which inputs the function being optimised, whether it is a minimisation or maximisation problem,
	// a random seed for reproduction purposes (with random behaviour when randomSeed = 0), and the boolean printInfo
	// which prints information as the optimisation proceeds. Using this constructor assumes this algorithm should 
	// have a maximum of 50*d samples, where d is the problem dimension.
	AmoebaSolver(Function* function, bool min = true, int randomSeed = 0, bool printInfo = false);

	// Constructor with the additional input maxEval, which specifies the maximum number of function samples.
	AmoebaSolver(Function* function, int maxEval, bool min = true, int randomSeed = 0, bool printInfo = false);

	virtual ~AmoebaSolver() override;

	// Main method which finds best possible point given budget; returns best point found.
	virtual VectorXd optimise() override;

	// Extra function so that the initial simplex and values are given, with heavy tune in mind
	VectorXd optimise(vector<VectorXd> simplexPoints, vector<double> simplexValues);

	int evalsUsed_;							// Total evaluations used so far.
	
};

// Genetic Algorithm (GA) solver implementation (see https://en.wikipedia.org/wiki/Genetic_algorithm).
class GASolver : public AuxSolver{
	public:

	// Constructor which inputs the function being optimised, whether it is a minimisation or maximisation problem,
	// a random seed for reproduction purposes (with random behaviour when randomSeed = 0), and the boolean printInfo
	// which prints information as the optimisation proceeds. Using this constructor assumes this algorithm should 
	// have a maximum of 10000 samples with a population of size 50 at each iteration.
	GASolver(Function* function, bool min = true, int randomSeed = 0, bool printInfo = false);

	// Constructor with the additional inputs maxEval and populationSize, which specifies the maximum number of function samples
	// and the population size at each iteration, respectively.
	GASolver(Function* function, int maxEval, int populationSize, bool min = true, int randomSeed = 0, bool printInfo = false);

	virtual ~GASolver() override;

	// Main method which finds best possible point given budget; returns best point found.
	virtual VectorXd optimise() override;

	// Extra function so that the initial population and values are given, with heavy tune in mind
	VectorXd optimise(vector<VectorXd> population, vector<double> populationVals);

	int populationSize_;					// Population at each iteration
	
};

// Dynamic Hill Climber solver implementation (see "Dynamic hill climbing" by Michael De La Maza and Deniz Yuret)
// A local optimiser which finds local minima and gets restarted as far as possible from previously evaluated points.
class DynamicHillClimberSolver : public AuxSolver{
	public:

	// Constructor which inputs the function being optimised, whether it is a minimisation or maximisation problem,
	// a random seed for reproduction purposes (with random behaviour when randomSeed = 0), and the boolean printInfo
	// which prints information as the optimisation proceeds. Using this constructor assumes this algorithm should 
	// have a maximum of 5000 samples.
	DynamicHillClimberSolver(Function* function, bool min = true, int randomSeed = 0, bool printInfo = false);

	// Constructor with the additional input maxEval, which specifies the maximum number of function samples.
	DynamicHillClimberSolver(Function* function, int maxEval, bool min = true, int randomSeed = 0, bool printInfo = false);

	virtual ~DynamicHillClimberSolver() override;

	// Main method which finds best possible point given budget; returns best point found.
	virtual VectorXd optimise() override;

	// Extra function so that the initial point is given, with heavy tune in mind
	VectorXd optimise(VectorXd startingPoint, double startingPointVal);

	// A function which initialises all the potential directions for future steps of the algorithm.
	// For every dimension of the problem, two potential directions are defined with a magnitude of 
	// half of the domain of that particular dimension.
	vector<VectorXd> initialiseDirections();
	
};

// Implementation of Toal's Heavy Tune for Kriging hyperparameter tuning (see "Kriging hyperparameter tuning strategies" by David JJ Toal, Neil W. Bressloff, and Andy J. Keane)
// Note that not all the details are given in the paper, only a rough outline. This implementation filled in the gaps where the information was not available in the paper, 
// mainly with the parameters of the hill climber and the genetic algorithm.
// Essentially this heavy tune solver first runs Amoeba to convergence twice, then runs a Genetic Algorithm, and finally a Dynamic Hill Climber. 
class HeavyTuneSolver : public AuxSolver{
	public:

	// Constructor which inputs the function being optimised, whether it is a minimisation or maximisation problem,
	// a random seed for reproduction purposes (with random behaviour when randomSeed = 0), and the boolean printInfo
	// which prints information as the optimisation proceeds. Using this constructor assumes this algorithm should 
	// have a maximum of 10000 samples, and the genetic algorithm will use a population of size 50 at each iteration.
	HeavyTuneSolver(Function* function, bool min = true, int randomSeed = 0, bool printInfo = false);

	// Constructor with the additional input maxEval, which specifies the maximum number of function samples and the population size of 
	// the GA algorithm. It uses at most 50d iterations of an amoeba solver twice, followed by an even split of the remaining evaluations
	// between GA and hill climber.
	HeavyTuneSolver(Function* function, int maxEval, int popSizeGA, bool min = true, int randomSeed = 0, bool printInfo = false);

	virtual ~HeavyTuneSolver() override;

	// Main method which finds best possible point given budget; returns best point found.
	virtual VectorXd optimise() override;

	int popSizeGA_;							// Population size of GA solver

};

// Custom solver inspired by Toal's Heavy Tuner. Does two rounds of Amoeba search, followed by 80% of the total remaining evaluations
// using Accelerated Random Search (ARS) and 20% of the total remaining evaluations on hill climbing. The idea is to explore the space
// with ARS and then conduct a final round of local optimisation. Seems to perform no worse (and at times better) than Heavy Tune.
class RandomHeavyTuneSolver : public AuxSolver{
	public:

	// Constructor which inputs the function being optimised, whether it is a minimisation or maximisation problem,
	// a random seed for reproduction purposes (with random behaviour when randomSeed = 0), and the boolean printInfo
	// which prints information as the optimisation proceeds. Using this constructor assumes this algorithm should 
	// have a maximum of 10000 samples, and ARS should use 10 spheres.
	RandomHeavyTuneSolver(Function* function, bool min = true, int randomSeed = 0, bool printInfo = false);

	// Constructor with the additional input maxEval, which specifies the maximum number of function samples, and numSphereARS which specifies
	// the number of spheres in the ARS algorithm. It uses at most 50d iterations of an amoeba solver twice, followed by a 80-20 split of the 
	// remaining evaluations between ARS and hill climber.
	RandomHeavyTuneSolver(Function* function, int maxEval, int numSphereARS, bool min = true, int randomSeed = 0, bool printInfo = false);

	virtual ~RandomHeavyTuneSolver() override;

	// Main method which finds best possible point given budget; returns best point found.
	virtual VectorXd optimise() override;

	int numSphereARS_;						// Number of spheres when running ARS solver
};

#endif