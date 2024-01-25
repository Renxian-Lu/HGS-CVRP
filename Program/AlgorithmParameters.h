//
// Created by chkwon on 3/23/22.
//

// This header file must be readable in C.

#ifndef ALGORITHMPARAMETERS_H
#define ALGORITHMPARAMETERS_H

struct AlgorithmParameters {
	int nbGranular;			// Granular search parameter, limits the number of moves in the RI local search
	int mu;					// Minimum population size
	int lambda;				// Number of solutions created before reaching the maximum population size (i.e., generation size)
	int nbElite;			// Number of elite individuals
	int nbClose;			// Number of closest solutions/individuals considered when calculating diversity contribution

	int nbIterPenaltyManagement;  // Number of iterations between penalty updates
	double targetFeasible;	      // Reference proportion for the number of feasible individuals, used for the adaptation of the penalty parameters
	double penaltyDecrease;	      // Multiplier used to decrease penalty parameters if there are sufficient feasible individuals
	double penaltyIncrease;	      // Multiplier used to increase penalty parameters if there are insufficient feasible individuals

	int seed;				// Random seed. Default value: 0
	int nbIter;				// Nb iterations without improvement until termination (or restart if a time limit is specified). Default value: 20,000 iterations
	int nbIterTraces;       // Number of iterations between traces display during HGS execution
	double timeLimit;		// CPU time limit until termination in seconds. Default value: 0 (i.e., inactive)
	int useSwapStar;		// Use SWAP* local search or not. Default value: 1. Only available when coordinates are provided.

	/* Add by Renxian Lu
	   Use crossover operator or not. Default value: "OX". 
	   1: Ordered crossover (OX)
	   2: Cycle Crossover (CX)
	   3: Partially Mapped Crossover (PMX)
	   4: Edge Recombination Crossover (ERX)
	   5: Heuristic Crossover (HX)
	   6: OX and PMX (each with same probability)
	   7: OX, PMX and CX (each with same probability)
	   8: OX, PMX, CX and HX (each with same probability)
	   9: All crossover operators (each with same probability)
	   10: AOS (Adaptive Operator Selection)
	   else: OX
	*/
	int useCrossover;
	int iterationAOS;       // Number of iterations (only for AOS)
	// End of addition
};

#ifdef __cplusplus
extern "C"
#endif
struct AlgorithmParameters default_algorithm_parameters();

#ifdef __cplusplus
void print_algorithm_parameters(const AlgorithmParameters & ap);
#endif

#endif //ALGORITHMPARAMETERS_H
