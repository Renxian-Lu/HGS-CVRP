/*MIT License

Copyright(c) 2020 Thibaut Vidal

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#ifndef GENETIC_H
#define GENETIC_H

#include "Population.h"
#include "Individual.h"

struct DoubleIndexValuePair {
    int index1;
    int value1;
    int index2;
    int value2;
};

class Genetic
{
public:

	Params & params;				// Problem parameters
	Split split;					// Split algorithm
	LocalSearch localSearch;		// Local Search structure
	Population population;			// Population (public for now to give access to the solutions, but should be be improved later on)
	Individual offspring;			// First individual to be used as input for the crossover
	typedef void (Genetic::*CrossoverFunction)(Individual&, const Individual&, const Individual&);
	std::vector<double> crossoverProbabilities = {0.2, 0.2, 0.2, 0.2, 0.2}; // Probability of each crossover operator


	// OX Crossover
	void crossoverOX(Individual & result, const Individual & parent1, const Individual & parent2);
	// CX Crossover
	void crossoverCX(Individual & result, const Individual & parent1, const Individual & parent2);
	// PMX Crossover
	void crossoverPMX(Individual & result, const Individual & parent1, const Individual & parent2);
	// Edge Recombination Crossover
	void crossoverER(Individual & result, const Individual & parent1, const Individual & parent2);
	// Heuristic Crossover
	void crossoverHX(Individual & result, const Individual & parent1, const Individual & parent2);
	// Find a mapping element in the parent2 for PMX
	int findElementInParent2(int start, int end, int index, std::vector<int>& subvectorP2, const Individual & parent1, const Individual & parent2, Individual & result);
	void crossoverSelection(Individual & result, const Individual & parent1, const Individual & parent2);
	int rwsSelection(Individual & result, const Individual & parent1, const Individual & parent2, int nbIter);

	// Running the genetic algorithm until maxIterNonProd consecutive iterations or a time limit
    void run() ;

	// Constructor
	Genetic(Params & params);

private:
    std::mt19937 randomEngine{std::random_device{}()};
	std::vector<CrossoverFunction> crossoverFunctions;
	int k = 0;
};

#endif
