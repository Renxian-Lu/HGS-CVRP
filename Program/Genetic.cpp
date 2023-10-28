#include "Genetic.h"
#include <unordered_map>
#include<algorithm>
#include <unordered_set>
#include <cstdlib>  // For rand() and srand()
#include <random>

void Genetic::run()
{	
	/* INITIAL POPULATION */
	population.generatePopulation();

	int nbIter;
	int nbIterNonProd = 1;
	if (params.verbose) std::cout << "----- STARTING GENETIC ALGORITHM" << std::endl;

	for (nbIter = 0 ; nbIterNonProd <= params.ap.nbIter && (params.ap.timeLimit == 0 || (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC < params.ap.timeLimit) ; nbIter++)
	{	
		/* SELECTION AND CROSSOVER */
		if(params.ap.useCrossover == 1)
			crossoverOX(offspring, population.getBinaryTournament(),population.getBinaryTournament());
		else if (params.ap.useCrossover == 2)
			crossoverCX(offspring, population.getBinaryTournament(),population.getBinaryTournament());
		else if (params.ap.useCrossover == 3)
			crossoverPMX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
		else if (params.ap.useCrossover == 4)
			crossoverER(offspring, population.getBinaryTournament(), population.getBinaryTournament());
		else if (params.ap.useCrossover == 5)
			crossoverHX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
		else if (params.ap.useCrossover == 6) 
		{
			// Choose one of the crossover randomly OX or PMX and each crossover has 50% chance
			// Seed the random number generator with the current iteration
    		srand(nbIter);
			// Generate a random number between 0 and 1 (50% chance)
    		int randomChoice = rand() % 2;

			if (randomChoice == 0) 
			{
				crossoverOX(offspring, population.getBinaryTournament(),population.getBinaryTournament()); // 50% chance: Perform OX
			} else {
				crossoverPMX(offspring, population.getBinaryTournament(), population.getBinaryTournament()); // 50% chance: Perform PMX
			}
		} 
		else if (params.ap.useCrossover == 7)
		{
			// Seed the random number generator with the current time
			std::mt19937 rng(nbIter);

			// Create a uniform distribution to choose between the methods
    		std::uniform_int_distribution<int> dist(0, 2);  // 0, 1, or 2 for the three methods

			// Generate a random number to choose the method
    		int selectedMethod = dist(rng);

			// Choose the method based on the random number
			switch (selectedMethod) {
				case 0:
					crossoverOX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
					break;
				case 1:
					crossoverPMX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
					break;
				case 2:
					crossoverCX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
					break;
			}
		}
		else if (params.ap.useCrossover == 8)
		{
			// Seed the random number generator with the current time
			std::mt19937 rng(nbIter);

			// Create a uniform distribution to choose between the methods
    		std::uniform_int_distribution<int> dist(0, 3);  // 0, 1, 2 or 3 for the three methods

			// Generate a random number to choose the method
    		int selectedMethod = dist(rng);

			// Choose the method based on the random number
			switch (selectedMethod) {
				case 0:
					crossoverOX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
					break;
				case 1:
					crossoverPMX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
					break;
				case 2:
					crossoverCX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
					break;
				case 3:
					crossoverHX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
					break;
			}
		}
		else if (params.ap.useCrossover == 9)
		{
			// Seed the random number generator with the current time
			std::mt19937 rng(nbIter);

			// Create a uniform distribution to choose between the methods
    		std::uniform_int_distribution<int> dist(0, 4);  // 0, 1, 2, 3 or 4 for the three methods

			// Generate a random number to choose the method
    		int selectedMethod = dist(rng);

			// Choose the method based on the random number
			switch (selectedMethod) {
				case 0:
					crossoverOX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
					break;
				case 1:
					crossoverPMX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
					break;
				case 2:
					crossoverCX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
					break;
				case 3:
					crossoverER(offspring, population.getBinaryTournament(), population.getBinaryTournament());
					break;
				case 4:
					crossoverHX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
					break;
			}
		}
		else if (params.ap.useCrossover == 10) 
		{	
			// std::cout << "----- k:" << k << std::endl;
			( this->*crossoverFunctions[0] )(offspring, population.getBinaryTournament(), population.getBinaryTournament());
			// crossoverSelection(offspring, population.getBinaryTournament(),population.getBinaryTournament());	
		} 
		else {
			crossoverOX(offspring, population.getBinaryTournament(),population.getBinaryTournament());
		}

		/* LOCAL SEARCH */
		localSearch.run(offspring, params.penaltyCapacity, params.penaltyDuration);
		bool isNewBest = population.addIndividual(offspring,true);
		if (!offspring.eval.isFeasible && params.ran()%2 == 0) // Repair half of the solutions in case of infeasibility
		{
			localSearch.run(offspring, params.penaltyCapacity*10., params.penaltyDuration*10.);
			if (offspring.eval.isFeasible) isNewBest = (population.addIndividual(offspring,false) || isNewBest);
		}

		/* TRACKING THE NUMBER OF ITERATIONS SINCE LAST SOLUTION IMPROVEMENT */
		if (isNewBest) nbIterNonProd = 1;
		else nbIterNonProd ++ ;

		/* DIVERSIFICATION, PENALTY MANAGEMENT AND TRACES */
		if (nbIter % params.ap.nbIterPenaltyManagement == 0) population.managePenalties();
		if (nbIter % params.ap.nbIterTraces == 0) population.printState(nbIter, nbIterNonProd);

		/* FOR TESTS INVOLVING SUCCESSIVE RUNS UNTIL A TIME LIMIT: WE RESET THE ALGORITHM/POPULATION EACH TIME maxIterNonProd IS ATTAINED*/
		if (params.ap.timeLimit != 0 && nbIterNonProd == params.ap.nbIter)
		{
			population.restart();
			nbIterNonProd = 1;
		}
	}
	if (params.verbose) std::cout << "----- GENETIC ALGORITHM FINISHED AFTER " << nbIter << " ITERATIONS. TIME SPENT: " << (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC << std::endl;
}

void Genetic::crossoverOX(Individual & result, const Individual & parent1, const Individual & parent2)
{
	// std::cout << "----- OX is used" << std::endl;

	// Frequency table to track the customers which have been already inserted
	std::vector <bool> freqClient = std::vector <bool> (params.nbClients + 1, false);

	// Picking the beginning and end of the crossover zone
	std::uniform_int_distribution<> distr(0, params.nbClients-1);
	int start = distr(params.ran);
	int end = distr(params.ran);

	// Avoid that start and end coincide by accident
	while (end == start) end = distr(params.ran);

	// Copy from start to end
	int j = start;
	while (j % params.nbClients != (end + 1) % params.nbClients)
	{
		result.chromT[j % params.nbClients] = parent1.chromT[j % params.nbClients];
		freqClient[result.chromT[j % params.nbClients]] = true;
		j++;
	}

	// Fill the remaining elements in the order given by the second parent
	for (int i = 1; i <= params.nbClients; i++)
	{
		int temp = parent2.chromT[(end + i) % params.nbClients];
		if (freqClient[temp] == false)
		{
			result.chromT[j % params.nbClients] = temp;
			j++;
		}
	}

	// Complete the individual with the Split algorithm
	split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverCX(Individual & result, const Individual & parent1, const Individual & parent2)
{
	// Frequency table to track the customers which have been already inserted
	std::vector <bool> freqClient = std::vector <bool> (params.nbClients + 1, false);


	std::unordered_map<int, DoubleIndexValuePair> HT1;
	for (int j = 0; j < params.nbClients; ++j) 
	{
		DoubleIndexValuePair tempPair;
		tempPair.value2 = parent2.chromT[j % params.nbClients];
		tempPair.index2 = j;
		tempPair.value1 = parent1.chromT[j % params.nbClients];
		HT1[parent2.chromT[j % params.nbClients]] = tempPair;
    }

	std::vector<std::vector<DoubleIndexValuePair>> Cycles;
	for (int j = 0; j < params.nbClients; ++j) {
        if (!freqClient[j]) {
            int cycleStart = j;
            std::vector<DoubleIndexValuePair> tempCycle;

            do {
                DoubleIndexValuePair tempPair = HT1[parent1.chromT[cycleStart]];
                tempCycle.push_back(tempPair);
                freqClient[tempPair.index2] = true;
                cycleStart = tempPair.index2;
            } while (cycleStart != j);

            Cycles.push_back(tempCycle);
        }
    }

	int counter = 0;
    for (const auto &C : Cycles) {
        for (const auto &tempPair : C) {
            if (counter % 2 == 0) {
                result.chromT[tempPair.index2] = tempPair.value1;
                // child2.alleles[tempPair.index2] = tempPair.value2;
            } else {
                result.chromT[tempPair.index2] = tempPair.value2;
                // child2.alleles[tempPair.index2] = tempPair.value1;
            }
        }
        ++counter;
    }

	// Complete the individual with the Split algorithm
	split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverPMX(Individual & result, const Individual & parent1, const Individual & parent2) 
{
	// Frequency table to track the customers which have been already inserted
	std::vector <bool> freqClient = std::vector <bool> (params.nbClients + 1, false);

    int size = params.nbClients;

	// Picking the beginning and end of the crossover zone
	std::uniform_int_distribution<> distr(0, size-1);
	int start = distr(params.ran);
	int end = distr(params.ran);

	// Avoid that start and end coincide by accident
	while (end == start) end = distr(params.ran);

	// Copy from start to end
	int j = start;
	while (j % params.nbClients != (end + 1) % params.nbClients)
	{
		result.chromT[j % params.nbClients] = parent1.chromT[j % params.nbClients];
		freqClient[result.chromT[j % params.nbClients]] = true;
		j++;
	}

	std::vector<int> subvectorP1;
	std::vector<int> subvectorP2;

	if(start > end)
	{
		// Get the subvector of 1. parent from start to end
		subvectorP1.assign(parent1.chromT.begin() + start, parent1.chromT.end());
		subvectorP1.insert(subvectorP1.end(), parent1.chromT.begin(), parent1.chromT.begin() + end + 1);
		// Get the subvector of 2. parent from start to end
		subvectorP2.assign(parent2.chromT.begin() + start, parent2.chromT.end());
		subvectorP2.insert(subvectorP2.end(), parent2.chromT.begin(), parent2.chromT.begin() + end + 1);
	} else {
		// Get the subvector of 1. parent from start to end
		subvectorP1.assign(parent1.chromT.begin() + start, parent1.chromT.begin() + end + 1);
		// Get the subvector of 2. parent from start to end
		subvectorP2.assign(parent2.chromT.begin() + start, parent2.chromT.begin() + end + 1);

	}

	int segmentLength = subvectorP1.size();

	// Write a for to loop the element in segment
	for (int i = 0; i < segmentLength; ++i) 
	{
		int index = start + i;
		int indexP2;

		// check if the element segment[i] is within parent2.chromT[start:end]
		if (!std::count(subvectorP1.begin(), subvectorP1.end(), parent2.chromT[index % params.nbClients]))
		{
			// Find the corresponding mapping element in parent2
			indexP2 = findElementInParent2(start, end, index, subvectorP2, parent1, parent2, result);
			result.chromT[indexP2] = parent2.chromT[index % params.nbClients];
			freqClient[parent2.chromT[index % params.nbClients]] = true;
		}
	}

	// Fill the remaining elements according the order by parent2
	for (int i = 0; i < params.nbClients; i++)
	{
		int temp = parent2.chromT[i];
		if (freqClient[temp] == false)
		{
			result.chromT[i] = temp;
		}
	}

	// Complete the individual with the Split algorithm
	split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverER(Individual & result, const Individual & parent1, const Individual & parent2) 
{
	// Create Neighbor List
	std::unordered_map<int, std::unordered_set<int>> neighborList;
	int alleleCount = params.nbClients;

	for (int j = 0; j < alleleCount; ++j) {
		neighborList[parent1.chromT[j]] = {};
	}

	for (int j = 0; j < alleleCount; ++j) {
		int prevIndex = (j - 1 + alleleCount) % alleleCount;
		int nextIndex = (j + 1) % alleleCount;

		neighborList[parent1.chromT[j]].insert(parent1.chromT[prevIndex]);
		neighborList[parent1.chromT[j]].insert(parent1.chromT[nextIndex]);

		neighborList[parent2.chromT[j]].insert(parent2.chromT[prevIndex]);
		neighborList[parent2.chromT[j]].insert(parent2.chromT[nextIndex]);
	}

	// Generate a random number 0 or 1 to select the parent
    int randomNum = std::rand() % 2;
    const Individual& selectedParent = (randomNum == 0) ? parent1 : parent2;

	int currentIndex = 0;
	int currentNode = selectedParent.chromT[0];

	std::unordered_set<int> visitedNodes = { currentNode };
	std::unordered_set<int> unvisitedNodes(selectedParent.chromT.begin(), selectedParent.chromT.end());
	unvisitedNodes.erase(currentNode);

	// Perform Edge Recombination Crossover
	while (currentIndex < alleleCount - 1) 
	{
		result.chromT[currentIndex] = currentNode;
		++currentIndex;

		visitedNodes.insert(currentNode);

		// Crossing out the current node from the neighbor lists
		for (auto& [node, neighbors] : neighborList) 
		{
			neighbors.erase(currentNode);
		}

		// Choose next node with the fewest neighbors
		int minNeighborCount = std::numeric_limits<int>::max();
		std::vector<int> nextNodes;

		for (int neighbor : neighborList[currentNode]) 
		{
			// Count the number of neighbors of the neighbor
			int neighborCount = neighborList[neighbor].size();

			// Check if the current neighbor has fewer neighbors than the current minimum
			if (neighborCount < minNeighborCount) {
				minNeighborCount = neighborCount;
				nextNodes.clear();
				nextNodes.push_back(neighbor);
			} else if (neighborCount == minNeighborCount) {
				nextNodes.push_back(neighbor);
			}
		}

		// If there are no nodes, choose a random unvisited node
		if (nextNodes.empty()) 
		{
			std::vector<int> unvisitedNodesVec(unvisitedNodes.begin(), unvisitedNodes.end());
			std::shuffle(unvisitedNodesVec.begin(), unvisitedNodesVec.end(), params.ran);
			currentNode = unvisitedNodesVec.front();
			unvisitedNodes.erase(currentNode);
		} else {
			// Choose a random node from the list of nodes with fewest neighbors
			std::shuffle(nextNodes.begin(), nextNodes.end(), params.ran);
			currentNode = nextNodes.front();
			unvisitedNodes.erase(currentNode);
		}

		visitedNodes.insert(currentNode);
		unvisitedNodes.erase(currentNode);
	}

	// Add the last node
	result.chromT[currentIndex] = currentNode;

	// Complete the individual with the Split algorithm
	split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverHX(Individual & result, const Individual & parent1, const Individual & parent2)
{	
	// Frequency table to track the customers which have been already inserted
	std::vector <bool> freqClient = std::vector <bool> (params.nbClients + 1, false);

	// Picking a random starting point exclude the depot
	std::uniform_int_distribution<> distr(1, params.nbClients);
	int start = distr(params.ran);
	int j = start;
	result.chromT[0] = j;
	freqClient[j] = true;

	int i = 1;
	do
	{
		// Get the position of node j in parent1 and parent2
		int indexP1 = std::distance(parent1.chromT.begin(), std::find(parent1.chromT.begin(), parent1.chromT.end(), j));
		int indexP2 = std::distance(parent2.chromT.begin(), std::find(parent2.chromT.begin(), parent2.chromT.end(), j));

		// Get the successor of node j in parent1 and parent2
		int successorP1 = (indexP1 == params.nbClients - 1) ? parent1.chromT[0] : parent1.chromT[indexP1+1];
		int successorP2 = (indexP2 == params.nbClients - 1) ? parent2.chromT[0] : parent2.chromT[indexP2+1];

		// Check if successorP1 and successorP2 are true in the freqClient
		bool isP1 = freqClient[successorP1];
		bool isP2 = freqClient[successorP2];

		// If both are available, then choose the one with the shorter distance edge
		if (!isP1 && !isP2)
		{
			double distanceP1 = params.timeCost[j][successorP1];
			double distanceP2 = params.timeCost[j][successorP2];
			if (distanceP1 < distanceP2)
			{
				result.chromT[i] = successorP1;
				freqClient[successorP1] = true;
				j = successorP1;
			} else {
				result.chromT[i] = successorP2;
				freqClient[successorP2] = true;
				j = successorP2;
			}
		} else if (!isP1) {
			result.chromT[i] = successorP1;
			freqClient[successorP1] = true;
			j = successorP1;
		} else if (!isP2) {
			result.chromT[i] = successorP2;
			freqClient[successorP2] = true;
			j = successorP2;
		} else {
			// If no available nodes from node j, then choose the other nodes with the shortest distance edge
			double minDistance = std::numeric_limits<double>::max();
			int minIndex = 0;
			for (int k = 1; k < params.nbClients+1; k++)
			{
				if (params.timeCost[j][k] < minDistance && !freqClient[k] && k != j)
				{
					minDistance = params.timeCost[j][k];
					minIndex = k;
				}
			}
			result.chromT[i] = minIndex;
			freqClient[minIndex] = true;
			j = minIndex;
		}
		i++;
	} while (i < params.nbClients);

	// Complete the individual with the Split algorithm
	split.generalSplit(result, parent1.eval.nbRoutes);
}

int Genetic::findElementInParent2(int start, int end, int index, std::vector<int>& subvectorP2, const Individual & parent1, const Individual & parent2, Individual & result)
{
	int elementP1 = parent1.chromT[index % params.nbClients];
	int elementP2 = parent2.chromT[index % params.nbClients];

	// Index of element with value elementP1 in parent2
	int indexP2 = std::distance(parent2.chromT.begin(), std::find(parent2.chromT.begin(), parent2.chromT.end(), elementP1));

	// Check if the corresponding mapping element in parent1 is in the subvector
	auto it = std::find(subvectorP2.begin(), subvectorP2.end(), elementP1);

	// If yes, then find further, if no, return the index of parent2
	if (it != subvectorP2.end()) {
		return findElementInParent2(start, end, indexP2, subvectorP2, parent1, parent2, result);
	} else {
		return indexP2;
	}
}

// Todo: Implement the crossover selection
void Genetic::crossoverSelection(Individual & result, const Individual & parent1, const Individual & parent2)
{

}

// Todo: roulette wheel selection
int Genetic::rwsSelection(Individual & result, const Individual & parent1, const Individual & parent2, int nbIter)
{
	int crossoverMethodIndex = 0;
	// Seed the random number generator with the current time
	std::mt19937 rng(nbIter);
	// Create a uniform real distribution in the range [0, 1]
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double r = dist(rng);
	// Calculate the cumulative probability
  	double cumulativeProb = 0;

	// Find the operator whose probability interval contains r
	for (const double& probability : crossoverProbabilities) {
		cumulativeProb += probability;
		if (cumulativeProb >= r) {
			// Get the index of the crossover method
			crossoverMethodIndex = std::distance(crossoverProbabilities.begin(), std::find(crossoverProbabilities.begin(), crossoverProbabilities.end(), probability));
			return crossoverMethodIndex;
		}
	}
	return crossoverProbabilities.back(); // Return the last crossover method if no one is found (this should not happen if probabilities are normalized)
}

Genetic::Genetic(Params & params) : 
	params(params), 
	split(params),
	localSearch(params),
	population(params,this->split,this->localSearch),
	offspring(params){
		// Initialize the crossover functions vector with your crossover methods
    crossoverFunctions = 
	{
        &Genetic::crossoverOX,
        &Genetic::crossoverCX,
        &Genetic::crossoverPMX,
        &Genetic::crossoverER,
        &Genetic::crossoverHX
    };
	}

