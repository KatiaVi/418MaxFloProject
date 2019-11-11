//
// Created by Ekaterina Villevald on 2019-11-02.
//

#include <cstdlib>
#include <ctime>
#include <cmath>
#include <numeric>
#include <iostream>
#include <queue>
#include <algorithm>
#include <random>

#include "genetic.h"


// Todo: (1) Malloc + free memory for all arrays
// Todo: (2) Address the other todos
// Todo: (3) Attempt to compile and run this


/* foundValidSolution
 * checks whether there exists a solution where the vertices are locally balanced and all vertices are blocked
 * updates BALANCE and BLOCKED
 */
int GeneticSequentialSolver::foundValidSolution(int num_vertices, int source, int sink, float** &capacities){

  //TODO: if only one solution left in the population that solution is returned (not sure if this is right)
  if (currentPopulationSize == 1){
    return 0;
  }
  for (int p=0; p < currentPopulationSize; p++){

    // Check Local Balance Property of solutions
    for (int i=0; i < num_vertices; i++) {
      BALANCE[p][i] = 1;
      if (i == sink)  continue;
      // vertex not locally balanced (input flow != output flow at some vertex)
      if (inputFlow[p][i] != outputFlow[p][i]){
        BALANCE[p][i] = 0;
      }
    }

    // Check that flow at source vertex is blocked
    // (all edges in the flow are saturated meaning fij = cij and vertex j is blocked)
    std::queue<int> vertexQ;
    int currentVertex = sink;
    BLOCKED[p][sink] = false; //TODO: is this the right value to initialize BLOCKED[p][sink] to?

    for (int i=0; i < num_vertices; i++){
      if (capacities[i][currentVertex] > 0) vertexQ.push(i);
    }

    while (not vertexQ.empty()){
      currentVertex = vertexQ.front();
      vertexQ.pop();
      BLOCKED[p][currentVertex] = true;

      for (int i=0; i < num_vertices; i++){
        if (capacities[currentVertex][i] > 0 and
        not (solutions[p][currentVertex][currentVertex] == capacities[i][currentVertex] or BLOCKED[p][currentVertex])){
            BLOCKED[p][currentVertex] = false;
        }
      }
      if (not BLOCKED[p][currentVertex]) break;
    }

    if (BLOCKED[p][source] and std::accumulate(BALANCE[p], BALANCE[p]+num_vertices, 0) == num_vertices) return p;
  }
  return -1;
}


/* initialize
 * creates the solutions array and initializes the flow randomly so that each edge
 * has a flow that is between 0 and the capacity of that edge
 */
void GeneticSequentialSolver::initialize(int num_vertices, float** &capacities){
  srand (static_cast <unsigned> (time(0)));

  for (int p = 0; p < currentPopulationSize; p++) {
    float **solution = new float*[num_vertices];
    bool *blockedP = new bool[num_vertices];
    int *balanceP = new int[num_vertices];

    for (int i = 0; i < num_vertices; i++) {

      solution[i] = new float[num_vertices];

      for (int j = 0; j < num_vertices; j++) {
        solution[i][j] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/capacities[i][j]));
      }
      blockedP[i] = false;
      balanceP[i] = 0;
    }
    excessFlow.push_back(new float[num_vertices]);
    inputFlow.push_back(new float[num_vertices]);
    outputFlow.push_back(new float[num_vertices]);
    totalFlow.push_back(new float[num_vertices]);

    BALANCE.push_back(balanceP);
    BLOCKED.push_back(blockedP);

    solutions.push_back(solution);
  }
}

/* evaluate
 * computes the excessFlow, inputFlow, outputFlow and totalFlow matrices for a generation
 * checks locally balanced condition and computes what flows are blocked if necessary
 * TODO: will we need to recalculate all of these for every individual for every generation
 * I'm thinking no since reproduction will only add 1 new individual, the others shouldn't be
 * changed? (think some more about this)
 */
void GeneticSequentialSolver::evaluate(int num_vertices, int sink, int source){
  for (int p=0; p < currentPopulationSize; p++) {
    for (int i = 0; i < num_vertices; i++) {
      for (int j=0; j < num_vertices; j++){
        inputFlow[p][i] += solutions[p][j][i];
        outputFlow[p][i] += (-solutions[p][i][j]);
      }
      excessFlow[p][i] = inputFlow[p][i] + outputFlow[p][i];
      totalFlow[p][i] = inputFlow[p][i] - outputFlow[p][i];
    }
    excessFlow[p][sink] = inputFlow[p][sink] + outputFlow[p][source];
    totalFlow[p][sink] = inputFlow[p][sink];
  }
}

/* reproduction
 * calculates the fitness scores (and probabilities)
 * for all solutions and chooses 2 solutions randomly with the calculated probabilities
 * do a crossover + mutate if necessary
 */
void GeneticSequentialSolver::reproduction(int num_vertices, float** capacities, int source, int sink){
  std::vector<float> fitnessScores;
  float fitnessSum = 0;

  for (int p=0; p < currentPopulationSize; p++){
    float sumOfFlow = 0.0;
    float sumOfCapacities = 0.0;
    for (int i=0; i < num_vertices; i++){
      sumOfFlow  += (std::accumulate(solutions[p][i], solutions[p][i] + num_vertices, 0.0));
      sumOfCapacities += (std::accumulate(capacities[i], capacities[i] + num_vertices, 0.0));
    }
    fitnessScores.push_back(
        (std::accumulate(BALANCE[p], BALANCE[p]+num_vertices, 0.0))/(num_vertices - 2.0f) -
        abs(std::accumulate(excessFlow[p], excessFlow[p]+num_vertices, 0.0))/sumOfFlow +
        sumOfFlow/sumOfCapacities);
    fitnessSum += fitnessScores[p];
  }

  //std::cout<<"Got all the fitness scores\n";

  std::random_device rd1;
  std::mt19937 rng1(rd1());
  std::discrete_distribution<> dist1(fitnessScores.begin(), fitnessScores.end());

  int individual1 = dist1(rng1);
  fitnessScores.erase(fitnessScores.begin() + individual1);

  std::random_device rd2;
  std::mt19937 rng2(rd2());
  std::discrete_distribution<> dist2(fitnessScores.begin(), fitnessScores.end());
  int individual2 = dist2(rng2);
  if (individual2 >= individual1){
    individual2 += 1;
  }
  //std::cout << "ind1: " << individual1 << "ind 2: " << individual2 << "\n";

  //std::cout << "generated random individuals for crossOver\n";


  float** newSolution = crossOver(individual1, individual2, num_vertices, source, sink, capacities);

 // std::cout << "finished crossOver\n";

  // delete the 2 old solutions and add new solution to generation
  // erase the old solution with higher index first so that the index of
  // the other old solutions doesn't change when you want to erase it
  (individual1 > individual2) ? solutions.erase(solutions.begin() + individual1) : solutions.erase(solutions.begin() + individual2);
  (individual1 > individual2) ? solutions.erase(solutions.begin() + individual2) : solutions.erase(solutions.begin() + individual1);
  solutions.push_back(newSolution);

  currentPopulationSize -= 1;
}

float cap(float value, float capacity){
  if (value < 0) return 0;
  if (value > capacity) return capacity;
  return value;
}

/* crossOver
 * compares the energies of all pairs of vertices in 2 solutions and updates the final solution
 * matrix with the dominant vertex, also can mutate with a probability determined in the mutation function
 * if the vertex has already been updates in the crossOver, then need to take average of old vertex
 * and the new vertex
 * TODO: since this step decreases the size of the population by 1 is there a way to ignore this
 * population in future steps (i.e. when we call re-evaluate after the reproduce step) (DONE: changed solutions
 * to a vector)
 * */
float** GeneticSequentialSolver::crossOver(int solution1, int solution2, int num_vertices, int source, int sink, float** capacities){
  float** newSolution = new float*[num_vertices];
  for (int i=0; i < num_vertices; i++){
    newSolution[i] = new float[num_vertices];
  }
//  std::cout<<"set newSolutionStuff1\n";


  char newEdgeAdded[num_vertices*num_vertices]; //Todo: maybe change to bitstring to make more efficient

  //TODO: paper says the source and sink vertices are not used during crossover. But what should they be set to then?
  newSolution[source][source] = 0.0f;
  newSolution[sink][sink] = 0.0f;

  //std::cout<<"set newSolutionStuff2\n";

  for (int i=0; i < num_vertices; i++){

    int inCapacity = 0;
    int outCapacity = 0;
    for (int j=0; j < num_vertices; j++){
      inCapacity += capacities[j][i];
      outCapacity += capacities[i][j];
    }

    float energy1 = k*abs(excessFlow[solution1][i]) +
        abs(std::min(inCapacity, outCapacity) - std::max(inputFlow[solution1][i], -1.0f*outputFlow[solution1][i]));
    float energy2 = k*abs(excessFlow[solution2][i]) +
        abs(std::min(inCapacity, outCapacity) - std::max(inputFlow[solution2][i], -1.0f*outputFlow[solution2][i]));

    int newSolutionIndex = (energy1 < energy2) ? solution1 : // energy level of solution1 is lower
                          ((energy2 < energy1) ? solution2 : // energy level of solution2 is lower
                          ((rand() > RAND_MAX/2) ? solution1 : solution2)); // randomly choose solution1 or solution2

    for (int j=0; j < num_vertices; j++){
      // add new edges (j->i and i->j) to the newSolution
      // possibly mutate the flow through the edge
      float newEdgeValueIJ = solutions[newSolutionIndex][i][j];
      float newEdgeValueJI = solutions[newSolutionIndex][j][i];

      float mutationAdjustment1 = mutation(i, j, energy1);
      float mutationAdjustment2= mutation(i, j, energy2);


      if (newEdgeAdded[i*num_vertices + j]) {
        newSolution[i][j] = cap((newSolution[i][j] + newEdgeValueIJ) / 2.0f + mutationAdjustment1, capacities[i][j]);
      }
      else{
        newSolution[i][j] = cap(newEdgeValueIJ + mutationAdjustment1, capacities[i][j]);
        newEdgeAdded[i*num_vertices + j] =  1;
      }


      if (newEdgeAdded[j*num_vertices + i]) {
        newSolution[j][i] = cap((newSolution[j][i] + newEdgeValueJI) / 2.0f + mutationAdjustment2, capacities[j][i]);
      }
      else{
        newSolution[j][i] = cap(newEdgeValueJI + mutationAdjustment2, capacities[j][i]);
        newEdgeAdded[j*num_vertices + i] =  1;
      }
    }
  }
  return newSolution;
}

/* mutation
 * determine probability of mutation for a vertex
 * use that probability to determine whether a mutation will occur or not
 * if it will mutate in the direction that will lower the energy level of the vertex
 * if it is above 0, and if it is at 0 mutate in either direction
 * */
float GeneticSequentialSolver::mutation(int individualIndex, int vertex, float energyLevel){
  float mutationProb = mutation_rate + sqrt(abs(excessFlow[individualIndex][vertex])/ totalFlow[individualIndex][vertex]);

  std::random_device rd;
  std::mt19937 gen1(rd());
  std::discrete_distribution<> dist{ mutationProb, 1 - mutationProb };
  int willMutate = dist(gen1);

  if (willMutate and energyLevel == 0){
    return ((rand() > RAND_MAX/2) ? mutation_increment : -1.0f * mutation_increment);
  }
  if (willMutate and energyLevel > 0){
    return -1.0f * mutation_increment;
    /* TODO: not sure about this, paper says vertex with an energy level above
     * zero mutation always occurs in the direction which will lower the energy level of the vertex
     * Why is this the case?
     */
  }
  return 0.0f;
}

/* solve
 * initializes the solutions and calls evaluate (if there is a solution where
 * everything is blocked and locally balanced) then output otherwise
 * run reproduction and evaluate again
 * */
void GeneticSequentialSolver::solve(MaxFlowInstance &input, MaxFlowSolution &output){

  MaxFlowSolution solution;

  initialize(input.inputGraph.num_vertices, input.inputGraph.capacities);
  printSolutions(input.inputGraph.num_vertices);

  evaluate(input.inputGraph.num_vertices, input.sink, input.source);
  int solutionIndex = foundValidSolution(
      input.inputGraph.num_vertices,
      input.source,
      input.sink,
      input.inputGraph.capacities);


  int ctr = 0;
  //Todo: could place a cap on the number of generations here to speed up algo
  while (solutionIndex < 0){
    reproduction(input.inputGraph.num_vertices, input.inputGraph.capacities, input.source, input.sink);
    if (ctr == 99) printSolutions(input.inputGraph.num_vertices);
    ctr += 1;

    evaluate(input.inputGraph.num_vertices, input.sink, input.source);

    //printSolutions(input.inputGraph.num_vertices);

    solutionIndex = foundValidSolution(
        input.inputGraph.num_vertices,
        input.source,
        input.sink,
        input.inputGraph.capacities);
  }

  solution.flow = solutions[solutionIndex];
  solution.maxFlow = outputFlow[solutionIndex][input.source];
  output = solution;
}

//
//void GeneticSequentialSolver::printTotalFlow(int num_vertices) {
//  for (int p=0; p < currentPopulationSize; p++){
//    std::cout << "solution #" << p << "\n";
//    for (int i=0; i < num_vertices; i++){
//      std::cout << solutions[p][i][j] << " ";
//    }
//  }
//}

void GeneticSequentialSolver::printSolutions(int num_vertices) {
  for (int p=0; p < currentPopulationSize; p++){
    std::cout << "solution #" << p << "\n";
    for (int i=0; i < num_vertices; i++){
      for (int j=0; j < num_vertices; j++){
        std::cout << solutions[p][i][j] << " ";
      }
      std::cout << "\n";
    }
  }
}
