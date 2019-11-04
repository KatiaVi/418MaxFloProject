//
// Created by Ekaterina Villevald on 2019-11-02.
//

#include <cstdlib>
#include <ctime>
#include <numeric>
#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>

#include "genetic.h"


// Todo: (1) Malloc + free memory for all arrays
// Todo: (2) Address the other todos
// Todo: (3) Attempt to compile and run this

/* solve
 * initializes the solutions and calls evaluate (if there is a solution where
 * everything is blocked and locally balanced) then output otherwise
 * run reproduction and evaluate again
 * */
MaxFlowSolution *solve(MaxFlowInstance *input){

  initialize(input->inputGraph.num_vertices, input->inputGraph.capacities);
  evaluate(input->inputGraph.num_vertices, input->sink, input->source);
  solutionIndex = foundValidSolution(
      input->inputGraph.num_vertices,
      input->source,
      input->sink,
      input->inputGraph.capacities);

  //Todo: could place a cap on the number of generations here to speed up algo
  while (solutionIndex < 0){
    reproduction(input->inputGraph.num_vertices, input->inputGraph.capacities);
    evaluate(input->inputGraph.num_vertices, input->sink, input->source);
    solutionIndex = foundValidSolution(
        input->inputGraph.num_vertices,
        input->source,
        input->sink,
        input->inputGraph.capacities);
  }
  return solutions[solutionIndex];
}

/* foundValidSolution
 * checks whether there exists a solution where the vertices are locally balanced and all vertices are blocked
 * updates BALANCE and BLOCKED
 */
int foundValidSolution(int num_vertices, int source, int sink, flow &capacities){


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
      currentVertex = vertexQ.pop();
      BLOCKED[p][currentVertex] = true;

      for (int i=0; i < num_vertices; i++){
        if (capacities[currentVertex][i] > 0 and
        not (solutions[p][currentVertex][currentVertex] == capacities[i][currentVertex] or BLOCKED[p][currentVertex])){
            BLOCKED[p][currentVertex] = false;
        }
      }
      if (not BLOCKED[p][currentVertex]) break;
    }

    if (BLOCKED[p][source] and std::accumulate(BALANCE[p], BALANCE[p]+num_vertices) == num_vertices) return p;
  }
  return -1;
}


/* initialize
 * creates the solutions array and initializes the flow randomly so that each edge
 * has a flow that is between 0 and the capacity of that edge
 */
void initialize(int num_vertices, flow &capacities){
  srand (static_cast <unsigned> (time(0)));

  for (int p=0; p < POPULATION_SIZE; p++) {
    for (int i = 0; i < num_vertices; i++) {
      for (int j = 0; j < num_vertices; j++) {
        solutions[p][i][j] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/capacities[i][j]));
      }
      BLOCKED[p][i] = false;
      BALANCE[p][i] = 0;
    }
  }
}

/* evaluate
 * computes the excessFlow, inputFlow, outputFlow and totalFlow matrices for a generation
 * checks locally balanced condition and computes what flows are blocked if necessary
 * TODO: will we need to recalculate all of these for every individual for every generation
 * I'm thinking no since reproduction will only add 1 new individual, the others shouldn't be
 * changed? (think some more about this)
 */
void evaluate(int num_vertices, int sink, int source){
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

int countBalance(bool a, bool b){ return  }

/* reproduction
 * calculates the fitness scores (and probabilities)
 * for all solutions and chooses 2 solutions randomly with the calculated probabilities
 * do a crossover + mutate if necessary
 */
void reproduction(int num_vertices, flow capacities){
  float fitnessScores[currentPopulationSize];
  float fitnessSum = 0;

  for (int p=0; p < currentPopulationSize; p++){
    float sumOfFlow = 0;
    float sumOfCapacities = 0;
    for (int i=0; i < num_vertices; i++){
      sumOfFlow += std::accumulate(solutions[p][i], solutions[p][i] + num_vertices, 0);
      sumOfCapacities += std::accumulate(capacities[p][i], capacities[p][i] + num_vertices, 0);
    }
    fitnessScores[p] =
        (std::accumulate(BALANCE[p], BALANCE[p]+num_vertices,0))/(num_vertices - 2.0f) -
        math.abs(std::accumulate(excessFlow[p], excessFlow[p]+num_vertices,0))/sumOfFlow +
        sumOfFlow/sumOfCapacities;
    fitnessSum += fitnessScores[p];
  }

  boost::random::random_device rng;
  boost::random::discrete_distribution<> dist(fitnessScores);

  int individual1 = dist(rng);
  int individual2 = dist(rng);
  while (individual1 == individual2) individual2 = dist(rng);

  flow newSolution = crossOver(individual1, individual2);

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
flow crossOver(int solution1, int solution2, int num_vertices, int source, int sink, flow capacities){
  flow newSolution;
  char newEdgeAdded[num_vertices*num_vertices]; //Todo: maybe change to bitstring to make more efficient

  //TODO: paper says the source and sink vertices are not used during crossover. But what should they be set to then?
  newSolution[source] = 0.0f;
  newSolution[sink] = 0.0f;


  for (int i=0; i < num_vertices; i++){

    int inCapacity = 0;
    int outCapacity = 0;
    for (int j=0; j < num_vertices; j++){
      inCapacity += capacities[j][i];
      outCapacity += capacities[i][j];
    }

    energy1 = k*math.abs(excessFlow[solution1][i]) +
        math.abs(math.min(inCapacity, outCapacity) - math.max(inputFlow[solution1][i], -1.0f*outputFlow[solution1][i]));
    energy2 = k*math.abs(excessFlow[solution2][i]) +
        math.abs(math.min(inCapacity, outCapacity) - math.max(inputFlow[solution2][i], -1.0f*outputFlow[solution2][i]));

    int newSolutionIndex = (energy1 < energy2) ? solution1 : // energy level of solution1 is lower
                          ((energy2 < energy1) ? solution2 : // energy level of solution2 is lower
                          ((rand() > RAND_MAX/2) ? solution1 : solution2)) // randomly choose solution1 or solution2

    for (int j=0; j < num_vertices; j++){
      // add new edges (j->i and i->j) to the newSolution
      // possibly mutate the flow through the edge
      float newEdgeValueIJ = solutions[newSolutionIndex][i][j];
      float newEdgeValueJI = solutions[newSolutionIndex][j][i]

      float mutationAdjustment1 = mutate(i, probOfMutation);
      float mutationAdjustment2= mutate(i, probOfMutation);


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
float mutation(int individualIndex, int vertex, float energyLevel){
  float mutationProb = x + math.sqrt(math.abs(excessFlow[individualIndex][vertex])/ totalFlow[individualIndex][i]);

  boost::mt19937 gen1;
  float probabilityMutation[] = { mutationProb, 1-mutationProb };
  boost::random::discrete_distribution<> dist(probabilityMutation);
  int willMutate = dist(gen1);
  if (willMutate and energyLevel == 0){
    return ((rand() > RAND_MAX/2) ? mutation_increment : -1.0f*mutation_increment);
  }
  if (willMutate and energyLevel > 0){
    return -1.0f*mutation_increment;
    /* TODO: not sure about this, paper says vertex with an energy level above
     * zero mutation always occurs in the direction which will lower the energy level of the vertex
     * Why is this the case?
     */
  }
  return 0.0f;


}