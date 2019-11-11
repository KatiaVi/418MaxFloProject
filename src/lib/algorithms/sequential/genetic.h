#include "../../world.h"
#include <vector>

//
// Created by Ekaterina Villevald on 2019-11-02.
//

// Notes: the paper used as a reference has integer capacities for the edges but
// since we want to solve a more generic problem where edges have non-int capacities
// we adjusted some steps of the algo (marked with ADJUSTED)

// a vertex i is locally balanced if its excess flow = 0 (or inflow at i + outflow at i = 0)
// a vertex i is blocked if (1) forall outgoing edges flow_ij = capacity_ij OR (2) vertex j is blocked

#ifndef INC_418MAXFLOPROJECT_GENETIC_H
#define INC_418MAXFLOPROJECT_GENETIC_H

#define POPULATION_SIZE 100;

class GeneticSequentialSolver{
 public:
  // look into shared, unique, reference stuff to make this more effecient
  void solve(MaxFlowInstance &input, MaxFlowSolution &output);

  /* Debug Helper Functions */
  void printSolutions(int num_vertices);

 private:
  float k = 2.0f; // factor for the balance of a vertex in its energy level calculation (cross over step)
  float mutation_rate = 0.02f;
  float mutation_increment = 1.0f;
  int currentPopulationSize = POPULATION_SIZE;

  std::vector<int*> BALANCE; // BALANCE[p][i] is 1 if vertex i in individual p is locally balanced (defined above) and 0 otws
  std::vector<bool*> BLOCKED; // BLOCKED[p][i] is true if vertex i in individual p has blocked flow (defined above) and 0 otws

  std::vector<float**> solutions;
  std::vector<float*> excessFlow;
  std::vector<float*> inputFlow;
  std::vector<float*> outputFlow;
  std::vector<float*> totalFlow;

  void initialize(int num_vertices, float** &capacities);
  void evaluate(int num_vertices, int sink, int source);
  void reproduction(int num_vertices, float** capacities, int source, int sink);
  float** crossOver(int solution1, int solution2, int num_vertices, int source, int sink, float** capacities);
  float mutation(int individualIndex, int vertex, float energyLevel);

  /* Helper functions for the above functions */
  int foundValidSolution(int num_vertices, int source, int sink, float** &capacities);

  };

#endif //INC_418MAXFLOPROJECT_GENETIC_H
