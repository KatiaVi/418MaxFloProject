#include "../../world.h"

//
// Created by Ekaterina Villevald on 2019-11-02.
//

// Notes: the paper used as a reference has integer capacities for the edges but
// since we want to solve a more generic problem where edges have non-int capacities
// we adjusted some steps of the algo (marked with ADJUSTED)

#ifndef INC_418MAXFLOPROJECT_GENETIC_H
#define INC_418MAXFLOPROJECT_GENETIC_H

#define POPULATION_SIZE 10;
typedef float** flow:

class GeneticSequentialSolver{
 public:
  // look into shared, unique, reference stuff to make this more effecient
  MaxFlowSolution *solve(MaxFlowInstance &input);
 private:
  float k = 2.0f; // factor for the balance of a vertex in its energy level calculation (cross over step)
  float mutation_rate = 0.02f;
  float truncation_rate = 0.2f;

  bool BALANCE[POPULATION_SIZE]; // BALANCE[k] is true if all vertices in individual k are locally balanced (defined above)
  bool* BLOCKED[POPULATION_SIZE]; // BLOCKED[k][i] is true if vertex i in individual k has blocked flow (defined above)

  flow solutions[POPULATION_SIZE];
  flow excessFlow[POPULATION_SIZE];
  flow inputFlow[POPULATION_SIZE];
  flow outputFlow[POPULATION_SIZE];
  flow totalFlow[POPULATION_SIZE];

  void initialize(int num_vertices, flow &capacities);
  void evaluate();
  void reproduction(int numBalancedVertices, flow* excessFlow);
  float** crossOver(flow &solution1, flow &solution2);
  float** mutation(flow &excessFlow, float totalFlow); // Note to self: assimilation check should happen in this function too?
};

#endif //INC_418MAXFLOPROJECT_GENETIC_H
