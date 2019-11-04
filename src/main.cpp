//
// Created by Ekaterina Villevald on 2019-11-04.
//

#include <iostream>
#include "lib/world.h"
#include "lib/algorithms/sequential/genetic.h"


int main ( int argc, char * argv[] )
{
  // try calling genetic algo solve function on a sample graph
  std::cout << argc << argv[0] << "Hello!\n";

  float capacities[6][6];
  Graph testGraph1;

  testGraph1.num_vertices = 6;
  testGraph1.num_edges = 8;
  testGraph1.capacities = capacities;

  testGraph1.capacities[0][1] = 4.0f;
  testGraph1.capacities[0][2] = 2.0f;
  testGraph1.capacities[1][3] = 3.0f;
  testGraph1.capacities[2][3] = 2.0f;
  testGraph1.capacities[2][4] = 3.0f;
  testGraph1.capacities[3][2] = 1.0f;
  testGraph1.capacities[3][5] = 2.0f;
  testGraph1.capacities[4][5] = 4.0f;

  int source = 0;
  int sink = 5;

  MaxFlowInstance inputInstance;
  inputInstance.inputGraph = testGraph1;
  inputInstance.sink = sink;
  inputInstance.source = source;

  MaxFlowSolution gSolution;

  GeneticSequentialSolver gSolver;
  gSolution = gSolver.solve(inputInstance);

  // gSolution.maxFlow should be 5
  // http://www.cs.cmu.edu/afs/cs/academic/class/15451-s16/www/lectures/lec12-flow1.pdf
  // (using this example graph)

  return 0;
}