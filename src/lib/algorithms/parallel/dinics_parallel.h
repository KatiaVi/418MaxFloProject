//
// Created by Ekaterina Villevald on 2019-11-09.
//

#include "../../world.h"
#include "../../timing.h" 
#include <vector>


class DinicsParallelSolver {
 public:
  // look into shared, unique, reference stuff to make this more effecient
  void solve(MaxFlowInstance *input, MaxFlowSolution *output);

 private:
  int num_vertices;
  int *levels; // levels of a node
  int **flows; // flow solution
  int **capacities;
  std::vector<std::vector<int>> edges; // index = vertex and the element is an array of vertices that vertex is adjacent to
  Timer t; 

  void initialize(MaxFlowInstance *input);
  bool BFS(int source, int sink);
  int sendFlow(int currentVertex, int flow, int sink, int *start);
};

