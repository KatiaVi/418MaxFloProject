//
// Created by Ekaterina Villevald on 2019-11-09.
//

#include "../../world.h"
#include "../../timing.h" 
#include <vector>

#ifndef INC_418MAXFLOPROJECT_DINICS_H
#define INC_418MAXFLOPROJECT_DINICS_H


class DinicsSequentialSolver{
 public:
  // look into shared, unique, reference stuff to make this more effecient
  void solve(MaxFlowInstance &input, MaxFlowSolution &output);

 private:
  int num_vertices;
  int *levels; // levels of a node
  float **flows; // flow solution
  float **capacities;
  std::vector<std::vector<int>> edges; // index = vertex and the element is an array of vertices that vertex is adjacent to
  Timer t; 

  void initialize(MaxFlowInstance &input);
  bool BFS(int source, int sink);
  float sendFlow(int currentVertex, int flow, int sink, int *start);

};


#endif //INC_418MAXFLOPROJECT_DINICS_H


