//
// Created by Ekaterina Villevald on 2019-11-02.
//

#ifndef INC_418MAXFLOPROJECT_MAXFLOWINSTANCE_H
#define INC_418MAXFLOPROJECT_MAXFLOWINSTANCE_H


struct Graph{
  int num_vertices;
  int num_edges;
  float** capacities; // num_vertices x num_vertices array where cij = flow capacity of edge i->j
};

class MaxFlowInstance{
 public:
  Graph inputGraph;
  int source; // id of source vertex
  int sink; // id of sink vertex
};

class MaxFlowSolution{
 public:
  float maxFlow;
  float** flow; // num_vertices x num_vertices array where cij = flow capacity of edge i->j
};

#endif //INC_418MAXFLOPROJECT_MAXFLOWINSTANCE_H
