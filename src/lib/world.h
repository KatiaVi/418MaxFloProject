//
// Created by Ekaterina Villevald on 2019-11-02.
//

#ifndef INC_418MAXFLOPROJECT_MAXFLOWINSTANCE_H
#define INC_418MAXFLOPROJECT_MAXFLOWINSTANCE_H

#include <vector>
#include <bits/stdc++.h>
#include<utility>


struct Graph{
  int num_vertices;
  int num_edges;
  int** capacities; // num_vertices x num_vertices array where cij = flow capacity of edge i->j
};

struct GraphSmall{
  int num_vertices;
  int num_edges;
  std::vector<int> *edges;
  std::vector<std::pair<int,int>>* capacities;
  // adjacency list where capacities[i] stores (j,c) if there's an
  // edge from i to j with capacity c
};

class MaxFlowInstance{
 public:
  Graph inputGraph;
  int source; // id of source vertex
  int sink; // id of sink vertex
};

class MaxFlowSolution{
 public:
  int maxFlow;
  int** flow; // num_vertices x num_vertices array where cij = flow capacity of edge i->j
};

class MaxFlowInstanceSmall{
 public:
  GraphSmall inputGraph;
  int source; // id of source vertex
  int sink; // id of sink vertex
};

class MaxFlowSolutionSmall{
 public:
  int maxFlow;
  std::vector<std::pair<int,int>>* flow;
  // adjacency list where flows[i] stores (j,f) if there's an
  // edge from i to j with flow of f through it
};

#endif //INC_418MAXFLOPROJECT_MAXFLOWINSTANCE_H
