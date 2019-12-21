//
// Created by Ekaterina Villevald on 2019-11-09.
//

#include "../../world.h"
#include "../../timing.h"
#include "tbb/concurrent_vector.h"
#include <vector>
#include <utility>

#define BAGSPLIT_CUTOFF 128

class DinicsParallelSolver {
 public:
  void solve(MaxFlowInstance *input, MaxFlowSolution *output);
  void smallSolve(MaxFlowInstanceSmall *input, MaxFlowSolutionSmall *output);


 private:
  int num_vertices;
  int *levels; // levels of a node
  int **flows; // flow solution
  std::vector<std::pair<int,int>>* flowsSmall;
  int **capacities;
  std::vector<std::pair<int,int>>* capacitiesSmall; // index = vertex and the element is an array of vertices that vertex is adjacent to
  std::vector<int>* edges; // index = vertex and the element is an array of vertices that vertex is adjacent to
  int currentLevel = 0;

  Timer t;
  Timer BFSTimer;
  Timer sendFlowTimer;
  double BFStime = 0;
  double sendFlowTime = 0;

  void initialize(MaxFlowInstance *input);
  bool BFS(int source, int sink);
  bool parallelBFS(int source, int sink);
  void processLevel(std::vector<int> &oldVertexQ, std::vector<int> &newVertexQ);
  int sendFlow(int currentVertex, int flow, int sink, int *start);

  void smallInitialize(MaxFlowInstanceSmall *input);
  bool smallBFS(int source, int sink);
  int smallSendFlow(int currentVertex, int flow, int sink, int *start);
  bool smallParallelBFS(int source, int sink);
  void smallProcessLevel(std::vector<int> &oldVertexQ, std::vector<int> &newVertexQ);

};

