//
// Created by Ekaterina Villevald on 2019-12-09.
//

#include "../../world.h"
#include "../../timing.h"
#include <atomic>
#include <set>
#include <queue>
#include <unordered_set>
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_set.h"

using namespace std;

class PushRelabelParallelSolverV2{
 public:
  void pushRelabel(MaxFlowInstance *input, MaxFlowSolution *output);
 private:
  int** flows; // x coordinate is the source, y coordinate is the sink, for edges

  int *excessPerVertex; // excess flow on each vertex also exists
  atomic_int *d; // the labels
  atomic_int *addedExcess;
  atomic_bool *isDiscovered;
  atomic_bool *updated;
  //int **discoveredVertices;
  int *copyOfLabels;
  int *copyOfExcess;
  //int **residual;
  //vector<std::pair<int,int>> *residual;
  vector<vector<std::pair<int,int>>> residual;

  vector<int> *reverseResiduals;
  //vector<tbb::concurrent_vector<int>> reverseResiduals;
  vector<tbb::concurrent_vector<int>> discoveredVertices;
  int *work;
  //  int *outDegrees;

  int *active; // replace with a queue
  Timer t;
  unordered_set<int> workingSet;
  // int *workingSet;
  void initialize(MaxFlowInstance *input);
  void preflow(MaxFlowInstance *input);
  int existsActiveNode(MaxFlowInstance *input);
  void globalRelabel(int numVertices, int source, int sink);
};


