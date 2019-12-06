#include <queue>
#include <cmath>
#include <iostream>
#include <limits.h>
#include <omp.h>

#include "dinics_parallel.h"


/* Note:
 * The code and implementation was taken from https://www.geeksforgeeks.org/dinics-algorithm-maximum-flow/
 * it was modified to fit our data structures for the maxflow input instance and the maxflow solution
 */

void DinicsParallelSolver::smallInitialize(MaxFlowInstanceSmall *input){
  num_vertices = input->inputGraph.num_vertices;
  levels = (int*)calloc(num_vertices, sizeof(int));
  flowsSmall = (std::vector<std::pair<int,int>>*)calloc(num_vertices, sizeof(std::vector<std::pair<int,int>>));
  capacitiesSmall = input->inputGraph.capacities;
  edges = input->inputGraph.edges;



  #pragma omp parallel for
  for (int i=0; i < num_vertices; i++){
    for (int j=0; j < capacitiesSmall[i].size(); j++){
      std::pair<int,int> initialFlow = std::make_pair(std::get<0>(capacitiesSmall[i][j]), 0);
      std::pair<int,int> initialFlowReverse = std::make_pair(i, 0);

      flowsSmall[i].push_back(initialFlow);
      flowsSmall[std::get<0>(capacitiesSmall[i][j])].push_back(initialFlowReverse);
    }
  }

  printf("Capacities:\n");
  for (int i=0; i < num_vertices; i++){
    for (int j=0; j < capacitiesSmall[i].size(); j++){
      printf("%d: (%d, %d)\n", i, std::get<0>(capacitiesSmall[i][j]), std::get<1>(capacitiesSmall[i][j]));
    }
    std::cout << "\n";
  }

  printf("Flows:\n");
  for (int i=0; i < num_vertices; i++){
    for (int j=0; j < flowsSmall[i].size(); j++){
      printf("%d: (%d, %d)\n", i, std::get<0>(flowsSmall[i][j]), std::get<1>(flowsSmall[i][j]));
    }
    std::cout << "\n";
  }

  printf("Edges:\n");
  for (int i=0; i < num_vertices; i++){
    for (int j=0; j < edges[i].size(); j++){
      printf("%d: %d\n", i, edges[i][j]);
    }
    std::cout << "\n";
  }
}

bool DinicsParallelSolver::smallBFS(int source, int sink){
  std::cout << "in smallBFS\n";
  BFSTimer.reset();
  #pragma parallel for private(num_vertices)
  for (int i = 0 ; i < num_vertices ; i++)
    levels[i] = -1;
  #pragma omp barrier

  levels[source] = 0;  // Level of source vertex

  // Create a queue, enqueue source vertex
  // and mark source vertex as visited here
  // level[] array works as visited array also.
  std::queue<int> vertexQ;
  vertexQ.push(source);

  while (!vertexQ.empty())
  {
    int current = vertexQ.front();
    vertexQ.pop();

    //TODO: Does initializing the reusable accesses to their own variables actually help with performance?
    std::vector<int> neighbors = edges[current];
    std::vector<std::pair<int,int>> flowNeighbors = flowsSmall[current];
    std::vector<std::pair<int,int>> capacitiesNeighbors = capacitiesSmall[current];
    int currentLevel = levels[current];


    #pragma parallel for firstprivate(flowNeighbors, capacitiesNeighbors)
    for (int i = 0; i < neighbors.size(); i++)
    {
      int child = neighbors[i];
      std::vector<std::pair<int,int>>::iterator flowToChild = std::find_if(flowNeighbors.begin(), flowNeighbors.end(),
                                     [child](std::pair<int,int> tup){ return (std::get<0>(tup) == child); });
//
//      std::cout << "(" <<
//        std::get<0>((*flowToChild)) <<
//            ", " << std::get<1>((*flowToChild)) << ")\n";
//
//      std::cout << "(" << std::get<0>(capacitiesNeighbors[i])  <<
//          ", " << std::get<1>(capacitiesNeighbors[i]) << ")\n";

      if (levels[child] < 0  && std::get<1>((*flowToChild)) < std::get<1>(capacitiesNeighbors[i]))
      {
        // Level of current vertex is,
        // level of parent + 1
        levels[child] = currentLevel + 1;
        #pragma omp critical
        {
          vertexQ.push(child);
        }
      }
    }
  }
  // IF we can not reach to the sink we
  // return false else true
  BFStime += BFSTimer.elapsed();
  std::cout << "finished smallBFS\n";
  return levels[sink] < 0 ? false : true;
}
int DinicsParallelSolver::smallSendFlow(int current, int flow, int sink, int *start){
//  std::cout << "in smallSendFlow with current: " << current << "\n";
// Sink reached
  if (current == sink) {
    return flow;
  }

  std::vector<std::pair<int,int>> tempCapacities = capacitiesSmall[current];
//  std::cout << "got capacities\n";

  std::vector<std::pair<int,int>> tempFlows = flowsSmall[current];
  std::vector<int> tempEdges = edges[current];
  int currentLevelPlusOne = levels[current] + 1;
  // Traverse all adjacent edges one -by - one.
//  std::cout << "got intermediate structures\n";

  for (  ; start[current] < num_vertices; start[current]++) {
    // Pick next edge from adjacency list of current
    __builtin_prefetch((const int*)(&tempEdges[start[current]+3]),0,0);
    printf("current = %d start[%d] = %d\n", current, current, start[current]);

//    printf("current = %d start[%d] = %d\n", current, start[current]);
    int i = start[current];
    std::vector<int>::iterator childIter = (std::find_if(tempEdges.begin(), tempEdges.end(),
        [i](int e){ return (e == i); }));
    if (childIter != tempEdges.end()) {
      int child = *childIter;
      printf("current = %d start[%d] = %d and child = %d\n", current, current, i, child);

      std::vector < std::pair < int, int >> ::iterator
      flowToChild = std::find_if(tempFlows.begin(), tempFlows.end(),
                                 [child](std::pair<int, int> tup) { return (std::get<0>(tup) == child); });

      std::vector < std::pair < int, int >> ::iterator
      flowFromChild = std::find_if(flowsSmall[child].begin(), flowsSmall[child].end(),
                                   [current](std::pair<int, int> tup) { return (std::get<0>(tup) == current); });

      std::vector < std::pair < int, int >> ::iterator
      capacityOnEdge = std::find_if(tempCapacities.begin(), tempCapacities.end(),
                                    [i](std::pair<int, int> tup) { return (std::get<0>(tup) == i); });

//      std::cout << "(" <<
//                std::get<0>((*flowToChild)) <<
//                ", " << std::get<1>((*flowToChild)) << ")\n";
//
//      std::cout << "(" <<
//                std::get<0>((*flowFromChild)) <<
//                ", " << std::get<1>((*flowFromChild)) << ")\n";
//
//      std::cout << "(" <<
//                std::get<0>((*capacityOnEdge)) <<
//                ", " << std::get<1>((*capacityOnEdge)) << ")\n";

      if (levels[child] == currentLevelPlusOne && std::get<1>(*flowToChild) < std::get<1>(*capacityOnEdge)) {
        // find minimum flow from current to sink
        int curr_flow = std::min(flow, std::get<1>(*capacityOnEdge) - std::get<1>(*flowToChild));
        int temp_flow = smallSendFlow(child, curr_flow, sink, start);
        printf("GOT temp_flow: %d\n", temp_flow);

        // flow is greater than zero
        if (temp_flow > 0) {
          std::vector<std::pair<int, int>>
              newFlows = {std::make_pair(std::get<0>(*flowToChild), std::get<1>(*flowToChild) + temp_flow),
                          std::make_pair(current, std::get<1>(*flowFromChild) - temp_flow)};

//          std::cout << "(" <<
//                    std::get<0>((*flowToChild)) <<
//                    ", " << std::get<1>((*flowToChild)) << ")\n";
//
//          std::cout << "(" <<
//                    std::get<0>((*flowFromChild)) <<
//                    ", " << std::get<1>((*flowFromChild)) << ")\n";

          std::iter_swap(flowToChild, newFlows.begin()); // does this change it in the actual flowsSmall structure
          std::iter_swap(flowFromChild, newFlows.end());

          return temp_flow;
        }
      }
    }
  }
  return 0;
}

void DinicsParallelSolver::smallSolve(MaxFlowInstanceSmall *input, MaxFlowSolutionSmall *output){
  std::cout << "in smallSolve\n";
  omp_set_num_threads(omp_get_max_threads());
  t.reset();
  smallInitialize(input);
  int totalFlow = 0;

  // Augment the flow while there is path
  // from source to sink
  while (smallBFS(input->source, input->sink))
  {
    // store how many edges are visited
    // from V { 0 to V }
    int *start = (int*)calloc(input->inputGraph.num_vertices + 1, sizeof(int)); // reset to all 0s instead

    sendFlowTimer.reset();
    int flow = smallSendFlow(input->source, INT_MAX, input->sink, start);
    sendFlowTime += sendFlowTimer.elapsed();


    printf("Flows:\n");
    for (int i=0; i < num_vertices; i++){
      for (int j=0; j < flowsSmall[i].size(); j++){
        printf("%d: (%d, %d)\n", i, std::get<0>(flowsSmall[i][j]), std::get<1>(flowsSmall[i][j]));
      }
      std::cout << "\n";
    }

    // while flow is not zero in graph from S to D
    while (flow > 0){
      // Add path flow to overall flow
      totalFlow += flow;

      sendFlowTimer.reset();
      flow = smallSendFlow(input->source, INT_MAX, input->sink, start);
      sendFlowTime += sendFlowTimer.elapsed();
    }
  }
  double time = t.elapsed();
  printf("Dinics Time: %6fs\n", time);
  printf("Time spent in BFS: %6fs\n", BFStime);
  printf("Time spent in sendFlow: %6fs\n", sendFlowTime);

  output->maxFlow = totalFlow;
  output->flow = flowsSmall;
}


void DinicsParallelSolver::initialize(MaxFlowInstance *input){
  num_vertices = input->inputGraph.num_vertices;
  levels = (int*)calloc(num_vertices, sizeof(int));
  flows = (int**)calloc(num_vertices, sizeof(int*));
  capacities = input->inputGraph.capacities;
  edges = (std::vector<int>*)calloc(num_vertices, sizeof(std::vector<int>));

  #pragma omp parallel for
  for (int i = 0; i < num_vertices; i++) {
    flows[i] = (int*)calloc(num_vertices, sizeof(int));
    //int* tempCapacities = capacities[i]; //Todo: Does this actually help?
    for (int j = 0; j < num_vertices; j++) {
      //__builtin_prefetch((const int*)(&capacities[i][j+250]),0,0);
      if(/*tempCapacities[j]*/ capacities[i][j]) {
          edges[i].push_back(j);
        }
    }
  }
}

int DinicsParallelSolver::sendFlow(int current, int flow, int sink, int start[])
{
  // Sink reached
  if (current == sink) {
    return flow;
  }

  int* tempCapacities = capacities[current];
  int* tempFlows = flows[current];
  std::vector<int> tempEdges = edges[current];
  int currentLevelPlusOne = levels[current] + 1;
  // Traverse all adjacent edges one -by - one.

  for (  ; start[current] < tempEdges.size(); start[current]++) {
    // Pick next edge from adjacency list of current
    __builtin_prefetch((const int*)(&tempEdges[start[current]+3]),0,0);

    int child = tempEdges[start[current]];
    if (levels[child] == currentLevelPlusOne && tempFlows[child] < tempCapacities[child])
    {
      // find minimum flow from current to sink
      int curr_flow = std::min(flow, tempCapacities[child] - tempFlows[child]);
      int temp_flow = sendFlow(child, curr_flow, sink, start);

      // flow is greater than zero
      if (temp_flow > 0)
      {
        // add flow  to current edge
        tempFlows[child] += temp_flow;

        // subtract flow from reverse edge
        // of current edge
        flows[child][current] -= temp_flow;
        return temp_flow;
      }
    }
  }
  return 0;
}


// Finds if more flow can be sent from s to t.
// Also assigns levels to nodes.
bool DinicsParallelSolver::BFS(int source, int sink)
{
  BFSTimer.reset();
  #pragma parallel for private(num_vertices)
  for (int i = 0 ; i < num_vertices ; i++)
    levels[i] = -1;
  #pragma omp barrier

  levels[source] = 0;  // Level of source vertex

  // Create a queue, enqueue source vertex
  // and mark source vertex as visited here
  // level[] array works as visited array also.
  std::queue<int> vertexQ;
  vertexQ.push(source);

  while (!vertexQ.empty())
  {
    int current = vertexQ.front();
    vertexQ.pop();

    //TODO: Does initializing the reusable accesses to their own variables actually help with performance?
    std::vector<int> neighbors = edges[current];
    int* flowNeighbors = flows[current];
    int* capacitiesNeighbors = capacities[current];
    int currentLevel = levels[current];

    #pragma parallel for firstprivate(flowNeighbors, capacitiesNeighbors)
    for (int i = 0; i < neighbors.size(); i++)
    {
      int child = neighbors[i];
      
      if (levels[child] < 0  && flowNeighbors[child] < capacitiesNeighbors[child])
      {
        // Level of current vertex is,
        // level of parent + 1
        levels[child] = currentLevel + 1;
        #pragma omp critical
        {
          vertexQ.push(child);
        }
      }
    }
  }
  // IF we can not reach to the sink we
  // return false else true
  BFStime += BFSTimer.elapsed();
  return levels[sink] < 0 ? false : true;
}

void DinicsParallelSolver::processLevel(std::vector<int> &oldVertexQ, std::vector<int> &newVertexQ){

    if (oldVertexQ.size() > BAGSPLIT_CUTOFF) {
      #pragma omp parallel
      #pragma omp single nowait
      {
        const int half_size = oldVertexQ.size() / 2;
        std::vector<int> A(oldVertexQ.begin(), oldVertexQ.begin() + half_size);
        std::vector<int> B(oldVertexQ.begin() + half_size, oldVertexQ.end());
        std::vector<int> newVertexQB;

        #pragma omp task shared(newVertexQB) untied
        {
          processLevel(B, newVertexQB);
        }
        processLevel(A, newVertexQ);

        #pragma omp taskwait
        newVertexQ.insert(newVertexQ.end(), newVertexQB.begin(), newVertexQB.end());
      }
    }
    else {
    for (int i = 0; i < oldVertexQ.size(); i++){
      int v = oldVertexQ[i];
      std::vector<int> tempVec = edges[v];
      int tempVecSize = tempVec.size();
      int* flowNeighbors = flows[v];
      int* capacitiesNeighbors = capacities[v];

      #pragma omp parallel for shared(newVertexQ) firstprivate(tempVec, flowNeighbors, capacitiesNeighbors)
      for (int j=0; j < tempVecSize; j++){
        int w = tempVec[j];
        if (levels[w] < 0 && flowNeighbors[w] < capacitiesNeighbors[w]) {
          levels[w] = currentLevel + 1;

          #pragma omp critical
          {
            newVertexQ.push_back(w);
          }
        }
      }
    }
  }
}

// Finds if more flow can be sent from s to t.
// Also assigns levels to nodes.
bool DinicsParallelSolver::parallelBFS(int source, int sink)
{
  BFSTimer.reset();

  #pragma omp parallel for
  for (int i = 0 ; i < num_vertices ; i++)
    levels[i] = -1;
  #pragma omp barrier


  levels[source] = 0;
  currentLevel = 0;
  std::vector<int> vertexQ;
  vertexQ.push_back(source);

  while (vertexQ.size() > 0){
    std::vector<int> vertexQNextLevel;
    processLevel(vertexQ, vertexQNextLevel);
    currentLevel += 1;
    vertexQ = vertexQNextLevel;
  }

  BFStime += BFSTimer.elapsed();
  return levels[sink] < 0 ? false : true;
}

void DinicsParallelSolver::solve(MaxFlowInstance *input, MaxFlowSolution *output){
  omp_set_num_threads(omp_get_max_threads());
  t.reset();
  initialize(input);
  int totalFlow = 0;

  // Augment the flow while there is path
  // from source to sink
  while (BFS(input->source, input->sink))
  {
    // store how many edges are visited
    // from V { 0 to V }
    int *start = (int*)calloc(input->inputGraph.num_vertices + 1, sizeof(int));

    sendFlowTimer.reset();
    int flow = sendFlow(input->source, INT_MAX, input->sink, start);
    sendFlowTime += sendFlowTimer.elapsed();

    // while flow is not zero in graph from S to D
    while (flow > 0){
      // Add path flow to overall flow
      totalFlow += flow;

      sendFlowTimer.reset();
      flow = sendFlow(input->source, INT_MAX, input->sink, start);
      sendFlowTime += sendFlowTimer.elapsed();
    }
  }
  double time = t.elapsed();
  printf("Dinics Time: %6fs\n", time);
  printf("Time spent in BFS: %6fs\n", BFStime);
  printf("Time spent in sendFlow: %6fs\n", sendFlowTime);

  output->maxFlow = totalFlow;
  output->flow = flows;
}