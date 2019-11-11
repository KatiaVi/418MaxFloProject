#include <queue>
#include <cmath>
#include <iostream>
#include "dinics.h"

/* Note:
 * The code and implementation was taken from https://www.geeksforgeeks.org/dinics-algorithm-maximum-flow/
 * it was modified to fit our data structures for the maxflow input instance and the maxflow solution
 */

void DinicsSequentialSolver::initialize(MaxFlowInstance &input){
  num_vertices = input.inputGraph.num_vertices;
  levels = new int[num_vertices];
  flows = new float*[num_vertices];
  capacities = input.inputGraph.capacities;


  for (int i = 0; i < num_vertices; i++) {
    flows[i] = new float[num_vertices];
    std::vector<int> adj;
    for (int j = 0; j < num_vertices; j++) {
      flows[i][j] = 0;
      if (capacities[i][j] > 0){
        adj.push_back(j);
      }
    }
    edges.push_back(adj);
  }
}

float DinicsSequentialSolver::sendFlow(int current, int flow, int sink, int start[])
{
  // Sink reached
  if (current == sink)
    return flow;

  // Traverse all adjacent edges one -by - one.
  for (  ; start[current] < edges[current].size(); start[current]++) {
    // Pick next edge from adjacency list of current
    int child = edges[current][start[current]];

    if (levels[child] == levels[current]+1 && flows[current][child] < capacities[current][child])
    {
      // find minimum flow from current to sink
      float curr_flow = fmin(flow, capacities[current][child] - flows[current][child]);

      float temp_flow = sendFlow(child, curr_flow, sink, start);

      // flow is greater than zero
      if (temp_flow > 0)
      {
        // add flow  to current edge
        flows[current][child] += temp_flow;

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
bool DinicsSequentialSolver::BFS(int source, int sink)
{
  for (int i = 0 ; i < num_vertices ; i++)
    levels[i] = -1;

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

    for (int i = 0; i < edges[current].size(); i++)
    {
      int child = edges[current][i];
      if (levels[child] < 0  && flows[current][child] < capacities[current][child])
      {
        // Level of current vertex is,
        // level of parent + 1
        levels[child] = levels[current] + 1;

        vertexQ.push(child);
      }
    }
  }

  // IF we can not reach to the sink we
  // return false else true
  return levels[sink] < 0 ? false : true;
}

void DinicsSequentialSolver::solve(MaxFlowInstance &input, MaxFlowSolution &output){

  initialize(input);

  int totalFlow = 0;

  // Augment the flow while there is path
  // from source to sink
  while (BFS(input.source, input.sink))
  {
    // store how many edges are visited
    // from V { 0 to V }
    int *start = new int[input.inputGraph.num_vertices + 1];
    int flow = sendFlow(input.source, INT_MAX, input.sink, start);

    // while flow is not zero in graph from S to D
    while (flow > 0){
      // Add path flow to overall flow
      totalFlow += flow;
      flow = sendFlow(input.source, INT_MAX, input.sink, start);
    }
  }

  output.maxFlow = totalFlow;
  output.flow = flows;
}