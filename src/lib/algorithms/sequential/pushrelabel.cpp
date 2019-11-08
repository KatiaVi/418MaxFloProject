//
// Created by Wynne Yao on 2019-11-04.
//

// @TODO: look at the FIFO version! 

#include <cstdlib>
#include <ctime>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <vector>

#include "pushrelabel.h"
#include "./world.h" 

using namespace std; 


MaxFlowSolution *solve(MaxFlowInstance *input){
  preflow(input); 
} 

void initialize(MaxFlowInstance *input) { 
  int numVertices = input->inputGraph.num_vertices;
  for (int i = 0; i < numVertices; i++) { 
    d[i] = 0; 
    active[i] = 0;
    excessPerVertex[i] = 0;  
    for (int j = 0; j < numVertices; j++) { 
      flows[i][j] = 0;
    } 
  }
}

// check that my preflow is correct 
void preflow(MaxFlowInstance *input) { 
  initialize(input); 

  int numVertices = input->inputGraph.num_vertices; 

  d[input->source] = numVertices; 

  for (int i = 0; i < numVertices; i++) { 
    d[i] = 0; 
    for (int j = 0; j < numVertices; j++) { 
      flows[i][j] = 0; // flow going from vertex i to vertex j is 0 
    }
  }

  // find all vertices adjacent to s 
  float **cap = input->inputGraph.capacities; 
  for (int i = 0; i < numVertices; i++) { 
    if (cap[input->source][i] != 0) { 
      // then it is an adjacent edge 
      flows[input->source][i] = cap[input->source][i]; 
      excessPerVertex[i] = flows[input->source][i];
      // add residual flow 
      flows[i][input->source] = -flows[input->source][i]; 
    }
  }
  // initialize active nodes to be all that have non zero excess
  
}

int existsActiveNode(MaxFlowInstance *input) { 
  int numVertices = input->inputGraph.num_vertices;
  for (int i = 0; i < numVertices; i++) { 
    if (excessPerVertex[i] > 0 && i != input->source && i != input->sink) { //@TODO: initialize active at some point? 
      return i; 
    }
  }
  return -1; 
}

bool isAdmissible(int u, int j) { 
  return (d[u] == d[j]+1); 
}

int findMinLabel(vector<int> outgoingEdges) { 
  int min = numeric_limits<int>::max(); 
  for (int i = 0; i < outgoingEdges.size(); i++) { 
    if (d[outgoingEdges[i]] < min) { 
      min = d[outgoingEdges[i]]; 
    }
  }
  return min; 
}

bool push(int numVertices, float **cap, int u) { 
  for (int v = 0; v < numVertices; v++) { 
      if (flows[u][v] != 0) { // it's a valid edge 
        if (flows[u][v] == cap[u][v]) { 
          continue; 
        }
        if (d[u] > d[v]) { // push if the height of the adjacent is smaller 
          // push flow = min(remaining flow on edge, excess flow)
          int flow = min(cap[u][v] - flows[u][v], excessPerVertex[u]); 
          
          // reduce excess flow for overflowing vertex 
          excessPerVertex[u] -= flow; 

          // increase excess flow for adjacent 
          excessPerVertex[v] += flow; 

          // add residual flow 
          flows[u][v] += flow; 
          flows[v][u] -= flow; //@TODO: i think this is what they're doing 
          return true; 
        }
      }
    }
    return false; 
}

void relabel(int numVertices, float **cap, int u) { 
  float minHeight = numeric_limits<int>::max(); 
  for (int i = 0; i < numVertices; i++) { 
    if (flows[u][i] != 0) { 
      if (flows[u][i] == cap[u][i]) { 
        continue; 
      }

      // update minimum height
      if (d[excessPerVertex[i]] < minHeight) { 
        minHeight = d[excessPerVertex[i]]; 
      }
    } 
  }
  d[u] = minHeight + 1; //@TODO: might be shifted in? 
}

void pushRelabel(MaxFlowInstance *input) { 
  preflow(input); 

  int numVertices = input->inputGraph.num_vertices;
  int u = existsActiveNode(input); 
  float **cap = input->inputGraph.capacities; 
  
  while (u != -1) { 
    // look at all outgoing edges of u 
    // vector<int> outgoingAdmissibleEdges; 
    // vector<int> outgoingEdges;  
    if (!push(numVertices, cap, u)) { 
      relabel(numVertices, cap, u); 
    } 
    /*for (int j = 0; j < numVertices; j++) { 
      if (flows[u][j] != 0 && isAdmissible(u, j)) { 
        outgoingAdmissibleEdges.push_back(j);
        // if flow is equal to capacity then can't push?
        if (flows[u][j] == cap[u][j]) { 
          continue; 
        }
        if (flows[u][j] < minOutFlow) { // find the min out flow (for separate step)
          minOutFlow = flows[u][j]; 
        }
      } else if (flows[u][j] != 0) { 
        outgoingEdges.push_back(j); 
      }
    } 
    if (outgoingAdmissibleEdges.size() > 0) { // do the push 
      float delta = min(minOutFlow, excessPerVertex[u]); 
      excessPerVertex[u] -= delta; 
      for (int i = 0; i < outgoingAdmissibleEdges.size(); i++) { 
        excessPerVertex[outgoingAdmissibleEdges[i]] += delta; 
      }
    } else { // do the relabel - @TODO: i dont think this is right 
      d[u] = 1+findMinLabel(outgoingEdges); 
    }*/ 
  }
}

