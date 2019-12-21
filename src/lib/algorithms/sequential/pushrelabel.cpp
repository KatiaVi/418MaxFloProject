//
// Created by Wynne Yao on 2019-11-04.
//

#include <cstdlib>
#include <ctime>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <vector>
#include <queue> 
#include <limits.h> 

#include "pushrelabel.h"

using namespace std; 

void PushRelabelSequentialSolver::initialize(MaxFlowInstance *input) { 
  int numVertices = input->inputGraph.num_vertices;
  d = (int*)calloc(numVertices, sizeof(int)); 
  excessPerVertex = (int*)calloc(numVertices, sizeof(int)); 
  flows = new int*[numVertices]; 
  for (int i = 0; i < numVertices; i++) { 
    flows[i] = new int[numVertices]; 
  }
}

void PushRelabelSequentialSolver::preflow(MaxFlowInstance *input) { 
  initialize(input); 
  int numVertices = input->inputGraph.num_vertices; 

  for (int i = 0; i < numVertices; i++) { 
    d[i] = 0; 
    for (int j = 0; j < numVertices; j++) { 
      flows[i][j] = 0; // flow going from vertex i to vertex j is 0 
    }
  }

  d[input->source] = numVertices; 
  // find all vertices adjacent to s 
  int **cap = input->inputGraph.capacities; 

  for (int i = 0; i < numVertices; i++) { 
    if (cap[input->source][i] != 0 && (input->source != i)) { 
      // then it is an adjacent edge 
      flows[input->source][i] = cap[input->source][i]; 
      excessPerVertex[i] = flows[input->source][i];
      // add residual flow 
      flows[i][input->source] = -flows[input->source][i]; 
    }
  }
  for (int i = 0; i < numVertices; i++) { 
    if (i != input->source && i != input->sink && excessPerVertex[i] > 0) { 
      activeQueue.push(i); 
    }
  }
  
  // initialize active nodes to be all that have non zero excess
}

bool PushRelabelSequentialSolver::push(int numVertices, int **cap, int u, int sink) { 
  for (int v = 0; v < numVertices; v++) { 
        
    if (d[u] == d[v]+1 && (cap[u][v] - flows[u][v] > 0)) { // push if the height of the adjacent is smaller
      // push flow = min(remaining flow on edge, excess flow)
      int flow = min(cap[u][v] - flows[u][v], excessPerVertex[u]); 
      // reduce excess flow for overflowing vertex 
      excessPerVertex[u] -= flow; 

      if (excessPerVertex[u] > 0) { 
        activeQueue.push(v);
      }

      // increase excess flow for adjacent 
    
      excessPerVertex[v] += flow; 
      if (excessPerVertex[v] > 0) { 
         activeQueue.push(v);
      }

      // add residual flow 
      flows[u][v] += flow; 
      flows[v][u] -= flow; 
      return true; 
    }
  }
  return false; 
}

void PushRelabelSequentialSolver::relabel(int numVertices, int **cap, int u) { 
  int minHeight = INT_MAX; 
  for (int v = 0; v < numVertices; v++) { 
    if ((cap[u][v] - flows[u][v]) > 0) { 
      minHeight = min(minHeight, d[v]); 
    } 
  }
  d[u] = minHeight + 1;
  activeQueue.push(u); 
}



void PushRelabelSequentialSolver::pushRelabel(MaxFlowInstance *input, MaxFlowSolution *output) { 
  t.reset();
  preflow(input);  
  int numVertices = input->inputGraph.num_vertices;

  int **cap = input->inputGraph.capacities; 

  while (!activeQueue.empty()) {
    int u = activeQueue.front(); 
    activeQueue.pop(); 

    while ((u != input->source) && (u != input->sink) && (excessPerVertex[u] > 0)) { 
      if (!push(numVertices, cap, u, input->sink)) { 
        relabel(numVertices, cap, u); 
      } 
    }
  }
  double time = t.elapsed(); 
  printf("Push-Relabel time: %6fs\n", time); 

  output->maxFlow = excessPerVertex[input->sink]; 
  output->flow = flows; 
   
}
