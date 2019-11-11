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
  active = (int*)calloc(numVertices, sizeof(int)); 
  excessPerVertex = (float*)calloc(numVertices, sizeof(float)); 
  flows = new float*[numVertices]; 
  for (int i = 0; i < numVertices; i++) { 
    flows[i] = new float[numVertices]; 
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
  float **cap = input->inputGraph.capacities; 

  for (int i = 0; i < numVertices; i++) { 
    if (cap[input->source][i] != 0 && (input->source != i)) { 
      // then it is an adjacent edge 
      flows[input->source][i] = cap[input->source][i]; 
      //@TODO: make sure i have floats for the flows 
      excessPerVertex[i] = flows[input->source][i];
      // add residual flow 
      flows[i][input->source] = -flows[input->source][i]; 
    }
  }
  
  // initialize active nodes to be all that have non zero excess
}

int PushRelabelSequentialSolver::existsActiveNode(MaxFlowInstance *input) { 
  int numVertices = input->inputGraph.num_vertices;
  for (int i = 0; i < numVertices; i++) { 
    if (excessPerVertex[i] > 0 && i != input->source && i != input->sink) { //@TODO: initialize active at some point?  
      return i; 
    }
  }
  return -1; 
}

bool PushRelabelSequentialSolver::push(int numVertices, float **cap, int u, int sink) { 
  
  for (int v = 0; v < numVertices; v++) { 
        
    if (d[u] > d[v] && (cap[u][v] - flows[u][v] > 0)) { // push if the height of the adjacent is smaller
      // push flow = min(remaining flow on edge, excess flow)
      float flow = min(cap[u][v] - flows[u][v], excessPerVertex[u]); //@TODO: bug: adding 0 here 
      // reduce excess flow for overflowing vertex 
      excessPerVertex[u] -= flow; 

      // increase excess flow for adjacent 
      excessPerVertex[v] += flow; 

      // add residual flow 
      flows[u][v] += flow; 
      flows[v][u] -= flow; 
      return true; 
    }
  }
    return false; 
}

void PushRelabelSequentialSolver::relabel(int numVertices, float **cap, int u) { 
  for (int i = 0; i < numVertices; i++) { 
    if (cap[u][i] - flows[u][i] != 0) { 
      assert(excessPerVertex[u] > 0 && d[u] <= d[i]); 
    }
  }
  
  int minHeight = INT_MAX; 
  for (int v = 0; v < numVertices; v++) { 
    if ((cap[u][v] - flows[u][v]) > 0) { 
      minHeight = min(minHeight, d[v]); 
    } 
  }
  d[u] = minHeight + 1;
}



void PushRelabelSequentialSolver::pushRelabel(MaxFlowInstance *input, MaxFlowSolution *output) { 
  preflow(input); 
  int numVertices = input->inputGraph.num_vertices;
   
  
  float **cap = input->inputGraph.capacities; 
  int u = existsActiveNode(input); 

  while (u != -1) {
    if ((u != input->source) && (u != input->sink) && (excessPerVertex[u] > 0)) { 
      //printf("u: %d\n", u);
      if (!push(numVertices, cap, u, input->sink)) { 
        //printf("relabel\n"); 
        relabel(numVertices, cap, u); 
      } 

      //testing code 
      /* for (int i = 0; i < numVertices; i++) { 
        for (int j = 0; j < numVertices; j++) { 
          if (flows[i][j] != 0) { 
            printf("flows[%d][%d]: %f\n", i, j, flows[i][j]); 
          }
        }
      }
      for (int i = 0; i < numVertices; i++) { 
        printf("d[%d]: %d\n", i, d[i]); 
      }
      for (int i = 0; i < numVertices; i++) { 
        printf("excessPerVertex[%d]: %f\n", i, excessPerVertex[i]); 
      } */ 
    }
    u = existsActiveNode(input); 
  }

  output->maxFlow = excessPerVertex[input->sink]; 
  output->flow = flows; 
   
}
