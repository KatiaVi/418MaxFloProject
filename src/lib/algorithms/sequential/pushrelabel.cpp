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

// check that my preflow is correct 
void preflow(MaxFlowInstance *input) { 
   
  Graph g = input->inputGraph; 
  int numVertices = input->inputGraph.num_vertices; 
  for (int i = 0; i < numVertices; i++) {  
    heights[i] = 0; 
  }

  heights[input->source] = numVertices; 

  for (int i = 0; i < numVertices; i++) { 
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
      excessPerVertex[i] = cap[input->source][i]; 
    }
  }
  // initialize active nodes to be all that have non zero excess
  
}

bool existsActiveNode(MaxFlowInstance *input) { 
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

void pushRelabel(MaxFlowInstance *input) { 
  int numVertices = input->inputGraph.num_vertices;
  int u = existsActiveNode(input); 
  float minOutFlow = numeric_limits<int>::max(); 
  while (u != -1) { 
    // look at all outgoing edges of u 
    vector<int> outgoingAdmissibleEdges; 
    vector<int> outgoingEdges;  
    for (int j = 0; j < numVertices; j++) { 
      if (flows[u][j] != 0 && isAdmissible(u, j)) { 
        outgoingAdmissibleEdges.push_back(j);
        if (flows[u][j] < minOutFlow) {
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
    }
  }
}

