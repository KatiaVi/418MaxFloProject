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
#include <atomic> 

#include "pushrelabel_parallel.h"

using namespace std; 

void PushRelabelParallelSolver::initialize(MaxFlowInstance *input) {
  int numVertices = input->inputGraph.num_vertices;
  d = (int*)calloc(numVertices, sizeof(int)); 
  active = (int*)calloc(numVertices, sizeof(int)); 
  excessPerVertex = (float*)calloc(numVertices, sizeof(float)); 
  copyOfExcess = (float*)calloc(numVertices, sizeof(float)); 
  addedPerVertex = (float*)calloc(numVertices, sizeof(float));
  isDiscovered = (atomic_flag*)malloc(numVertices*sizeof(atomic_flag));
  workingSet = (int*)calloc(numVertices, sizeof(int)); 
  discoveredVertices = new int*[numVertices]; 
  copyOfLabels = (int*)calloc(numVertices, sizeof(int)); 
  
  flows = new float*[numVertices]; 

//  #pragma omp parallel for (local copy of isDiscovered[i])
  for (int i = 0; i < numVertices; i++) { 
    flows[i] = new float[numVertices]; 
    discoveredVertices[i] = new int[numVertices]; 
    workingSet[i] = -1; // initialize all to -1 so that only the active ones get added
    isDiscovered[i].clear();
  }
}

void PushRelabelParallelSolver::preflow(MaxFlowInstance *input) {
  initialize(input); 
  int numVertices = input->inputGraph.num_vertices; 

  #pragma omp parallel for
  for (int i = 0; i < numVertices; i++) { 
    d[i] = 0; 
    for (int j = 0; j < numVertices; j++) { 
      flows[i][j] = 0; // flow going from vertex i to vertex j is 0
      discoveredVertices[i][j] = -1;  
    }
  }

  d[input->source] = numVertices; 
  // find all vertices adjacent to s 
  float **cap = input->inputGraph.capacities; 

//  #pragma omp parallel for
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
//  #pragma omp parallel for
  for (int i = 0; i < numVertices; i++) { 
    if (i != input->source && i != input->sink && excessPerVertex[i] > 0) { 
      workingSet[i] = i; 
    }
  }
  // initialize active nodes to be all that have non zero excess
}

int PushRelabelParallelSolver::existsActiveNode(MaxFlowInstance *input) {
  int numVertices = input->inputGraph.num_vertices;
  for (int i = 0; i < numVertices; i++) { 
    if (excessPerVertex[i] > 0 && i != input->source && i != input->sink) { //@TODO: initialize active at some point?  
      return i; 
    }
  }
  return -1; 
}

bool PushRelabelParallelSolver::push(int numVertices, float **cap, int u, int sink) {
  
  for (int v = 0; v < numVertices; v++) { 
        
    if (d[u] == d[v]+1 && (cap[u][v] - flows[u][v] > 0)) { // push if the height of the adjacent is smaller
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

void PushRelabelParallelSolver::relabel(int numVertices, float **cap, int u) {
  
  int minHeight = INT_MAX; 
  for (int v = 0; v < numVertices; v++) { 
    if ((cap[u][v] - flows[u][v]) > 0) { 
      minHeight = min(minHeight, d[v]); 
    } 
  }
  d[u] = minHeight + 1;
}


void PushRelabelParallelSolver::pushRelabel(MaxFlowInstance *input, MaxFlowSolution *output) {
  t.reset();
  preflow(input); 
  int numVertices = input->inputGraph.num_vertices;
  
  float **cap = input->inputGraph.capacities; 
  // int u = existsActiveNode(input); 
  while (true) {
    std::cout << ""
    //@TODO: make working set a vector to check if empty 
    bool isEmpty = true; 
    for (int i = 0; i < numVertices; i++) { 
      if (workingSet[i] != -1) { 
        isEmpty = false; 
      }
    }
    if (isEmpty) { 
      break; 
    }

//    #pragma omp parallel for (local copy of isDiscovered[i])
    for (int i = 0 ; i < numVertices; i++) { 
      int v = i; 
      if (workingSet[i] != -1) { // checks if the vertex is in the working set 
        for (int j = 0 ; j < numVertices; j++) { 
          discoveredVertices[v][j] = 0; // reinitialize - @TODO: slow? 
        }
        copyOfLabels[v] = d[v]; 
        // copyOfExcess[v] = excessPerVertex[v];
        int e = excessPerVertex[v]; // copy of excess 
        while (e > 0) { 
          int newLabel = numVertices; 
          bool skipped = false;

 //         #pragma omp parallel for (local copy of isDiscovered[w])
          for (int w = 0; w < numVertices; w++) { 
            if (cap[v][w] > 0) { 
              if (e == 0) { 
                break; // has already pushed out all excess flow 
              }
              bool admissible = (d[v] == d[w] +1); 
              if (excessPerVertex[w]) { 
                bool win = (d[v] == d[w]+1) || (d[v] < d[w]-1) || (d[v] == d[w] and v < w); 
                if (admissible && !win) { 
                  skipped = true; 
                  continue; // continue to the next edge 
                }  
              }
              if (admissible && (cap[v][w] - flows[v][w] > 0)) { 
                float delta = min(cap[v][w] - flows[v][w], excessPerVertex[v]); 
                // add residual flow 
                flows[v][w] += delta; 
                flows[v][w] -= delta; 

                e -= delta; 
                // atomic fetch-and-add
                addedPerVertex[w] += delta; 
                if (w != input->sink && isDiscovered[w].test_and_set()) {
                  discoveredVertices[v][w] = 1; // @TODO: make discoveredVertices a vector?
                }
              } 
              if (cap[v][w] - flows[v][w] > 0 && d[w] >= copyOfLabels[v]) { 
                newLabel = min(newLabel, d[w]+1); 
              }
            }
          }
          if (e == 0 || skipped) { 
            break; 
          }
          copyOfLabels[v] = newLabel; 
          if (copyOfLabels[v] == numVertices) { 
            break; 
          }
        }
        addedPerVertex[v] = e - excessPerVertex[v]; 
        if (e and isDiscovered[v].test_and_set()) { 
          discoveredVertices[v][v] = 1; 
        }
      }
    }

    // the update of everything
//    #pragma omp parallel for (local copy of isDiscovered[i])
    for (int i = 0; i < numVertices; i++) { 
      if (workingSet[i] != -1) { 
        d[i] = copyOfLabels[i]; 
        excessPerVertex[i] += addedPerVertex[i]; 
        addedPerVertex[i] = 0; 
        isDiscovered[i].clear(); 
        // the concat of the newly discovered vertices to the working set 
        for (int j = 0; j < numVertices; j++) { 
          if (discoveredVertices[i][j] != -1 && d[i] < numVertices) { 
            workingSet[j] = j; 
          }
        }
      }
    }

//    #pragma omp parallel for
    for (int i = 0; i < numVertices; i++) { 
      if (workingSet[i] != -1) { 
        excessPerVertex[i] += addedPerVertex[i];
        addedPerVertex[i] = 0; 
        isDiscovered[i].clear(); 
      }
    }
  } 
  // old code 
  // while (u != -1) {
  //   if ((u != input->source) && (u != input->sink) && (excessPerVertex[u] > 0)) { 
  //     //printf("u: %d\n", u);
  //     if (!push(numVertices, cap, u, input->sink)) { 
  //       //printf("relabel\n"); 
  //       relabel(numVertices, cap, u); 
  //     } 
  //   }
  //   u = existsActiveNode(input); 
  // }
  double time = t.elapsed(); 
  printf("Push-Relabel Parallel time: %6fms\n", time);

  output->maxFlow = excessPerVertex[input->sink]; 
  output->flow = flows; 
   
}
