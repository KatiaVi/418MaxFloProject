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

#include "pushrelabel_parallel1.h"

using namespace std; 

void PushRelabelParallelSolver::initialize(MaxFlowInstance *input) { 
  int numVertices = input->inputGraph.num_vertices;
  d = (int*)calloc(numVertices, sizeof(int)); 
  dCopies = (int*)calloc(numVertices, sizeof(int)); 
  active = (int*)calloc(numVertices, sizeof(int)); 
  excessPerVertex = (int*)calloc(numVertices, sizeof(int)); 
  excessChanges = (atomic_int*)malloc(numVertices*sizeof(atomic_int)); 
  flows = new int*[numVertices]; 
  #pragma omp parallel for
  for (int i = 0; i < numVertices; i++) { 
    flows[i] = new int[numVertices]; 
    excessChanges[i] = ATOMIC_VAR_INIT(0); 
  }
}

void PushRelabelParallelSolver::preflow(MaxFlowInstance *input) { 
  initialize(input); 
  int numVertices = input->inputGraph.num_vertices; 

  #pragma omp parallel for 
  for (int i = 0; i < numVertices; i++) { 
    d[i] = 0; 
    dCopies[i] = 0; 
    for (int j = 0; j < numVertices; j++) { 
      flows[i][j] = 0; // flow going from vertex i to vertex j is 0 
    }
  }
  

  d[input->source] = numVertices; 
  dCopies[input->source] = numVertices; 
  // find all vertices adjacent to s 
  int **cap = input->inputGraph.capacities; 
  #pragma omp barrier

  #pragma omp parallel for 
  for (int i = 0; i < numVertices; i++) { 
    if (cap[input->source][i] != 0 && (input->source != i)) { 
      // then it is an adjacent edge 
      flows[input->source][i] = cap[input->source][i]; 
      //@TODO: make sure i have ints for the flows 
      excessPerVertex[i] = flows[input->source][i];
      excessChanges[i] = flows[input->source][i];
      // add residual flow 
      flows[i][input->source] = -flows[input->source][i]; 
    }
  }
  #pragma omp barrier
  // initialize active nodes to be all that have non zero excess
}

int PushRelabelParallelSolver::existsActiveNode(int numVertices, int source, int sink) { 
  for (int i = 0; i < numVertices; i++) { 
    if (excessPerVertex[i] > 0 && i != source && i != sink) { //@TODO: initialize active at some point?  
      return i; 
    }
  }
  return -1; 
}

// finds active nodes and in parallel, on each active node, pushes on their admissible edges (sequential)
void PushRelabelParallelSolver::pushOnActiveNodes(MaxFlowInstance *input) {
  int numVertices = input->inputGraph.num_vertices;
  int **cap = input->inputGraph.capacities;
  //@TODO: will this work in the way i think it will? 
  #pragma omp parallel for 
  for (int i = 0; i < numVertices; i++) { 
    if (excessPerVertex[i] > 0 && i != input->source && i != input->sink) { //@TODO: initialize active at some point?  
       // loop through its edges and see if they are admissible 
       for (int j = 0; j < numVertices; j++) { 
         if (!push(numVertices, cap, i, input->sink)) { 
           relabel(numVertices, cap, i); 
         }
       }
    }
  }
  // put a barrier - @TODO: is there an implicit barrier at the end of functions anyway? 
  #pragma omp barrier 
}

bool PushRelabelParallelSolver::push(int numVertices, int **cap, int u, int sink) { 
  
  for (int v = 0; v < numVertices; v++) { 
        
    if (d[u] == d[v]+1 && (cap[u][v] - flows[u][v] > 0)) { // push if the height of the adjacent is smaller
      // push flow = min(remaining flow on edge, excess flow)
      int flow = min(cap[u][v] - flows[u][v], excessPerVertex[u]); //@TODO: bug: adding 0 here 
      // reduce excess flow for overflowing vertex 
      atomic_fetch_add(&excessChanges[u], -flow); 
      atomic_fetch_add(&excessChanges[v], flow); 
      // excessChanges[u] -= flow; 
      // excessChanges[v] += flow; 

      // add residual flow 
      flows[u][v] += flow; 
      flows[v][u] -= flow; 
      return true; 
    }
  }
  return false; 
}

void PushRelabelParallelSolver::relabel(int numVertices, int **cap, int u) { 
  int minHeight = INT_MAX; 
  #pragma omp parallel for 
  for (int v = 0; v < numVertices; v++) { 
    if ((cap[u][v] - flows[u][v]) > 0) { 
      minHeight = min(minHeight, d[v]); //@TODO: I think it is d not dCopies here 
    } 
  }
  dCopies[u] = minHeight + 1;
  #pragma omp barrier 
}

void PushRelabelParallelSolver::updateLabelsAndExcess(int numVertices, int source) { // where to call this? 
  #pragma omp parallel for 
  for (int i = 0; i < numVertices; i++) { 
    d[i] = dCopies[i]; 
    excessPerVertex[i] += excessChanges[i]; //@TODO: I dont think that you need to know what is in 
    // the working set if you clear the excessChanges and dCopies on each round 
    // but not sure, can try both for measuring performance 
  }
  #pragma omp barrier //@TODO: not sure if this is needed after a omp parallel for 
  #pragma omp parallel for 
  for (int i = 0; i < numVertices; i++) { 
    dCopies[i] = 0; 
    excessChanges[i] = 0; 
  } // must reinitialize excessChanges as what i initialized initially? or just to 0?  
  dCopies[source] = numVertices;
  #pragma omp barrier
}


void PushRelabelParallelSolver::pushRelabel(MaxFlowInstance *input, MaxFlowSolution *output) { 
  t.reset();
  preflow(input); 
  int numVertices = input->inputGraph.num_vertices;

  int **cap = input->inputGraph.capacities;
  int u = existsActiveNode(numVertices, input->source, input->sink); 
  while (u != -1) { //@TODO: could replace with something that actually gets the active set? so that dont have to loop through so many times 
    pushOnActiveNodes(input);
    updateLabelsAndExcess(numVertices, input->source); 
    u = existsActiveNode(numVertices, input->source, input->sink); 
    printf("u: %d\n", u); 
  }
}


// void PushRelabelParallelSolver::pushRelabel(MaxFlowInstance *input, MaxFlowSolution *output) { 
//   t.reset();
//   preflow(input); 
//   int numVertices = input->inputGraph.num_vertices;

//   int **cap = input->inputGraph.capacities; 
//   // int u = existsActiveNode(input); 

//   while (!activeQueue.empty()) {
//     int u = activeQueue.front(); 
//     activeQueue.pop(); 

//     while ((u != input->source) && (u != input->sink) && (excessPerVertex[u] > 0)) { 
//       //printf("u: %d\n", u);
//       if (!push(numVertices, cap, u, input->sink)) { 
//         //printf("relabel\n"); 
//         relabel(numVertices, cap, u); 
//       } 

//       //testing code 
//       /* for (int i = 0; i < numVertices; i++) { 
//         for (int j = 0; j < numVertices; j++) { 
//           if (flows[i][j] != 0) { 
//             printf("flows[%d][%d]: %f\n", i, j, flows[i][j]); 
//           }
//         }
//       }
//       for (int i = 0; i < numVertices; i++) { 
//         printf("d[%d]: %d\n", i, d[i]); 
//       }
//       for (int i = 0; i < numVertices; i++) { 
//         printf("excessPerVertex[%d]: %f\n", i, excessPerVertex[i]); 
//       } */ 
//     }
//     // u = existsActiveNode(input); 
//   }
//   double time = t.elapsed(); 
//   printf("Push-Relabel time: %6fs\n", time); 

//   output->maxFlow = excessPerVertex[input->sink]; 
//   output->flow = flows; 
   
// }
