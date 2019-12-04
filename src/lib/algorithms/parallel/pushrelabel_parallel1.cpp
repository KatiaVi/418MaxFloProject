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
  // #pragma omp parallel for
  for (int i = 0; i < numVertices; i++) { 

void PushRelabelParallelSolver::preflow(MaxFlowInstance *input) { 
  initialize(input); 
  int numVertices = input->inputGraph.num_vertices; 

  // #pragma omp parallel for 
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
  // #pragma omp barrier

  // #pragma omp parallel for 
  for (int i = 0; i < numVertices; i++) { 
    if (cap[input->source][i] != 0 && (input->source != i)) { 
      // then it is an adjacent edge 
      flows[input->source][i] = cap[input->source][i]; 
      //@TODO: make sure i have ints for the flows 
      excessPerVertex[i] = flows[input->source][i];
      // excessChanges[i] = flows[input->source][i];
      // add residual flow 
      flows[i][input->source] = -flows[input->source][i]; 
    }
  }
  // #pragma omp barrier
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



void PushRelabelParallelSolver::push(int numVertices, int **cap, int u, int sink) { 
  printf("u in push: %d\n", u); 
  dCopies[u] = d[u]; 
  int e = excessPerVertex[u]; 
  while (e > 0) {
    int newLabel = numVertices;
    for (int v = 0; v < numVertices; v++) { 
      // printf("u: %d, v: %d, d[u]: %d, d[v]: %d, cap[u][v]: %d, flows[u][v]: %d\n", u, v, d[u], d[v], cap[u][v], flows[u][v]); 
      if (d[u] == d[v]+1 && (cap[u][v] - flows[u][v] > 0)) { // push if the height of the adjacent is smaller
        // push flow = min(remaining flow on edge, excess flow)
        int flow = min(cap[u][v] - flows[u][v], excessPerVertex[u]); //@TODO: bug: adding 0 here 
        // reduce excess flow for overflowing vertex 
        e -= flow; 
        //atomic_fetch_add(&excessChanges[u], -flow); 
        atomic_fetch_add(&excessChanges[v], flow); 
        
        // why do you need a separate step for updating the excess and labels? - labels I get, but not really the excesses 
        // excessChanges[u] -= flow; 
        // excessChanges[v] += flow; 

        // add residual flow 
        flows[u][v] += flow; 
        flows[v][u] -= flow; 
        // return true; 
      } 
      if ((cap[u][v] - flows[u][v] > 0) && d[v] >= dCopies[u]) {
        newLabel = min(newLabel, d[v]+1); 
      }
    }
    if (e == 0) { // if pushed all the flow, then you won't need to relabel (might increase but won't assign it)
      break; //@TODO: when is excessPerVertex[i] set to 0? 
    }
    dCopies[u] = newLabel; 
    if (dCopies[u] == numVertices) { //@TODO: not sure why 
      break; 
    }
  }
  
  // needsRelabeling[u] = true; 
  // return false; 
}

void PushRelabelParallelSolver::updateLabelsAndExcess(int numVertices, int source) { // where to call this? 
  // #pragma omp parallel for 
  for (int i = 0; i < numVertices; i++) { 
    d[i] = dCopies[i]; 
    excessPerVertex[i] += excessChanges[i]; //@TODO: maybe needs to be an atomic fetch and add? maybe not idk 
  }
  // #pragma omp barrier //@TODO: not sure if this is needed after a omp parallel for 

  // #pragma omp parallel for 
  for (int i = 0; i < numVertices; i++) { 
    dCopies[i] = 0; 
    excessChanges[i] = 0; 
  } // @TODO: must reinitialize excessChanges as what i initialized initially? or just to 0?  
  dCopies[source] = numVertices;
  // #pragma omp barrier
  // for (int i = 0; i < numVertices; i++) { 
  //   printf("d[%d]: %d\n", i, d[i]); 
  //   printf("excessPerVertex[%d]: %d\n", i, excessPerVertex[i]); 
  // }
}

// note that you cannot do a parallel for inside a parallel for 
void PushRelabelParallelSolver::relabel(int numVertices, int **cap, int u) { 
  int minHeight = INT_MAX; 

  for (int v = 0; v < numVertices; v++) { 
    if ((cap[u][v] - flows[u][v]) > 0) { 
      printf("dCopies[%d]: %d\n", v, dCopies[v]); 
      minHeight = min(minHeight, dCopies[v]); //@TODO: I think this might be dCopies[v] instead? 
       
    } 
  }
  dCopies[u] = minHeight + 1;
  printf("setting dCopies[%d] = %d\n", u, dCopies[u]);
}

// finds active nodes and in parallel, on each active node, pushes on their admissible edges (sequential)
void PushRelabelParallelSolver::pushOnActiveNodes(MaxFlowInstance *input) {
  int numVertices = input->inputGraph.num_vertices;
  int **cap = input->inputGraph.capacities;
  // #pragma omp parallel for 
  for (int i = 0; i < numVertices; i++) { 
    if (excessPerVertex[i] > 0 && i != input->source && i != input->sink) { //@TODO: initialize active at some point?  
      // loop through its edges and see if they are admissible 
      //  printf("i: %d\n", i); 
      push(numVertices, cap, i, input->sink); 
      
    }
  }
  // #pragma omp barrier 
  // after doing all the pushes, do all the relabels 
  // #pragma omp parallel for 
 
  // #pragma omp barrier 
  // print all the labels and the excesses after this step 

  
}

void PushRelabelParallelSolver::pushRelabel(MaxFlowInstance *input, MaxFlowSolution *output) { 
  t.reset();
  preflow(input); 
  int numVertices = input->inputGraph.num_vertices;

  int **cap = input->inputGraph.capacities;
  int u = existsActiveNode(numVertices, input->source, input->sink); 
  while (u != -1) { //@TODO: could replace with something that actually gets the active set? so that dont have to loop through so many times 
    pushOnActiveNodes(input); //safe i think 
    printf("active u: %d\n", u); 
    printf("excessPerVertex[1] after push: %d\n", excessPerVertex[1]); 
    // apparently there always exists and active node? 
    updateLabelsAndExcess(numVertices, input->source); 
    u = existsActiveNode(numVertices, input->source, input->sink); 
  }
  double time = t.elapsed(); 
  printf("Push-Relabel time: %6fs\n", time); 

  output->maxFlow = excessPerVertex[input->sink]; 
  output->flow = flows; 
}



