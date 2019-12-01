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
  excessPerVertex = (int*)calloc(numVertices, sizeof(int)); 
  copyOfExcess = (int*)calloc(numVertices, sizeof(int)); 
  addedPerVertex = (int*)calloc(numVertices, sizeof(int));
  isDiscovered = (atomic_flag*)malloc(numVertices*sizeof(atomic_flag));
  workingSet = (int*)calloc(numVertices, sizeof(int)); 
  discoveredVertices = new int*[numVertices]; 
  copyOfLabels = (int*)calloc(numVertices, sizeof(int)); 
  
  flows = new int*[numVertices]; 

//  #pragma omp parallel for (local copy of isDiscovered[i])
  for (int i = 0; i < numVertices; i++) { 
    flows[i] = new int[numVertices]; 
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
  copyOfLabels[input->source] = numVertices; 
  // find all vertices adjacent to s 
  int **cap = input->inputGraph.capacities; 

//  #pragma omp parallel for
  for (int i = 0; i < numVertices; i++) { 
    if (cap[input->source][i] != 0 && (input->source != i)) { 
      // then it is an adjacent edge 
      flows[input->source][i] = cap[input->source][i]; 
      //@TODO: make sure i have ints for the flows 
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

bool PushRelabelParallelSolver::push(int numVertices, int **cap, int u, int sink) {
  
  for (int v = 0; v < numVertices; v++) { 
        
    if (d[u] == d[v]+1 && (cap[u][v] - flows[u][v] > 0)) { // push if the height of the adjacent is smaller
      // push flow = min(remaining flow on edge, excess flow)
      int flow = min(cap[u][v] - flows[u][v], excessPerVertex[u]); //@TODO: bug: adding 0 here 
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

void PushRelabelParallelSolver::relabel(int numVertices, int **cap, int u) {
  
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
  
  int **cap = input->inputGraph.capacities; 
  // int u = existsActiveNode(input);
  while (true) {
    for (int k=0; k < numVertices; k++){
      std::cout << d[k] << " ";
    }
    std::cout << "\n";
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
      printf("v: %d\n", v);

      if (workingSet[i] != -1) { // checks if the vertex is in the working set 
        printf("new v from working set: %d\n", v); 
        for (int j = 0 ; j < numVertices; j++) { 
          discoveredVertices[v][j] = -1; // reinitialize - @TODO: slow? 
        }
        copyOfLabels[v] = d[v]; 
        // copyOfExcess[v] = excessPerVertex[v];
        int e = excessPerVertex[v]; // copy of excess 
        while (e > 0) { 
          int newLabel = numVertices; 
          bool skipped = false;

          // #pragma omp parallel for (local copy of isDiscovered[w])
          for (int w = 0; w < numVertices; w++) {
            printf("testing edge (%d, %d)\n", v, w); 
            if (cap[v][w] > 0 && cap[v][w] - flows[v][w] > 0) { // it's a residual edge 
              if (e == 0) {
                printf("excess[%d] is 0\n", v); 
                //workingSet[v] = -1;
                break; // has already pushed out all excess flow 
              }
              bool admissible = (copyOfLabels[v] == d[w] + 1);
              printf("edge from %d to %d is admissible: %d\n", v, w, admissible); 
              // std::cout << admissible << "<-- admissible\n";
              if (excessPerVertex[w]) {
                bool win = (d[v] == d[w]+1) || (d[v] < d[w]-1) || (d[v] == d[w] and v < w);
                if (admissible && !win) { 
                  skipped = true; 
                  continue; // continue to the next edge 
                }  
              }

              if (admissible && (cap[v][w] - flows[v][w] > 0)) {
                int delta = min(cap[v][w] - flows[v][w], excessPerVertex[v]); 
                printf("pushing %f flow from %d to %d\n", delta, v, w); 
                // add residual flow 
                flows[v][w] += delta; 
                flows[w][v] -= delta; 

                e -= delta; 
                printf("new excess flow for %d: %f\n", v, e); 
                //@TODO: atomic fetch-and-add
                addedPerVertex[w] += delta;  // MAKE ATOMIC
                
                if (w != input->sink && (isDiscovered[w]).test_and_set()) { //@TODO ???
                  printf("HERE\n"); 
                  discoveredVertices[v][w] = 1; // @TODO: make discoveredVertices a vector?
                }
                printf("AFTER PUSH - e: %f\n", v, e);
                printf("AFTER PUSH - addedPerVertex[%d]: %f\n", w, addedPerVertex[w]);
              }
              // std::cout << d[w] << "<- d[w]\n";
              if (cap[v][w] - flows[v][w] > 0 && d[w] >= copyOfLabels[v]) {
                // NOT HERE

                newLabel = min(newLabel, d[w]+1);
                // printf("new label for %d is: %d\n", v, newLabel); 
                // std::cout << "new label: " << newLabel << d[w] << "<- d[w] " << w << "<- w " << v << "<-- v\n";
              }
            }
          }

          if (e == 0 || skipped) {
            /* if (e == 0) { 
              workingSet[v] = -1; // we may be removing it in the wrong place since the update depends on things still in working set 
            }*/ 
            break;
          }
          // std::cout << newLabel << " newLabel after for loop\n";
          printf("new label for %d is: %d\n", v, newLabel);
          copyOfLabels[v] = newLabel; 
          if (copyOfLabels[v] == numVertices) {
            break;
          }
        }
        addedPerVertex[v] = e - excessPerVertex[v]; // @TODO: oh this is the amount added to v?
        printf("addedPerVertex[%d]: %f\n", v, addedPerVertex[v]);
        printf("excess:%f\n", e);
        if (e > 0 && isDiscovered[v].test_and_set()) {
          discoveredVertices[v][v] = 1; 
        }
      }
    }

    printf("updating everything\n"); 
    // the update of everything
    // #pragma omp parallel for (local copy of isDiscovered[i])
    for (int i = 0; i < numVertices; i++) { 
      if (workingSet[i] != -1) { 
        d[i] = copyOfLabels[i]; // 2 was removed from the working set too early to have its copy of labels be updated 
        printf("d[%d]: %d\n", i, d[i]); 
        excessPerVertex[i] += addedPerVertex[i]; 
        printf("addedPerVertex[%d]: %f\n", i, addedPerVertex[i]); 
        addedPerVertex[i] = 0; 
        isDiscovered[i].clear(); 
      }
      if (workingSet[i] != -1 && d[i] < numVertices) { 
        printf("here\n"); 
        // the concat of the newly discovered vertices to the working set 
        for (int j = 0; j < numVertices; j++) { 
          printf("%d has these discoveredVertices[%d][%d]: %d\n", i, i, j, discoveredVertices[i][j]); 
          if (discoveredVertices[i][j] != -1) {
            workingSet[j] = j; 
          }
        }
      } 
      if (excessPerVertex[i] == 0) { 
        workingSet[i] = -1; 
      } 
      // } else { // take things out of the working set 
      //   workingSet[i] = -1; //@TODO: not sure why they aren't just recalculating the workingSet by looking at the excesses (why the discovered vertices thing?)
      // }
    }
    for (int k = 0; k < numVertices; k++) { 
      printf("excessPerVertex[%d]: %f\n", k, excessPerVertex[k]); 
    }
    
    
    

//    #pragma omp parallel for
    for (int i = 0; i < numVertices; i++) { 
      if (workingSet[i] != -1) { 
        excessPerVertex[i] += addedPerVertex[i];
        printf("excess: %f\n", excessPerVertex[i]);
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
