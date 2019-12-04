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
#include <set> 
#include <queue> 

#include "pushrelabel_parallel.h"

using namespace std; 

void PushRelabelParallelSolver::initialize(MaxFlowInstance *input) {
  int numVertices = input->inputGraph.num_vertices;
  d = (int*)calloc(numVertices, sizeof(int)); 
  active = (int*)calloc(numVertices, sizeof(int)); 
  excessPerVertex = (int*)calloc(numVertices, sizeof(int)); 
  copyOfExcess = (int*)calloc(numVertices, sizeof(int)); 
  addedPerVertex = (int*)calloc(numVertices, sizeof(int));
  isDiscovered = (atomic_bool*)malloc(numVertices*sizeof(atomic_bool));
  workingSet.clear(); 
  // workingSet = (int*)calloc(numVertices, sizeof(int)); 
  discoveredVertices = new int*[numVertices]; 
  copyOfLabels = (int*)calloc(numVertices, sizeof(int)); 
  work = (int*)calloc(numVertices, sizeof(int)); 
  
  flows = new int*[numVertices]; 
  residual = new int*[numVertices]; // every time I update flows I should be updating residuals 

//  #pragma omp parallel for (local copy of isDiscovered[i])
  for (int i = 0; i < numVertices; i++) { 
    flows[i] = new int[numVertices]; 
    residual[i] = new int[numVertices]; 
    discoveredVertices[i] = new int[numVertices]; 
    // workingSet[i] = -1; // initialize all to -1 so that only the active ones get added
    isDiscovered[i] = ATOMIC_FLAG_INIT;// isDiscovered[i].clear();
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
      residual[input->source][i] = cap[input->source][i] - flows[input->source][i]; //@START
      //@TODO: make sure i have ints for the flows 
      excessPerVertex[i] = flows[input->source][i];
      // add residual flow 
      flows[i][input->source] = -flows[input->source][i]; 
      residual[i][input->source] = cap[i][input->source] - flows[i][input->source]; 
    }
  }
//  #pragma omp parallel for
  for (int i = 0; i < numVertices; i++) { 
    if (i != input->source && excessPerVertex[i] > 0) { //@TODO: took out the && i != input->sink 
      // workingSet[i] = i; 
      workingSet.insert(i); 
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

void PushRelabelParallelSolver::globalRelabel(int numVertices, int source, int sink) { 
  for (int i = 0; i < numVertices; i++) { 
    d[i] = numVertices; 
  }
  d[sink] = 0; 
  vector<int> q; 
  q.push_back(sink); 
  while (!q.empty()) { 
    for (int v : q) { 
      for (int j = 0; j < numVertices; j++) { 
        discoveredVertices[v][j] = -1; 
      }
      for (int w = 0; w < numVertices; w++) { 
        if (w != source && residual[w][v] > 0) { 
          //@TODO: make this atomic compare and swap 
          if (w != sink && d[w] == numVertices) { 
            d[w] = d[v]+1;
            discoveredVertices[v][w] = 1; 
          }
        }
      }
    }
    // done with parallel prefix sum 
    vector<int> newq; 
    for (int v : q) { 
      for (int j = 0; j < numVertices; j++) { 
        if (discoveredVertices[v][j] != -1) { 
          newq.push_back(j); 
        }
      }
    }
    q.swap(newq); 
  }
}


void PushRelabelParallelSolver::pushRelabel(MaxFlowInstance *input, MaxFlowSolution *output) {
  t.reset();
  preflow(input); 
  int numVertices = input->inputGraph.num_vertices;
  int numEdges = input->inputGraph.num_edges;
  
  int **cap = input->inputGraph.capacities; 
  int workSinceLastGR = INT_MAX; 
  float freq = 0.5; 
  int a = 6; 
  int beta = 12; 
  // int u = existsActiveNode(input);
  while (true) {

    //@TODO: make working set a vector to check if empty 
    // bool isEmpty = true; 
    // for (int i = 0; i < numVertices; i++) { 
    //   if (workingSet[i] != -1) { 
    //     printf("i is in the working set: %d\n", i); 
    //     isEmpty = false; 
    //   }
    // }
    // if (isEmpty) { 
    //   break; 
    // }

    // update redisual capacities
    std::cout << "Residuals:\n";
    for (int i=0; i < numVertices; i++){
      for (int j=0; j < numVertices; j++){
        residual[i][j] = cap[i][j] - flows[i][j];
        // if (residual[i][j] > 0) {
        //   printf("(%d,%d): %d\n", i, j, residual[i][j]);
        // }
      }
      // std::cout << "\n";s
    }
    // std::cout << "\n";
    // std::cout << "Excess:\n";
    // for (int j=0; j < numVertices; j++){
    //   std::cout << excessPerVertex[j] << " ";
    // }
    // std::cout << "\n";

    if (freq * workSinceLastGR > a * numVertices + numEdges) { 
      printf("does a global relabel\n"); 
      workSinceLastGR = 0; 
      globalRelabel(numVertices, input->source, input->sink); 
      for (int i = 0; i < numVertices; i++) { 
        copyOfLabels[i] = d[i];
      }
      
      // @TODO: parallel array comprehension using map/filter 
      set<int> newWorkingSet; 
      for (int v : workingSet) { 
        if (d[v] < numVertices) { 
          newWorkingSet.insert(v); 
        }
      }
      workingSet.swap(newWorkingSet); 
    }

    if (workingSet.empty()) { 
      break; 
    }

    


    std::cout << "Working set:\n";
    for (int v: workingSet){
      std::cout << v << " ";
    }
    std::cout << "\n";
//    #pragma omp parallel for (local copy of isDiscovered[i])
    for (int v : workingSet) { 
      printf("new v from working set: %d\n", v);

        for (int j = 0 ; j < numVertices; j++) { 
          discoveredVertices[v][j] = -1; // reinitialize @TODO: 
          // @TODO: make discoveredVertices into an array of sets otherwise will have way too big of a matrix for larger test cases, also then fill in with not just 1 but the actual vertex index   
        }
        copyOfLabels[v] = d[v]; 
        // copyOfExcess[v] = excessPerVertex[v];
        int e = excessPerVertex[v]; // copy of excess 
        printf("excess for %d: %d\n", v, e); 
        work[v] = 0; 
        while (e > 0) { 
          int newLabel = numVertices; 
          bool skipped = false;

          int numEdgesScanned = 0; 
          // #pragma omp parallel for (local copy of isDiscovered[w])
          for (int w = 0; w < numVertices; w++) {
            printf("testing edge (%d, %d)\n", v, w); 
            // printf("cap[%d][%d]: %d\n", v, w, cap[v][w]); 
            if (residual[v][w] /*cap[v][w]- flows[v][w]*/ > 0) { // it's a residual edge // cap[v][w] > 0 && - maybe is trying to push back onto start?
              if (e == 0) {
                printf("excess[%d] is 0\n", v); 
                //workingSet[v] = -1;
                break; // has already pushed out all excess flow 
              }
              numEdgesScanned += 1; 
              bool admissible = (copyOfLabels[v] == d[w] + 1);
              printf("edge from %d to %d is admissible: %d\n", v, w, admissible); 
              // if (excessPerVertex[w]) {
              //   bool win = (d[v] == d[w]+1) || (d[v] < d[w]-1) || (d[v] == d[w] and v < w);
              //   if (admissible && !win) { 
              //     skipped = true; 
              //     continue; // continue to the next edge 
              //   }  
              // }

              if (admissible && (residual[v][w] /* cap[v][w]- flows[v][w]*/ > 0)) {
                
                int delta = min(residual[v][w] /* cap[v][w] - flows[v][w]*/, e);
                printf("pushing %d flow from %d to %d\n", delta, v, w); 
                // add residual flow 
                flows[v][w] += delta; 
                flows[w][v] -= delta; 

                e -= delta; 
                printf("new excess flow for %d: %d\n", v, e); 
                //@TODO: atomic fetch-and-add
                addedPerVertex[w] += delta;  // MAKE ATOMIC
                if (w != input->sink) { // && (isDiscovered[w]).exchange(true)) { //@TODO ??? why dont u just add it to the set 
                  // printf("discovers %d\n", w); 
                  discoveredVertices[v][w] = 1; // @TODO: make discoveredVertices a vector?
                }
                printf("AFTER PUSH - e: %d\n", e);
                printf("AFTER PUSH - addedPerVertex[%d]: %d\n", w, addedPerVertex[w]);
                //printf("AFTER PUSH - (cap[%d][%d]-flows[%d][%d]): %d\n", v, w, v, w, (cap[v][w] - flows[v][w]));
              }
              printf("d[%d] = %d \n", w, d[w]);

              if (residual[v][w] /* cap[v][w]- flows[v][w]*/ > 0 && d[w] >= copyOfLabels[v]) {
                newLabel = min(newLabel, d[w]+1);
                // printf("new label for %d is: %d\n", v, newLabel);
              }
            }
          }

          if (e == 0 || skipped) {
            break;
          }
          printf("new label for %d is: %d\n", v, newLabel);
          copyOfLabels[v] = newLabel; 
          work[v] += numEdgesScanned + beta; 
          if (copyOfLabels[v] == numVertices) {
            break;
          }
        }
        printf("addedPerVertex[%d]: %d\n", v, addedPerVertex[v]);
        printf("e=%d and excessPerVertex[%d]: %d\n", e, v, excessPerVertex[v]);

        addedPerVertex[v] += (e - excessPerVertex[v]); //@TODO: problematic line - because doesn't take into account the changes that have already been made to the excess - overwrites addedPerVertex if it was discovered by another thing in the workingSet
        printf("addedPerVertex[%d]: %d\n", v, addedPerVertex[v]);
        std::cout << "AddedExcess:\n";
        for (int j=0; j < numVertices; j++){
          std::cout << addedPerVertex[j] << " ";
        }
        std::cout << "\n";

        // the line below doesnt really make sense - oh would make sense if it was skipped - @TODO: add back later 
        // could be if broke out at this line: if (copyOfLabels[v] == numVertices)
        if (e > 0) { // this line is wrong because e will always be = 0 after the while loop
           discoveredVertices[v][v] = 1;
        }
    }

    // for (int i = 0 ; i < numVertices; i++) {
    //   int v = i;

    //   if (workingSet[i] != -1) { // checks if the vertex is in the working set 
        
        
    //     
    //   }
    // }

    printf("updating everything\n"); 
    // the update of everything
    // #pragma omp parallel for (local copy of isDiscovered[i])
    for (int i = 0; i < numVertices; i++) { 
      // if (workingSet[i] != -1) { 
      d[i] = copyOfLabels[i]; // 2 was removed from the working set too early to have its copy of labels be updated 
      // printf("d[%d]: %d\n", i, d[i]); 
      excessPerVertex[i] += addedPerVertex[i]; 
      // printf("addedPerVertex[%d]: %f\n", i, addedPerVertex[i]); 
      addedPerVertex[i] = 0; 
      
      isDiscovered[i] = ATOMIC_FLAG_INIT; // isDiscovered[i].clear(); 
      // }
    }

    std::cout << "New D:\n";
    for (int j=0; j < numVertices; j++){
      std::cout << d[j] << " ";
    }
    std::cout << "\n";

    for (int v : workingSet) { 
      workSinceLastGR += work[v]; 
    }
    
    // could create a copy of the working set and work off of that 
    // make working set an actual set in c++ 
    // create new working set 
    set<int> workingSetNew; 
    for (int i : workingSet) {
      //if (d[i] < numVertices) {
        for (int j = 0; j < numVertices; j++) { 
          printf("%d has these discoveredVertices[%d][%d]: %d\n", i, i, j, discoveredVertices[i][j]); 
          if (discoveredVertices[i][j] != -1 and d[j] < numVertices) { // && excessPerVertex[j] > 0
            workingSetNew.insert(j); 
          }
        }
      //}
    }
    workingSet.swap(workingSetNew);

    //@TODO: to debug: print out the things in workingSetNew 
    // for (int i = 0; i < numVertices; i++) { 
    //   if (workingSet[i] != -1 && d[i] < numVertices) { //&& excessPerVertex[i] > 0 makes it stop running forever 
    //     // the concat of the newly discovered vertices to the working set 
        
    //   } 
    // }
    
    for (int i : workingSet) { 
      printf("in working set: %d\n", i); 
    }
    // break; 
//    #pragma omp parallel for
    for (int i : workingSet) { 
      // if (workingSet[i] != -1) { 
        excessPerVertex[i] += addedPerVertex[i]; // @TODO: the discovered one might have had stuff added to it 
        // printf("excess: %f\n", excessPerVertex[i]);
        addedPerVertex[i] = 0; 
        isDiscovered[i] = ATOMIC_FLAG_INIT; //isDiscovered[i].clear(); 
      // }
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
