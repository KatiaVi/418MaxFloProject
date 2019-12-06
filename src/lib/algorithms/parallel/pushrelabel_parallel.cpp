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
  d = (atomic_int*)calloc(numVertices, sizeof(atomic_int)); 
  active = (int*)calloc(numVertices, sizeof(int)); 
  excessPerVertex = (int*)calloc(numVertices, sizeof(int)); 
  copyOfExcess = (int*)calloc(numVertices, sizeof(int)); 
  addedExcess = (int*)calloc(numVertices, sizeof(int));
  isDiscovered = (atomic_bool*)malloc(numVertices*sizeof(atomic_bool));
  workingSet.clear(); 
  
  copyOfLabels = (int*)calloc(numVertices, sizeof(int)); 
  work = (int*)calloc(numVertices, sizeof(int)); 
  
  flows = (int**)calloc(numVertices, sizeof(int*));
  discoveredVertices = (int**)calloc(numVertices, sizeof(int*));
  residual = (int**)calloc(numVertices, sizeof(int*));

  #pragma omp parallel for
  for (int i = 0; i < numVertices; i++) { 
    flows[i] = (int*)calloc(numVertices, sizeof(int));
    discoveredVertices[i] = (int*)calloc(numVertices, sizeof(int));
    residual[i] = (int*)calloc(numVertices, sizeof(int));
    isDiscovered[i] = ATOMIC_FLAG_INIT;
  }
}

void PushRelabelParallelSolver::preflow(MaxFlowInstance *input) {
  initialize(input); 
  int numVertices = input->inputGraph.num_vertices; 

  d[input->source] = numVertices; 
  copyOfLabels[input->source] = numVertices; 
  // find all vertices adjacent to s 
  int **cap = input->inputGraph.capacities; 

  #pragma omp parallel for 
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
    for (int j = 0; j < numVertices; j++) { 
      residual[i][j] = cap[i][j] - flows[i][j]; 
    }
  }

  // don't put a parallel for here! 
  for (int i = 0; i < numVertices; i++) { 
    if (i != input->source && excessPerVertex[i] > 0) { //@TODO: took out the && i != input->sink 
      // workingSet[i] = i; 
      workingSet.insert(i); 
    }
  }
}

void PushRelabelParallelSolver::globalRelabel(int numVertices, int source, int sink) { 
  #pragma omp parallel for 
  for (int i = 0; i < numVertices; i++) { 
    d[i] = numVertices; 
  }
  d[sink] = 0; 
  vector<int> q; 
  q.push_back(sink); 
  while (!q.empty()) {  
    #pragma omp parallel for  
    for (int i = 0; i < q.size(); i++) { 
      int v = q[i];  
      for (int j = 0; j < numVertices; j++) { // clearing this is a little slow 
        discoveredVertices[v][j] = 0; 
      }
      for (int w = 0; w < numVertices; w++) { 
        if (w != source && residual[w][v] > 0) { //@TODO: spends a lot of time checking this line in global relabel for all the vertices 
          //@TODO!: make this atomic compare and swap 
          int tmp = numVertices; 
          if (w != sink && d[w].compare_exchange_strong(tmp, d[v]+1)) { 
            discoveredVertices[v][w] = 1; 
          }
          // if (w != sink && d[w] == numVertices) { 
          //   d[w] = d[v]+1;
          //   discoveredVertices[v][w] = 1; 
          // }
        }
      }
    }
    //@TODO: done with parallel prefix sum 
    vector<int> newq; 
    for (int v : q) { 
      for (int j = 0; j < numVertices; j++) { 
        if (discoveredVertices[v][j]) { 
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
  
  while (true) {
    // std::cout << "Residuals:\n";
    // #pragma omp parallel for
    // for (int i=0; i < numVertices; i++){
    //   for (int j=0; j < numVertices; j++){
    //     residual[i][j] = cap[i][j] - flows[i][j]; // this is where a lot of time is - hmm 
    //     // if (residual[i][j] > 0) {
    //     //   printf("(%d,%d): %d\n", i, j, residual[i][j]);
    //     // }
    //   }
    //   // std::cout << "\n";s
    // }
    // std::cout << "\n";
    // std::cout << "Excess:\n";
    // for (int j=0; j < numVertices; j++){
    //   std::cout << excessPerVertex[j] << " ";
    // }
    // std::cout << "\n";

    if (freq * workSinceLastGR > a * numVertices + numEdges) { 
      // printf("does a global relabel\n"); 
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
    
    // std::cout << "Working set:\n";
    // for (int v: workingSet){
    //   std::cout << v << " ";
    // }
    // std::cout << "\n"; 
    #pragma omp parallel 
    { 
      #pragma omp single 
      { 
        for (auto iter = workingSet.begin(); iter != workingSet.end(); iter++) { 
      // printf("new v from working set: %d\n", v);
          #pragma omp task 
          { 
            int v = *iter; 
            for (int j = 0 ; j < numVertices; j++) { 
              discoveredVertices[v][j] = 0; // reinitialize @TODO: 
              // @TODO: make discoveredVertices into an array of sets otherwise will have way too big of a matrix for larger test cases, also then fill in with not just 1 but the actual vertex index   
            }
            copyOfLabels[v] = d[v]; 
            // copyOfExcess[v] = excessPerVertex[v];
            int e = excessPerVertex[v]; // copy of excess 
            // printf("excess for %d: %d\n", v, e); 
            work[v] = 0; 
            while (e > 0) { 
              int newLabel = numVertices; 
              bool skipped = false;

              int numEdgesScanned = 0; 
              // #pragma omp parallel for (local copy of isDiscovered[w])
              for (int w = 0; w < numVertices; w++) {
                // printf("testing edge (%d, %d)\n", v, w); 
                // printf("cap[%d][%d]: %d\n", v, w, cap[v][w]); 
                if (residual[v][w] /*cap[v][w]- flows[v][w]*/ > 0) { // it's a residual edge // cap[v][w] > 0 && - maybe is trying to push back onto start?
                  if (e == 0) {
                    // printf("excess[%d] is 0\n", v); 
                    //workingSet[v] = -1;
                    break; // has already pushed out all excess flow 
                  }
                  numEdgesScanned += 1; 
                  bool admissible = (copyOfLabels[v] == d[w] + 1);
                  // printf("edge from %d to %d is admissible: %d\n", v, w, admissible); 
                  // if (excessPerVertex[w]) {
                  //   bool win = (d[v] == d[w]+1) || (d[v] < d[w]-1) || (d[v] == d[w] and v < w);
                  //   if (admissible && !win) { 
                  //     skipped = true; 
                  //     continue; // continue to the next edge 
                  //   }  
                  // }

                  if (admissible && (residual[v][w] /* cap[v][w]- flows[v][w]*/ > 0)) {
                    
                    int delta = min(residual[v][w] /* cap[v][w] - flows[v][w]*/, e);
                    // printf("pushing %d flow from %d to %d\n", delta, v, w); 
                    // add residual flow 
                    flows[v][w] += delta; 
                    flows[w][v] -= delta; 

                    e -= delta; 
                    // printf("new excess flow for %d: %d\n", v, e); 
                    //@TODO: atomic fetch-and-add
                    addedExcess[w] += delta;  // MAKE ATOMIC
                    if (w != input->sink) { // && (isDiscovered[w]).exchange(true)) { //@TODO ??? why dont u just add it to the set 
                      // printf("discovers %d\n", w); 
                      discoveredVertices[v][w] = 1; // @TODO: make discoveredVertices a vector?
                    }
                    // printf("AFTER PUSH - e: %d\n", e);
                    // printf("AFTER PUSH - addedExcess[%d]: %d\n", w, addedExcess[w]);
                    //printf("AFTER PUSH - (cap[%d][%d]-flows[%d][%d]): %d\n", v, w, v, w, (cap[v][w] - flows[v][w]));
                  }
                  // printf("d[%d] = %d \n", w, d[w]);

                  if (residual[v][w] /* cap[v][w]- flows[v][w]*/ > 0 && d[w] >= copyOfLabels[v]) {
                    newLabel = min(newLabel, d[w]+1);
                    // printf("new label for %d is: %d\n", v, newLabel);
                  }
                }
              }

              if (e == 0 || skipped) {
                break;
              }
              // printf("new label for %d is: %d\n", v, newLabel);
              copyOfLabels[v] = newLabel; 
              work[v] += numEdgesScanned + beta; 
              if (copyOfLabels[v] == numVertices) {
                break;
              }
            }
            // printf("addedExcess[%d]: %d\n", v, addedExcess[v]);
            // printf("e=%d and excessPerVertex[%d]: %d\n", e, v, excessPerVertex[v]);

            addedExcess[v] += (e - excessPerVertex[v]); //@TODO: problematic line - because doesn't take into account the changes that have already been made to the excess - overwrites addedExcess if it was discovered by another thing in the workingSet
            // printf("addedExcess[%d]: %d\n", v, addedExcess[v]);
            // std::cout << "AddedExcess:\n";
            // for (int j=0; j < numVertices; j++){
            //   std::cout << addedExcess[j] << " ";
            // }
            // std::cout << "\n";

            // the line below doesnt really make sense - oh would make sense if it was skipped - @TODO: add back later 
            // could be if broke out at this line: if (copyOfLabels[v] == numVertices)
            if (e > 0) { 
              discoveredVertices[v][v] = 1;
            }
          }
        }
      }
    }
    

    // for (int i = 0 ; i < numVertices; i++) {
    //   int v = i;

    //   if (workingSet[i] != -1) { // checks if the vertex is in the working set 
        
        
    //     
    //   }
    // }

    // printf("updating everything\n"); 
    // the update of everything
    #pragma omp parallel for  
    for (int i = 0; i < numVertices; i++) { 
      // if (workingSet[i] != -1) { 
      d[i] = copyOfLabels[i]; // 2 was removed from the working set too early to have its copy of labels be updated 
      // printf("d[%d]: %d\n", i, d[i]); 
      excessPerVertex[i] += addedExcess[i]; 
      // printf("addedExcess[%d]: %f\n", i, addedExcess[i]); 
      addedExcess[i] = 0; 
      
      isDiscovered[i] = ATOMIC_FLAG_INIT; // isDiscovered[i].clear(); 
      // }
    }

    // std::cout << "New D:\n";
    // for (int j=0; j < numVertices; j++){
    //   std::cout << d[j] << " ";
    // }
    // std::cout << "\n";
    
    // for (int v : workingSet) { 
    //   workSinceLastGR += work[v]; 
    // }
    
    // could create a copy of the working set and work off of that 
    // make working set an actual set in c++ 
    // create new working set 
    set<int> workingSetNew; 
    // this loop was slowwww 
    for (int i : workingSet) {
      //if (d[i] < numVertices) {
        workSinceLastGR += work[i]; // combined in here 
        for (int j = 0; j < numVertices; j++) { 
          //@TODO: put the update to residual capacities here? - 2nd line is needed
          residual[i][j] = cap[i][j] - flows[i][j]; 
          residual[j][i] = cap[j][i] - flows[j][i]; // do this in the initializer 
          // printf("%d has these discoveredVertices[%d][%d]: %d\n", i, i, j, discoveredVertices[i][j]); 
          if (discoveredVertices[i][j] and d[j] < numVertices) { // && excessPerVertex[j] > 0
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
    
    /*for (int i : workingSet) { 
      printf("in working set: %d\n", i); 
    }*/ 
    // break; 
    // have some kind of for here
    #pragma omp parallel  
    { 
      #pragma omp single 
      { 
        for (auto iter = workingSet.begin(); iter != workingSet.end(); iter++) { 
          #pragma omp task 
          { 
            int i = *iter; 
            excessPerVertex[i] += addedExcess[i]; // @TODO: the discovered one might have had stuff added to it 
            // printf("excess: %f\n", excessPerVertex[i]);
            addedExcess[i] = 0; 
            isDiscovered[i] = ATOMIC_FLAG_INIT;
          }
        } 
      }
    }
    // for (int i : workingSet) { 
    //   // if (workingSet[i] != -1) { 
    //     excessPerVertex[i] += addedExcess[i]; // @TODO: the discovered one might have had stuff added to it 
    //     // printf("excess: %f\n", excessPerVertex[i]);
    //     addedExcess[i] = 0; 
    //     isDiscovered[i] = ATOMIC_FLAG_INIT; //isDiscovered[i].clear(); 
    //   // }
    // }
  } 

  double time = t.elapsed(); 
  printf("Push-Relabel Parallel time: %6fs\n", time);

  output->maxFlow = excessPerVertex[input->sink]; 
  output->flow = flows; 
   
}
