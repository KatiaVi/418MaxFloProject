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
#include <unordered_set>
#include <iostream>

#include "pushrelabel_parallel.h"
#include "tbb/concurrent_vector.h"

using namespace std; 

void PushRelabelParallelSolver::initialize(MaxFlowInstance *input) {
  int numVertices = input->inputGraph.num_vertices;
  d = (atomic_int*)calloc(numVertices, sizeof(atomic_int)); 
  active = (int*)calloc(numVertices, sizeof(int)); 
  excessPerVertex = (int*)calloc(numVertices, sizeof(int)); 
  copyOfExcess = (int*)calloc(numVertices, sizeof(int));
  addedExcess = (atomic_int*)calloc(numVertices, sizeof(atomic_int));
  workingSet.clear(); 
  
  copyOfLabels = (int*)calloc(numVertices, sizeof(int)); 
  work = (int*)calloc(numVertices, sizeof(int)); 
  
  flows = (int**)calloc(numVertices, sizeof(int*));
  reverseResiduals = (vector<int> *)calloc(numVertices, sizeof(vector<int>));
  residual = (int**)calloc(numVertices, sizeof(int*));

  #pragma omp parallel for
  for (int i = 0; i < numVertices; i++) { 
    flows[i] = (int*)calloc(numVertices, sizeof(int));
    residual[i] = (int*)calloc(numVertices, sizeof(int));
  }
}

void PushRelabelParallelSolver::preflow(MaxFlowInstance *input) {
  initialize(input); 
  int numVertices = input->inputGraph.num_vertices; 

  d[input->source] = numVertices; 
  copyOfLabels[input->source] = numVertices; 
  int **cap = input->inputGraph.capacities;

  tbb::concurrent_vector<int> cv = {};
  std::vector<std::pair<int,int>> vec = {};

  #pragma omp parallel for 
  for (int i = 0; i < numVertices; i++) { 
    if (cap[input->source][i] != 0 && (input->source != i)) {
      // then it is an adjacent edge 
      flows[input->source][i] = cap[input->source][i];
      residual[input->source][i] = cap[input->source][i] - flows[input->source][i]; 

      excessPerVertex[i] = flows[input->source][i];
      // add residual flow 
      flows[i][input->source] = -flows[input->source][i];
    }
    for (int j = 0; j < numVertices; j++) {
      if (cap[i][j] - flows[i][j] > 0) {
          residual[i][j] = cap[i][j] - flows[i][j];
      }
    }

  }

  for (int i = 0; i < numVertices; i++) {
    discoveredVertices.push_back(cv);
    if (i != input->source && excessPerVertex[i] > 0) { 
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
  tbb::concurrent_vector<int> q;
  q.push_back(sink);

  #pragma omp parallel for
  for (int v=0; v < numVertices; v++) {
    for (int w = 0; w < numVertices; w++) {
      if (residual[v][w] > 0) {
        #pragma omp critical
        {
          reverseResiduals[w].push_back(v);
        }
      }
    }
  }


  while (!q.empty()) {  
    #pragma omp parallel for
    for (int i = 0; i < q.size(); i++) { 
      int v = q[i];
      discoveredVertices[v].clear();

      for (int i = 0; i < reverseResiduals[v].size(); i++){
        int w = reverseResiduals[v][i];
        int tmp = numVertices;
        if (d[w].compare_exchange_strong(tmp, d[v]+1)) {
          discoveredVertices[v].push_back(w);
        }
      }
    }

    tbb::concurrent_vector<int> newq;
    for (int v : q) { 
      #pragma omp parallel for
      for (int j=0; j < discoveredVertices[v].size(); j++){
        newq.push_back(discoveredVertices[v][j]);
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

    if (freq * workSinceLastGR > a * numVertices + numEdges) { 
      workSinceLastGR = 0;
      globalRelabel(numVertices, input->source, input->sink);
      #pragma omp parallel for
      for (int i = 0; i < numVertices; i++) { 
        copyOfLabels[i] = d[i];
      }
     
      unordered_set<int> newWorkingSet; 
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

    #pragma omp parallel 
    { 
      #pragma omp single
      { 
        for (auto iter = workingSet.begin(); iter != workingSet.end(); iter++) { 
          #pragma omp task untied
          {
            int v = *iter;
            discoveredVertices[v].clear();
            copyOfLabels[v] = d[v]; 
            int e = excessPerVertex[v]; // copy of excess
            work[v] = 0;
            while (e > 0) { 
              int newLabel = numVertices; 
              bool skipped = false;

              int numEdgesScanned = 0; 
              
              for (int w = 0; w < numVertices; w++) {

                if (residual[v][w] > 0) {
                  if (e == 0) {
                    break; // has already pushed out all excess flow
                  }
                    numEdgesScanned += 1;
                    bool admissible = (copyOfLabels[v] == d[w] + 1);
                    if (admissible) {
                      int delta = min(residual[v][w], e);

                      // add residual flow
                      flows[v][w] += delta;
                      flows[w][v] -= delta;

                      e -= delta;
                      atomic_fetch_add(&addedExcess[w], delta);

                      if (w != input->sink) {
                        discoveredVertices[v].push_back(w);
                      }
                    }
                    if (d[w] >= copyOfLabels[v]) {
                      newLabel = min(newLabel, d[w] + 1);
                    }
                }
              }

              if (e == 0 || skipped) {
                break;
              }

              copyOfLabels[v] = newLabel;
              work[v] += numEdgesScanned + beta;
              if (copyOfLabels[v] == numVertices) {
                break;
              }
            }

            addedExcess[v] += (e - excessPerVertex[v]);
            if (e > 0) {
              discoveredVertices[v].push_back(v);
            }
          }
        }
      }
      #pragma omp taskwait 
    }

    #pragma omp parallel for  
    for (int i = 0; i < numVertices; i++) { 
      d[i] = copyOfLabels[i];
      excessPerVertex[i] += addedExcess[i];
      addedExcess[i] = 0;
      reverseResiduals[i].clear();
    }

    unordered_set<int> workingSetNew; 
    for (int i : workingSet) {
      workSinceLastGR += work[i]; 
      for (int k = 0; k < discoveredVertices[i].size(); k++){
        int j = discoveredVertices[i][k];
        if(d[j] < numVertices) {

          residual[i][j] = cap[i][j] - flows[i][j];
          residual[j][i] = cap[j][i] - flows[j][i];


          workingSetNew.insert(j); 
        }
      }
      residual[i][input->sink] = cap[i][input->sink] - flows[i][input->sink]; 
      residual[input->sink][i] = cap[input->sink][i] - flows[input->sink][i];
    }
    workingSet.swap(workingSetNew);


    #pragma omp parallel  
    { 
      #pragma omp single 
      { 
        for (auto iter = workingSet.begin(); iter != workingSet.end(); iter++) { 
          #pragma omp task 
          { 
            int i = *iter; 
            excessPerVertex[i] += addedExcess[i]; 
            addedExcess[i] = 0;
          }
        } 
      }
    }
  } 

  double time = t.elapsed(); 
  printf("Push-Relabel Parallel time: %6fs\n", time);

  output->maxFlow = excessPerVertex[input->sink]; 
  output->flow = flows; 
   
}
