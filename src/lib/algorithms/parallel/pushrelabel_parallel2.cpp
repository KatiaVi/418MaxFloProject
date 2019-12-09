//
// Created by Ekaterina Villevald on 2019-12-09.
//

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

#include "pushrelabel_parallel2.h"
#include "tbb/concurrent_vector.h"

using namespace std;

void PushRelabelParallelSolverV2::initialize(MaxFlowInstance *input) {
  int numVertices = input->inputGraph.num_vertices;
  d = (atomic_int*)calloc(numVertices, sizeof(atomic_int));
  active = (int*)calloc(numVertices, sizeof(int));
  excessPerVertex = (int*)calloc(numVertices, sizeof(int));
  copyOfExcess = (int*)calloc(numVertices, sizeof(int));
  // outDegrees = (int*)calloc(numVertices, sizeof(int));
  addedExcess = (atomic_int*)calloc(numVertices, sizeof(atomic_int));
  workingSet.clear();

  copyOfLabels = (int*)calloc(numVertices, sizeof(int));
  work = (int*)calloc(numVertices, sizeof(int));

  flows = (int**)calloc(numVertices, sizeof(int*));
  reverseResiduals = (vector<int> *)calloc(numVertices, sizeof(vector<int>));
//  residual = (int**)calloc(numVertices, sizeof(int*));
//  residual = (vector<std::pair<int,int>> *)calloc(numVertices, sizeof(vector<std::pair<int,int>>));



#pragma omp parallel for
  for (int i = 0; i < numVertices; i++) {
    flows[i] = (int*)calloc(numVertices, sizeof(int));
//    residual[i] = (int*)calloc(numVertices, sizeof(int));
  }
}

void PushRelabelParallelSolverV2::preflow(MaxFlowInstance *input) {
  initialize(input);
  int numVertices = input->inputGraph.num_vertices;

  d[input->source] = numVertices;
  copyOfLabels[input->source] = numVertices;
  // find all vertices adjacent to s
  int **cap = input->inputGraph.capacities;

  tbb::concurrent_vector<int> cv = {};
  std::vector<std::pair<int,int>> vec = {};

#pragma omp parallel for
  for (int i = 0; i < numVertices; i++) {
    if (cap[input->source][i] != 0 && (input->source != i)) {
      // then it is an adjacent edge
      flows[input->source][i] = cap[input->source][i];
      //   residual[input->source][i] = cap[input->source][i] - flows[input->source][i]; // -> moved this to another for loop so that this one could be parallelized

      excessPerVertex[i] = flows[input->source][i];
      // add residual flow
      flows[i][input->source] = -flows[input->source][i];
      // residual[i].push_back(make_pair(input->source, cap[i][input->source] - flows[i][input->source]));
      // residual[i][input->source] = cap[i][input->source] - flows[i][input->source];

      //  residual[i].push_back(make_pair(input->source, cap[i][input->source] - flows[i][input->source]))
    }
//    std::cout << "yo\n";

  }
//  std::cout << "finished preflow1\n";


  // don't put a parallel for here!
//  tbb::concurrent_vector<int> cv = {};
  for (int i = 0; i < numVertices; i++) {
    residual.push_back(vec);
    discoveredVertices.push_back(cv);

    for (int j = 0; j < numVertices; j++) {
      if (cap[i][j] - flows[i][j] > 0) {
        residual[i].push_back(make_pair(j, cap[i][j] - flows[i][j]));
        //  residual[i][j] = cap[i][j] - flows[i][j];
      }
    }

    if (i != input->source && excessPerVertex[i] > 0) { //@TODO: took out the && i != input->sink
      workingSet.insert(i);
    }
  }
  std::cout << "finished preflow\n";
}

void PushRelabelParallelSolverV2::globalRelabel(int numVertices, int source, int sink) {
#pragma omp parallel for
  for (int i = 0; i < numVertices; i++) {
    d[i] = numVertices;
  }
  d[sink] = 0;
  tbb::concurrent_vector<int> q;
  q.push_back(sink);

#pragma omp parallel for
  for (int v=0; v < numVertices; v++) {
    for (int k = 0; k < residual[v].size(); k++) {
      int w = residual[v][k].first;
      if (residual[v][k].second <= 0) {
        //    std::cout << "residual is negative\n";
      }
#pragma omp critical
      {
        reverseResiduals[w].push_back(v);
      }
      //   }
    }
  }
  std::cout<< "set reverseResiduals\n";


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

    //@TODO: done with parallel prefix sum
    tbb::concurrent_vector<int> newq;
    for (int v : q) {
      // dont put a parallel for here
#pragma omp parallel for
      for (int j=0; j < discoveredVertices[v].size(); j++){
        newq.push_back(discoveredVertices[v][j]);
      }
    }
    q.swap(newq);
  }
}


void PushRelabelParallelSolverV2::pushRelabel(MaxFlowInstance *input, MaxFlowSolution *output) {
  t.reset();
  preflow(input);
//  std::cout << "preflow done\n";



  int numVertices = input->inputGraph.num_vertices;
  int numEdges = input->inputGraph.num_edges;

//  for(int i=0; i < numVertices; i++){
//    std::cout << "vertex: " << i << ":\n";
//    for(int j=0; j < residual[i].size(); j++){
//      std::cout << residual[i][j].first << " ";
//    }
//    std::cout << "\n";
//  }

  int **cap = input->inputGraph.capacities;
  int workSinceLastGR = INT_MAX;
  float freq = 0.5;
  int a = 6;
  int beta = 12;

  while (true) {

    if (freq * workSinceLastGR > a * numVertices + numEdges) {
      workSinceLastGR = 0;
      globalRelabel(numVertices, input->source, input->sink);
//      std::cout << "relabel!\n";
#pragma omp parallel for
      for (int i = 0; i < numVertices; i++) {
        copyOfLabels[i] = d[i];
      }

      // @TODO: parallel array comprehension using map/filter
      unordered_set<int> newWorkingSet;
      for (int v : workingSet) {
        if (d[v] < numVertices) {
          newWorkingSet.insert(v);
        }
      }
      workingSet.swap(newWorkingSet);
    }

//    for(int i=0; i < numVertices; i++){
//      std::cout << "vertex: " << i << ":\n";
//      for(int j=0; j < residual[i].size(); j++){
//        std::cout << residual[i][j].first << " ";
//      }
//      std::cout << "\n";
//    }

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
              for (int k = 0; k < residual[v].size(); k++) {
                int w = residual[v][k].first;

                //    if (residual[v][w] > 0) {
                if (e == 0) {
                  break; // has already pushed out all excess flow
                }
                numEdgesScanned += 1;
                bool admissible = (copyOfLabels[v] == d[w] + 1);
                if (admissible) {

                  // int delta = min(residual[v][w], e);
                  int delta = min(residual[v][k].second, e);
//                      std::cout << "made a delta\n";


                  // add residual flow
                  flows[v][w] += delta;
                  flows[w][v] -= delta;

                  e -= delta;
                  atomic_fetch_add(&addedExcess[w], delta);

                  if (w != input->sink) {
                    //@TODO ??? why dont u just add it to the set
                    discoveredVertices[v].push_back(w);//[w] = 1; // @TODO: make discoveredVertices a vector?
                  }
                }
                if (d[w] >= copyOfLabels[v]) {
                  newLabel = min(newLabel, d[w] + 1);
                }
                //    }
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

            addedExcess[v] += (e - excessPerVertex[v]); //@TODO: problematic line - because doesn't take into account the changes that have already been made to the excess - overwrites addedExcess if it was discovered by another thing in the workingSet

            // the line below doesnt really make sense - oh would make sense if it was skipped - @TODO: add back later
            // could be if broke out at this line: if (copyOfLabels[v] == numVertices)
            if (e > 0) {
              discoveredVertices[v].push_back(v);//[v] = 1;
            }
          }
        }
      }
#pragma omp taskwait
    }

    // printf("updating everything\n");
    // the update of everything
#pragma omp parallel for
    for (int i = 0; i < numVertices; i++) {
      d[i] = copyOfLabels[i];
      excessPerVertex[i] += addedExcess[i];
      addedExcess[i] = 0;
      reverseResiduals[i].clear();
    }

    //  std::cout << "GOT HERE\n";

    // could create a copy of the working set and work off of that
    // make working set an actual set in c++
    // create new working set
    unordered_set<int> workingSetNew;
    // this loop was slowwww
    for (int i : workingSet) {
      workSinceLastGR += work[i]; // combined in here
      for (int k = 0; k < discoveredVertices[i].size(); k++){
        int j = discoveredVertices[i][k];
        if(d[j] < numVertices) {
          std::vector<std::pair<int,int>>::iterator residualij = std::find_if(residual[i].begin(), residual[i].end(),
                                                                              [j](std::pair<int,int> tup){ return (tup.first== j); });

          std::vector<std::pair<int,int>>::iterator residualji = std::find_if(residual[j].begin(), residual[j].end(),
                                                                              [i](std::pair<int,int> tup){ return (tup.first== i); });

          std::vector<std::pair<int, int>>
              newResiduals= {std::make_pair(j, cap[i][j] - flows[i][j]),
                             std::make_pair(i, cap[j][i] - flows[j][i])};

          //  std::cout << "old size: " << residual[i].size() << "\n";

          if (residualij != residual[i].end() && newResiduals.front().second > 0){
            (*residualij).swap(newResiduals.front());
          }
          else if (residualij != residual[i].end() && newResiduals.front().second <= 0){
            residual[i].erase(residualij);
            //    std::cout << "new size: " << residual[i].size() << "\n";
          }
          else if (newResiduals.front().second > 0){
            residual[i].push_back(newResiduals.front());
            //    std::cout << "new size: " << residual[i].size() << "\n";

          }

          if (residualji != residual[j].end() && newResiduals.back().second > 0){
            (*residualji).swap(newResiduals.back());
          }
          else if (residualji != residual[j].end() && newResiduals.back().second <= 0){
            residual[j].erase(residualji);
          }
          else if (newResiduals.back().second > 0){
            residual[j].push_back(newResiduals.back());
          }

//          (*residualij).swap(newResiduals.front());
//          (*residualji).swap(newResiduals.back());

//          printf("(%d,%d)\n", (*residualij).first, (*residualij).second);
//          printf("(%d,%d)\n", (*residualji).first, (*residualji).second);

//         std::cout << "swapped these residuals\n";

//          residual[i][j] = cap[i][j] - flows[i][j];
//          residual[j][i] = cap[j][i] - flows[j][i];


          workingSetNew.insert(j);
        }
      }
      int t = input->sink;
//      printf("sink: %d\n", t);


      std::vector<std::pair<int,int>>::iterator residualisink = std::find_if(residual[i].begin(), residual[i].end(),
                                                                             [t](std::pair<int,int> tup){ return (tup.first== t); });
      std::vector<std::pair<int,int>>::iterator residualsinki = std::find_if(residual[t].begin(), residual[t].end(),
                                                                             [i](std::pair<int,int> tup){ return (tup.first== i); });

      std::vector<std::pair<int, int>>
          newResiduals = {std::make_pair(t, cap[i][t] - flows[i][t]),
                          std::make_pair(i, cap[t][i] - flows[t][i])};

      if (residualisink != residual[i].end() && newResiduals.front().second > 0) {
        (*residualisink).swap(newResiduals.front());
      }
      else if(residualisink != residual[i].end() && newResiduals.front().second <= 0){
        residual[i].erase(residualisink);
      }
      else if (newResiduals.front().second > 0){
        residual[i].push_back(newResiduals.front());
      }

      if (residualsinki != residual[t].end() && newResiduals.back().second > 0) {
        (*residualsinki).swap(newResiduals.back());
      }
      else if(residualsinki != residual[t].end() && newResiduals.back().second <= 0){
        residual[t].erase(residualsinki);
      }
      else if (newResiduals.back().second > 0){
        residual[t].push_back(newResiduals.back());
      }




//      residual[i][input->sink] = cap[i][input->sink] - flows[i][input->sink]; //@TODO: maybe can do this access in a more cache friendly way
//      residual[input->sink][i] = cap[input->sink][i] - flows[input->sink][i];
    }
    workingSet.swap(workingSetNew);

//    for(int i=0; i < numVertices; i++){
//      std::cout << i << ":\n";
//      for(int j=0; j < residual[i].size(); j++){
//        std::cout << residual[i][j].second << " ";
//      }
//      std::cout << "\n";
//    }



    //@TODO: to debug: print out the things in workingSetNew
    // for (int i = 0; i < numVertices; i++) {
    //   if (workingSet[i] != -1 && d[i] < numVertices) { //&& excessPerVertex[i] > 0 makes it stop running forever
    //     // the concat of the newly discovered vertices to the working set

    //   }
    // }

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

