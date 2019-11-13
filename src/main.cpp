//
// Created by Ekaterina Villevald on 2019-11-04.
//

#include <iostream>
#include "lib/world.h"
#include "lib/algorithms/sequential/genetic.h"
#include "lib/algorithms/sequential/pushrelabel.h"
#include "lib/algorithms/sequential/dinics.h"

void generateCapacities(int numVertices, int numEdges, float **capacities){ 
  capacities = new float*[numVertices]; 
  for (int i = 0; i < numVertices; i++) { 
    capacities[i] = new float[numVertices]; 
  }
  for (int i = 0; i < numVertices; i++) { 
    for (int j = 0; j < numVertices; j++) { 
      capacities[i][j] = 0.0f; 
    }
  }
}

int main ( int argc, char * argv[] )
{
  // try calling genetic algo solve function on a sample graph
  // std::cout << argc << argv[0] << "Hello!\n";

  // maxFlow should be 5
  // taken from: http://www.cs.cmu.edu/afs/cs/academic/class/15451-s16/www/lectures/lec12-flow1.pdf
  float **capacities1; 
  int numVertices = 6; 
  int numEdges = 8; 
  capacities1 = new float*[numVertices];
  for (int i = 0; i < numVertices; i++) {
    capacities1[i] = new float[numVertices];
  }
  for (int i = 0; i < numVertices; i++) { 
    for (int j = 0; j < numVertices; j++) { 
      capacities1[i][j] = 0.0f; 
    }
  }
  
  Graph testGraph1;
  testGraph1.num_vertices = 6;
  testGraph1.num_edges = 8;
  testGraph1.capacities = capacities1;

  std::cout << "creating test graph\n";

  testGraph1.capacities[0][1] = 4.0f;
  testGraph1.capacities[0][2] = 2.0f;
  testGraph1.capacities[1][3] = 3.0f;
  testGraph1.capacities[2][3] = 2.0f;
  testGraph1.capacities[2][4] = 3.0f;
  testGraph1.capacities[3][2] = 1.0f;
  testGraph1.capacities[3][5] = 2.0f;
  testGraph1.capacities[4][5] = 4.0f;

  for (int i=0; i < 6; i++){
    for (int j=0; j < 6; j++){
      std::cout << testGraph1.capacities[i][j] << " ";
    }
    std::cout << "\n";
  }

  std::cout << "finished creating test graph\n";

  // should be 6 
  Graph testGraph2; 
  float **capacities2; 
  int numVertices2 = 4; 
  int numEdges2 = 5; 
  capacities2 = new float*[numVertices2]; 
  for (int i = 0; i < numVertices2; i++) { 
    capacities2[i] = new float[numVertices2]; 
  }
  for (int i = 0; i < numVertices2; i++) { 
    for (int j = 0; j < numVertices2; j++) { 
      capacities2[i][j] = 0.0f; 
    }
  }
  testGraph2.num_vertices = numVertices2; 
  testGraph2.num_edges = numEdges2; 
  testGraph2.capacities = capacities2; 

  testGraph2.capacities[0][1] = 2; 
  testGraph2.capacities[0][2] = 4; 
  testGraph2.capacities[1][2] = 3; 
  testGraph2.capacities[1][3] = 1; 
  testGraph2.capacities[2][3] = 5; 

  // test graph from wikipedia: https://en.wikipedia.org/wiki/Push%E2%80%93relabel_maximum_flow_algorithm
  // should be 14 
  Graph testGraph3; 
  float **capacities3; 
  int numVertices3 = 6; 
  int numEdges3 = 8; 
  capacities3 = new float*[numVertices3]; 
  for (int i = 0; i < numVertices3; i++) { 
    capacities3[i] = new float[numVertices3]; 
  }
  for (int i = 0; i < numVertices3; i++) { 
    for (int j = 0; j < numVertices3; j++) { 
      capacities3[i][j] = 0.0f; 
    }
  }
  testGraph3.num_vertices = numVertices3; 
  testGraph3.num_edges = numEdges3; 
  testGraph3.capacities = capacities3; 

  testGraph3.capacities[0][1] = 15; 
  testGraph3.capacities[0][3] = 4; 
  testGraph3.capacities[1][2] = 12; 
  testGraph3.capacities[2][5] = 7; 
  testGraph3.capacities[2][3] = 3; 
  testGraph3.capacities[3][4] = 10; 
  testGraph3.capacities[4][1] = 5; 
  testGraph3.capacities[4][5] = 10; 

  MaxFlowInstance inputInstance;
  inputInstance.inputGraph = testGraph1;
  inputInstance.sink = 5;
  inputInstance.source = 0;

  MaxFlowInstance inputInstance2;
  inputInstance2.inputGraph = testGraph2;
  inputInstance2.sink = 3;
  inputInstance2.source = 0;

  MaxFlowInstance inputInstance3;
  inputInstance3.inputGraph = testGraph3;
  inputInstance3.sink = numVertices3-1;
  inputInstance3.source = 0;

  MaxFlowSolution gSolution;
  GeneticSequentialSolver gSolver;
  gSolver.solve(inputInstance, gSolution);
  std::cout << "genetic maxflow is " << gSolution.maxFlow << "\n";

  PushRelabelSequentialSolver prSolver; 
  MaxFlowSolution prSolution; 
  //prSolver.pushRelabel(&inputInstance, &prSolution); 
  prSolver.pushRelabel(&inputInstance3, &prSolution); 
  assert(prSolution.maxFlow == 12); 
  printf("push relabel maxflow: %f\n", prSolution.maxFlow);

  MaxFlowSolution dSolution;
  DinicsSequentialSolver dSolver;
  dSolver.solve(inputInstance3, dSolution);
  printf("dinic maxflow: %f\n", dSolution.maxFlow);


  

  return 0;
}