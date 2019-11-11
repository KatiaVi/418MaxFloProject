//
// Created by Ekaterina Villevald on 2019-11-04.
//

#include <iostream>
#include "lib/world.h"
#include "lib/algorithms/sequential/genetic.h"
<<<<<<< Updated upstream
#include "lib/algorithms/sequential/pushrelabel.h"
=======
#include "lib/algorithms/sequential/dinics.h"
>>>>>>> Stashed changes

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
  //generateCapacities(6, 8, capacities1);
  
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


  int source = 0;
  int sink = 5;

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

  MaxFlowInstance inputInstance;
  inputInstance.inputGraph = testGraph1;
  inputInstance.sink = sink;
  inputInstance.source = source;

  MaxFlowInstance inputInstance2;
  inputInstance2.inputGraph = testGraph2;
  inputInstance2.sink = 3;
  inputInstance2.source = 0;


  MaxFlowSolution gSolution;
  GeneticSequentialSolver gSolver;
  gSolver.solve(inputInstance, gSolution);

  std::cout << "genetic maxflow is " << gSolution.maxFlow << "\n";
  gSolver.printSolutions(inputInstance.inputGraph.num_vertices);
  /*GeneticSequentialSolver gSolver;
  gSolution = gSolver.solve(inputInstance);*/ 
  PushRelabelSequentialSolver prSolver; 
  MaxFlowSolution prSolution; 
  //prSolver.pushRelabel(&inputInstance, &prSolution); 
  prSolver.pushRelabel(&inputInstance, &prSolution); 
  printf("push relabel maxflow: %f\n", prSolution.maxFlow);

  std::cout << "dinic maxflow is " << gSolution.maxFlow << "\n";

//  MaxFlowSolution dSolution;
//  DinicsSequentialSolver dSolver;
//  dSolver.solve(inputInstance, dSolution);
//
//  std::cout << "maxflow is " << dSolution.maxFlow << "\n";
  //gSolver.printSolutions(inputInstance.inputGraph.num_vertices);

  // gSolution.maxFlow should be 5
  // http://www.cs.cmu.edu/afs/cs/academic/class/15451-s16/www/lectures/lec12-flow1.pdf
  // (using this example graph)

  return 0;
}