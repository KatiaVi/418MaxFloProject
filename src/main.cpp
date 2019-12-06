//
// Created by Ekaterina Villevald on 2019-11-04.
//

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>

#include "lib/world.h"
#include "lib/algorithms/sequential/pushrelabel.h"
#include "lib/algorithms/sequential/dinics.h"
#include "lib/algorithms/parallel/pushrelabel_parallel.h"
#include "lib/algorithms/parallel/dinics_parallel.h"


typedef std::vector<std::string> stringvec;

void read_directory(const char* name, stringvec& v)
{
  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir (name)) != NULL) {
    /* print all the files and directories within directory */
    while ((ent = readdir (dir)) != NULL) {
      v.push_back(ent->d_name);
    }
    closedir (dir);
  } else {
    /* could not open directory */
    perror ("Can't open directory");
  }
}


/* TODO: Remove this deadcode (not being used anywhere) */
void generateCapacities(int numVertices, int numEdges, int **capacities){ 
  capacities = new int*[numVertices]; 
  for (int i = 0; i < numVertices; i++) { 
    capacities[i] = new int[numVertices]; 
  }
  for (int i = 0; i < numVertices; i++) { 
    for (int j = 0; j < numVertices; j++) { 
      capacities[i][j] = 0.0f; 
    }
  }
}

int main ( int argc, char * argv[] ) {

  // Get all the txt file Test Cases
  stringvec allFiles;
  stringvec testFiles;
  read_directory("tests", allFiles);
  for (int i = 0; i < allFiles.size(); i++){
    if (allFiles[i].find("delaunay_n14.txt") != string::npos){ //@TODO: change this to run specific test cases
      testFiles.push_back(allFiles[i]);
    }
  }

  // Initialize Instances of Solvers
  PushRelabelSequentialSolver prSolver;
  PushRelabelParallelSolver parallelPrSolver;
  DinicsSequentialSolver dSolver;
  DinicsParallelSolver parallelDSolver;

  for (std::string testFileName : testFiles) {

    // Initialize Structures for Test Case
    MaxFlowInstance inputInstance;
    MaxFlowInstanceSmall inputInstanceSmall;
    Graph testGraph1;
    GraphSmall testSmallGraph;
    int **capacities1;
    std::vector<std::pair<int,int>>* capacities2;
    std::vector<int>* edges2;
    MaxFlowSolution prSolution;
    MaxFlowSolution parallelPrSolution;
    MaxFlowSolution dSolution;
    MaxFlowSolution parallelDSolution;
    MaxFlowSolutionSmall parallelDSolutionSmall;



    // Open and Read File containing Test Case
    const string &fileName = "tests/" + testFileName;
    std::ifstream file;
    file.open(fileName);

    if (!file.is_open()) {
      std::cout << "not open file\n";
    }

    if (file.is_open()) {

      std::string line;
      while (std::getline(file, line)) {

        if (line[0] == 'p') {
          std::vector<std::string> results;
          std::istringstream iss(line);
          for(std::string s; iss >> s; )
            results.push_back(s);

          testGraph1.num_vertices = std::stoi(results[2]);
          testGraph1.num_edges = std::stoi(results[3]);

          testSmallGraph.num_vertices = std::stoi(results[2]);
          testSmallGraph.num_edges = std::stoi(results[3]);

          capacities1 = new int *[testGraph1.num_vertices];
          capacities2 = (std::vector<std::pair<int,int>>*)calloc(
              testSmallGraph.num_vertices, sizeof(std::vector<std::pair<int,int>>));
          edges2 = (std::vector<int>*)calloc(testSmallGraph.num_vertices, sizeof(std::vector<int>));

          for (int i = 0; i < testGraph1.num_vertices; i++) {
            capacities1[i] = new int[testGraph1.num_vertices];
          }
          for (int i = 0; i < testGraph1.num_vertices; i++) {
            for (int j = 0; j < testGraph1.num_vertices; j++) {
              capacities1[i][j] = 0.0f;
            }
          }
          testGraph1.capacities = capacities1;
          testSmallGraph.capacities = capacities2;
          testSmallGraph.edges = edges2;
        }

        if (line[0] == 'n') {
          std::vector<std::string> results;
          std::istringstream iss(line);
          for(std::string s; iss >> s; )
            results.push_back(s);

          if (results[2] == "s") {
            inputInstance.source = std::stoi(results[1]) - 1;
            inputInstanceSmall.source = std::stoi(results[1]) - 1;
          } else if (results[2] == "t") {
            inputInstance.sink = std::stoi(results[1]) - 1;
            inputInstanceSmall.sink = std::stoi(results[1]) - 1;
          }
        }

        if (line[0] == 'a') {
          std::vector<std::string> results;
          std::istringstream iss(line);
          for(std::string s; iss >> s; )
            results.push_back(s);

          testGraph1.capacities[std::stoi(results[1]) - 1][std::stoi(results[2]) - 1] = std::stoi(results[3]);

          testSmallGraph.capacities[std::stoi(results[1]) - 1]
          .push_back(pair<int,int>(std::stoi(results[2]) - 1, std::stoi(results[3])));
          testSmallGraph.edges[std::stoi(results[1]) - 1].push_back(std::stoi(results[2]) - 1);
        }
      }
      file.close();
    }

    inputInstance.inputGraph = testGraph1;
    inputInstanceSmall.inputGraph = testSmallGraph;
    std::cout << "Test Case: " << testFileName << "\n";
    printf("----------------------------\n");
    std::cout << "Input Instance Info:\n";
    std::cout << inputInstance.sink << " <- sink\n";
    std::cout << inputInstance.source << " <- source\n";
    std::cout << inputInstance.inputGraph.num_vertices << " <- numVertices\n";
    std::cout << inputInstance.inputGraph.num_edges << " <- numEdges\n";
    printf("\n");

    std::cout << "Small Input Instance Info:\n";
    std::cout << inputInstanceSmall.sink << " <- sink\n";
    std::cout << inputInstanceSmall.source << " <- source\n";
    std::cout << inputInstanceSmall.inputGraph.num_vertices << " <- numVertices\n";
    std::cout << inputInstanceSmall.inputGraph.num_edges << " <- numEdges\n";
    printf("\n");

    // Solve Maxflow with solvers
    std::cout << "Results Info:\n";

    if (argc == 1) {
      prSolver.pushRelabel(&inputInstance, &prSolution);
      std::cout << "push relabel maxflow: " << prSolution.maxFlow << "\n";

      parallelPrSolver.pushRelabel(&inputInstance, &parallelPrSolution);
      std::cout << "parallel push relabel maxflow: " << parallelPrSolution.maxFlow << "\n";

      dSolver.solve(&inputInstance, &dSolution);
      std::cout << "dinic maxflow: " << dSolution.maxFlow << "\n";

      parallelDSolver.solve(&inputInstance, &parallelDSolution);
      std::cout << "parallel dinic maxflow: " << parallelDSolution.maxFlow << "\n";
    }
    else if (strcmp(argv[1], "-seq") == 0){
      if (strcmp(argv[2], "dinics") == 0){
        dSolver.solve(&inputInstance, &dSolution);
        std::cout << "dinic maxflow: " << dSolution.maxFlow << "\n";
      }
      if (strcmp(argv[2], "pushrelabel") == 0){
        prSolver.pushRelabel(&inputInstance, &prSolution);
        std::cout << "push relabel maxflow: " << prSolution.maxFlow << "\n";
      }
    }

    else if (strcmp(argv[1], "-par") == 0){
      if (strcmp(argv[2], "dinics") == 0){
        //parallelDSolver.solve(&inputInstance, &parallelDSolution);
        parallelDSolver.smallSolve(&inputInstanceSmall, &parallelDSolutionSmall);
        std::cout << "parallel dinic maxflow: " << parallelDSolution.maxFlow << "\n";
      }
      if (strcmp(argv[2], "pushrelabel") == 0){
        parallelPrSolver.pushRelabel(&inputInstance, &parallelPrSolution);
        std::cout << "parallel push relabel maxflow: " << parallelPrSolution.maxFlow << "\n";
      }
    }
    std::cout << "\n";

  }

  return 0;
}