//
// Created by Ekaterina Villevald on 2019-11-04.
//

#include <iostream>
#include <boost/filesystem.hpp>


#include "lib/world.h"
#include "lib/algorithms/sequential/genetic.h"
#include "lib/algorithms/sequential/pushrelabel.h"
#include "lib/algorithms/sequential/dinics.h"

#include <fstream>
#include <boost/algorithm/string.hpp>

typedef std::vector<std::string> stringvec;
struct path_leaf_string
{
  std::string operator()(const boost::filesystem::directory_entry& entry) const
  {
    return entry.path().leaf().string();
  }
};

void read_directory(const std::string& name, stringvec& v)
{
  boost::filesystem::path p(name);
  boost::filesystem::directory_iterator start(p);
  boost::filesystem::directory_iterator end;
  std::transform(start, end, std::back_inserter(v), path_leaf_string());
}


/* TODO: Remove this deadcode (not being used anywhere) */
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

  // Get all the txt file Test Cases
  stringvec allFiles;
  stringvec testFiles;
  read_directory("tests", allFiles);
  for (int i = 0; i < allFiles.size(); i++){
    if (allFiles[i].find(".txt") != string::npos){
      testFiles.push_back(allFiles[i]);
    }
  }

  // Initialize Instances of Solvers
  PushRelabelSequentialSolver prSolver;
  DinicsSequentialSolver dSolver;

  for (std::string testFileName : testFiles) {

    // Initialize Structures for Test Case
    MaxFlowInstance inputInstance;
    Graph testGraph1;
    float **capacities1;
    MaxFlowSolution prSolution;
    MaxFlowSolution dSolution;


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
          boost::split(results, line, [](char c) { return c == ' '; });

          testGraph1.num_vertices = std::stoi(results[2]);
          testGraph1.num_edges = std::stoi(results[3]);

          capacities1 = new float *[testGraph1.num_vertices];
          for (int i = 0; i < testGraph1.num_vertices; i++) {
            capacities1[i] = new float[testGraph1.num_vertices];
          }
          for (int i = 0; i < testGraph1.num_vertices; i++) {
            for (int j = 0; j < testGraph1.num_vertices; j++) {
              capacities1[i][j] = 0.0f;
            }
          }
          testGraph1.capacities = capacities1;
        }

        if (line[0] == 'n') {
          std::vector<std::string> results;
          boost::split(results, line, [](char c) { return c == ' '; });

          if (results[2] == "s") {
            inputInstance.source = std::stoi(results[1]) - 1;
          } else if (results[2] == "t") {
            inputInstance.sink = std::stoi(results[1]) - 1;
          }

        }

        if (line[0] == 'a') {
          std::vector<std::string> results;
          boost::split(results, line, [](char c) { return c == ' '; });

          testGraph1.capacities[std::stoi(results[1]) - 1][std::stoi(results[2]) - 1] = std::stof(results[3]);
        }
      }
      file.close();
    }
    inputInstance.inputGraph = testGraph1;
    std::cout << "Test Case: " << testFileName << "\n";
    printf("----------------------------\n");
    std::cout << "Input Instance Info:\n";
    std::cout << inputInstance.sink << " <- sink\n";
    std::cout << inputInstance.source << " <- source\n";
    std::cout << inputInstance.inputGraph.num_vertices << " <- numVertices\n";
    std::cout << inputInstance.inputGraph.num_edges << " <- numEdges\n";
    printf("\n");

    // Solve Maxflow with solvers
    std::cout << "Results Info:\n";
    prSolver.pushRelabel(&inputInstance, &prSolution);
    dSolver.solve(&inputInstance, &dSolution);

    printf("push relabel maxflow: %f\n", prSolution.maxFlow);
    printf("dinic maxflow: %f\n", dSolution.maxFlow);
    printf("\n");

  }

  return 0;
}