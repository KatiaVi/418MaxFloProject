#include "../../world.h"
#include <queue> 
#include "../../timing.h"
#include <atomic>
using namespace std; 

class PushRelabelParallelSolver{ 
  public: 
    void pushRelabel(MaxFlowInstance *input, MaxFlowSolution *output); 
  private:
    int** flows; // x coordinate is the source, y coordinate is the sink, for edges 
    int *excessPerVertex; // excess flow on each vertex also exists 
    atomic_int* excessChanges; 
    int *d; // the labels
    int *dCopies;  
    int *active; // replace with a queue 
    queue<int> activeQueue; 
    int totalFlow; 
    Timer t; 

    void pushOnActiveNodes(MaxFlowInstance *input); 
    void updateLabelsAndExcess(int numVertices, int source); 
    void initialize(MaxFlowInstance *input);
    void preflow(MaxFlowInstance *input);
    bool push(int numVertices, int **cap, int u, int sink);
    void relabel(int numVertices, int **cap, int u);
    int existsActiveNode(int numVertices, int source, int sink);
};

