#include "../../world.h"
#include <queue> 
#include "../../timing.h"
using namespace std; 

class PushRelabelSequentialSolver{ 
  public: 
    void pushRelabel(MaxFlowInstance *input, MaxFlowSolution *output); 
  private:
    int** flows; // x coordinate is the source, y coordinate is the sink, for edges 
    int *excessPerVertex; // excess flow on each vertex also exists 
    int *d; // the labels 
    queue<int> activeQueue; 
    int totalFlow; 
    Timer t; 

    void initialize(MaxFlowInstance *input);
    void preflow(MaxFlowInstance *input);
    bool push(int numVertices, int **cap, int u, int sink);
    void relabel(int numVertices, int **cap, int u);
    int existsActiveNode(MaxFlowInstance *input);
};

