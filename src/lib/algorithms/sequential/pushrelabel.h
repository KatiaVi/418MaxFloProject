#include "../../world.h"
#include <queue> 
using namespace std; 

class PushRelabelSequentialSolver{ 
  public: 
    void pushRelabel(MaxFlowInstance *input, MaxFlowSolution *output); 
  private:
    float** flows; // x coordinate is the source, y coordinate is the sink, for edges 
    float *excessPerVertex; // excess flow on each vertex also exists 
    int *d; // the labels 
    int *active; // replace with a queue 
    queue<int> activeQueue; 
    float totalFlow; 

    void initialize(MaxFlowInstance *input);
    void preflow(MaxFlowInstance *input);
    bool push(int numVertices, float **cap, int u, int sink);
    void relabel(int numVertices, float **cap, int u);
    int existsActiveNode(MaxFlowInstance *input);
};

