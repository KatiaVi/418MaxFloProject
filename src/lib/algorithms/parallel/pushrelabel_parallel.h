#include "../../world.h"
#include "../../timing.h"
#include <atomic>
using namespace std; 

class PushRelabelParallelSolver{
  public: 
    void pushRelabel(MaxFlowInstance *input, MaxFlowSolution *output); 
  private:
    float** flows; // x coordinate is the source, y coordinate is the sink, for edges 
    
    float *excessPerVertex; // excess flow on each vertex also exists 
    int *d; // the labels 
    float *addedPerVertex; // added for prsn 
    atomic_flag *isDiscovered; // added for prsn
    int **discoveredVertices; 
    int *copyOfLabels; 
    float *copyOfExcess; 

    int *active; // replace with a queue 
    Timer t; 
    int *workingSet;  
    void initialize(MaxFlowInstance *input);
    void preflow(MaxFlowInstance *input);
    bool push(int numVertices, float **cap, int u, int sink);
    void relabel(int numVertices, float **cap, int u);
    int existsActiveNode(MaxFlowInstance *input);
};

