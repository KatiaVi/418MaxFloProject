struct Graph{
  int num_vertices;
  int num_edges;
  float** capacities; // num_vertices x num_vertices array where cij = flow capacity of edge i->j
};

int* heights; // heights of vertices 
float** flows; // x coordinate is the source, y coordinate is the sink, for edges 
float** excess; // excess flow on each edge 

// excess flow on each vertex also exists 
float *excessPerVertex;
int *d; // the labels 

int *active; 