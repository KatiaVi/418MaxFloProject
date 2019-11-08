float** flows; // x coordinate is the source, y coordinate is the sink, for edges 

// excess flow on each vertex also exists 
float *excessPerVertex;
int *d; // the labels 

int *active; 