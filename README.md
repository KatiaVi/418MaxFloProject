# 418MaxFloProject
Parallel and Sequential implementations of Dinic's, Push-Relabel and Genetic algorithms for solving the graph maxflow problem.

### How to Build Project:
```
mkdir build && cd build
cmake .. && make
cd ..
```

### How to run (runs main.cpp):

In the directory where you created the build directory:
```
./maxflow
```

Todo (for Katia):
Sequential Genetic Algorithm:
4. Create a small test case to run the genetic sequential code on
5. (Discuss with Wynne) Re-evaluate the design choices made (class types, functions, private/public variables)
6. Try to hook up a larger test case (from our actual workload)
7. Run Sequential algorithm for larger test case + debug any issues
8. Create tester/modify Makefile? to add timer to evaluate speed of algo
8. Tune the sequential algorithm for performance

