# 418MaxFloProject
Parallel and Sequential implementations of Dinic's, Push-Relabel and Genetic algorithms for solving the graph maxflow problem.

### How to Build Project:
```
mkdir build && cd build
cmake .. && make
cd ..
```

### How to run (runs main.cpp):

#### To Run All Maxflow Algorithms
In the directory where you created the build directory:
```
./maxflow
```

#### To Run a Specific Maxflow Algorithm
In the directory where you created the build directory:
(i.e. sequentual dinics)
```
./maxflow -seq dinics
```
(i.e. parallel push relabel)
```
./maxflow -par pushrelabel
```




