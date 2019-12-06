set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED  ON)
set(CMAKE_BUILD_TYPE "RELEASE")
#set(CMAKE_BUILD_TYPE "DEBUG")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

set(CXX_COVERAGE_COMPILE_FLAGS "-g -O2 -std=c++11 -fvisibility=hidden -lpthread -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_COVERAGE_COMPILE_FLAGS}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#add_library(LibModule main.cpp)
#target_link_libraries(LibModule lib/)

add_executable(maxflow main.cpp
lib/algorithms/sequential/pushrelabel.cpp lib/algorithms/sequential/pushrelabel.h
lib/algorithms/parallel/pushrelabel_parallel.cpp lib/algorithms/parallel/pushrelabel_parallel.h
lib/algorithms/sequential/dinics.cpp lib/algorithms/sequential/dinics.h
lib/algorithms/parallel/dinics_parallel.cpp lib/algorithms/parallel/dinics_parallel.h
)
# Problem here may not be able to compile the other parallel algo because reach limit on the number
# of executables attached to the maxflow command?
# lib/algorithms/sequential/dinics.cpp lib/algorithms/sequential/dinics.h
# lib/algorithms/parallel/dinics_parallel.cpp lib/algorithms/parallel/dinics_parallel.h


target_include_directories(maxflow PRIVATE lib)


