set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE "RELEASE")
#set(CMAKE_BUILD_TYPE "DEBUG")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

set(CXX_COVERAGE_COMPILE_FLAGS "-g -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_COVERAGE_COMPILE_FLAGS}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#add_library(LibModule main.cpp)
#target_link_libraries(LibModule lib/)

add_executable(maxflow main.cpp lib/algorithms/sequential/genetic.cpp lib/algorithms/sequential/genetic.h lib/algorithms/sequential/pushrelabel.cpp lib/algorithms/sequential/pushrelabel.h)
#target_link_libraries(maxflow LibModule)

target_include_directories(maxflow PRIVATE lib)


