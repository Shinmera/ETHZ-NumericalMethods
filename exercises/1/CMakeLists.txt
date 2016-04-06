project(ProblemSheet1)
cmake_minimum_required(VERSION 2.6.2)

# Tested with g++ and clang
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall -O3 -Wextra -Wno-deprecated-declarations -pedantic" )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    # Use -O3 for strong optimization. We use stdc++11. Compile with Wall to get all warnings.
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall -O3 -Wextra -Wno-deprecated-declarations -pedantic" )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

# If your compiler already finds and includes Eigen you can safely remove the following two lines
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

### Problem Sheet 1
# make kron: for problem 2
set(SRCS kron.cpp)
add_executable(kron ${SRCS}) 

# make multAmin: for problem 3
set(SRCS multAmin.cpp)
add_executable(multAmin ${SRCS}) 

# make multAmin: for problem 4
set(SRCS matPow.cpp)
add_executable(matPow ${SRCS}) 

# make arrowmatvec: for problem 1
set(SRCS arrowmatvec.cpp)
add_executable(arrowmatvec ${SRCS}) 

# make arrowmatvec2: for problem 1
set(SRCS arrowmatvec2.cpp)
add_executable(arrowmatvec2 ${SRCS}) 
