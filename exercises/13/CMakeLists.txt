project(ProblemSheet7)
cmake_minimum_required(VERSION 2.6.2)

# Tested with g++ and clang
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall -O3 -Wextra -Wno-deprecated-declarations -pedantic" )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    # Use -O3 for strong optimization. We use stdc++11. Compile with Wall to get all warnings.
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall -Wextra -Wno-deprecated-declarations -pedantic" )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

# If your compiler already finds and includes Eigen you can safely remove the following two lines
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

### Problem Sheet 12
# make rk3prey: problem 1
set(SRCS rk3prey.cpp)
add_executable(rk3prey ${SRCS})

# make system: problem 2
set(SRCS system.cpp)
add_executable(system ${SRCS})

# make rk3prey: problem 3
set(SRCS taylorprey.cpp)
add_executable(taylorprey ${SRCS})

### Problem Sheet 13
# make matrix_ode: problem 3
set(SRCS matrix_ode.cpp)
add_executable(matrix_ode ${SRCS})

# make stabrk: problem 3
set(SRCS stabrk.cpp)
add_executable(stabrk ${SRCS})

# make LV: problem 1
set(SRCS LV.cpp)
add_executable(LV ${SRCS})
