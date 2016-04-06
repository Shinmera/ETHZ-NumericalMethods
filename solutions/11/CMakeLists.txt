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

### Problem Sheet 10
# make quadinf: problem 2
set(SRCS quadinf.cpp)
add_executable(quadinf ${SRCS})

# make piecewise_hermite_interpolation: problem 1
set(SRCS legendre.cpp)
add_executable(legendre ${SRCS})

### Problem Sheet 11
# make laserquad: problem 1 (compile and link)
set(SRCS quadsingint.cpp gauleg.cpp)
add_executable(quadsingint ${SRCS})
# Equivalent to:
# g++ [compiler opts.] -c quadsingint.cpp
# g++ [compiler opts.] -c laserquad.cpp
# g++ [compiler opts.]  quadsingint.o gauleg.o -o quadsingint

# make laserquad: problem 5 (compile and link)
set(SRCS laserquad.cpp gauleg.cpp)
add_executable(laserquad ${SRCS})


set(SRCS quadU.cpp)
add_executable(quadU ${SRCS})
