cmake_minimum_required(VERSION 3.22)

project(GuerrieriZ)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# compilo e collego l'esempio alla libreria
add_executable(test test.cpp)