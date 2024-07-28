cmake_minimum_required(VERSION 3.24)
project(CompMod)

set(CMAKE_CXX_STANDARD 17)

add_executable(CompMod main.cpp solver_fvm.cpp solver_fvm.h verifier.cpp verifier.h)
