cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
link_directories(${CMAKE_SOURCE_DIR}/ext/lib/)

project(SPOPT CXX)

add_library(problem_data STATIC ./src/ProblemData.cpp)
target_include_directories(problem_data PUBLIC ./include ./ext/include/)
target_link_libraries(problem_data libyaml-cpp.a)

add_library(solver STATIC ./src/Solver.cpp)
target_include_directories(solver PUBLIC ./include ./ext/include/)
target_link_libraries(solver libyaml-cpp.a libaa.a blas lapack)

add_library(dual_lagrangian_solver STATIC ./src/DualLagrangianSolver.cpp)
target_include_directories(dual_lagrangian_solver PUBLIC ./include ./ext/include)
target_link_libraries(dual_lagrangian_solver solver libyaml-cpp.a)

add_library(hsde_solver STATIC ./src/HSDESolver.cpp)
target_include_directories(hsde_solver PUBLIC ./include ./ext/include)
target_link_libraries(hsde_solver solver libyaml-cpp.a)

add_executable(SPOPT ./src/SPOPT.cpp)
target_include_directories(SPOPT PRIVATE ./include ./ext/include)
target_link_libraries(SPOPT problem_data hsde_solver dual_lagrangian_solver solver libyaml-cpp.a libboost_program_options-mt.a)
