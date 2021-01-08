cmake_minimum_required(VERSION 3.1)

add_compile_options(-std=c++17 -O3)
link_directories(${CMAKE_SOURCE_DIR}/ext/lib/ ${CMAKE_SOURCE_DIR}/ext/lib/matlab/ ${CMAKE_SOURCE_DIR}/ext/lib/mosek/ ${CMAKE_SOURCE_DIR}/ext/lib/scs/)

project(SPOPT CXX)

add_library(problem_data STATIC ./src/ProblemData.cpp)
if (APPLE)
    target_compile_definitions(problem_data PRIVATE __BUILD_WITH_MATLAB__)
endif()
target_include_directories(problem_data PUBLIC ./include ./ext/include/)
if (APPLE)
    target_link_libraries(problem_data libyaml-cpp.a libmat.dylib libmex.dylib libmx.dylib)
else()
    target_link_libraries(problem_data libyaml-cpp.a)
endif()

add_library(solver STATIC ./src/Solver.cpp)
target_include_directories(solver PUBLIC ./include ./ext/include/)
target_link_libraries(solver libyaml-cpp.a libaa.a blas lapack)

add_library(dual_lagrangian_solver STATIC ./src/DualLagrangianSolver.cpp)
target_include_directories(dual_lagrangian_solver PUBLIC ./include ./ext/include)
target_link_libraries(dual_lagrangian_solver solver libyaml-cpp.a)

add_library(hsde_solver STATIC ./src/HSDESolver.cpp)
target_include_directories(hsde_solver PUBLIC ./include ./ext/include)
target_link_libraries(hsde_solver solver libyaml-cpp.a)

add_library(mosek_solver STATIC ./src/MOSEKSolver.cpp)
target_include_directories(mosek_solver PUBLIC ./include ./ext/include ./ext/include/mosek)
if(APPLE)
    target_link_libraries(mosek_solver solver libmosek64.9.2.dylib libfusion64.9.2.dylib)
else()
    target_link_libraries(mosek_solver solver libmosek64.so.9.2 libfusion64.so.9.2)
endif()

add_library(scs_solver STATIC ./src/SCSSolver.cpp)
target_include_directories(scs_solver PUBLIC ./include/ ./ext/include/scs ./ext/include/linsys ./ext/include)
if(APPLE)
    target_link_libraries(scs_solver solver libscsdir.a libscsdir.dylib libscsindir.a libscsindir.dylib)
else()
    target_link_libraries(scs_solver solver libscsdir.a libscsdir.so libscsindir.a libscsindir.so)
endif()

add_executable(SPOPT ./src/SPOPT.cpp)
target_include_directories(SPOPT PRIVATE ./include ./ext/include)
target_link_libraries(SPOPT problem_data hsde_solver dual_lagrangian_solver mosek_solver scs_solver solver libyaml-cpp.a libboost_program_options-mt.a)
