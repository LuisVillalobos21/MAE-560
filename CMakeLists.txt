cmake_minimum_required(VERSION 3.10)

# Set your project name
project(MAE_560_Project_1)

# Add C++ standard (e.g., C++11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define path to third-party libraries
set(EIGEN3_INCLUDE_DIR C:\\Users\\luis2\\eigen-3.4.0\\eigen-3.4.0)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")

# Include directories
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/MAE_560_Project_1/common) 

# Add your source files
add_executable(
    Exp_Euler
    MAE_560_Project_1/Solvers/Heat_Eqn_Explicit_Euler/main.cpp
    MAE_560_Project_1/common/Mesh.cpp
    MAE_560_Project_1/common/SolutionVector.cpp
    MAE_560_Project_1/common/BoundaryConditions.cpp
    MAE_560_Project_1/common/Discretization.cpp
    MAE_560_Project_1/common/TimeIntegration.cpp
 "MAE_560_Project_1/Common/FormLinearSolveEqn.cpp")

add_executable(
    CN
    MAE_560_Project_1/Solvers/Heat_Eqn_CN/main.cpp
    MAE_560_Project_1/common/Mesh.cpp
    MAE_560_Project_1/common/SolutionVector.cpp
    MAE_560_Project_1/common/BoundaryConditions.cpp
    MAE_560_Project_1/common/Discretization.cpp
    MAE_560_Project_1/common/TimeIntegration.cpp
 "MAE_560_Project_1/Common/FormLinearSolveEqn.cpp")
