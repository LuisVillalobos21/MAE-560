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

# Add your source files (replace main.cpp with your source files)
add_executable(MAE_560_Project_1 MAE_560_Project_1/main.cpp)
