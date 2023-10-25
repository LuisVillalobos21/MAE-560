# CMakeList.txt : CMake project for MAE_560_Project_1, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (MAE_560_Project_1 "MAE_560_Project_1.cpp" "MAE_560_Project_1.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET MAE_560_Project_1 PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
