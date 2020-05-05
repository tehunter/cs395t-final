# CMakeList.txt : CMake project for TH_CS395T_final, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
add_executable (TH_CS395T_final "TH_CS395T_final.cpp" "TH_CS395T_final.h")

# TODO: Add tests and install targets if needed.

target_link_libraries(TH_CS395T_final PRIVATE Eigen3::Eigen)

