# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

project ("TH_CS395T_final")

include(ExternalProject)

### LOAD PACKAGES
# Threads
find_package(Threads REQUIRED)

# Eigen 3
find_package(Eigen3 REQUIRED NO_MODULE)

# Include sub-projects.
add_subdirectory ("TH_CS395T_final")
