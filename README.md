# CS395T Final Project

## General Build Instruction

This is a standard cmake project, with a few git submodules.
I used Visual Studio woth CMake GUI on Windows to compile and build. To build this project in Linux, run:
```
git clone --recursive https://github.com/tehunter/cs395t-final psim
cd TH_CS395T_final
mkdir build
cd build
cmake ..
make -j `nproc`
```

Ensure that CMake can find your Eigen and Libigl libraries.
Also check that the PACKAGE_VISUALIZER variable is on to compile the visualizer.
Without the visualizer, you can still run tests if the PACKAGE_TESTS variable is on.

## Brief introduction to the code structure

This repository organizes the source code in the following manner:
* `TH_CS395T_final/`: projectsource code
* `TH_CS395T_final/lib/`: core physical simulation code
* `TH_CS395T_final/vis/`: visualization of the current physical system, with
    [libigl](https://libigl.github.io/)
* `tests`: `gtest` unit tests
* `third-party/`: git submodules
