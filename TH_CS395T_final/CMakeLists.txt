﻿# CMakeList.txt : CMake project for TH_CS395T_final, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

option(PACKAGE_VISUALIZER              "Build a real-time visualizer"            ON)

# Add source to this project's executable.
# SIMULATION
add_library (FluidSimulation STATIC "lib/EulerState.h" "lib/EulerState.cpp" "lib/EulerSimulation.h" "lib/EulerSimulation.cpp")
if (MSVC)
    set_property(TARGET FluidSimulation PROPERTY
    MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
target_link_libraries(FluidSimulation PUBLIC Eigen3::Eigen)
target_include_directories(FluidSimulation PUBLIC ${GTEST_INCLUDE_DIRS})

# EXECUTABLE
add_executable (TH_CS395T_final "TH_CS395T_final.cpp" "TH_CS395T_final.h")
if (MSVC)
    set_property(TARGET TH_CS395T_final PROPERTY
    MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

# VISUALIZER
if (PACKAGE_VISUALIZER)
    find_package(LIBIGL REQUIRED)

    add_library (FluidVisualizer STATIC "vis/ParticleVisualizer.h" "vis/ParticleVisualizer.cpp" "vis/FluidScenarios.h" "vis/FluidScenarios.cpp")
    if (MSVC)
        set_property(TARGET FluidVisualizer PROPERTY
        MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()
    target_link_libraries(FluidVisualizer PUBLIC Eigen3::Eigen)

    if(SPDLOG_FOUND)
        target_link_libraries(FluidVisualizer PUBLIC spdlog::spdlog)
    else()
        target_include_directories(FluidVisualizer PUBLIC ${CMAKE_SOURCE_DIR}/third-party/spdlog/include)
    endif()

    target_link_libraries(FluidVisualizer PUBLIC igl::core igl::opengl_glfw igl::opengl_glfw_imgui)
    target_link_libraries(TH_CS395T_final PRIVATE igl::core igl::opengl_glfw)
    target_link_libraries(TH_CS395T_final PRIVATE FluidVisualizer)
	# add_subdirectory(sln)
endif()

target_link_libraries(TH_CS395T_final PRIVATE FluidSimulation)

# TODO: Add tests and install targets if needed.

target_link_libraries(TH_CS395T_final PRIVATE Eigen3::Eigen)

if(SPDLOG_FOUND)
    target_link_libraries(FluidSimulation PUBLIC spdlog::spdlog)
else()
    target_include_directories(FluidSimulation PUBLIC ${CMAKE_SOURCE_DIR}/third-party/spdlog/include)
endif()


