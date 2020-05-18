# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project ("TH_CS395T_final")

option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR      "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL             "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO           "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK             "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE           "Use Embree"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGUI"          ON)
option(LIBIGL_WITH_LIM              "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB           "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK            "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL           "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW      "Use GLFW"           ON)
option(LIBIGL_WITH_PNG              "Use PNG"            ON)
option(LIBIGL_WITH_PYTHON           "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN           "Use Tetgen"         ON)
option(LIBIGL_WITH_TRIANGLE         "Use Triangle"       OFF)
option(LIBIGL_WITH_VIEWER           "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_XML              "Use XML"            OFF)

option(CMAKE_USE_WIN32_THREADS_INIT "using WIN32 threads" ON)

if (MSVC)
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
endif()

include(ExternalProject)

### LOAD PACKAGES
# Threads
find_package(Threads REQUIRED)

# Eigen 3
find_package(Eigen3 REQUIRED NO_MODULE)

# LIBIGL
#include(${CMAKE_SOURCE_DIR}/cmake/igl_viewer.cmake)



if(NOT TARGET spdlog)
    # Stand-alone build
    find_package(spdlog CONFIG)
endif()



# TESTING - Added
option(PACKAGE_TESTS "Build the tests" OFF)
option(gtest_disable_pthreads "Disable uses of pthreads in gtest." ON)
if(PACKAGE_TESTS)
	add_subdirectory("${PROJECT_SOURCE_DIR}/third-party/googletest" "third-party/googletest")

	mark_as_advanced(
			BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
			gmock_build_tests gtest_build_samples gtest_build_tests
			gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
	)

	set_target_properties(gtest PROPERTIES FOLDER extern)
	set_target_properties(gtest_main PROPERTIES FOLDER extern)
	set_target_properties(gmock PROPERTIES FOLDER extern)
	set_target_properties(gmock_main PROPERTIES FOLDER extern)
endif()

# Include sub-projects.
add_subdirectory ("TH_CS395T_final")

include(GoogleTest)
if(PACKAGE_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()