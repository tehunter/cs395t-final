# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

project ("TH_CS395T_final")

if (MSVC)
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
endif()

include(ExternalProject)

### LOAD PACKAGES
# Threads
find_package(Threads REQUIRED)

# Eigen 3
find_package(Eigen3 REQUIRED NO_MODULE)


if(NOT TARGET spdlog)
    # Stand-alone build
    find_package(spdlog CONFIG)
endif()



# TESTING - Added
option(PACKAGE_TESTS "Build the tests" OFF)
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