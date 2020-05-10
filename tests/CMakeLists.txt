macro(package_add_test TESTNAME)
    # create an exectuable in which the tests will be stored
    add_executable(${TESTNAME} ${ARGN})
    if (MSVC)
        set_property(TARGET ${TESTNAME} PROPERTY
        MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()

    # link the Google test infrastructure, mocking library, and a default main fuction to
    # the test executable.  Remove g_test_main if writing your own main function.

    #target_include_directories(${TESTNAME} lib)
    #target_link_directories(${TESTNAME} lib)
    
    target_link_libraries(${TESTNAME} gtest gmock gtest_main)
    target_link_libraries(${TESTNAME} FluidSimulation)
    if(SPDLOG_FOUND)
        target_link_libraries(${TESTNAME} PRIVATE spdlog::spdlog)
    else()
        target_include_directories(${TESTNAME} PRIVATE ${CMAKE_SOURCE_DIR}/third-party/spdlog/include)
    endif()


    # gtest_discover_tests replaces gtest_add_tests,
    # see https://cmake.org/cmake/help/v3.10/module/GoogleTest.html for more options to pass to it
    gtest_discover_tests(${TESTNAME}
        # set a working directory so your project root so that you can find test data via paths relative to the project root
        WORKING_DIRECTORY ${PROJECT_DIR}
        PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
    )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)

    add_test(MyTests ${TESTNAME})
endmacro()

package_add_test(EulerState test_EulerState.cpp)  
package_add_test(EulerSimulation test_EulerSimulation.cpp)  