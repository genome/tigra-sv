cmake_minimum_required(VERSION 2.8)

project(example)

set(CMAKE_INSTALL_PREFIX "/usr")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/build-common/cmake")
include(TestHelper)
include(VersionHelper)
include(ProjectHelper)

include_directories("src/lib")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE release CACHE STRING 
        "Options: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
    message("No CMAKE_BUILD_TYPE specified, defaulting to release")
endif ()

if (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_C_FLAGS "-std=c99")
    set(CMAKE_CXX_FLAGS "-Wall -ansi -pedantic")
endif ()


# unit tests
find_package(GTest)
if(GTEST_FOUND)
    message("Google Test framework found, building unit tests")
    enable_testing(true)
    include_directories(src)
    add_projects(test/lib)
else()
    message("Google Test framework not found, no tests will be built")
    message("GTest is available at http://code.google.com/p/googletest/")
    message("Ubuntu users can likely sudo apt-get install gtest")
endif()

# main project
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
add_projects(src/lib)
add_projects(src/exe)

include(PackageInfo.cmake)
