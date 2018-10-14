cmake_minimum_required(VERSION 3.9)
project( polyhedron_gravity CXX )

# set(CMAKE_CXX_FLAGS "-lpthread -pthread")
# set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
# Add -Weffc++ for all sorts of good errors to fix
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3")

################################################################################
# GLOBAL SETUP
################################################################################
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Find OpenMP
set(ENV{OMP_NUM_THREADS} 4)
set(ENV{OMP_NESTED} TRUE)
find_package(OpenMP REQUIRED)

# CGAL Setup
set(ENV{CGAL_DIR} /usr/local/Cellar/cgal/)
find_package(CGAL REQUIRED COMPONENTS)
if (NOT CGAL_FOUND)
message(STATUS "This project requires the CGAL Library and will not compile")
return()
endif()
# include helper file
include( ${CGAL_USE_FILE} )
message(STATUS "CGAL Lib: ${CGAL_LIBRARIES}")

# Eigen setup
find_package(Eigen3 REQUIRED)
if (NOT Eigen3_FOUND)
message(STATUS "This project requires Eigen3 and will not be compiled.")
return()
endif()
include(${EIGEN3_USE_FILE})

# find_package(Threads REQUIRED)

# libigl setup
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           ON)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          OFF)
option(LIBIGL_WITH_PNG               "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  OFF)
option(LIBIGL_WITH_XML               "Use XML"            OFF)

find_package(LIBIGL QUIET)

################################################################################
# GOOGLE TESTING
################################################################################
add_subdirectory(${PROJECT_SOURCE_DIR}/extern/googletest/googletest)
ENABLE_TESTING()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

################################################################################
# BUILD TARGETS
################################################################################
include_directories(${PROJECT_SOURCE_DIR}/include)
add_executable(main src/main.cpp src/cube.cpp)
# target_link_libraries(main
################################################################################
# Testing example
################################################################################
ENABLE_TESTING()

set(test_all_src 
    tests/test_all.cpp
    tests/test_polyhedron.cpp
    tests/test_mesh.cpp
    tests/test_loader.cpp
    tests/test_libigl.cpp
    tests/test_cgal.cpp
    src/cgal.cpp
    src/mesh.cpp
    src/loader.cpp
    src/polyhedron.cpp)

add_executable(test_all ${test_all_src})
target_link_libraries(test_all gtest igl::core ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
add_test(NAME test_all COMMAND test_all
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
