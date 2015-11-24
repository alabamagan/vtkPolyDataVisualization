cmake_minimum_required(VERSION 2.8)
project(vtkVisualization)
set(VTK_DIR "/opt/VTK-6.3-GPU/lib/cmake/vtk-6.3")
#set(INCLUDE_DIR "/opt/VTK-6.3-GPU/include")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11")

set(CMAKE_BUILD_TYPE "release")

set(SOURCE_FILES
    main.cpp)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(vtkVisualization ${SOURCE_FILES})
target_link_libraries(vtkVisualization ${VTK_LIBRARIES})
