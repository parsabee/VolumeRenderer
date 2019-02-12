cmake_minimum_required(VERSION 3.10)

PROJECT(volume_renderer)
SET(VTK_DIR /Users/hank/Hartree/VTK/install/lib/cmake/vtk-8.1)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(volume_renderer vol_ren.cpp)

SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(volume_renderer ${VTK_LIBRARIES})
else()
  target_link_libraries(volume_renderer vtkHybrid)
endif()


