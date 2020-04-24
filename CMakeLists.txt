project(Glints)
cmake_minimum_required(VERSION 3.9)
set(CMAKE_CXX_STANDARD 11)


add_executable(glints src/pathtracer.cpp
src/brdf.cpp
src/Intersection.cpp
src/mesh.cpp
src/position_normal_bvh.cpp
src/RayVec.cpp
src/sphere.cpp
src/triangle.cpp
src/SrgbTransform.cpp
src/Texture.cpp)

target_include_directories(glints PUBLIC src)

find_package(GSL REQUIRED)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
target_link_libraries(glints PUBLIC OpenMP::OpenMP_CXX GSL::gsl GSL::gslcblas)
else()
target_link_libraries(glints PUBLIC  GSL::gsl GSL::gslcblas)
endif()
