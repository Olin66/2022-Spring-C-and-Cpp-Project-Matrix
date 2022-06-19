cmake_minimum_required(VERSION 3.10.2)
project(matrix)

set(CMAKE_CXX_STANDARD 17)
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
add_executable(matrix main.cpp test1.cpp test2.cpp test3.cpp)
target_link_libraries(matrix ${OpenCV_LIBS})