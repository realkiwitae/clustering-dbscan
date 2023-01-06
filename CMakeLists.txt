cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
project(dbscan-clustering) #name of your project

include_directories(./)
include_directories(./headers)
add_executable(dbscan main.cpp ./utils/cvutils.cpp ./dbscan.cpp)
find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )
target_compile_options(dbscan PRIVATE -Wall -ggdb3)
target_link_libraries(dbscan rt ${OpenCV_LIBS})