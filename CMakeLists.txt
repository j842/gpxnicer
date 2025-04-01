cmake_minimum_required(VERSION 3.10)
project(gpxnicer)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(CURL REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(OpenCV REQUIRED)

# Include directories
include_directories(
    ${CURL_INCLUDE_DIRS}
    ${LIBXML2_INCLUDE_DIR}
    ${OpenCV_INCLUDE_DIRS}
)

# Add executable
add_executable(gpxnicer src/gpxnicer.cpp)

# Link libraries
target_link_libraries(gpxnicer
    ${CURL_LIBRARIES}
    ${LIBXML2_LIBRARIES}
    ${OpenCV_LIBS}
)
