cmake_minimum_required(VERSION 3.0.0)
project(ephemeris_fitting VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
    "./third_party/eigen-3.4.0/"
    "./third_party/nlohmann_json-3.11.2/"
    "./src/"
)

file(GLOB SRC_LIST "src/*.cc")

add_executable(ephemeris_fitting ${SRC_LIST})

# Test
add_executable(test_gnsstime test/test_gnsstime.cc src/gnsstime.cc)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)