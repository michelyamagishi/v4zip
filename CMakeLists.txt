cmake_minimum_required(VERSION 3.14)
project(v4zip VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -fsanitize=address,undefined")
endif()

# Include directories
include_directories(include)

# Main executable
add_executable(v4zip src/main.cpp)

# Tests
option(BUILD_TESTS "Build unit tests" ON)


# Install
install(TARGETS v4zip DESTINATION bin)
install(DIRECTORY include/v4zip DESTINATION include)
