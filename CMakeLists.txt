cmake_minimum_required(VERSION 3.16)
project(v4zip VERSION 2.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT UNIX)
    message(FATAL_ERROR "The v4zip command-line program requires a POSIX system")
endif()
if(NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(FATAL_ERROR "v4zip requires a 64-bit build")
endif()

# Dependencies
find_path(ZSTD_INCLUDE_DIR NAMES zstd.h)
find_library(ZSTD_LIBRARY NAMES zstd)
if(NOT ZSTD_INCLUDE_DIR OR NOT ZSTD_LIBRARY)
    message(FATAL_ERROR "Install the libzstd development headers and library")
endif()

# Main executable
add_executable(v4zip
    src/main.cpp
    src/dna.cpp
    src/archive.cpp
    src/archive2.cpp
    src/side_streams.cpp
)

# Include directories and libraries
target_include_directories(v4zip PRIVATE include "${ZSTD_INCLUDE_DIR}")
target_link_libraries(v4zip PRIVATE "${ZSTD_LIBRARY}")
target_compile_definitions(v4zip PRIVATE V4ZIP_VERSION="${PROJECT_VERSION}")

# Compiler flags; Release optimization is supplied by CMake.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(v4zip PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Install
install(TARGETS v4zip RUNTIME DESTINATION bin)
