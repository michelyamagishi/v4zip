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

if(BUILD_TESTS)
    enable_testing()

    add_executable(test_v4group tests/test_v4group.cpp)
    add_executable(test_arithmetic_coder tests/test_arithmetic_coder.cpp)
    add_executable(test_roundtrip tests/test_roundtrip.cpp)
    add_executable(test_fastq_io tests/test_fastq_io.cpp)
    add_executable(test_quality_model tests/test_quality_model.cpp)
    add_executable(test_fastq_roundtrip tests/test_fastq_roundtrip.cpp)
    add_executable(test_robustness tests/test_robustness.cpp)

    add_test(NAME V4Group COMMAND test_v4group)
    add_test(NAME ArithmeticCoder COMMAND test_arithmetic_coder)
    add_test(NAME Roundtrip COMMAND test_roundtrip)
    add_test(NAME FastqIO COMMAND test_fastq_io)
    add_test(NAME QualityModel COMMAND test_quality_model)
    add_test(NAME FastqRoundtrip COMMAND test_fastq_roundtrip)
    add_test(NAME Robustness COMMAND test_robustness)
endif()

# Install
install(TARGETS v4zip DESTINATION bin)
install(DIRECTORY include/v4zip DESTINATION include)
