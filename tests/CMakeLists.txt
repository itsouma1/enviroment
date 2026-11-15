cmake_minimum_required(VERSION 3.16)
project(climate_tests CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(COVERAGE "Instrument the build for gcov coverage" OFF)
option(BUILD_FUZZERS "Build libFuzzer targets (requires clang)" OFF)

if(COVERAGE)
  add_compile_options(--coverage -O0 -g)
  add_link_options(--coverage)
endif()

file(GLOB CLIMATE_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../src/*.cpp")
list(FILTER CLIMATE_SOURCES EXCLUDE REGEX ".*/main\\.cpp$")
add_library(climate_core STATIC ${CLIMATE_SOURCES})
target_include_directories(climate_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../src")

enable_testing()

# 1. Self-checking C++ suite (the same binary `make test` runs).
add_executable(unit_suite test_all.cpp)
target_link_libraries(unit_suite PRIVATE climate_core)
add_test(NAME unit COMMAND unit_suite)

# 2. Floating-point differential / determinism harness.
add_executable(fp_diff fp_diff.cpp)
target_link_libraries(fp_diff PRIVATE climate_core)
add_test(NAME fp_diff COMMAND fp_diff)

# 3. GoogleTest, when the dependency is present.
set(GTEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../.tools/deps/googletest")
if(EXISTS "${GTEST_DIR}/CMakeLists.txt")
  add_subdirectory("${GTEST_DIR}" googletest EXCLUDE_FROM_ALL)
  add_executable(gtest_suite gtest_suite.cpp)
  target_link_libraries(gtest_suite PRIVATE climate_core gtest gtest_main)
  add_test(NAME google_test COMMAND gtest_suite)
else()
  message(STATUS "googletest not vendored: run tools/bootstrap_tests.sh to enable it")
endif()

# 4. RapidCheck property tests, when the dependency is present.
set(RC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../.tools/deps/rapidcheck")
if(EXISTS "${RC_DIR}/CMakeLists.txt")
  set(RC_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
  add_subdirectory("${RC_DIR}" rapidcheck EXCLUDE_FROM_ALL)
  add_executable(rapidcheck_suite rapidcheck_suite.cpp)
  target_link_libraries(rapidcheck_suite PRIVATE climate_core rapidcheck)
  add_test(NAME rapidcheck COMMAND rapidcheck_suite)
else()
  message(STATUS "rapidcheck not vendored: run tools/bootstrap_tests.sh to enable it")
endif()

# 5. libFuzzer targets. Clang-only; GCC has no -fsanitize=fuzzer.
if(BUILD_FUZZERS)
  add_executable(fuzz_lab fuzz_lab.cpp)
  target_link_libraries(fuzz_lab PRIVATE climate_core)
  target_compile_options(fuzz_lab PRIVATE -fsanitize=fuzzer,address)
  target_link_options(fuzz_lab PRIVATE -fsanitize=fuzzer,address)
  add_test(NAME fuzz_lab_smoke COMMAND fuzz_lab -runs=2000)
endif()

# 6. Coverage report target.
if(COVERAGE)
  add_custom_target(coverage
    COMMAND gcovr --root "${CMAKE_CURRENT_SOURCE_DIR}/.." --filter ".*/src/.*"
            --print-summary --txt coverage.txt
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    COMMENT "gcovr coverage report -> coverage.txt")
endif()
