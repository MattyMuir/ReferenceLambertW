# GLOBAL - ReferenceLambertW

# === Includes ===
cmake_policy(SET CMP0069 NEW)
include(cmake/macros.cmake)

# === Configuration ===
cmake_minimum_required(VERSION 3.21)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
set(WARNINGS OFF CACHE BOOL "Enables clang-tidy and warnings on clang-like compilers")
if (WARNINGS)
	add_compile_options("-Wall")
    add_compile_options("-Wextra")
    add_compile_options("-Wpedantic")
    enable_clangtidy()
endif()

set(REFERENCEW_NO_ITERATOR_DEBUGGING OFF CACHE BOOL "Disable iterator debugging")
if (REFERENCEW_NO_ITERATOR_DEBUGGING)
    add_compile_definitions(_HAS_ITERATOR_DEBUGGING=0)
endif()

set(REFERENCEW_MSVC_STATIC_RUNTIME, OFF CACHE BOOL "Use static C runtime library when taregting MSVC ABI")
if (REFERENCEW_MSVC_STATIC_RUNTIME)
    use_static_msvc_crt()
endif()

# === Projects ===
project("ReferenceLambertW")
add_subdirectory("src")
add_subdirectory("examples")
add_subdirectory("bench")

# === Tests ===
enable_testing()
set(CTEST_OUTPUT_ON_FAILURE ON)

add_subdirectory("tests")

add_test(NAME FloatW0 COMMAND tests 0)
add_test(NAME FloatWm1 COMMAND tests 1)
add_test(NAME DoubleW0 COMMAND tests 2)
add_test(NAME DoubleWm1 COMMAND tests 3)
add_test(NAME FloatW0Exhaustive COMMAND tests 4)
add_test(NAME FloatWm1Exhaustive COMMAND tests 5)