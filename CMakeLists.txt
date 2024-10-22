# GLOBAL - ReferenceLambertW

# === Includes ===
include(macros.cmake)

# === Configuration ===
cmake_minimum_required(VERSION 3.21)

# Options
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