include(CheckIPOSupported)

macro(enable_ipo project)
    check_ipo_supported(RESULT supported OUTPUT error)
    if (supported)
        message(STATUS "IPO enabled for project ${project}")
        set_property(TARGET ${project} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(STATUS "IPO not supported: ${error}")
    endif()
endmacro()

macro(set_arch project)
    message(STATUS "Enabled extended architecture for project ${project}")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        target_compile_options(${project} PUBLIC -march=native)
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        target_compile_options(${project} PUBLIC /arch:AVX2)
    endif()
endmacro()

macro(enable_mp)
    cmake_host_system_information(RESULT CORES QUERY NUMBER_OF_LOGICAL_CORES)
    set(CMAKE_BUILD_PARALLEL_LEVEL ${CORES})
endmacro()

macro(enable_clangtidy)
    find_program(CLANG_TIDY clang-tidy)
    set(CLANG_TIDY_COMMAND "${CLANG_TIDY}")
    if(CLANG_TIDY)
        set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_COMMAND})
        set(CMAKE_C_CLANG_TIDY   ${CLANG_TIDY_COMMAND})
    endif()
endmacro()

macro(use_static_msvc_crt)
    string(REPLACE "-D_DLL" "" CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
    string(REPLACE "-D_DLL" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
    string(REPLACE "-D_DLL" "" CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
    string(REPLACE "-D_DLL" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
    set(MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endmacro()