# LIBRARY PROJECT - ReferenceLambertW

# == Includes ==
include(GNUInstallDirs)

# === Create Library ===
add_library(ReferenceLambertW  "ReferenceW.cpp" "Interval.h" )

# === Libraries ===
find_package(PkgConfig)
pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)
target_link_libraries(ReferenceLambertW PRIVATE PkgConfig::mpfr)

# === Feature Enables ===
if (REFERENCEW_MSVC_STATIC_RUNTIME)
    set_property(TARGET ReferenceLambertW PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
set_property(TARGET ReferenceLambertW PROPERTY CXX_STANDARD 20)

# === Install ===
# paths for binaries and headers
install(TARGETS ${PROJECT_NAME}
    EXPORT "${PROJECT_NAME}Config"
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # bin
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME} # include/SomeLibrary
)

message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

# Export includes
install(FILES "../include/ReferenceW.h" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY "../src"
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        FILES_MATCHING
        PATTERN "*.h")

# CMake config name, namespace and path
install(
    EXPORT "${PROJECT_NAME}Config"
    FILE "${PROJECT_NAME}Config.cmake"
    NAMESPACE ReferenceLambertW::
    DESTINATION cmake
)