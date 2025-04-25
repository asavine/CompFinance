cmake_minimum_required(VERSION 3.30)

project(AAD_book VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD              20)
set(CMAKE_CXX_STANDARD_REQUIRED     ON)
set(CMAKE_CXX_EXTENSIONS            OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS   ON)

if(WIN32)
    message(STATUS "Windows detected")
    message(STATUS "Toolset: ${CMAKE_GENERATOR_TOOLSET}")
    set(dlloader_include_dir_platform ${CMAKE_SOURCE_DIR}/dl_loader/windows)
    add_compile_definitions(WIN_EXPORT)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG /openmp /Oi")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_RELEASE} -D_DEBUG /Od")
endif()

if(UNIX)
    message(STATUS "Unix system detected")
    set(dlloader_include_dir_platform ${CMAKE_SOURCE_DIR}/dl_loader/linux)
endif()

set(dlloader_include_dir ${CMAKE_SOURCE_DIR}/dl_loader)


MESSAGE(STATUS "the source directory is: ${CMAKE_SOURCE_DIR}")
MESSAGE(STATUS "the binary directory is: ${CMAKE_BINARY_DIR}")

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

if(ENABLE_WARNINGS)
    include(Warnings)
endif()

add_subdirectory(libraries)
