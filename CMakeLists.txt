# ===============================================================================
# Set the minimum required version of CMake                                    ||
# ===============================================================================
#set cmake version
cmake_minimum_required(VERSION 3.8)
#get version from file version.txt
file(STRINGS "${CMAKE_CURRENT_LIST_DIR}/version.txt" projectVersion)
#set project name and version
project(safetics_wasm VERSION ${projectVersion})
message("${PROJECT_NAME} Project version: ${projectVersion}")

# ===============================================================================
# Global internal Configuration                                                ||
# The output directory is set to the build directory                           ||
# ===============================================================================
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

## Manual configuration
# ================================================================================
# Emscripten Configuration                                                      ||
# Emcc is the compiler for WebAssembly (C target is emcc/ C++ target is em++)   ||
# ================================================================================
# if (DEFINED EMSCRIPTEN)
#     set(CMAKE_C_COMPILER "emcc")
#     set(CMAKE_CXX_COMPILER "em++")
# endif ()

# ================================================================================
# Require Module                                                                ||
# 1. CMake module path                                                          ||
# ================================================================================
# find_package (Eigen REQUIRED NO_MODULE)
# include_directories(/opt/homebrew/Cellar/eigen/3.4.0_1)
# ================================================================================
# Dependencies (Todo)                                                           ||
# 1. external libraries                                                         ||
# ================================================================================
# if (DEFINED EMSCRIPTEN)
#     if ("${CMAKE_SYSTEM}" MATCHES "Linux" OR "${CMAKE_SYSTEM}" MATCHES "Darwin-*")
#         find_package(OpenGL REQUIRED)
#     endif ()
# else ()
#     if (THREEPP_USE_EXTERNAL_GLFW)
#         find_package(glfw3 CONFIG REQUIRED)
#         if (NOT TARGET "glfw::glfw" AND TARGET "glfw")
#             add_library(glfw::glfw ALIAS glfw)
#         endif ()
#     endif ()
# endif ()

# ================================================================================
# Include directories                                                           ||
# ================================================================================
file(GLOB_RECURSE SOURCES "src/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_FLAGS "-I src/external/Eigen --no-entry -s EXPORT_ES6=1 -s WASM=1 -s ENVIRONMENT='web' -s EXPORT_NAME='safatics' -s USE_ES6_IMPORT_META=0 -s ASSERTIONS -O3 -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"cwrap\"]'"
)
