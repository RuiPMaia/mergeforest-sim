cmake_minimum_required(VERSION 3.21)

# This template attempts to be "fetch_content"-able
# so that it works well with tools like CPM or other
# manual dependency management

# Only set the cxx_standard if it is not set by someone else
if (NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
endif()

# strongly encouraged to enable this globally to avoid conflicts between
# -Wpedantic being enabled and -std=c++20 and -std=gnu++20 for example
# when compiling with PCH enabled
set(CMAKE_CXX_EXTENSIONS OFF)

# Set the project name and language
project(
  mergeforest_sim
  VERSION 1.0.0
  DESCRIPTION ""
  HOMEPAGE_URL "https://github.com/RuiPMaia/mergeforest-sim"
  LANGUAGES CXX C)

include(cmake/PreventInSourceBuilds.cmake)
include(ProjectOptions.cmake)


mergeforest_sim_setup_options()

mergeforest_sim_global_options()
include(Dependencies.cmake)
mergeforest_sim_setup_dependencies()

mergeforest_sim_local_options()

target_compile_features(mergeforest_sim_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

add_library(mergeforest_sim::mergeforest_sim_options ALIAS mergeforest_sim_options)
add_library(mergeforest_sim::mergeforest_sim_warnings ALIAS mergeforest_sim_warnings)

# Adding the src:
# add_subdirectory(mergeforest-sim)

file(GLOB src_files mergeforest-sim/*.cpp mergeforest-sim/mergeforest/*.cpp mergeforest-sim/gamma/*.cpp)

add_executable(mergeforest_sim ${src_files})
set_target_properties(mergeforest_sim PROPERTIES OUTPUT_NAME "mergeforest-sim")

target_link_libraries(
  mergeforest_sim
  PRIVATE mergeforest_sim::mergeforest_sim_options
          mergeforest_sim::mergeforest_sim_warnings)

target_link_system_libraries(	
  mergeforest_sim
  PRIVATE CLI11::CLI11
          fmt::fmt
          spdlog::spdlog
          toml11::toml11)

target_include_directories(
  mergeforest_sim
  PRIVATE "${CMAKE_SOURCE_DIR}")
