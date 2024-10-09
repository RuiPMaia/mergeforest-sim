include(cmake/SystemLink.cmake)
include(CMakeDependentOption)
include(CheckCXXCompilerFlag)


macro(mergeforest_sim_supports_sanitizers)
  if((CMAKE_CXX_COMPILER_ID MATCHES ".*Clang.*" OR CMAKE_CXX_COMPILER_ID MATCHES ".*GNU.*") AND NOT WIN32)
    set(SUPPORTS_UBSAN ON)
  else()
    set(SUPPORTS_UBSAN OFF)
  endif()

  if((CMAKE_CXX_COMPILER_ID MATCHES ".*Clang.*" OR CMAKE_CXX_COMPILER_ID MATCHES ".*GNU.*") AND WIN32)
    set(SUPPORTS_ASAN OFF)
  else()
    set(SUPPORTS_ASAN ON)
  endif()
endmacro()

macro(mergeforest_sim_setup_options)
  option(mergeforest_sim_ENABLE_HARDENING "Enable hardening" OFF)
  option(mergeforest_sim_ENABLE_COVERAGE "Enable coverage reporting" OFF)
  cmake_dependent_option(
    mergeforest_sim_ENABLE_GLOBAL_HARDENING
    "Attempt to push hardening options to built dependencies"
    OFF
    mergeforest_sim_ENABLE_HARDENING
    OFF)

  mergeforest_sim_supports_sanitizers()

  option(mergeforest_sim_ENABLE_IPO "Enable IPO/LTO" OFF)
  option(mergeforest_sim_WARNINGS_AS_ERRORS "Treat Warnings As Errors" ON)
  option(mergeforest_sim_ENABLE_USER_LINKER "Enable user-selected linker" OFF)
  option(mergeforest_sim_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" OFF)
  option(mergeforest_sim_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
  option(mergeforest_sim_ENABLE_SANITIZER_UNDEFINED "Enable undefined sanitizer" OFF)
  option(mergeforest_sim_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
  option(mergeforest_sim_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)
  option(mergeforest_sim_ENABLE_UNITY_BUILD "Enable unity builds" OFF)
  option(mergeforest_sim_ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)
  option(mergeforest_sim_ENABLE_CPPCHECK "Enable cpp-check analysis" OFF)
  option(mergeforest_sim_ENABLE_PCH "Enable precompiled headers" OFF)
  option(mergeforest_sim_ENABLE_CACHE "Enable ccache" OFF)

  if(NOT PROJECT_IS_TOP_LEVEL)
    mark_as_advanced(
      mergeforest_sim_ENABLE_IPO
      mergeforest_sim_WARNINGS_AS_ERRORS
      mergeforest_sim_ENABLE_USER_LINKER
      mergeforest_sim_ENABLE_SANITIZER_ADDRESS
      mergeforest_sim_ENABLE_SANITIZER_LEAK
      mergeforest_sim_ENABLE_SANITIZER_UNDEFINED
      mergeforest_sim_ENABLE_SANITIZER_THREAD
      mergeforest_sim_ENABLE_SANITIZER_MEMORY
      mergeforest_sim_ENABLE_UNITY_BUILD
      mergeforest_sim_ENABLE_CLANG_TIDY
      mergeforest_sim_ENABLE_CPPCHECK
      mergeforest_sim_ENABLE_COVERAGE
      mergeforest_sim_ENABLE_PCH
      mergeforest_sim_ENABLE_CACHE)
  endif()
  
endmacro()

macro(mergeforest_sim_global_options)
  if(mergeforest_sim_ENABLE_IPO)
    include(cmake/InterproceduralOptimization.cmake)
    mergeforest_sim_enable_ipo()
  endif()

  mergeforest_sim_supports_sanitizers()

  if(mergeforest_sim_ENABLE_HARDENING AND mergeforest_sim_ENABLE_GLOBAL_HARDENING)
    include(cmake/Hardening.cmake)
    if(NOT SUPPORTS_UBSAN 
       OR mergeforest_sim_ENABLE_SANITIZER_UNDEFINED
       OR mergeforest_sim_ENABLE_SANITIZER_ADDRESS
       OR mergeforest_sim_ENABLE_SANITIZER_THREAD
       OR mergeforest_sim_ENABLE_SANITIZER_LEAK)
      set(ENABLE_UBSAN_MINIMAL_RUNTIME FALSE)
    else()
      set(ENABLE_UBSAN_MINIMAL_RUNTIME TRUE)
    endif()
    message("${mergeforest_sim_ENABLE_HARDENING} ${ENABLE_UBSAN_MINIMAL_RUNTIME} ${mergeforest_sim_ENABLE_SANITIZER_UNDEFINED}")
    mergeforest_sim_enable_hardening(mergeforest_sim_options ON ${ENABLE_UBSAN_MINIMAL_RUNTIME})
  endif()
endmacro()

macro(mergeforest_sim_local_options)
  if(PROJECT_IS_TOP_LEVEL)
    include(cmake/StandardProjectSettings.cmake)
  endif()

  add_library(mergeforest_sim_warnings INTERFACE)
  add_library(mergeforest_sim_options INTERFACE)

  include(cmake/CompilerWarnings.cmake)
  mergeforest_sim_set_project_warnings(
    mergeforest_sim_warnings
    ${mergeforest_sim_WARNINGS_AS_ERRORS}
    ""
    ""
    ""
    "")

  if(mergeforest_sim_ENABLE_USER_LINKER)
    include(cmake/Linker.cmake)
    configure_linker(mergeforest_sim_options)
  endif()

  include(cmake/Sanitizers.cmake)
  mergeforest_sim_enable_sanitizers(
    mergeforest_sim_options
    ${mergeforest_sim_ENABLE_SANITIZER_ADDRESS}
    ${mergeforest_sim_ENABLE_SANITIZER_LEAK}
    ${mergeforest_sim_ENABLE_SANITIZER_UNDEFINED}
    ${mergeforest_sim_ENABLE_SANITIZER_THREAD}
    ${mergeforest_sim_ENABLE_SANITIZER_MEMORY})

  set_target_properties(mergeforest_sim_options PROPERTIES UNITY_BUILD ${mergeforest_sim_ENABLE_UNITY_BUILD})

  if(mergeforest_sim_ENABLE_PCH)
    target_precompile_headers(
      mergeforest_sim_options
      INTERFACE
      <vector>
      <string>
      <utility>)
  endif()

  if(mergeforest_sim_ENABLE_CACHE)
    include(cmake/Cache.cmake)
    mergeforest_sim_enable_cache()
  endif()

  include(cmake/StaticAnalyzers.cmake)
  if(mergeforest_sim_ENABLE_CLANG_TIDY)
    mergeforest_sim_enable_clang_tidy(mergeforest_sim_options ${mergeforest_sim_WARNINGS_AS_ERRORS})
  endif()

  if(mergeforest_sim_ENABLE_CPPCHECK)
    mergeforest_sim_enable_cppcheck(${mergeforest_sim_WARNINGS_AS_ERRORS} "" # override cppcheck options
    )
  endif()

  if(mergeforest_sim_ENABLE_COVERAGE)
    include(cmake/Tests.cmake)
    mergeforest_sim_enable_coverage(mergeforest_sim_options)
  endif()

  if(mergeforest_sim_WARNINGS_AS_ERRORS)
    check_cxx_compiler_flag("-Wl,--fatal-warnings" LINKER_FATAL_WARNINGS)
    if(LINKER_FATAL_WARNINGS)
      # This is not working consistently, so disabling for now
      # target_link_options(mergeforest_sim_options INTERFACE -Wl,--fatal-warnings)
    endif()
  endif()

  if(mergeforest_sim_ENABLE_HARDENING AND NOT mergeforest_sim_ENABLE_GLOBAL_HARDENING)
    include(cmake/Hardening.cmake)
    if(NOT SUPPORTS_UBSAN 
       OR mergeforest_sim_ENABLE_SANITIZER_UNDEFINED
       OR mergeforest_sim_ENABLE_SANITIZER_ADDRESS
       OR mergeforest_sim_ENABLE_SANITIZER_THREAD
       OR mergeforest_sim_ENABLE_SANITIZER_LEAK)
      set(ENABLE_UBSAN_MINIMAL_RUNTIME FALSE)
    else()
      set(ENABLE_UBSAN_MINIMAL_RUNTIME TRUE)
    endif()
    mergeforest_sim_enable_hardening(mergeforest_sim_options OFF ${ENABLE_UBSAN_MINIMAL_RUNTIME})
  endif()

endmacro()
