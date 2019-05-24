cmake_minimum_required(VERSION 2.8.11)

project(GPStudio)

set(DEBUG_MODE 0)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

#add_definitions(
#    -Wfatal-errors
#    )


if(MSVC)
  add_definitions ("-DBOOST_ALL_NO_LIB -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS -DNOMINMAX")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj /EHsc /fp:fast /wd4800 /wd4521 /wd4251 /wd4275 /wd4305 /wd4355 /wd4819 /wd4250 /wd4267 /wd4005 /wd4819")
  # no global code optimization to speed up compile
  #SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /GL /wd4800")
  #SET(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /LTCG")
  #SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")

  # /MANIFEST:NO") # please, don't disable manifest generation, otherwise crash at start for vs2008
else(MSVC)
  if (DEBUG_MODE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -ftemplate-depth-1024 -O0 -std=c++11")
  else (DEBUG_MODE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -ftemplate-depth-1024 -O0 -std=c++11")
  endif()
endif()

# the gcc-4.2.1 coming with MacOS X is not compatible with the OpenMP pragmas we use, so disabling OpenMP for it
if((NOT APPLE) OR (NOT CMAKE_COMPILER_IS_GNUCXX) OR (GCC_VERSION VERSION_GREATER 4.2.1) OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
  find_package(OpenMP)
endif()

if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  message (STATUS "Found OpenMP")
  if(MSVC90 OR MSVC10)
    if(MSVC90)
      set(OPENMP_DLL VCOMP90)
    elseif(MSVC10)
      set(OPENMP_DLL VCOMP100)
    endif(MSVC90)
    set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} /DELAYLOAD:${OPENMP_DLL}D.dll")
    set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /DELAYLOAD:${OPENMP_DLL}.dll")
  endif(MSVC90 OR MSVC10)
else(OPENMP_FOUND)
  message (STATUS "Not found OpenMP")
endif()

set(GPStudio_OUTPUT_LIB_DIR ${PROJECT_BINARY_DIR}/lib)
set(GPStudio_OUTPUT_BIN_DIR ${PROJECT_BINARY_DIR}/bin)
make_directory(${GPStudio_OUTPUT_LIB_DIR})
make_directory(${GPStudio_OUTPUT_BIN_DIR})
if(WIN32)
  foreach(config ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${config} CONFIG)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${CONFIG} "${GPStudio_OUTPUT_LIB_DIR}")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG} "${GPStudio_OUTPUT_BIN_DIR}")
    # ---[ Windows requires DLLs (shared libraries) to be installed in the same directory as executables
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${CONFIG} "${GPStudio_OUTPUT_BIN_DIR}")
  endforeach(config)
else(WIN32)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${GPStudio_OUTPUT_LIB_DIR}")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${GPStudio_OUTPUT_BIN_DIR}")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${GPStudio_OUTPUT_LIB_DIR}")
endif(WIN32)

add_subdirectory(proto)
add_subdirectory(GPStudio)
