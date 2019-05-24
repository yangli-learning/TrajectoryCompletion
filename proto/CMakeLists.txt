project(Proto)

#	Google Protocol Buffers
find_package(Protobuf REQUIRED)

# As of Ubuntu 14.04 protoc is no longer a part of libprotobuf-dev package and should be installed
# separately as in: sudo apt-get install protobuf-compiler
if(PROTOBUF_PROTOC_EXECUTABLE)
	message(STATUS "Found PROTOBUF Compiler: ${PROTOBUF_PROTOC_EXECUTABLE}")
else()
	message(FATAL_ERROR "Could not find PROTOBUF Compiler")
endif()

include_directories(${PROTOBUF_INCLUDE_DIR})
file(GLOB ProtoFiles "${CMAKE_CURRENT_SOURCE_DIR}/*.proto")
PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders ${ProtoFiles})

SOURCE_GROUP("Resources" FILES ${ProtoFiles})

set(lib_name proto)
add_library(${lib_name} ${ProtoFiles} ${ProtoSources} ${ProtoHeaders})
target_link_libraries(${lib_name} ${PROTOBUF_LIBRARIES})

if(WIN32 AND MSVC)
  set_target_properties(${lib_name} PROPERTIES LINK_FLAGS_RELEASE /OPT:REF)
elseif(CMAKE_SYSTEMname STREQUAL "Darwin")
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set_target_properties(${lib_name} PROPERTIES LINK_FLAGS -Wl)
  endif()
elseif(__COMPILER_PATHSCALE)
  set_target_properties(${lib_name} PROPERTIES LINK_FLAGS -mp)
else()
  set_target_properties(${lib_name} PROPERTIES LINK_FLAGS -Wl)
endif()

set_target_properties(${lib_name} PROPERTIES DEBUG_POSTFIX _debug)
set_target_properties(${lib_name} PROPERTIES RELEASE_POSTFIX _release)

# Create proto include directory
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/proto)

# Copy proto headers to /proto
foreach(header ${ProtoHeaders})
    ADD_CUSTOM_COMMAND(TARGET proto
        COMMAND cmake -E copy ${header}
        ${CMAKE_SOURCE_DIR}/proto/
        DEPENDS ${header}
        )
endforeach(header)
