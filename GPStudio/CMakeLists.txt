set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(Boost_USE_STATIC_LIBS  ON)
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF)


option(Boost_NO_SYSTEM_PATHS "BoostPath" TRUE) 
set(Boost_CUSTOM_PREFIX "/usr/share/" CACHE STRING "Boost library prefix")
if (Boost_NO_SYSTEM_PATHS)
  set(BOOST_ROOT ${Boost_CUSTOM_PREFIX})
  set(BOOST_INCLUDE_DIRS "${BOOST_ROOT}/include")
  set(BOOST_LIBRARY_DIRS "${BOOST_ROOT}/lib")
  message(STATUS "Boost include dirs:" ${BOOST_INCLUDE_DIRS})
  message(STATUS "Boost library dirs:" ${BOOST_LIBRARY_DIRS})
endif ()


find_package(Boost REQUIRED timer chrono exception system filesystem )
include_directories(${Boost_INCLUDE_DIRS})


if (CV_LIB_DIR)
	link_directories(${CV_LIB_DIR})
endif()
if (CV_INCLUDE_DIR)
	link_directories(${CV_INCLUDE_DIR})
endif()

find_package(OpenCV REQUIRED opencv_core opencv_imgproc opencv_flann)
if(3 VERSION_EQUAL ${OpenCV_VERSION_MAJOR})
  find_package(OpenCV REQUIRED opencv_imgcodecs)
else()
  find_package(OpenCV REQUIRED opencv_highgui)
endif()
include_directories(${OpenCV_INCLUDE_DIRS})


option(GFlags_NO_SYSTEM_PATHS "GFlagsPath" FALSE)
if (GFlags_NO_SYSTEM_PATHS)
set(GFLAGS_INCLUDE_DIR /usr/include)
set(GFLAGS_LIBS /usr/lib/x86_64-linux-gnu)
set(GFLAGS_STATIC_LIB "${GFLAGS_LIBS}/libgflags.a")
else()
find_package(gflags REQUIRED)    
endif()

#find_package(CUDA 5.5 REQUIRED)
#include_directories(${CUDA_INCLUDE_DIRS})

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
find_package(Qt5Core REQUIRED)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5OpenGL REQUIRED)
find_package(Qt5Concurrent REQUIRED)
#find_package(QCustomPlot REQUIRED) #qcustomplot
#include_directories(${QCustomPlot_INCLUDE_DIR})
#message(STATUS "QCustomPlot include dir:" ${QCustomPlot_INCLUDE_DIR})
#message(STATUS "QCustomPlot libraries:" ${QCustomPlot_LIBRARIES})

option(OSG_NO_SYSTEM_PATHS "OSGPath" true) 
if (OSG_NO_SYSTEM_PATHS)
  set(OSG_LIB_ROOT /usr/lib/x86_64-linux-gnu) # yang
  #set(OSG_LIB_ROOT /usr/lib)   # yangyan
  set(OPENSCENEGRAPH_LIBRARIES ${OSG_LIB_ROOT}/libosg.so ${OSG_LIB_ROOT}/libOpenThreads.so ${OSG_LIB_ROOT}/libosgDB.so ${OSG_LIB_ROOT}/libosgGA.so ${OSG_LIB_ROOT}/libosgManipulator.so ${OSG_LIB_ROOT}/libosgQt.so ${OSG_LIB_ROOT}/libosgUtil.so ${OSG_LIB_ROOT}/libosgViewer.so ${OSG_LIB_ROOT}/libosgWidget.so ${OSG_LIB_ROOT}/libosgText.so )
 
  set(OPENSCENEGRAPH_INCLUDE_DIRS /usr/include)

  message(STATUS "OSG include dir: " ${OPENSCENEGRAPH_INCLUDE_DIRS})
  message(STATUS "OSG lib dir: " ${OPENSCENEGRAPH_LIBRARIES})
else()
find_package(OpenSceneGraph REQUIRED osgViewer osgText osgDB osgGA osgQt osgManipulator osgUtil) 
endif()
include_directories(${OPENSCENEGRAPH_INCLUDE_DIRS})

find_package(GeoGraphicLib REQUIRED)
include_directories(${GEOGRAPHICLIB_INCLUDE_DIR})
add_definitions(-DGEOGRAPHICLIB_SHARED_LIB=0)

#	Google Protocol Buffers
find_package(Protobuf REQUIRED)
include_directories(${PROTOBUF_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/proto)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/impl)

# Set Qt files and resources here
set(uis       resources/main_window.ui)
set(resources  resources/main_window.qrc)

# Qt stuff
QT5_WRAP_UI(ui_srcs ${uis})
QT5_ADD_RESOURCES(resource_srcs ${resources})

# Put the ui in the windows project file
IF (${CMAKE_BUILD_TOOL} MATCHES "msdev")
  SET (srcs ${srcs} ${uis})
ENDIF (${CMAKE_BUILD_TOOL} MATCHES "msdev")
IF (${CMAKE_BUILD_TOOL} MATCHES "devenv")
  SET (srcs ${srcs} ${uis})
ENDIF (${CMAKE_BUILD_TOOL} MATCHES "devenv")

set(PARAMETER_FILENAME "${CMAKE_CURRENT_SOURCE_DIR}/parameters.ini" CACHE FILEPATH "Parameters filename")

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
                "${CMAKE_CURRENT_BINARY_DIR}/config.h"
                )

set(config      ${CMAKE_CURRENT_SOURCE_DIR}/parameters.ini)

set(incs    include/main_window.h
            include/adapter_widget.h
            include/osg_viewer_widget.h
            include/scene_widget.h
            include/axis_indicator.h
            include/stateset_manipulator.h
            include/threaded_painter.h
            include/toggle_handler.h
            include/update_callback.h
            include/point_intersector.h
            include/common.h
            include/color_map.h
            include/pick_handler.h
            include/renderable.h
            include/osg_utility.h
            include/information.h
            include/force_update_visitor.h
            include/task_dispatcher.h
            include/parameters.h
            include/singleton.h
            include/gps_trajectories.h
            include/bbox_2d.h
            include/l_shape_iterator.h
            include/l1skeleton.h 
            include/skeleton_graph.h
            include/dgraph.h
            include/voting_map_renderer.h
            include/l1skeleton_graph.h
            include/junction_network.h
            include/junction_graph.h
            include/process_graph.h
            include/viterbi.h
	    include/trajectory_completion.h
	    include/guiless.h
include/qcustomplot.h
            )

set(srcs    src/main.cpp
            src/main_window.cpp
            src/adapter_widget.cpp
            src/osg_viewer_widget.cpp
            src/scene_widget.cpp
            src/axis_indicator.cpp
            src/stateset_manipulator.cpp
            src/threaded_painter.cpp
            src/toggle_handler.cpp
            src/update_callback.cpp
            src/point_intersector.cpp
            src/color_map.cpp
            src/pick_handler.cpp
            src/renderable.cpp
            src/osg_utility.cpp
            src/information.cpp
            src/common.cpp
            src/force_update_visitor.cpp
            src/parameters.cpp
            src/task_dispatcher.cpp
            src/gps_trajectories.cpp
            src/bbox_2d.cpp
            src/l_shape_iterator.cpp
	    src/l1skeleton.cpp  
	    src/skeleton_graph.cpp
            src/dgraph.cpp
            src/voting_map_renderer.cpp
            src/l1skeleton_graph.cpp
            src/junction_network.cpp
            src/junction_graph.cpp
	    src/process_graph.cpp
            src/viterbi.cpp
	    src/trajectory_completion.cpp
	    src/guiless.cpp
src/qcustomplot.cpp
           )

set(incs_impl   )

# Organize files
SOURCE_GROUP("Resources" FILES ${uis} ${resources} ${config})
SOURCE_GROUP("Generated" FILES ${ui_srcs} ${moc_srcs} ${resource_srcs} ${CMAKE_CURRENT_BINARY_DIR}/config.h)
SET_SOURCE_FILES_PROPERTIES(${srcs} PROPERTIES OBJECT_DEPENDS "${ui_srcs}")

set(exe_name GPStudio)
add_executable(${exe_name} ${ui_srcs} ${moc_srcs} ${resource_srcs} ${incs} ${srcs} ${incs_impl} ${config})

target_link_libraries(${exe_name} proto gflags 
                      ${Boost_TIMER_LIBRARY} ${Boost_CHRONO_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_EXCEPTION_LIBRARY} ${Boost_GEOMETRY_LIBRARY} ${Boost_GRAPH_LIBRARY}
                      ${OpenCV_LIBS}
                      ${OPENSCENEGRAPH_LIBRARIES}
                      ${PROTOBUF_LIBRARIES}
                      ${GEOGRAPHICLIB_LIBRARIES}
                      Qt5::Core Qt5::Widgets Qt5::OpenGL Qt5::Concurrent  )

if(WIN32 AND MSVC)
  # to workaround osgDB-std::ofstream conflict
  set_target_properties(${exe_name} PROPERTIES LINK_FLAGS /FORCE:MULTIPLE)
  set_target_properties(${exe_name} PROPERTIES LINK_FLAGS_RELEASE /OPT:REF)
elseif(CMAKE_SYSTEMname STREQUAL "Darwin")
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set_target_properties(${exe_name} PROPERTIES LINK_FLAGS -Wl)
  endif()
elseif(__COMPILER_PATHSCALE)
  set_target_properties(${exe_name} PROPERTIES LINK_FLAGS -mp)
else()
  set_target_properties(${exe_name} PROPERTIES LINK_FLAGS -Wl,--as-needed)
endif()

set_target_properties(${exe_name} PROPERTIES COMPILE_DEFINITIONS "UNICODE")

set_target_properties(${exe_name} PROPERTIES DEBUG_POSTFIX _debug)
set_target_properties(${exe_name} PROPERTIES RELEASE_POSTFIX _release)
