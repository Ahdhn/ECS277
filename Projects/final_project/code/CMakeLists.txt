cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME "openvdb_ao")
set(PROJECT_NAME_TUTORIAL "openvdb_tutorial")

project(${PROJECT_NAME})

set (BOOST_ROOT "C:/dev/boost_1_65_1_x64")
set (BOOST_LIBRARYDIR "C:/dev/boost_1_65_1/lib64-msvc-14.0")
set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}; ${BOOST_ROOT}")

set(CMAKE_MODULE_PATH "C:/projects/openvdb/deps/vfxcmake/cmake")
set (ILMBase_ROOT "c:/ilmbase")
set (OpenEXR_ROOT "c:/openexr")
set (TBB_USE_DEBUG_BUILD 0)
set (TBB_ROOT_DIR "c:/projects/openvdb/deps/tbb44_20160526oss")
set (BLOSC_ROOT "c:/blosc")
set (ZLIB_ROOT "C:/zlib")
set (ZLIB_INCLUDE_DIR "C:/zlib/include")
set (ZLIB_LIBRARY "C:/zlib/lib/zlibstatic.lib")
set (GLEW_LOCATION "C:/projects/openvdb/deps/glew-1.10.0")
set (GLFW_LOCATION "C:/projects/openvdb/deps/glfw-3.2.bin.WIN64")
set (CPPUNIT_INCLUDE_DIR "c:/cppunit/include")
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install")

set (OPENVDB_LOCATION "C:/projects/openvdb/install")

include_directories(${CMAKE_CURRENT_SOURCE_DIR} )

# used by other sub projects
set(OPENVDB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set(Boost_DEBUG OFF)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.59.0 REQUIRED COMPONENTS iostreams system thread filesystem)

if(Boost_FOUND)
  message(STATUS "BOOST FOUND: " ${Boost_FOUND})
  message(STATUS "BOOST VERSION: " ${Boost_VERSION})
endif()

find_package(ILMBase REQUIRED)

if(ILMBASE_FOUND)
  message(STATUS "ILMBASE FOUND: " ${ILMBASE_FOUND})
  message(STATUS "ILMBASE VERSION: " ${ILMBASE_VERSION})
endif()

find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy)

if(TBB_FOUND)
  message(STATUS "TBB FOUND: " ${TBB_FOUND})
  message(STATUS "TBB VERSION: " ${TBB_VERSION})
endif()

find_package(ZLIB REQUIRED)

if(ZLIB_FOUND)
  message(STATUS "ZLIB FOUND: " ${ZLIB_FOUND})
  message(STATUS "ZLIB VERSION: " ${ZLIB_VERSION_STRING})
endif()

find_package(Blosc)

if(BLOSC_FOUND)
  message(STATUS "BLOSC FOUND: " ${BLOSC_FOUND})
  message(STATUS "BLOSC VERSION: " ${BLOSC_VERSION})
endif()

find_package(OpenVDB)

if(OpenVDB_FOUND)
  message(STATUS "OpenVDB FOUND: " ${OpenVDB_FOUND})
  message(STATUS "OpenVDB VERSION Major: " ${OpenVDB_MAJOR_VERSION})
  message(STATUS "OpenVDB VERSION Minor: " ${OpenVDB_MINOR_VERSION})  
endif()

find_package(OpenGL)
find_package(GLEW)

if(GLEW_FOUND)
  message(STATUS "GLEW FOUND: " ${GLEW_FOUND})  
endif()

find_package(GLFW)

if(GLFW_FOUND)
  message(STATUS "GLFW FOUND: " ${GLFW_FOUND})
  message(STATUS "GLFW VERSION: " ${GLFW_VERSION})
endif()





# collect source files
file(GLOB H_FILES ${PROJECT_NAME}/*.h
                  ${PROJECT_NAME}/viewer/*.h)
file(GLOB CXX_FILES ${PROJECT_NAME}/*.cc
                    ${PROJECT_NAME}/viewer/*.cc)
set(SOURCE_FILES ${H_FILES} ${CXX_FILES})

if(WIN32)
  add_definitions(-DNOMINMAX -DOPENVDB_STATICLIB
      -DOPENVDB_PRIVATE -DHALF_EXPORTS
      -DOPENVDB_DLL -DBOOST_ALL_NO_LIB
      -DOPENVDB_USE_BLOSC)
endif()

include_directories(
    ${ILMBASE_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${TBB_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
    ${BLOSC_INCLUDE_DIRS}
	${OpenVDB_INCLUDE_DIR}
	${GLFW_INCLUDE_DIR}
	${GLEW_INCLUDE_DIR})

add_definitions(-D_USE_MATH_DEFINES -DOPENVDB_USE_GLFW_3 -DGLEW_STATIC) 


set( INPUT_DIR "data" )
get_filename_component(INPUT_DIR ${INPUT_DIR} ABSOLUTE)
add_definitions( -DINPUT_DIR=${INPUT_DIR})

#########################################################
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME}
  ${Boost_LIBRARIES}
  ${ZLIB_LIBRARIES}
  optimized ${ILMBASE_Half_LIBRARY}
  ${TBB_LIBRARIES}
  optimized ${BLOSC_STATIC_LIBRARY}
  ${OpenVDB_LIBRARIES}
  ${GLFW_LIBRARIES}
  ${GLEW_LIBRARY}
  ${OPENGL_LIBRARIES})  
if(WIN32)
  set_target_properties(${PROJECT_NAME}
      PROPERTIES
      COMPILE_FLAGS "/bigobj")
endif()

















##########################################################
add_executable(${PROJECT_NAME_TUTORIAL} ${PROJECT_NAME_TUTORIAL}/openvdb_tutorial.cc)
target_link_libraries(${PROJECT_NAME_TUTORIAL}
  ${Boost_LIBRARIES}
  ${ZLIB_LIBRARIES}
  optimized ${ILMBASE_Half_LIBRARY}
  ${TBB_LIBRARIES}
  optimized ${BLOSC_STATIC_LIBRARY}
  ${OpenVDB_LIBRARIES})  
if(WIN32)
  set_target_properties(${PROJECT_NAME_TUTORIAL}
      PROPERTIES
      COMPILE_FLAGS "/bigobj")
endif()

