#RayCasting CMakeLists
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

set(PROJECT_NAME "VoroMultiRes")
project(${PROJECT_NAME} LANGUAGES CXX CUDA) 

#default build type is Release
#set(CMAKE_BUILD_TYPE Debug) #for dev only, TODO remove this 
if (CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE Release)
endif ()

########## Libraries/flags Starts Here ######################
find_package(CUDA REQUIRED)
include_directories("${CUDA_INCLUDE_DIRS}")
#set (GENCODE " -gencode arch=compute_35,code=sm_35")#K40
set (GENCODE " -gencode arch=compute_61,code=sm_61")#GTX 1050 TITAN Xp
#set (GENCODE " -gencode arch=compute_52,code=sm_52")#GTX 970
#set (GENCODE " -gencode arch=compute_70,code=sm_70")#TITAN V

message(GENCODE = "${GENCODE}")

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -lineinfo; -std=c++11; -expt-extended-lambda; -O3; -use_fast_math; -rdc=true; --ptxas-options=-v; -warn-spills; ${GENCODE})

string(APPEND CMAKE_CUDA_FLAGS ${GENCODE})#pass info to linker 

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -Wall -std=c++11")
########## Libraries/flags Ends Here ######################

########## inc/libs/exe/features Starts Here ######################
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(SRC_DIRECTORY "./" )
#convert relative path to absolute path
get_filename_component(SRC_DIRECTORY ${SRC_DIRECTORY} ABSOLUTE)
message( "SRC_DIRECTORY= " ${SRC_DIRECTORY})
set( MY_INCLUDES  "./" ${CMAKE_CURRENT_SOURCE_DIR}
					   ${SRC_DIRECTORY}/src
					   "../external/eigen-git-mirror")
					   
include_directories(${MY_INCLUDES})

add_definitions( -DINPUT_DIR=${SRC_DIRECTORY}/data/)
add_definitions( -DOUTPUT_DIR=${SRC_DIRECTORY}/viz/)

file(GLOB MY_SRC ${SRC_DIRECTORY}/src/*.h
				 ${SRC_DIRECTORY}/src/*.hxx
				 ${SRC_DIRECTORY}/src/*.hpp
				 ${SRC_DIRECTORY}/src/*.cu
				 ${SRC_DIRECTORY}/src/*.cuh
				 ${SRC_DIRECTORY}/src/*.cpp)

#create target exe
if (MSVC)
	add_executable(${PROJECT_NAME} "${MY_SRC}")
else()
	CUDA_ADD_EXECUTABLE(${PROJECT_NAME} "${MY_SRC}")
endif()

##add features to the library -> require it to compile with cuda c++11
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

##enable position-indep code 
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE  ON)

##enable the separate compilation on both device and host sides 
set_target_properties(${PROJECT_NAME} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(${PROJECT_NAME} ${CUDA_LIBRARIES} cuda)
########## inc/libs/exe/features Ends Here ######################