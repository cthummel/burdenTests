cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tools)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tools/tmp)

if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif ()

INCLUDE_DIRECTORIES(externals, src)

ADD_SUBDIRECTORY(externals)
ADD_SUBDIRECTORY(src)

