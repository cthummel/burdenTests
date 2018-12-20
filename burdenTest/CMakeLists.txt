cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
elseif (OPENMP_NOTFOUND)
    
endif()
if(APPLE)
	#set(OPENMP_LIBRARIES "/usr/local/Cellar/llvm/7.0.0/lib")
	#set(OPENMP_INCLUDES "/usr/local/Cellar/llvm/7.0.0/include")

	if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        	set(OpenMP_C "${CMAKE_C_COMPILER}")
        	set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp ")
        	set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp7")
        	set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        	set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        	set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    	endif()

    	if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      		set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
      		set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp ")
      		set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp7")
      		set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      		set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      		set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    	endif()

endif()


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tools)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tools/tmp)

if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif ()

INCLUDE_DIRECTORIES(externals, src)

ADD_SUBDIRECTORY(externals)
ADD_SUBDIRECTORY(src)

