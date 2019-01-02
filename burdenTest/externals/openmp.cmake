SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(OPENMP_PROJECT openmp_project CACHE INTERNAL "openmp_project name")
SET(OPENMP_DIR ${CMAKE_BINARY_DIR}/externals/openmp CACHE INTERNAL "openmp project directory")
ExternalProject_Add(${OPENMP_PROJECT}
        GIT_REPOSITORY https://github.com/llvm-mirror/openmp.git
        GIT_TAG master
        CONFIGURE_COMMAND  
        BUILD_COMMAND  cmake -B${OPENMP_DIR}/src/openmp_project/build -H${OPENMP_DIR}/src/openmp_project && make omp
        INSTALL_COMMAND 
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 0
        PREFIX ${OPENMP_DIR}
        CMAKE_CACHE_ARGS
                -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}


)

ExternalProject_Get_Property(${OPENMP_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${OPENMP_PROJECT} BINARY_DIR)

add_custom_command(
    OUTPUT ${BINARY_DIR}
    COMMAND make
    WORKING_DIRECTORY ${SOURCE_DIR}
)


MESSAGE("OPENMP_DIR: ${OPENMP_DIR}")
MESSAGE("OPENMP_BINARY_DIR: ${BINARY_DIR}")
MESSAGE("OPENMP_SRC_DIR: ${SOURCE_DIR}")

SET(OPENMP_LIB ${BINARY_DIR}/runtime/src/libgomp.dylib ${BINARY_DIR}/runtime/src/libomp.dylib ${BINARY_DIR}/runtime/src/libiomp5.dylib CACHE INTERNAL "OPENMP Library")
SET(OPENMP_INCLUDE ${BINARY_DIR}/runtime/src/ CACHE INTERNAL "OPENMP Include")