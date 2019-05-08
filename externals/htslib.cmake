SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(HTSLIB_PROJECT htslib_project CACHE INTERNAL "htslib_project name")
SET(HTSLIB_DIR ${CMAKE_BINARY_DIR}/externals/htslib CACHE INTERNAL "htslib project directory")
ExternalProject_Add(${HTSLIB_PROJECT}
        GIT_REPOSITORY https://github.com/samtools/htslib.git
        GIT_TAG master
        CONFIGURE_COMMAND autoreconf && ./configure --prefix=${CMAKE_SOURCE_DIR}/externals
        BUILD_COMMAND make
        INSTALL_COMMAND 
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        PREFIX ${HTSLIB_DIR}
        CMAKE_CACHE_ARGS
                -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}


)

ExternalProject_Get_Property(${HTSLIB_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${HTSLIB_PROJECT} BINARY_DIR)

MESSAGE("BINARY_DIR: ${BINARY_DIR}")
MESSAGE("SRC_DIR: ${SOURCE_DIR}")


SET(HTSLIB_LIB ${BINARY_DIR}/libgsl.a CACHE INTERNAL "HTSLIB Library")
SET(HTSLIB_INCLUDE ${SOURCE_DIR} CACHE INTERNAL "HTSLIB Include")