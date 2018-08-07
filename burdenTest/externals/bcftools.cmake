SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(BCFTOOLS_PROJECT bcftools_project CACHE INTERNAL "bcftools_project name")
SET(BCFTOOLS_DIR ${CMAKE_BINARY_DIR}/externals/bcftools CACHE INTERNAL "bcftools project directory")
ExternalProject_Add(${BCFTOOLS_PROJECT}
        GIT_REPOSITORY https://github.com/samtools/bcftools.git
        GIT_TAG master
        CONFIGURE_COMMAND autoreconf && ./configure --enable-libgsl
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        PREFIX ${BCFTOOLS_DIR}
        CMAKE_CACHE_ARGS
                -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}


)

ExternalProject_Get_Property(${BCFTOOLS_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${BCFTOOLS_PROJECT} BINARY_DIR)

MESSAGE("BINARY_DIR: ${BINARY_DIR}")

MESSAGE("SRC_DIR: ${SOURCE_DIR}")


#SET(HTSLIB_LIB ${BINARY_DIR}/libgsl.a CACHE INTERNAL "BCFTOOLS Library")
#SET(HTSLIB_INCLUDE ${SOURCE_DIR} CACHE INTERNAL "BCFTOOLS Include")