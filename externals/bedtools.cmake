SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(BEDTOOLS_PROJECT bedtools_project CACHE INTERNAL "bedtools_project name")
SET(BEDTOOLS_DIR ${CMAKE_BINARY_DIR}/externals/bedtools CACHE INTERNAL "bedtools project directory")
ExternalProject_Add(${BEDTOOLS_PROJECT}
        GIT_REPOSITORY https://github.com/arq5x/bedtools2.git
        GIT_TAG master
        CONFIGURE_COMMAND 
        BUILD_COMMAND make
        INSTALL_COMMAND 
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 0
        PREFIX ${BEDTOOLS_DIR}
        CMAKE_CACHE_ARGS
                -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}


)

ExternalProject_Get_Property(${BEDTOOLS_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${BEDTOOLS_PROJECT} BINARY_DIR)

add_custom_command(
    OUTPUT ${BINARY_DIR}
    COMMAND make
    WORKING_DIRECTORY ${SOURCE_DIR}
)



MESSAGE("BINARY_DIR: ${BINARY_DIR}")
MESSAGE("SRC_DIR: ${SOURCE_DIR}")