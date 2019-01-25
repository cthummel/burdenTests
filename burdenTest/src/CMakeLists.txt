#INCLUDE_DIRECTORIES(
#	${GSL_INCLUDE}
#	#Other third party functions here
#)

SET(MAIN_SRC main.cpp)

INCLUDE_DIRECTORIES(
	${PROJECT_SOURCE_DIR}
	${GSL_INCLUDE}
)

ADD_EXECUTABLE(test ${MAIN_SRC})

ADD_LIBRARY(burden_test_core STATIC
			     wsbt.cpp
			     input.cpp
			     output.cpp
			     cast.cpp
			     skat.cpp
			     geneInfo.cpp
			    # skato.cpp
			     genericBurdenTest.cpp
)

SET(CORE_LIB burden_test_core CACHE STRING "THE CORE LIBRARY")

TARGET_LINK_LIBRARIES(test 
	   ${CORE_LIB}
	   ${GSL_LIB}
	   ${OPENMP_LIB}
)

add_dependencies(test  ${BURDEN_TEST_EXTERNAL_PROJECT})





