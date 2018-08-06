if("master" STREQUAL "")
  message(FATAL_ERROR "Tag for git checkout should not be empty.")
endif()

set(run 0)

if("/uufs/chpc.utah.edu/common/home/u0261133/burdenTests/burdenTest/externals/externals/gsl/src/gsl_project-stamp/gsl_project-gitinfo.txt" IS_NEWER_THAN "/uufs/chpc.utah.edu/common/home/u0261133/burdenTests/burdenTest/externals/externals/gsl/src/gsl_project-stamp/gsl_project-gitclone-lastrun.txt")
  set(run 1)
endif()

if(NOT run)
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: '/uufs/chpc.utah.edu/common/home/u0261133/burdenTests/burdenTest/externals/externals/gsl/src/gsl_project-stamp/gsl_project-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E remove_directory "/uufs/chpc.utah.edu/common/home/u0261133/burdenTests/burdenTest/externals/externals/gsl/src/gsl_project"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/uufs/chpc.utah.edu/common/home/u0261133/burdenTests/burdenTest/externals/externals/gsl/src/gsl_project'")
endif()

# try the clone 3 times incase there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git" clone "https://github.com/ampl/gsl.git" "gsl_project"
    WORKING_DIRECTORY "/uufs/chpc.utah.edu/common/home/u0261133/burdenTests/burdenTest/externals/externals/gsl/src"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/ampl/gsl.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git" checkout master
  WORKING_DIRECTORY "/uufs/chpc.utah.edu/common/home/u0261133/burdenTests/burdenTest/externals/externals/gsl/src/gsl_project"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'master'")
endif()

execute_process(
  COMMAND "/usr/bin/git" submodule init
  WORKING_DIRECTORY "/uufs/chpc.utah.edu/common/home/u0261133/burdenTests/burdenTest/externals/externals/gsl/src/gsl_project"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to init submodules in: '/uufs/chpc.utah.edu/common/home/u0261133/burdenTests/burdenTest/externals/externals/gsl/src/gsl_project'")
endif()

execute_process(
  COMMAND "/usr/bin/git" submodule update --recursive
  WORKING_DIRECTORY "/uufs/chpc.utah.edu/common/home/u0261133/burdenTests/burdenTest/externals/externals/gsl/src/gsl_project"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/uufs/chpc.utah.edu/common/home/u0261133/burdenTests/burdenTest/externals/externals/gsl/src/gsl_project'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "/uufs/chpc.utah.edu/common/home/u0261133/burdenTests/burdenTest/externals/externals/gsl/src/gsl_project-stamp/gsl_project-gitinfo.txt"
    "/uufs/chpc.utah.edu/common/home/u0261133/burdenTests/burdenTest/externals/externals/gsl/src/gsl_project-stamp/gsl_project-gitclone-lastrun.txt"
  WORKING_DIRECTORY "/uufs/chpc.utah.edu/common/home/u0261133/burdenTests/burdenTest/externals/externals/gsl/src/gsl_project"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/uufs/chpc.utah.edu/common/home/u0261133/burdenTests/burdenTest/externals/externals/gsl/src/gsl_project-stamp/gsl_project-gitclone-lastrun.txt'")
endif()

