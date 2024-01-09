IF(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
	MESSAGE(FATAL_ERROR ": In-source builds are not allowed!")
ENDIF()

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.10 FATAL_ERROR)

SET(PROJECT_NAME Final)
PROJECT(${PROJECT_NAME})
ENABLE_LANGUAGE(Fortran)

IF(release_type)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
ELSE()
	SET(CMAKE_BUILD_TYPE Debug CACHE STRING "" FORCE)
ENDIF()
MESSAGE(====== BUILD_TYPE: ${CMAKE_BUILD_TYPE})

IF(MSVC_IDE)
	SET(Fortran_FLAGS          " /stand:f03 /threads")
	SET(Fortran_FLAGS_DEBUG    " /debug:full /warn:all /check:pointer /check:bounds /check:uninit /traceback")
	SET(Fortran_FLAGS_RELEASE  "/03 ")
	
ELSEIF("${CMAKE_CXX_COMPLIER_ID}" STREQUAL "GNU")
	SET(Fortran_FLAGS          " -std=f2008 -mthreads")
	SET(Fortran_FLAGS_DEBUG    " -g -wall -fcheck=pointer,bounds -fbacktrace")
	SET(Fortran_FLAGS_RELEASE  "-Ofast")

ENDIF()

ADD_SUBDIRECTORY(src)

ENABLE_TESTING()
ADD_SUBDIRECTORY(unitTest)

FILE(COPY doc DESTINATION ${CMAKE_CURRENT_BINARY_DIR})