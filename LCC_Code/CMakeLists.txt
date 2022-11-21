# @file SBFem/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.23.0)

project (H1_GD_HdivH1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_definitions(-D_AUTODIFF)

IF(APPLE)
	# Mac OS X specific code
	add_definitions(-DMACOSX)
	set(MACOSX MACOSX)
	FIND_LIBRARY(ACCELERATE_LIB Accelerate)
	link_libraries(${ACCELERATE_LIB})
	#   SET(OperatingSystem "Mac OS X")
ENDIF(APPLE)

IF (WIN32)
	# Prevents timespec redefinition problem with Visual Studio 2015
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
    if (NOT CMAKE_GENERATOR MATCHES "Ninja")
    	set(CMAKE_INSTALL_PREFIX "C:/Program Files/PZ" CACHE PATH "Installation path, prepended onto install directories." FORCE) 
	endif()
ELSE (WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")#this flag is for preventing errors when casting from double to float
	if ( NOT APPLE AND NOT CMAKE_BUILD_TYPE )
		set (CMAKE_BUILD_TYPE Release CACHE STRING "Options: Release Debug None MinSizeRel RelWithDebInfo" FORCE) 
	endif()
ENDIF()

include_directories(${PZ_INCLUDE_DIRS})

# Location to search files or libraries
set(SEARCH_DIRS ENV PATH /usr/local/opt/ /usr/local/opt/qt5/lib/cmake/ /usr/include /usr/local/include /usr/local/include/boost /usr/local/include/log4cxx /usr/local/include/tbb /usr/lib /usr/lib64 /usr/local/lib ./. ./.. ./../libs ./../libs/fad ./../libs/metis ./../libs/metis/include ./../libs/metis/lib ./../externallibs ./../externallibs/fad ./../externallibs/pthread ./../externallibs/metis ./../externallibs/pthread/include ./../externallibs/pthread/lib ./..metis ./../metis/include ./../metis/lib ./externallibs/lib ./externallibs/include ./../externallibs/lib ./../externallibs/include ./externallibs/pzlib/lib ./externallibs/pzlib/include ./../externallibs/pzlib/lib ./../externallibs/pzlib/include)


#-----
#enabling to use LOG4CXX library
option(USING_LOG4CXX "Whether the LOG4CXX library will be linked in" OFF)
if(USING_LOG4CXX)
	#Adding a variable to hold this definition
	set(LOG4CXX LOG4CXX)
	find_path(LOG4CXX_INCLUDE log4cxx PATHS ${SEARCH_DIRS} ${PROJECT_SOURCE_DIR}/../log4cxx ${PROJECT_SOURCE_DIR}/../log4cxx_src/include ${PROJECT_SOURCE_DIR}/../externallibs/log4cxx_src ${PROJECT_SOURCE_DIR}/../externallibs/log4cxx_src/include)
	find_library(LOG4CXX_LIB_RELEASE
			NAMES liblog4cxx.dylib log4cxx.lib liblog4cxx.so liblog4cxx.a
			PATHS ${SEARCH_DIRS} ./../log4cxx_src ./../log4cxx_src/lib ./../externallibs/log4cxx_src ./../externallibs/log4cxx_src/lib ./../externallibs/lib ./../externallibs/lib/Release
			)
	find_library(LOG4CXX_LIB_DEBUG
			NAMES liblog4cxx.dylib log4cxx.lib liblog4cxx.so liblog4cxx.a       # Change name here if debug lib has different name.
			PATHS ${SEARCH_DIRS} ./../log4cxx_src ./../log4cxx_src/lib ./../externallibs/log4cxx_src ./../externallibs/log4cxx_src/lib ./../externallibs/lib ./../externallibs/lib/Debug
			)
	set(LOG4CXX_LIB
			debug ${LOG4CXX_LIB_DEBUG}
			optimized ${LOG4CXX_LIB_RELEASE}
			)
	if(LOG4CXX_INCLUDE-NOTFOUND)
		set (LOG4CXX_INCLUDE "" CACHE PATH "Directory where log4cxx header files can be found")
	else()
		include_directories(${LOG4CXX_INCLUDE})
	endif()

	if(LOG4CXX_LIB-NOTFOUND)
		set (LOG4CXX_LIB "" CACHE PATH "Directory where the log4cxx library can be found")
	endif()
endif(USING_LOG4CXX)



#enabling pz library
find_path(PZ_INCLUDEINI NAMES pz_config.h PATHS ${SEARCH_DIRS})
set(PZ_INCLUDE "${PZ_INCLUDEINI}" "${PZ_INCLUDEINI}/Analysis" "${PZ_INCLUDEINI}/Common" "${PZ_INCLUDEINI}/Frontal" "${PZ_INCLUDEINI}/External" "${PZ_INCLUDEINI}/External/FAD" "${PZ_INCLUDEINI}/External/FAD/fad" "${PZ_INCLUDEINI}/External/FAD/TinyFadET" "${PZ_INCLUDEINI}/External/FAD/TinyFad" "${PZ_INCLUDEINI}/Geom" "${PZ_INCLUDEINI}/Integral" "${PZ_INCLUDEINI}/LinearSolvers" "${PZ_INCLUDEINI}/Material" "${PZ_INCLUDEINI}/Material/REAL" "${PZ_INCLUDEINI}/Material/Complex" "${PZ_INCLUDEINI}/Matrix" "${PZ_INCLUDEINI}/Mesh" "${PZ_INCLUDEINI}/Multigrid" "${PZ_INCLUDEINI}/Optimization" "${PZ_INCLUDEINI}/Post" "${PZ_INCLUDEINI}/Pre" "${PZ_INCLUDEINI}/Refine" "${PZ_INCLUDEINI}/Save" "${PZ_INCLUDEINI}/Shape" "${PZ_INCLUDEINI}/StrMatrix" "${PZ_INCLUDEINI}/SubStruct" "${PZ_INCLUDEINI}/Topology" "${PZ_INCLUDEINI}/Util")
find_library(PZ_LIB NAMES libpz.lib pz.lib libpz.a PATHS ${SEARCH_DIRS}  ${PROJECT_SOURCE_DIR}/../pz/lib)
if(PZ_INCLUDE-NOTFOUND)
	set(PZ_INCLUDE "" CACHE PATH "Directory where pz.h can be found")
endif()
if(PZ_LIB-NOTFOUND)
	set(PZ_LIBexit "" CACHE PATH "Directory where the pz library can be found")
endif()
include_directories(${PZ_INCLUDE})

#enabling pthread library
find_path(PTHREAD_INCLUDE NAMES PThread.h pthread.h PATHS ${SEARCH_DIRS} D:/Projetos/externallibs/include )
find_library(PTHREAD_LIB NAMES libpthread.dylib pthread.lib pthreadVC2.lib libPThread.so libPThread.a libpthread.so libpthread.a PATHS ${SEARCH_DIRS} D:/Projetos/externallibs/lib  ${PROJECT_SOURCE_DIR}/../pthread/lib)
if(PTHREAD_INCLUDE-NOTFOUND)
	set(PTHREAD_INCLUDE "" CACHE PATH "Directory where pthread.h can be found")
endif()
if(PTHREAD_LIB-NOTFOUND)
	set(PTHREAD_LIBexit "" CACHE PATH "Directory where the pthread library can be found")
endif()
include_directories(${PTHREAD_INCLUDE})

if(WIN32)
	#add a define saying it is a VC compiler
	set(VC "#Will add a define VC on pz_config.h")
	add_definitions(-DVC)
	#define _USE_MATH_DEFINES for example: M_PI constant definitions.
	add_definitions(-D_USE_MATH_DEFINES)
	#disabling VC warnings
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	#disabling All warnings
	#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W0")
	#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W0")
	#disabling STL warning for std::_Vector_alloc when LOG4CXX is enabled
	if (USING_LOG4CXX)
		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /wd4251")
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /wd4251")
	endif()
	#define use of pthread static lib.
	#enabling /bigobj
	add_definitions("/bigobj")

	FOREACH(FLAG_TYPE EXE MODULE SHARED)
		STRING (REPLACE "INCREMENTAL:YES" "INCREMENTAL" FLAG_TMP "${CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG}")
		STRING (REPLACE "INCREMENTAL:NO" "INCREMENTAL" FLAG_TMP ${FLAG_TMP})
		STRING (REPLACE "INCREMENTAL" "INCREMENTAL:NO" FLAG_TMP ${FLAG_TMP})
		STRING (REPLACE "/EDITANDCONTINUE" "" FLAG_TMP ${FLAG_TMP})
		SET(CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG "${FLAG_TMP}" CACHE STRING "Overriding default debug ${FLAG_TYPE} linker flags." FORCE)
		MARK_AS_ADVANCED (CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG)
	ENDFOREACH ()

endif()


add_subdirectory(LCCCommon)
add_subdirectory(LCCMaterials)
add_subdirectory(LCCAnalysis)
add_subdirectory(LCCAdaptive)
add_subdirectory(LCCEstimator)

include_directories(LCCCommon)
include_directories(LCCMaterials)
include_directories(LCCAnalysis)
include_directories(LCCAdaptive)
include_directories(LCCEstimator)

add_subdirectory(LCC_Projects)
