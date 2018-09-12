
########################################################################
cmake_minimum_required (VERSION 3.5.1)
project (Absolve)

########################################################################
# verbose 
set (CMAKE_VERBOSE_MAKEFILE ON)

########################################################################
# set (CMAKE_BUILD_TYPE "Debug")
set (CMAKE_BUILD_TYPE "Release")

########################################################################
include_directories(
./dep/hmmer-3.1b2/src
./dep/hmmer-3.1b2/easel
./dep/bwa-0.7.12
)

link_directories(
./dep/hmmer-3.1b2/src
./dep/hmmer-3.1b2/easel
./dep/bwa-0.7.12
)

########################################################################

if(UNIX AND NOT APPLE)
	find_library(MLIB m)
	if(MLIB STREQUAL "MLIB-NOTFOUND")
		message(FATAL_ERROR "could not find math library/libm, please install")
	endif()
endif()

find_library(ZLIB z)
if(ZLIB STREQUAL "ZLIB-NOTFOUND")
	message(FATAL_ERROR "could not find zlib library, please install")
endif()

find_library(PTHREADLIB pthread)
if(PTHREADLIB STREQUAL "PTHREADLIB-NOTFOUND")
	message(FATAL_ERROR "could not find posix thread/pthread library, please install")
endif()

find_library(BOOST_SYSTEM boost_system)
if(BOOST_SYSTEM STREQUAL "BOOST_SYSTEM-NOTFOUND")
	message(FATAL_ERROR "could not find boost_system library, please install")
endif()

find_library(BOOST_FILESYSTEM boost_filesystem)
if(BOOST_FILESYSTEM STREQUAL "BOOST_FILESYSTEM-NOTFOUND")
	message(FATAL_ERROR "could not find boost_filesystem library, please install")
endif()

find_library(BOOST_IOSTREAMS boost_iostreams)
if(BOOST_IOSTREAMS STREQUAL "BOOST_IOSTREAMS-NOTFOUND")
	message(FATAL_ERROR "could not find boost_iostreams library, please install")
endif()

find_library(BOOST_REGEX boost_regex)
if(BOOST_REGEX STREQUAL "BOOST_REGEX-NOTFOUND")
	message(FATAL_ERROR "could not find boost_regex library, please install")
endif()




########################################################################
add_executable(absolve_rel 
	src/aligner.cpp
	src/bwaaligner.cpp
	src/bwareader.cpp
	src/delimfile.cpp
	src/fasta.cpp
	src/ighmm.cpp
	src/main.cpp
	src/pipeline.cpp
	src/result.cpp
	src/runstats.cpp
	src/sswaligner.cpp
	src/ssw_cpp.cpp
	src/translate.cpp
	src/util.cpp
	src/bwamem_wrapper.c
	src/ssw.c
)

set_target_properties(absolve_rel PROPERTIES LINKER_LANGUAGE CXX)

if(UNIX AND NOT APPLE)
	target_link_libraries(absolve_rel 
		m
		z
		pthread
		boost_system
		boost_filesystem
		boost_iostreams
		boost_regex
		hmmer
		easel
		bwa
	)
elseif(UNIX AND APPLE)
	target_link_libraries(absolve_rel 
		z
		pthread
		boost_system
		boost_filesystem
		boost_iostreams
		boost_regex
		hmmer
		easel
		bwa
	)
else()
	message(FATAL_ERROR "only Linux and Darwin(MacOSX) are supported at this time")
endif()

