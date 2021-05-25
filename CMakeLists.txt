cmake_minimum_required(VERSION 3.15.0)

project(descur LANGUAGES C CXX Fortran)

# Start Project and note the language used. The BLAS and LAPACK libraries for ifort require a working C and C++ compiler.
project (DESCUR Fortran C CXX)

set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "ON")

# Set a directories to build all binary files.
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)
set (LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Add extra compiler flags to various builds.
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${EXTRA_RELEASE_FLAGS}")
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${EXTRA_DEBUG_FLAGS}")
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_FLAGS}")

# Fortran specific settings. The first setting tells the compiler to use the C preprocessor.
# The second places a common directory to place all of the module files.
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

find_package (BLAS REQUIRED)
if (NOT ${BLA_VENDOR} MATCHES "^Intel")
    find_package (LAPACK REQUIRED)
endif()

set (descur_sources "")
add_subdirectory(src)

add_library (descur STATIC ${descur_sources})
set_target_properties (descur PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/descur)
target_include_directories(descur PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY}/descur)

if (BLAS_FOUND)
	target_link_libraries (descur PUBLIC ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
	target_link_libraries (descur PUBLIC ${LAPACK_LIBRARIES})
endif ()

get_target_property(DESCUR_LINK_LIB descur LINK_LIBRARIES)
message(STATUS "descur linked libraries are ${DESCUR_LINK_LIB}")

# Define an executable and link all libraries.
add_executable (xdescur ${CMAKE_CURRENT_SOURCE_DIR}/src/descur.f)
target_link_libraries (xdescur descur)

