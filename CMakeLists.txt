# The name of our project is "HELLO". CMakeLists files in this project
# can
# refer to the root source directory of the project as
# ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (LUX)
SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_BUILD_TYPE Release)
ADD_DEFINITIONS(
	    -std=c++0x
			-O2
)
# Recurse into the "Hello" and "Demo" subdirectories. This does not
# actually
# cause another cmake executable to run. The same process will walk
# through
# the project's entire directory structure.
add_subdirectory (src) 
add_subdirectory (include)
