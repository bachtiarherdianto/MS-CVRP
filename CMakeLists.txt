cmake_minimum_required(VERSION 3.16)
project(MS-CVRP)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_VERBOSE_MAKEFILE ON)

option(ENABLE_VERBOSE "Enable verbose output" ON)
option(ENABLE_LOG_OUTPUT "Enable detailed log optimization" OFF)
option(ENABLE_GUIDANCE "Enable proposed guidance" ON)

message("-- Build options")

if(ENABLE_VERBOSE)
    message("--- Verbose output ENABLED")
    add_definitions(-DVERBOSE)
else()
    message("--- Verbose output DISABLED")
endif()

if(ENABLE_LOG_OUTPUT)
    message("--- Detailed log output ENABLED")
    add_definitions(-DLOG_OUTPUT)
else()
    message("--- Detailed log output DISABLED")
endif()

if(ENABLE_GUIDANCE)
    message("--- Feature-based guidance mechanism ENABLED")
    add_definitions(-DGUIDANCE)
else()
    message("--- Feature-based guidance mechanism DISABLED")
endif()

add_executable(MS-CVRP
	main.cpp
	Solver.cpp
	base.cpp
)