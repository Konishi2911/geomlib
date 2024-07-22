cmake_minimum_required(VERSION 3.20)
project(geomlib VERSION 0.10.2 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

include_directories(
	${CMAKE_SOURCE_DIR}/include
	${CMAKE_SOURCE_DIR}/third_party/lalib/include
	${CMAKE_SOURCE_DIR}/third_party/mathlib/include
)

find_package(GTest)
if (${GTest_FOUND}) 
	enable_testing()
	add_subdirectory(test)
endif()

include(third_party/lalib/cmake/backend.cmake)