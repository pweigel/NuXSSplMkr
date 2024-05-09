cmake_minimum_required (VERSION 3.5)

option(CMAKE_BUILD_TYPE "" "Release")
set(_FLAGS "-O2 -Wall -fPIC")
set(_FLAGS_DEBUG "-g -O0 -Wall -fPIC")
set(_FLAGS_RELEASE "-O2 -Wall -fPIC -s")
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(_FLAGS "${_FLAGS} -stdlib=libc++")
    set(_FLAGS_DEBUG "${_FLAGS_DEBUG} -stdlib=libc++")
    set(_FLAGS_RELEASE "${_FLAGS_RELEASE} -stdlib=libc++")
endif()
set(CMAKE_CXX_FLAGS ${_FLAGS})
set(CMAKE_CXX_FLAGS_DEBUG ${_FLAGS_DEBUG})
set(CMAKE_CXX_FLAGS_RELEASE ${_FLAGS_RELEASE})

# the project name
project(NuXSSplMkr)

# activating some global properties for the project
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(PkgConfig REQUIRED)
find_package(photospline REQUIRED)
if (photospline)
#########################
endif(photospline)
find_package(Boost REQUIRED COMPONENTS filesystem program_options)
find_package(GSL REQUIRED)
find_package(APFEL REQUIRED)

LIST (APPEND SOURCES
    ${PROJECT_SOURCE_DIR}/src/configuration.cpp
    ${PROJECT_SOURCE_DIR}/src/cross_section.cpp
    ${PROJECT_SOURCE_DIR}/src/phase_space.cpp
    ${PROJECT_SOURCE_DIR}/src/physconst.cpp
    ${PROJECT_SOURCE_DIR}/src/structure_function.cpp
    ${PROJECT_SOURCE_DIR}/src/tools.cpp
)

add_library(NuXSSplMkr ${SOURCES})
set_property(TARGET NuXSSplMkr PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(NuXSSplMkr)

include_directories(${PROJECT_SOURCE_DIR}/include)

if (Boost_FOUND)
    include_directories(${BOOST_INCLUDE_DIR})
else (Boost_FOUND)
    message(FATAL_ERROR "Boost not found!")
endif (Boost_FOUND)

# # LHAPDF
find_program(LHAPDF_CONFIG lhapdf-config REQUIRED)
if (LHAPDF_CONFIG)
    execute_process(COMMAND ${LHAPDF_CONFIG} --cflags OUTPUT_VARIABLE LHAPDF_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${LHAPDF_CONFIG} --libs OUTPUT_VARIABLE LHAPDF_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(LHAPDF_CXX_FLAGS ${LHAPDF_CXX_FLAGS} CACHE STRING INTERNAL)
    set(LHAPDF_LIBRARIES ${LHAPDF_LIBRARIES} CACHE STRING INTERNAL)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LHAPDF_CXX_FLAGS}")
    target_link_libraries(NuXSSplMkr ${LHAPDF_LIBRARIES})
else(LHAPDF_CONFIG)
    message(FATAL_ERROR "LHAPDF not found!")
endif(LHAPDF_CONFIG)

if (GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIR})
else (GSL_FOUND)
    message(FATAL_ERROR "GSL not found!")
endif (GSL_FOUND)

find_program(APFEL_CONFIG apfel-config REQUIRED)
if (APFEL_CONFIG)
    execute_process(COMMAND ${APFEL_CONFIG} --cppflags OUTPUT_VARIABLE APFEL_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${APFEL_CONFIG} --ldflags OUTPUT_VARIABLE APFEL_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(APFEL_CXX_FLAGS ${APFEL_CXX_FLAGS} CACHE STRING INTERNAL)
    set(APFEL_LIBRARIES ${APFEL_LIBRARIES} CACHE STRING INTERNAL)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${APFEL_CXX_FLAGS}")
    target_link_libraries(NuXSSplMkr ${APFEL_LIBRARIES})
else (APFEL_CONFIG)
    message(FATAL_ERROR "APFEL not found!")
endif (APFEL_CONFIG)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -Wunused")

# add_subdirectory(src)
# add_subdirectory(mains)

install(TARGETS NuXSSplMkr RUNTIME DESTINATION bin DESTINATION lib)
install(DIRECTORY include/NuXSSplMkr DESTINATION test_include)