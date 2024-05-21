cmake_minimum_required (VERSION 3.5)

option(CMAKE_BUILD_TYPE "" "Release")
# set(CMAKE_FIND_DEBUG_MODE 1)
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
#find_package(photospline REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem program_options)
find_package(GSL REQUIRED)
find_package(APFEL REQUIRED)

set(SOURCES
    ${PROJECT_SOURCE_DIR}/src/configuration.cpp
    ${PROJECT_SOURCE_DIR}/src/cross_section.cpp
    ${PROJECT_SOURCE_DIR}/src/phase_space.cpp
    ${PROJECT_SOURCE_DIR}/src/physconst.cpp
    ${PROJECT_SOURCE_DIR}/src/spline_maker.cpp
    ${PROJECT_SOURCE_DIR}/src/structure_function.cpp
    ${PROJECT_SOURCE_DIR}/src/tools.cpp
)

# Create the library
add_library(NuXSSplMkr ${SOURCES})

# Ensure the library is position independent
set_property(TARGET NuXSSplMkr PROPERTY POSITION_INDEPENDENT_CODE ON)

# Link required libraries to NuXSSplMkr
target_link_libraries(NuXSSplMkr ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${LHAPDF_LIBRARIES} ${APFEL_LIBRARIES})

include_directories(${PROJECT_SOURCE_DIR}/include)

if (Boost_FOUND)
    include_directories(${BOOST_INCLUDE_DIR})
else (Boost_FOUND)
    message(FATAL_ERROR "Boost not found!")
endif (Boost_FOUND)

# LHAPDF
find_program(LHAPDF_CONFIG lhapdf-config REQUIRED)
if (LHAPDF_CONFIG)
    execute_process(COMMAND ${LHAPDF_CONFIG} --cflags OUTPUT_VARIABLE LHAPDF_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${LHAPDF_CONFIG} --libs OUTPUT_VARIABLE LHAPDF_LIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)

    exec_program(${LHAPDF_CONFIG}
      ARGS --libdir
      OUTPUT_VARIABLE LHAPDF_LIBDIR
    )
    find_library(LHAPDF_LIBRARY NAMES LHAPDF PATHS ${LHAPDF_LIBDIR})
    set(LHAPDF_LIBRARY ${LHAPDF_LIBRARY} CACHE STRING INTERNAL)

    set(LHAPDF_CXX_FLAGS ${LHAPDF_CXX_FLAGS} CACHE STRING INTERNAL)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LHAPDF_CXX_FLAGS}")
    target_link_libraries(NuXSSplMkr ${LHAPDF_LIBRARY})
else(LHAPDF_CONFIG)
    message(FATAL_ERROR "LHAPDF not found!")
endif(LHAPDF_CONFIG)

if (GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIR})
else (GSL_FOUND)
    message(FATAL_ERROR "GSL not found!")
endif (GSL_FOUND)

find_program(PHOTOSPLINE_CONFIG photospline-config REQUIRED)
if (PHOTOSPLINE_CONFIG)
    execute_process(COMMAND ${PHOTOSPLINE_CONFIG} --cflags OUTPUT_VARIABLE PHOTOSPLINE_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${PHOTOSPLINE_CONFIG} --libs OUTPUT_VARIABLE PHOTOSPLINE_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(PHOTOSPLINE_CXX_FLAGS ${PHOTOSPLINE_CXX_FLAGS} CACHE STRING INTERNAL)
    set(PHOTOSPLINE_LIBRARIES ${PHOTOSPLINE_LIBRARIES} CACHE STRING INTERNAL)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PHOTOSPLINE_CXX_FLAGS}")
    target_link_libraries(NuXSSplMkr ${PHOTOSPLINE_LIBRARIES})
else (PHOTOSPLINE_CONFIG)
    message(FATAL_ERROR "photospline-config not found!")
endif (PHOTOSPLINE_CONFIG)

#find_program(APFEL_CONFIG apfel-config REQUIRED)
#if (APFEL_CONFIG)
#    execute_process(COMMAND ${APFEL_CONFIG} --cppflags OUTPUT_VARIABLE APFEL_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
#    execute_process(COMMAND ${APFEL_CONFIG} --ldflags OUTPUT_VARIABLE APFEL_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)
#    set(APFEL_CXX_FLAGS ${APFEL_CXX_FLAGS} CACHE STRING INTERNAL)
#    set(APFEL_LIBRARIES ${APFEL_LIBRARIES} CACHE STRING INTERNAL)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${APFEL_CXX_FLAGS}")
#    target_link_libraries(NuXSSplMkr ${APFEL_LIBRARIES})
#else (APFEL_CONFIG)
#    message(FATAL_ERROR "APFEL not found!")
#endif (APFEL_CONFIG)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -Wunused")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

add_subdirectory(mains)

install(TARGETS NuXSSplMkr RUNTIME DESTINATION bin LIBRARY DESTINATION lib)
