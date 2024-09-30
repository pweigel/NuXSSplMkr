LIST (APPEND MAIN_SOURCES
    ${PROJECT_SOURCE_DIR}/mains/make_structure_functions.cpp
    ${PROJECT_SOURCE_DIR}/mains/create_lhapdf_structure_functions.cpp
    ${PROJECT_SOURCE_DIR}/mains/grid2spline.cpp
    ${PROJECT_SOURCE_DIR}/mains/calculate_dsdxdy.cpp
    ${PROJECT_SOURCE_DIR}/mains/calculate_dsdxdy_small.cpp
    ${PROJECT_SOURCE_DIR}/mains/calculate_nusquids_dsdy.cpp
    ${PROJECT_SOURCE_DIR}/mains/calculate_nusquids_sigma.cpp
    ${PROJECT_SOURCE_DIR}/mains/calculate_dsdy.cpp
    ${PROJECT_SOURCE_DIR}/mains/calculate_total_xs.cpp
    ${PROJECT_SOURCE_DIR}/mains/evolve_pdf.cpp
)

foreach(SOURCE ${MAIN_SOURCES})
    # Get the file name without the extension
    get_filename_component(EXE_NAME ${SOURCE} NAME_WE)

    # Add executable
    add_executable(${EXE_NAME} ${SOURCE})

    # Link the executable with the NuXSSplMkr library and other required libraries
    target_link_libraries(${EXE_NAME} PRIVATE ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${LHAPDF_LIBRARIES} ${PHOTOSPLINE_LIBRARIES} ${APFEL_LIBRARIES} NuXSSplMkr)

    # Include directories
    target_include_directories(${EXE_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include ${LHAPDF_INCLUDE_DIR} ${APFEL_INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${GSL_INCLUDE_DIR})

    #set_target_properties(${EXE_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

endforeach()