LIST (APPEND EXECUTABLES
    ${PROJECT_SOURCE_DIR}/mains/make_sf_splines.cpp
    ${PROJECT_SOURCE_DIR}/mains/make_sf_grid.cpp
    ${PROJECT_SOURCE_DIR}/mains/make_all_sf_splines.cpp
    ${PROJECT_SOURCE_DIR}/mains/calculate_dsdxdy.cpp
    ${PROJECT_SOURCE_DIR}/mains/calculate_dsdy.cpp
    ${PROJECT_SOURCE_DIR}/mains/calculate_xs.cpp
)

foreach(SOURCE_PATH ${EXECUTABLES})
    get_filename_component(EXECUTABLE_NAME ${SOURCE_PATH} NAME_WE)
    add_executable(${EXECUTABLE_NAME} ${SOURCE_PATH})
    target_include_directories(${EXECUTABLE_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)
    target_link_libraries(${EXECUTABLE_NAME} 
        PUBLIC 
        NuXSSplMkr 
        ${LHAPDF_LIBRARIES} 
        ${APFEL_LIBRARIES} 
        ${GSL_LIBRARIES} 
        ${BOOST_LIBRARIES}
    )

    install(TARGETS ${EXECUTABLE_NAME}
        RUNTIME DESTINATION bin
    )

endforeach()

